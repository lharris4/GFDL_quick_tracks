#!/usr/bin/env python3
import os
import sys
import numpy as np
from datetime import datetime

#This file reads TRACKER .dat output

# =====================================================================
# CONFIGURATION SETTINGS
# =====================================================================
# Absolute directory path pointing to your high-res native model runs
NATIVE_HEAD = "/archive/Lucas.Harris/SHiELD/202407/20130101.00Z.C3072.xs24v2.amip/history/"
FILE_SUFFIX = "_C3072_11520x5760.fre.nc"

SEARCH_RADIUS_DEG = 2.0  # Spatial bounding box around the coarse center

try:
    from netCDF4 import Dataset, date2index
except ImportError:
    print("[CRITICAL] Python package 'netCDF4' is missing. Run: pip install netCDF4")
    sys.exit(1)

# Handle Raw Tracker Input file via arguments
if len(sys.argv) < 2:
    print("Usage: python3 refine_tracks.py <tracker_output.dat>")
    print("Example: python3 refine_tracks.py 2013.try.coarse.dat")
    sys.exit(1)

input_file = sys.argv[1]
if input_file.endswith(".dat"):
    output_file = input_file[:-4] + ".refine.dat"
else:
    output_file = input_file + ".refine.dat"

if not os.path.exists(input_file):
    print(f"Error: Master tracker file '{input_file}' not found.")
    sys.exit(1)

# =====================================================================
# TIMELINE BUILDER SUBROUTINES
# =====================================================================
print("Scanning history segments to build active interval timelines...")
segment_datetimes = []
if os.path.exists(NATIVE_HEAD):
    for folder in os.listdir(NATIVE_HEAD):
        if len(folder) == 10 and folder.isdigit():
            try:
                segment_datetimes.append(datetime.strptime(folder, "%Y%m%d%H"))
            except ValueError:
                pass
segment_datetimes.sort()

def find_correct_segment(track_time, segments):
    if not segments:
        return None
    if track_time <= segments[0]:
        return segments[0].strftime("%Y%m%d%H")
    for i in range(len(segments) - 1):
        if segments[i] < track_time <= segments[i+1]:
            return segments[i].strftime("%Y%m%d%H")
    return segments[-1].strftime("%Y%m%d%H")

# =====================================================================
# CORE ITERATION & PARSING STAGE (Tracker Raw Output Format)
# =====================================================================
print(f"Processing raw tracker data from: {input_file}")
print(f"Refined raw output will be written to: {output_file}")
print(f"Looking for files in {NATIVE_HEAD}")

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line_idx, line in enumerate(infile, 1):
        parts = line.split()
        
        # Raw lines must contain at least the 10 core metrics
        if len(parts) < 10:
            if parts:  # Log unexpected formatting lines if they appear
                print(f"[WARNING] Skipping poorly formatted line {line_idx}: {line.strip()}")
            outfile.write(line)
            continue
            
        # Parse tracking row variables based on 1-indexed specification mapping:
        # 0=LineNo, 1=Year, 2=Month, 3=Day, 4=Hour, 5=StormIdx, 6=Lon, 7=Lat, 8=MinPres, 9=MaxWind
        try:
            year, month, day, hour = int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])
            coarse_lon, coarse_lat = float(parts[6]), float(parts[7])
            wind = float(parts[9])
        except ValueError:
            print(f"[WARNING] Skipping unparseable text line {line_idx}: {line.strip()}")
            outfile.write(line)
            continue

        track_time = datetime(year, month, day, hour)
        segment_dir = find_correct_segment(track_time, segment_datetimes)
        
        # Fallback to coarse tracker value if file isn't found
        refined_wind = wind
        
        if segment_dir is None:
            print(f"[ERROR] Native Resolution input directory could not be resolved for line {line_idx} ({track_time.strftime('%Y-%m-%d %H:00')})")
        else:
            u_file = os.path.join(NATIVE_HEAD, segment_dir, f"uas{FILE_SUFFIX}")
            v_file = os.path.join(NATIVE_HEAD, segment_dir, f"vas{FILE_SUFFIX}")
            
            if not os.path.exists(u_file):
                print(f"[ERROR] Missing high-resolution variable input file on disk: {u_file}")
            elif not os.path.exists(v_file):
                print(f"[ERROR] Missing high-resolution variable input file on disk: {v_file}")
            else:
                try:
                    with Dataset(u_file, 'r') as nc_u, Dataset(v_file, 'r') as nc_v:
                        native_lons = np.mod(nc_u.variables['grid_xt'][:], 360.0)
                        native_lats = nc_u.variables['grid_yt'][:]
                        target_lon = np.mod(coarse_lon, 360.0)
                        
                        lon_idx = np.where(np.abs(native_lons - target_lon) <= SEARCH_RADIUS_DEG)[0]
                        lat_idx = np.where(np.abs(native_lats - coarse_lat) <= SEARCH_RADIUS_DEG)[0]
                        
                        if len(lon_idx) > 0 and len(lat_idx) > 0:
                            lon_start, lon_end = lon_idx[0], lon_idx[-1] + 1
                            lat_start, lat_end = lat_idx[0], lat_idx[-1] + 1
                            
                            time_var = nc_u.variables['time']
                            time_idx = date2index(track_time, time_var, select='nearest')
                            
                            u_patch = nc_u.variables['uas'][time_idx, lat_start:lat_end, lon_start:lon_end]
                            v_patch = nc_v.variables['vas'][time_idx, lat_start:lat_end, lon_start:lon_end]
                            
                            wind_speed_patch = np.sqrt(u_patch**2 + v_patch**2)
                            refined_wind = np.nanmax(wind_speed_patch)

                except Exception as e:
                    print(f"[ERROR] Failed to query high-resolution variables at line {line_idx}. Error: {e}")

        # Update ONLY the 10th column (Maximum Wind, 0-indexed position 9)
        parts[9] = f"{refined_wind:.2f}"
        
        # Create a universally clean, decimal-aligned row format.
        # This mirrors the tracker spacing while ensuring your sorter script
        # has zero trouble reading the data.
        formatted_line = (
            f"{int(parts[0]):10d}"       # 1. Line Number
            f"{int(parts[1]):8d}"        # 2. Year
            f"{int(parts[2]):4d}"        # 3. Month
            f"{int(parts[3]):4d}"        # 4. Day
            f"{int(parts[4]):4d}"        # 5. Hour
            f"{int(parts[5]):8d}"        # 6. Storm Index ID
            f"{float(parts[6]):14.2f}"   # 7. Center Lon
            f"{float(parts[7]):14.2f}"   # 8. Center Lat
            f"{float(parts[8]):14.2f}"   # 9. Minimum Pressure
            f"{float(parts[9]):14.2f}"   # 10. Maximum Wind Speed (Refined!)
        )
        
        # Dynamically append any remaining miscellaneous diagnostics 
        # (columns 11 to 16+) so nothing is lost, keeping them decimal aligned.
        for extra_part in parts[10:]:
            try:
                # Try formatting as a clean decimal float
                formatted_line += f"{float(extra_part):14.2f}"
            except ValueError:
                # Fallback safely if it's the scientific E-notation or integer string
                formatted_line += f"{extra_part:>14s}"
                
        outfile.write(formatted_line + "\n")

print("Processing Complete.")
