#!/usr/bin/env python3
import os
import sys
import numpy as np
import cartopy.io.shapereader as shpreader
from shapely.geometry import MultiLineString, LineString

def generate_coast_file_via_reader(output_filename='coast.dat', resolution='110m'):
    print(f"--- DEBUG STEP 1: Requesting shapefile download path via Cartopy core ---")
    try:
        # Request the native download and grab the absolute string path on your machine
        shpfilename = shpreader.natural_earth(
            resolution=resolution, category='physical', name='coastline'
        )
        print(f"[SUCCESS] Shapefile located/downloaded successfully.")
        print(f"Target location on disk: {shpfilename}")
    except Exception as e:
        print(f"[CRITICAL ERROR] Failed to acquire or download the Natural Earth file.")
        print(f"Details: {e}")
        sys.exit(1)

    print(f"\n--- DEBUG STEP 2: Instantiating the raw shapefile reader layer ---")
    if not os.path.exists(shpfilename):
        print(f"[ERROR] Path string returned by Cartopy doesn't exist on disk!")
        sys.exit(1)
        
    reader = shpreader.Reader(shpfilename)
    # Grab records to inspect if geometry is actually populated
    records = list(reader.records())
    print(f"Total features identified inside shapefile: {len(records)}")
    
    if len(records) == 0:
        print("[ERROR] Shapefile is empty or corrupted. Please try deleting the cache directory.")
        sys.exit(1)

    print(f"\n--- DEBUG STEP 3: Parsing polygon coordinate vectors ---")
    output_lines = []
    
    for record in records:
        geom = record.geometry
        if geom is None:
            continue

        # FIXED: Natural Earth 'coastline' stores paths as LineStrings, not Polygons!
        if isinstance(geom, LineString):
            lines = [geom]
        elif isinstance(geom, MultiLineString):
            lines = list(geom.geoms)
        else:
            continue
            
        for line in lines:
            # Extract raw coordinate paths from the linestring boundary
            coords = np.array(line.coords)
            
            # Map longitudes cleanly to your [0, 360] framework
            lons = np.mod(coords[:, 0], 360.0)  
            lats = coords[:, 1]
            
            for lon, lat in zip(lons, lats):
                output_lines.append([lon, lat])
            
            # Append a break of NaNs between independent strokes so Matplotlib lifts the pen
            output_lines.append([np.nan, np.nan])
            
    print(f"Total lines prepared for data generation: {len(output_lines)}")

    if len(output_lines) > 0:
        output_array = np.array(output_lines)
        np.savetxt(output_filename, output_array, fmt='%10.4f', header='Longitude   Latitude')
        print(f"\n[SUCCESS] '{output_filename}' generated successfully with {len(output_array)} rows.")
    else:
        print("\n[FAIL] No line coordinates could be mapped from the shapefile features.")

if __name__ == '__main__':
    generate_coast_file_via_reader('coast.dat', resolution='110m')
