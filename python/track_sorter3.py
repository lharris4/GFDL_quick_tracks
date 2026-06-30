#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# 1. Handle the Filename via Command-Line Arguments
if len(sys.argv) < 2:
    print("Error: Missing input filename argument.")
    print("Usage: python3 track_analyzer.py <input_filename.txt>")
    sys.exit(1)

filename = sys.argv[1]
storm_file = filename + '.storm'
hurricane_file = filename + '.hurricane'

tlen = 12
TClen = 6
tm_c = 252

hurThresh = 32.5
TSThresh = 17.5
equatorward_of_lat = 40
maskres = 0.25

# Geographic Basin Configurations
zooms = (
    {'Name': 'World', 'ShortName': 'world', 'corners': (0, -80, 360, 80)}, 
    {'Name': 'Northern Hemisphere', 'ShortName': 'NH', 'corners': (0, 0, 360, 60)},
    {'Name': 'Western Pacific', 'ShortName': 'WPac', 'corners': (100, 0, 180, 60)},
    {'Name': 'Eastern Pacific', 'ShortName': 'EPac', 'corners': (220, 0, 270, 40)},
    {'Name': 'North Atlantic', 'ShortName': 'NAtl', 'corners': (250, 0, 370, 60)},
)

fstorms, fTCsum, fHTsum = {}, {}, {}
nstorm, nTC, nhur = {}, {}, {}
stormcount = {'cy': np.zeros(len(zooms), dtype=np.int32),
              'TC': np.zeros(len(zooms), dtype=np.int32),
              'hur': np.zeros(len(zooms), dtype=np.int32)}

sumfmt = '%05d++++++++++++++++++++%05d   %06d    %12.4f  %4d %4d\n'

for zoom in zooms:
    basin = zoom['ShortName']
    fstorms[basin] = open(f"{filename}.allstorms.{basin}.txt", 'w')
    fTCsum[basin]  = open(f"{filename}.TS.{basin}.txt", 'w')
    fHTsum[basin]  = open(f"{filename}.C15w.{basin}.txt", 'w')
    nstorm[basin], nTC[basin], nhur[basin] = 0, 0, 0

# Open Core Storm and Hurricane Master File Handles for writing
master_hur_fid = open(hurricane_file, 'w')

# Histogram allocations for Density Map Layers
hist_nx, hist_ny = 72, 45
Hcy = np.zeros((hist_nx, hist_ny))
Htc = np.zeros((hist_nx, hist_ny))
Hhur = np.zeros((hist_nx, hist_ny))
genesis = np.zeros((hist_nx, hist_ny))

def get_matching_basins(lon_pt, lat_pt):
    matching = []
    lon_eval = lon_pt + 360.0 if lon_pt < 0 else lon_pt
    for i, zoom in enumerate(zooms):
        c = zoom['corners']
        if (c[0] <= lon_eval <= c[2]) and (c[1] <= lat_pt <= c[3]):
            matching.append((i, zoom['ShortName']))
    return matching

# Load Raw Tracking Matrix Data safely
if not os.path.exists(filename):
    print(f"Error: Input dataset file '{filename}' does not exist on disk.")
    sys.exit(1)

q = np.loadtxt(filename)

# Global Track Map Workspace Layout
fig, ax = plt.subplots(figsize=(12, 6))
ax.set_xlim([0, 360])
ax.set_ylim([-60, 60])
ax.grid(True, linestyle=':', alpha=0.6)

i = 0
while True:
    i += 1
    trackind = np.where(q[:, 5] == i)[0]
    if len(trackind) == 0:
        break
        
    if np.isnan(q[trackind[0], 7]):
        continue
        
    dtdays = 0.25
        
    if len(trackind) >= tlen:
        ax.plot(q[trackind, 6], q[trackind, 7], color='cyan', linewidth=0.5, alpha=0.7)
        
        basins_initial = get_matching_basins(q[trackind[0], 6], q[trackind[0], 7])
        for idx, name in basins_initial:
            stormcount['cy'][idx] += 1
            nstorm[name] += 1
            sumstr = sumfmt % (nstorm[name], nstorm[name], i, np.max(q[trackind, 9]), np.min(q[trackind, 8]), len(trackind))
            fstorms[name].write(sumstr)
            for nt in trackind:
                fstorms[name].write(f"{int(q[nt, 1]):4d} {int(q[nt, 2]):2d} {int(q[nt, 3]):2d} {int(q[nt, 4]):2d}   "
                                    f"{q[nt, 6]:6.2f}   {q[nt, 7]:6.2f}   {q[nt, 8]:6.2f}   {q[nt, 9]:6.2f}\n")

    TC = np.where((q[trackind, 9] >= TSThresh) & (q[trackind, 15] >= tm_c))[0]
    hur = np.where((q[trackind, 9] >= hurThresh) & (q[trackind, 15] >= tm_c))[0]
    wc = np.where(q[trackind, 15] >= tm_c)[0]
    
    if len(wc) < 8 or len(trackind) < tlen:
        continue

    found = 0
    if len(TC) >= TClen:
        lon_raw = q[trackind, 6].copy()
        df = np.append(0, np.diff(lon_raw))
        side = np.cumprod(np.sign(180 - np.abs(df))) < 0
        
        if np.any(side):
            sg = -np.sign(df[side][0])
            lon2 = lon_raw.copy()
            lon3 = lon_raw.copy()
            lon3[side] += 360 * sg
            lon2[~side] -= 360 * sg
            ax.plot(lon2[TC], q[trackind[TC], 7], 'k-', linewidth=2)
            ax.plot(lon3[TC], q[trackind[TC], 7], 'k-', linewidth=2)
        else:
            ax.plot(lon_raw[TC], q[trackind[TC], 7], 'k-', linewidth=2)
            
        ax.plot(lon_raw[TC[0]], q[trackind[TC[0]], 7], 'ko', markersize=5)
        found = 1

    if found > 0:
        g_lon, g_lat = q[trackind[TC[0]], 6], q[trackind[TC[0]], 7]
        active_basins = get_matching_basins(g_lon, g_lat)
        
        for idx, name in active_basins:
            stormcount['TC'][idx] += 1
            nTC[name] += 1
            sumstr = sumfmt % (nTC[name], nTC[name], i, np.max(q[trackind[TC], 9]), np.min(q[trackind[TC], 8]), len(TC))
            fTCsum[name].write(sumstr)
            for nt in trackind[TC]:
                fTCsum[name].write(f"{int(q[nt, 1]):4d} {int(q[nt, 2]):2d} {int(q[nt, 3]):2d} {int(q[nt, 4]):2d}   "
                                   f"{q[nt, 6]:6.2f}   {q[nt, 7]:6.2f}   {q[nt, 8]:6.2f}   {q[nt, 9]:6.2f}\n")

        # 2. Correctly populate and write the Master Hurricane Out File to disk
        if len(hur) > 0 and np.abs(g_lat) <= equatorward_of_lat:
            ax.plot(q[trackind[hur], 6], q[trackind[hur], 7], 'r-', linewidth=2)
            
            # Write structured metadata headers directly to master hurricane log
            master_hur_fid.write(f"STORM ID: {i:05d} | MaxWind: {np.max(q[trackind[hur], 9]):6.2f} m/s\n")
            
            for nt in trackind[hur]:
                row_str = (f"{int(q[nt, 1]):4d} {int(q[nt, 2]):2d} {int(q[nt, 3]):2d} {int(q[nt, 4]):2d}   "
                           f"{q[nt, 6]:6.2f}   {q[nt, 7]:6.2f}   {q[nt, 8]:6.2f}   {q[nt, 9]:6.2f}\n")
                master_hur_fid.write(row_str)
            
            for idx, name in active_basins:
                stormcount['hur'][idx] += 1
                nhur[name] += 1
                sumstr = sumfmt % (nhur[name], nhur[name], i, np.max(q[trackind[hur], 9]), np.min(q[trackind[hur], 8]), len(hur))
                fHTsum[name].write(sumstr)
                for nt in trackind[hur]:
                    fHTsum[name].write(row_str)

        # 3. Correctly write individual localized storm tracking data files to disk
        st_file_path = f"{storm_file}.{i}"
        with open(st_file_path, 'w') as st_fid:
            for nt in trackind:
                st_fid.write(f"{int(q[nt, 1]):4d} {int(q[nt, 2]):2d} {int(q[nt, 3]):2d} {int(q[nt, 4]):2d}   "
                             f"{q[nt, 6]:6.2f}   {q[nt, 7]:6.2f}   {q[nt, 8]:6.2f}   {q[nt, 9]:6.2f}\n")
                             
        ax.text(q[trackind[0], 6] + 2, q[trackind[0], 7], str(i), fontsize=8)

        # Compute Density matrix indices
        Hd, xedges, yedges = np.histogram2d(q[trackind, 6], q[trackind, 7], bins=[hist_nx, hist_ny], range=[[0, 360], [-90, 90]])
        Hcy += Hd * dtdays
        Htc += Hd * dtdays
        if len(hur) > 0:
            Hhur += Hd * dtdays
        
        Hg, _, _ = np.histogram2d([g_lon], [g_lat], bins=[hist_nx, hist_ny], range=[[0, 360], [-90, 90]])
        genesis += Hg

plt.savefig(f"{filename}.tracks_global.png", dpi=150, bbox_inches='tight')

# Density Map Export Function
def plot_density_map(H, title, filename_out, cmap):
    fig_d, ax_d = plt.subplots(figsize=(10, 5))
    X, Y = np.meshgrid(xedges, yedges)
    mesh = ax_d.pcolormesh(X, Y, H.T, cmap=cmap, shading='flat')
    ax_d.set_xlim([0, 360])
    ax_d.set_ylim([-60, 60])
    plt.colorbar(mesh, ax=ax_d)
    ax_d.set_title(title)
    plt.savefig(filename_out, dpi=150, bbox_inches='tight')
    plt.close(fig_d)

plot_density_map(Hcy, 'Cyclone Storm-Days Density', f"{filename}.density_cy.png", 'Greens')
plot_density_map(Htc, 'Tropical Cyclone Storm-Days Density', f"{filename}.density_TC.png", 'Greys')
plot_density_map(Hhur, 'Hurricane Storm-Days Density', f"{filename}.density_hur.png", 'Reds')
plot_density_map(genesis, 'Tropical Genesis Events', f"{filename}.genesis.png", 'Purples')

# Close open log summaries
master_hur_fid.close()
for f_dict in [fstorms, fTCsum, fHTsum]:
    for f in f_dict.values():
        f.close()
plt.close(fig)

print("Analysis Execution Complete.")
