#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# =====================================================================
# 1. SETUP & CONFIGURATION
# =====================================================================
if len(sys.argv) < 2:
    print("Error: Missing input filename argument.")
    print("Usage: python3 track_analyzer.py <input_filename.txt>")
    sys.exit(1)

filename = sys.argv[1]
storm_file = filename + '.storm'
hurricane_file = filename + '.hurricane'

tlen, TClen, tm_c = 12, 6, 252
hurThresh, TSThresh, equatorward_of_lat = 32.5, 17.5, 40
maskres = 0.25

zooms = (
    {'Name': 'World', 'ShortName': 'world', 'corners': (0, -80, 360, 80)}, 
    {'Name': 'Northern Hemisphere', 'ShortName': 'NH', 'corners': (0, 0, 360, 60)},
    {'Name': 'Western Pacific', 'ShortName': 'WPac', 'corners': (100, 0, 180, 60)},
    {'Name': 'Eastern Pacific', 'ShortName': 'EPac', 'corners': (220, 0, 270, 40)},
    {'Name': 'North Atlantic', 'ShortName': 'NAtl', 'corners': (250, 0, 370, 60)},
)

# Dynamic file handles and counters initialized cleanly via dictionary comprehensions
fstorms = {z['ShortName']: open(f"{filename}.allstorms.{z['ShortName']}.txt", 'w') for z in zooms}
fTCsum  = {z['ShortName']: open(f"{filename}.TS.{z['ShortName']}.txt", 'w') for z in zooms}
fHTsum  = {z['ShortName']: open(f"{filename}.C15w.{z['ShortName']}.txt", 'w') for z in zooms}

nstorm = {z['ShortName']: 0 for z in zooms}
nTC    = {z['ShortName']: 0 for z in zooms}
nhur   = {z['ShortName']: 0 for z in zooms}

stormcount = {
    'cy': np.zeros(len(zooms), dtype=np.int32),
    'TC': np.zeros(len(zooms), dtype=np.int32),
    'hur': np.zeros(len(zooms), dtype=np.int32)
}

hist_nx, hist_ny = 72, 45
Hcy = np.zeros((hist_nx, hist_ny))
Htc = np.zeros((hist_nx, hist_ny))
Hhur = np.zeros((hist_nx, hist_ny))
genesis = np.zeros((hist_nx, hist_ny))

# =====================================================================
# 2. HELPER SUBROUTINES (CONSOLIDATED LOGIC)
# =====================================================================

def get_matching_basins(lon_pt, lat_pt):
    matching = []
    lon_eval = lon_pt + 360.0 if lon_pt < 0 else lon_pt
    for i, zoom in enumerate(zooms):
        c = zoom['corners']
        if (c[0] <= lon_eval <= c[2]) and (c[1] <= lat_pt <= c[3]):
            matching.append((i, zoom['ShortName']))
    return matching

def write_storm_data(file_handle, count, storm_id, max_w, min_p, track_indices):
    """Centralized formatting block for header lines and tracking step records."""
    sumfmt = '%05d++++++++++++++++++++%05d   %06d    %12.4f  %4d %4d\n'
    file_handle.write(sumfmt % (count, count, storm_id, max_w, min_p, len(track_indices)))
    for nt in track_indices:
        file_handle.write(f"{int(q[nt, 1]):4d} {int(q[nt, 2]):2d} {int(q[nt, 3]):2d} {int(q[nt, 4]):2d}   "
                          f"{q[nt, 6]:6.2f}   {q[nt, 7]:6.2f}   {q[nt, 8]:6.2f}   {q[nt, 9]:6.2f}\n")

def log_to_basins(basin_matches, count_dict, file_dict, key, storm_id, max_w, min_p, indices):
    """Dynamically updates statistics counts and logs output lines for matching basins."""
    for idx, name in basin_matches:
        stormcount[key][idx] += 1
        count_dict[name] += 1
        write_storm_data(file_dict[name], count_dict[name], storm_id, max_w, min_p, indices)

def plot_wrapped_lines(ax, lons, lats, fmt, linewidth=1, alpha=1.0, zorder=2):
    lons_wrapped = np.mod(lons, 360.0)
    all_breaks = np.unique(np.concatenate([
        np.where(np.abs(np.diff(lons_wrapped)) > 180.0)[0],
        np.where(np.isnan(lons_wrapped) | np.isnan(lats))[0]
    ]))
    
    start_idx = 0
    for brk in all_breaks:
        sub_lon, sub_lat = lons_wrapped[start_idx:brk + 1], lats[start_idx:brk + 1]
        valid = ~np.isnan(sub_lon) & ~np.isnan(sub_lat)
        if np.any(valid):
            ax.plot(sub_lon[valid], sub_lat[valid], fmt, linewidth=linewidth, alpha=alpha, zorder=zorder)
        start_idx = brk + 1
        
    if start_idx < len(lons_wrapped):
        valid = ~np.isnan(lons_wrapped[start_idx:]) & ~np.isnan(lats[start_idx:])
        if np.any(valid):
            ax.plot(lons_wrapped[start_idx:][valid], lats[start_idx:][valid], fmt, linewidth=linewidth, alpha=alpha, zorder=zorder)

# =====================================================================
# 3. CORE FILE PROCESSING & VISUALIZATION
# =====================================================================
if not os.path.exists(filename):
    print(f"Error: Input dataset file '{filename}' does not exist on disk.")
    sys.exit(1)

q = np.loadtxt(filename)

date_text = "Analysis Period: Data Timeframe Unspecified"
if len(q) > 0:
    date_text = f"Analysis Period: {int(q[0,1]):04d}-{int(q[0,2]):02d}-{int(q[0,3]):02d} to {int(q[-1,1]):04d}-{int(q[-1,2]):02d}-{int(q[-1,3]):02d}"

fig, ax = plt.subplots(figsize=(12, 6))
ax.set_xlim([0, 360])
ax.set_ylim([-60, 60])
ax.set_aspect('equal', adjustable='box')
ax.set_xticks(np.arange(0, 361, 60))
ax.set_yticks(np.arange(-60, 61, 20))
ax.grid(True, linestyle=':', alpha=0.5, color='gray')
ax.text(0.0, -.12, date_text, transform=ax.transAxes, fontsize=9, color='dimgray', ha='left', va='top')

if os.path.exists('coast.dat'):
    coast = np.loadtxt('coast.dat')
    plot_wrapped_lines(ax, coast[:, 0], coast[:, 1], fmt='k-', linewidth=0.5, alpha=1.0, zorder=1)

master_hur_fid = open(hurricane_file, 'w')

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
        plot_wrapped_lines(ax, q[trackind, 6], q[trackind, 7], fmt='c-', linewidth=0.5, alpha=0.7, zorder=2)
        basins_initial = get_matching_basins(q[trackind[0], 6], q[trackind[0], 7])
        log_to_basins(basins_initial, nstorm, fstorms, 'cy', i, np.max(q[trackind, 9]), np.min(q[trackind, 8]), trackind)

    TC = np.where((q[trackind, 9] >= TSThresh) & (q[trackind, 15] >= tm_c))[0]
    hur = np.where((q[trackind, 9] >= hurThresh) & (q[trackind, 15] >= tm_c))[0]
    wc = np.where(q[trackind, 15] >= tm_c)[0]
    
    if len(wc) < 8 or len(trackind) < tlen:
        continue

    found = 0
    if len(TC) >= TClen:
        plot_wrapped_lines(ax, q[trackind[TC], 6], q[trackind[TC], 7], fmt='k-', linewidth=2, zorder=3)
        ax.plot(np.mod(q[trackind[TC], 6], 360.0)[0], q[trackind[TC[0]], 7], 'ko', markersize=5, zorder=4)
        found = 1

    if found > 0:
        g_lon, g_lat = q[trackind[TC[0]], 6], q[trackind[TC[0]], 7]
        active_basins = get_matching_basins(g_lon, g_lat)
        
        # Log to Tropical Storm basin files
        log_to_basins(active_basins, nTC, fTCsum, 'TC', i, np.max(q[trackind[TC], 9]), np.min(q[trackind[TC], 8]), trackind[TC])

        # Log to Hurricane basin files & Master Hurricane log
        if len(hur) > 0 and np.abs(g_lat) <= equatorward_of_lat:
            plot_wrapped_lines(ax, q[trackind[hur], 6], q[trackind[hur], 7], fmt='r-', linewidth=2, zorder=3)
            ax.text(np.mod(q[trackind[hur], 6], 360.0)[0] + 2, q[trackind[hur], 7][0], str(i), fontsize=8, zorder=5)
            
            master_hur_fid.write(f"STORM ID: {i:05d} | MaxWind: {np.max(q[trackind[hur], 9]):6.2f} m/s\n")
            log_to_basins(active_basins, nhur, fHTsum, 'hur', i, np.max(q[trackind[hur], 9]), np.min(q[trackind[hur], 8]), trackind[hur])
        else:
            ax.text(np.mod(q[trackind, 6], 360.0)[0] + 2, q[trackind[0], 7], str(i), fontsize=8, zorder=5)

        # Write out individual localized storm profile record to disk
        with open(f"{storm_file}.{i}", 'w') as st_fid:
            for nt in trackind:
                st_fid.write(f"{int(q[nt, 1]):4d} {int(q[nt, 2]):2d} {int(q[nt, 3]):2d} {int(q[nt, 4]):2d}   "
                             f"{q[nt, 6]:6.2f}   {q[nt, 7]:6.2f}   {q[nt, 8]:6.2f}   {q[nt, 9]:6.2f}\n")

        Hd, xedges, yedges = np.histogram2d(q[trackind, 6], q[trackind, 7], bins=[hist_nx, hist_ny], range=[[0, 360], [-90, 90]])
        Hcy += Hd * dtdays
        Htc += Hd * dtdays
        if len(hur) > 0:
            Hhur += Hd * dtdays
        
        Hg, _, _ = np.histogram2d([g_lon], [g_lat], bins=[hist_nx, hist_ny], range=[[0, 360], [-90, 90]])
        genesis += Hg

plt.savefig(f"{filename}.tracks_global.png", dpi=150, bbox_inches='tight')

# =====================================================================
# 4. EXPORT MAP DENSITY LAYERS
# =====================================================================
def plot_density_map(H, title, filename_out, cmap):
    fig_d, ax_d = plt.subplots(figsize=(12, 6))
    X, Y = np.meshgrid(xedges, yedges)
    mesh = ax_d.pcolormesh(X, Y, H.T, cmap=cmap, shading='flat', zorder=1)
    
    ax_d.set_xlim([0, 360])
    ax_d.set_ylim([-60, 60])
    ax.set_aspect('equal', adjustable='box')
    ax_d.set_xticks(np.arange(0, 361, 60))
    ax_d.set_yticks(np.arange(-60, 61, 20))
    ax_d.grid(True, linestyle=':', alpha=0.5, color='gray')
    
    if os.path.exists('coast.dat'):
        plot_wrapped_lines(ax_d, coast[:, 0], coast[:, 1], fmt='k-', linewidth=0.5, alpha=1.0, zorder=2)
        
    ax_d.text(0.0, -0.12, date_text, transform=ax_d.transAxes, fontsize=9, color='dimgray', ha='left', va='top')
    plt.colorbar(mesh, ax=ax_d, pad=0.02, shrink=0.8)
    ax_d.set_title(title, pad=15)
    plt.savefig(filename_out, dpi=150, bbox_inches='tight')
    plt.close(fig_d)

plot_density_map(Hcy, 'Cyclone Storm-Days Density', f"{filename}.density_cy.png", 'Greens')
plot_density_map(Htc, 'Tropical Cyclone Storm-Days Density', f"{filename}.density_TC.png", 'Greys')
plot_density_map(Hhur, 'Hurricane Storm-Days Density', f"{filename}.density_hur.png", 'Reds')
plot_density_map(genesis, 'Tropical Cyclone Genesis Events', f"{filename}.genesis.png", 'Purples')

master_hur_fid.close()
for f_dict in [fstorms, fTCsum, fHTsum]:
    for f in f_dict.values():
        f.close()
plt.close(fig)

print("Analysis Execution Complete.")
