# Quick Tracks: The simple cyclone tracker

Quick Tracks (QT) is intended to be a simple and efficient means of objectively identifying tropical and extratropical cyclones in gridded atmospheric data, typically either model or reanalysis output. QT is especially useful for identifying large numbers of cyclones in long datasets to produce climatologies or statistical seasonal/annual forecasts. QT is used heavily within GFDL to track TCs in the [SHiELD models](www.gfdl.noaa.gov/shield) or in [SPEAR](www.gfdl.noaa.gov); and previously in HiRAM and FLOR. A description of the method is given in the appendix of [Harris, Lin, and Tu (JClim, 2016)](https://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-15-0389.1). 

In this context, "Quick" can be more considered to refer to the speed at which the tracker can be set up, not necessarily execution speed.

**Caveats:** QT is not intended to track weak, disorganized systems over their entire lifetime. For this purpose a more comprehensive, rigorous tracker is recommended, such as Tim Marchok's GFDL Vortex Tracker. As with any objective cyclone identification method QT may not correctly identify 100% of all cyclones, nor will all identified disturbances necessarily be coherent cyclones. QT may also not be quick enough.

**NOTE** GFDL can only provide minimal support for this tracker.

## Compilation and installation

**Requirements:** A Fortran 90 compiler and the NetCDF and UDUnits libraries are required to compile and run the tracker executable. The post-processing scripts require Python 2.6 with the NumPy and Pandas libraries. A source of NetCDF formatted input data is needed to run QT; examples include the GFDL models or the MERRA reanalysis, both of which are freely available for download.

To compile (using the Intel Fortran compiler in this example):

          `ifort -O2 src/tracker.F90 -o bin/track.exe -L/usr/local/netcdf4/lib -L/usr/local/hdf5/lib -L/opt/intel/icc/11.1.073/lib/intel64 -L/usr/local/udunits-1.12.11/lib/ -I/usr/local/udunits-1.12.11/include/ -lnetcdf -lnetcdff -lhdf5 -lhdf5_hl -lz -limf -ludunits -lm`

Your command line may vary based on your choice of compiler and system configuration. To compile with missing-value checking enabled, add `-DCHECK_MISSING`.

See the file 'compile' for more compile commands.

## Running QT

Several scripts have been supplied to be run with sample output from the GFDL High Resolution Atmosphere Model (HiRAM). The script can run on both lat-lon global domains and on native-grid nests, at any time frequency, although the elapsed time between outputs should be no longer than six hours. Add the path to the directory containing the input files, and run the scripts to run QT on the data. 

To run on a different dataset, edit the entries in the nlist namelist in the runscript to those in your input file.

## Post-processing QT

After running the tracker a text file output/####.txt will be created. The file contains the tracked cyclone data. The script `track_sorter_z1y.py` will then filter the tracks to only count cyclones which satisfy a series of requirements, the 'long-lived tropical cyclones' of Chen and Lin (2012,2013), as well as counting the number of hurricane-strength cyclones. The script then plots the tracks, in a number of pre-defined basins, produces density plots of the tracks, and produces a text file containing a summarized report of every identified cyclone.

Run the script using:

`python track_sorter_z1y.py output/####.txt'

If you used the runscript to run the tracker then this script has been run automatically.

If you want only certain dates:

`python scripts/track_sorter_z1y.py -s '2005-08-01' -e '2005-08-31' output/####.txt'

(Either the start date or ending date can be omitted.)

The "storms" summary file can be then used as input to the intercomparison scripts in the scripts/ directory. In particular trackstat.py will make comparisons to IBTrACS observed tropical cyclone data.

## QT input options

There is a host of command line and namelist features that can be applied.

## Neat features

Tracker:

- Land/sea mask
- Non-uniform grids
- Restricted track domain

## Tracker verification statistics

TBA.

## Limitations

The track sorter uses the Pandas libraries, which was a real pain to use, and results in undecipherable code. 

## Technical tracker description

**Credits:** QT was written by Lucas Harris. QT was heavily influenced by the comprehensive GFDL tracker by Tim Marchok, and the TStorm tracker of Joe Sirutis and Ming Zhao. Parts of the tracker code base are originally from Shian-Jiann Lin. Further testing and feedback by Jan-Huey Chen and ChiaYing Tu (Academia Sinica).  All personnel are from GFDL unless otherwise denoted.

Please send questions and comments to Lucas.Harris@noaa.gov .
