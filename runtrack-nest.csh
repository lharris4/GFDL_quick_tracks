#!/bin/bash --norc
#PBS -N runtrack-nest
#PBS -l walltime=36:00:00
#PBS -j oe
#PBS -q batch
#PBS -o /home/lmh/tmp/

source $HOME/bin/fms/sys/environ.bashrc

module load ifort netcdf gcp

cd /home/lmh/research/seasonal/quick_tracks/trunk/
#set -x
set -e

head=/archive/Lucas.Harris/HiRAM_subs/
if [[ -z ${runname+x} ]] #check if $runname is not set
then
    if [[ $# > 0 ]]
    then
	runname=$1
    else
	runname=C384n3-GHS_fore_nh_0313_100701
    fi
fi
#runname=
#runname=C384n3-GHS_subfore_nh_0313_100701

andir=${head}/${runname}/Analysis/track/
mkdir -p $andir/
mkdir -p $andir/stdout
mkdir -p $andir/input
mkdir -p $andir/output

cd ${head}/${runname}/

files=`ls --color=never -1 history/????????/????????.atmos_4xdaily.nest02.nc`
dmget $files
nfiles=`echo $files | wc -w`
files=`ls --color=never -1 $files  | sed 's/^\(.*\)$/"\1"/' | tr '\n' ','`

#Corrected grid_spec for GHS --- lmh 29 jan 16
#Note that for a different grid this gridfile will need to be changed
gridfile=/archive/Lucas.Harris/model/input/grid/C384-nest-GHS-lnoref/grid_spec.nest02.nc
inputname="${andir}/input/${runname}-nest.nml" 

touch $inputname
cat > $inputname <<EOF
&nlist
    infile = $files
    outfile = '${andir}/output/${runname}-nest.dat',
    latlon_grid = .false.
    periodicx = .false.
    grid_file = "$gridfile"
    ncontours = 1,
    warm_core_check = .true.
    cint_slp = 2.
&end
&flist
    name_lon_2d = 'grid_lont' ! revised --- lmh 29 jan 16
    name_lat_2d = 'grid_latt' ! revised --- lmh 29 jan 16
&end
EOF

/home/Lucas.Harris/research/seasonal/quick_tracks/trunk/bin/track.exe  $inputname  2>&1 | tee ${andir}/stdout/runtrack-${runname}-nest.out 

#set +x
module purge
module load gcc netcdf/4.2 python/2.7.3

python /home/lmh/research/seasonal/quick_tracks/trunk/python/track_sorter_z1y.py -h 29.5  ${andir}/output/${runname}-nest.dat &
#\rm ${andir}/output/${runname}-nest.dat*{EPac,NH,world,WPac}*

wait

