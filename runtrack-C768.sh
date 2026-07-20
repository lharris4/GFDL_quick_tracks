#!/bin/bash --norc
#PBS -N runtrack
#PBS -l walltime=36:00:00
#PBS -j oe
#PBS -q batch
#PBS -o /home/lmh/tmp/

#source $HOME/bin/fms/sys/environ.bashrc

module load netcdf-fortran udunits

cd /home/lmh/research/seasonal/quick_tracks/trunk/
set -x #echo
set -e

head=/archive/Lucas.Harris/SHiELD/202407/
if [[ -z ${runname+x} ]] #check if $runname is not set
then
    if [[ $# > 0 ]]
    then
	runname=$1
    else
	runname=20220730.00Z.C768.C768_rt2022.diagsA
    fi
fi

andir=${head}/${runname}/Analysis/track/
mkdir -p $andir/
mkdir -p $andir/stdout
mkdir -p $andir/input
mkdir -p $andir/output

cd ${head}/${runname}/

#The files:


for var in "uas" "vas" "psl" "TMP500_300"
do
    files=`ls history/??????????/${var}_C768_3072x1536.fre.nc`
    dmget $files &
done
wait
filelist=$(ls -1 history/??????????/uas_C768_3072x1536.fre.nc  | tr '\n' ' ' | sed 's/uas/\#/')
	   

inputname="${andir}/input/${runname}-global.nml" 

touch $inputname
#infile needs to be a list
cat > $inputname <<EOF
&nlist
   infile = '$filelist', 
   outfile = '$andir/output/NH.dat',
   lat_s = 5
   lat_e = 50
   lon_s = 0
   lon_e = 360
   ncontours = 1,
   warm_core_check = .T.,
   vort_thresh = 0.
   cint_slp = 2.
   one_variable_per_file = .true.
   dist_threshold = 250.e3
   latlon_grid = .T.

&end

&flist
   name_u_ref = 'uas'
   name_v_ref = 'vas'
   name_SLP = 'psl'
   name_TM = 'TMP500_300'
   name_lon_2d = 'grid_xt' 
   name_lat_2d = 'grid_yt' 
&end
EOF

#THE MAGIC LINE avoids mysterious segfaults
ulimit -s unlimited

#/home/Lucas.Harris/research/seasonal/quick_tracks/trunk/bin/track_checkmissing.exe  $inputname 2>&1 | tee ${andir}/stdout/runtrack-${runname}.out 
#/home/Lucas.Harris/research/seasonal/quick_tracks/trunk/bin/track-debug.exe  $inputname 2>&1 | tee ${andir}/stdout/runtrack-${runname}.out 
/home/Lucas.Harris/research/seasonal/quick_tracks/trunk/bin/track.exe  $inputname 2>&1 | tee ${andir}/stdout/runtrack-${runname}.out || exit

