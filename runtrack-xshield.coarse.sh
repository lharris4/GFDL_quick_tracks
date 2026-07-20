#!/bin/bash --norc
#PBS -N runtrack
#PBS -l walltime=36:00:00
#PBS -j oe
#PBS -q batch
#PBS -o /home/lmh/tmp/

#source $HOME/bin/fms/sys/environ.bashrc

module load netcdf-fortran udunits

cd /home/lmh/research/seasonal/quick_tracks/trunk/
#set -x #echo
set -e

head=/archive/Lucas.Harris/SHiELD/202407/
if [[ -z ${runname+x} ]] #check if $runname is not set
then
    if [[ $# > 0 ]]
    then
	runname=$1
    else
	runname=20130101.00Z.C3072.xs24v2.amip
    fi
fi
year=2013

andir=${head}/${runname}/Analysis/track/
mkdir -p $andir/
mkdir -p $andir/stdout
mkdir -p $andir/input
mkdir -p $andir/output

cd ${head}/${runname}/

#The files:


for var in "us_coarse" "vs_coarse" "psl" "TMP500_300"
do
    files=`ls history/${year}??????/${var}_C3072_1440x720.fre.nc`
    dmget $files &
done
wait
#filelist=$(ls -1 history/${year}??????/us_coarse_C3072_1440x720.fre.nc  | tr '\n' ',' | sed 's/us_coarse/\#/g')
#Expand the glob into a standard bash array
files=( history/${year}??????/us_coarse_C3072_1440x720.fre.nc )
# Perform string substitution on all array items: replace 'us' with '#'
# The syntax "${array[@]/search/replace}" modifies every element in place.
substituted_files=( "${files[@]/us_coarse/\#}" )
# Format into a comma-separated list of quoted strings for the namelist
filelist=$(printf "'%s'\n" "${substituted_files[@]}" | paste -sd, -)
echo "$filelist"

inputname="${andir}/input/${runname}-global.nml" 

touch $inputname
#infile needs to be a list
cat > $inputname <<EOF
&nlist
   infile = $filelist, !no comma at end
   outfile = '$andir/output/${year}.coarse.dat',
   ncontours = 1,
   lat_s = -50.
   lat_e = 50.
   max_cyrad = 1200.
   warm_core_check = .T.,
   vort_thresh = 0.
   cint_slp = 2.
   one_variable_per_file = .true.
   dist_threshold = 250.e3
   latlon_grid = .T.

&end

&flist
   name_u_ref = 'us_coarse'
   name_v_ref = 'vs_coarse'
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

