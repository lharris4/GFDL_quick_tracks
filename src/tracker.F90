! Compiling on analysis:
!    module load ifort
!    module load netcdf
!    # Full speed run:
!    ifort -O2 src/tracker.F90 -o bin/track.exe -L/usr/local/netcdf4/lib -L/usr/local/hdf5/lib -L/opt/intel/icc/11.1.073/lib/intel64 -L/usr/local/udunits-1.12.11/lib/ -I/usr/local/udunits-1.12.11/include/ -lnetcdf -lnetcdff -lhdf5 -lhdf5_hl -lz -limf -ludunits -lm
!    # Debug run:
!    ifort -O2 src/tracker.F90 -o bin/track.exe -L/usr/local/netcdf4/lib -L/usr/local/hdf5/lib -L/opt/intel/icc/11.1.073/lib/intel64 -L/usr/local/udunits-1.12.11/lib/ -I/usr/local/udunits-1.12.11/include/ -lnetcdf -lnetcdff -lhdf5 -lhdf5_hl -lz -limf -ludunits -lm -traceback -check bounds

!Proposed to-do items:
! 1. Better track connection algorithm, that should be independent of the numerical ordering of storms -- DONE: added 14feb13
! 2. Better warm core detection, along the lines of Bob Hart's
!    algorithm that uses relative temperatures instead of an absolute threshold -- DONE 15mar12
! 3. Faster flood fill algorithm (possible?)
! 4. More accurate trajectory calculation:
!    a. Use 3d vector wind then project onto sphere (as in fv_update_phys.F90 in the dycore)
!    b. Higher-order trajectory computation (ie 2nd order Runge Kutta)
! 5. Land/sea mask, particularly useful for hurricanes -- DONE mar12 (no longer works as of 2014?)
! 6. Output metadata with ASCII output?
! 7. Output NetCDF file showing regions identified as certain cyclones? -- DONE 13oct11
! 8. Write out a proper x- and y- coordinate so GrADS can read in the .nc output. -- DONE dec14
! 9. Run tracker on individual tile, or on a nested grid -- DONE mar13
!10. Should be able to detect lon, lat, time automatically for a x-y-t variable
!11. Vorticity tracking, implemented by Ethan Coffel; need to merge.
!12. Determine number of files from infile, and disregard nfiles -- DONE oct13
!13. Julian calendar fix for dates prior to adoption of Gregorian calendar (need a more permanent fix)
!14. Each variable in a different file. -- DONE early 2014
!15. Compute actual distance from origin in flood_fill instead of counting cells -- DONE 11jun14
!16. Do the same for the warm-core offset -- DONE mid 2014, debugged/tested 2jan15
!17. Total rainfall computation -- DONE 2dec14
!18. OpenMP?

! 7may13: Added code to check for missing values. This however does slow down the tracker about 10%.
! 10oct13: started cleaning up code in preparation for public release.

     program main
      implicit         none

#include <netcdf.inc>

      integer, parameter :: MAXSTORM = 600
      integer, parameter :: NUMFIELDS = 13
      real, parameter :: RADIUS = 6371.e3 !m
      integer, parameter :: MAXFILES = 400
      real, parameter :: GRAV = 9.80
      real, parameter :: rdgas = 287.04
      real, parameter :: rvgas = 461.50
      integer, parameter :: MISSING_VALUE = -999

      integer:: im, jm, km
      integer:: nt
      integer:: im_avg, jm_avg
      integer:: im_t, jm_t
      integer:: im_r, jm_r
      integer:: im_g, jm_g
      integer:: ntimes, ntimes_avg
      real*4, allocatable:: lon(:,:), lat(:,:), area(:,:), lev(:), ilev(:), peln(:)
      real*4, allocatable :: timeval(:) !*4 or *8 ?
      real*4, allocatable :: timeval_avg(:,:) 
      real*8 :: time8
      real*4, allocatable:: slp(:,:), tm(:,:), vort(:,:), wk(:,:), slp_s(:,:)
      real*4, allocatable :: zsurf(:,:), T(:,:,:), precip(:,:)
      real*4, allocatable:: u_ref(:,:), v_ref(:,:)     ! 10-m winds
      real*4, allocatable:: u_steer(:,:), v_steer(:,:)     ! steering 500 mb winds
      integer, allocatable:: id(:,:), id_out(:,:), centers(:,:), id_warm(:,:)
      real*4, allocatable :: topo(:,:), land(:,:)
      logical, allocatable :: landtopomask_wind(:,:)
      integer, allocatable :: regions(:,:)
      real, allocatable :: regions_wind(:,:), regions_snow(:,:)

      integer, allocatable :: warm_core_center(:,:)

      real*4, allocatable :: lon_in(:), lat_in(:)

      !Storms found at a particular time
      integer, allocatable :: allstorms(:,:), allstorms_m1(:,:), allstorms_m2(:,:)
      real, allocatable :: allstormval(:,:), allstormval_m1(:,:), allstormval_m2(:,:)
      logical, allocatable :: claimed(:)
      integer :: stormids(MAXSTORM)

      integer:: fid, fidavg, fidtopo, fidrgn, fidgrid
      integer:: i, j, n, ii, jj

      real*4:: speed
      real:: wind, dist
      integer:: dt_hr = MISSING_VALUE
      integer:: year=MISSING_VALUE, month=MISSING_VALUE, day=MISSING_VALUE, hour=MISSING_VALUE, tmin
      real :: tsec
      integer:: year_prev=2010, month_prev=12, day_prev=1, hour_prev=0, tmin_prev
      real :: tsec_prev
      integer:: year_avg, month_avg, day_avg, hour_avg, tmin_avg
      real :: tsec_avg
      integer:: nstorm, nstorm_m1, nstorm_m2, nstormtot, nprev
      integer:: trackcount=1, nfile
      real :: distmin
      integer :: minstorm, prevstorm, np2, dum(2)
      real :: extp(2)
      real :: tm_c = MISSING_VALUE
      real*4::pi
      real :: zvir, ginv
      character*50 :: region_format

#ifdef CHECK_MISSING
      real*4 :: slp_fillvalue = -1.e10
      real*4 :: vort_fillvalue = -1.e10
      real*4 :: tm_fillvalue = -1.e10
      real*4 :: u_fillvalue = 1.e15
      real*4 :: v_fillvalue = 1.e15
#endif

      character*400 :: grid_file = ''
      logical :: debug = .false.

      character*400 :: fname_tmp

      logical :: use_grid_file = .false.

      integer :: max_cycells_warm_offset
      real :: xdist, rad


      real    :: dx_avg, dphi_avg
      real    :: rdx_avg, rdphi_avg

      integer :: status, varid
      character(len=40) :: units, lonname, latname, dumunit

!!!!! nlist options: tracker runtime options
      character*400 ::  infile(MAXFILES)               !! Input file names; assumed to all have the same x and y
                                                       !! dimensions and the same dimension and variable names,
                                                       !! but time dimension may differ.

      character*400 ::  infile_avg(MAXFILES)           !! Input file names for time-averaged fields; format must
                                                       !! be identical to the corresponding files in infile

      character*400 ::  outfile                        !! Name of output ASCII format file

      integer       ::  nfiles=-1                       !! Depreciated, do not use

      real          ::  dist_threshold = 750.e3        !! Maximum distance, in meters, allowed between storm centers
                                                       !! identified at successive times (after performing an
                                                       !! extrapolation on the position of the cyclone center at the
                                                       !! earlier time) to be combined into a track.
      
      logical       :: extrap_traj = .true.            !! Whether to perform a extrapolation of cyclones using an
                                                       !! earlier timestep (instead of the steering-level winds in
                                                       !! u_steer, v_steer; steering winds are always used on the
                                                       !! first time analyzed by the tracker) to give a first
                                                       !! estimate of where a storm center should be at the following
                                                       !! time; used in track connection process. Requires
                                                       !! extrap_adv_traj = .true.
      
      logical       :: periodicx  = .true.             !! Whether domain is assumed periodic in the x direction. If
                                                       !! using a global grid or a cartesian doubly-periodic domain
                                                       !! this should be set; if the input domain is only a subset
                                                       !! of a periodic domain (such as only over a portion of the
                                                       !! world; ie. a limited-area domain or one face of a cubed-
                                                       !! sphere grid) then set to .false.

      logical       :: latlon_grid = .true.            !! Whether the input is defined on a regular latitude-longitude
                                                       !! grid. If this is false the tracker will look for the latitudes
                                                       !! and longitudes of each grid cell in the variables name_lat_2d
                                                       !! and name_lon_2d; if use_grid_file = .true. it will look for this
                                                       !! file in grid_file.

      logical       :: extrap_adv_traj = .true.        !! Whether to perform a extrapolation of cyclones using an earlier 
                                                       !! timestep; if extrap_traj = .false. then only the steering winds
                                                       !! are used 

      logical       :: avg_adv_winds = .true.          !! ?? Need to sort out advective options

      real          :: cint_slp = 4.                   !! Contour interval used to identify a cyclone center

      integer       :: ncontours = 1                   !! Number of closed contours needed to identify a cyclone center

      logical       :: warm_core_check = .false.       !! Whether to also find upper-tropospheric temperature extrema (see
                                                       !! name_tm in nlist) near the centers of identified cyclones.

      real          :: vort_thresh = 0.15e-3           !! Minimum cyclonic (hemisphere-dependent) vorticity magnitude needed
                                                       !! in a cyclone region to identify a cyclone center.

      real*4        :: maxSLPthreshold = 1021.         !! Maximum admissible central SLP. This value should be set to 0. if
                                                       !! tracking anomalous SLP.

      real          :: lon_s = 0., lon_e = 360.
      real          :: lat_s = -90., lat_e = 90.       !! Domain over which to restrict search for cyclone centers (in
                                                       !! degrees). This may depend on your input file's domain; some
                                                       !! data uses a longitude range of [-180,180] instead of [0,360].

      character*400 ::  idfile = ''                    !! Creates an output NetCDF file that contains identified cyclone
                                                       !! regions and cyclone centers, grouped by track number. The
                                                       !! resulting files may be used; if blank no file will be created.

      logical       :: do_sort = .true.                !! When enabled, sorts storms from lowest to highest central SLP
                                                       !! before doing search for storms to connect into a track. This is
                                                       !! to ensure that strong cyclones will be connected first to
                                                       !! strong cyclones.

      logical       :: save_untracked_storms = .false. !! Whether to write to the output file storms that are identified
                                                       !! but are not connected into a track.

      real          :: max_cyrad = 3000.               !! Maximum radius for cyclone size, in km. (Converted to meters below)

      real          :: max_warm_offset = 500.          !! Maximum distance warm core allowed from cyclone center, in km.
                                                       !! (Converted to meters below)

      real          :: max_warm_cyrad = 880.           !! Maximum radius for warm-core contours, in km. 
                                                       !! (Converted to meters below)

      logical       :: julian_calendar_fix = .false.   !! Subtracts 13 days (assumes 20th/21st century) from date. This is
                                                       !! useful for coupled simulations that count days from 0001, since the
                                                       !! GFDL model computes the ACTUAL Gregorian date after 4 October 1582,
                                                       !! but does not give the actual number of days if they are being counted
                                                       !! from before 1582. The UDUNITS library assumes that it is being given
                                                       !! the Julian day (the actual number of days) if days are counted from
                                                       !! before 1582.

      logical       :: one_variable_per_file = .false. !! Option for each variable to be in its own file. Input filename must
                                                       !! have the form 'input.#.nc' (one hash) in which the variable name
                                                       !! (given by the flist variables below) replaces the hash-mark. Input
                                                       !! files must otherwise be identical; in particular the dimids must be 
                                                       !! the same.       

      real          :: precip_max_radius = 500.        !! Maximum radius from center for counting precipitation as part of a storm.
                                                       !! km, converted to meters below. Requires infile_avg be set.

      !Land mask files: experimental features, may be removed
      character*400 :: topo_file = '' !! Filename containing surface topography
      logical       :: use_wind_landmask =   .false.
      logical       :: use_wind_oceanmask =   .false.
      logical       :: useregionmask = .false.
      character*400 :: regionmask = ''

      !Other experimental features that may not be implemented
      logical       :: extrap_rk = .false. !Not implemented
      logical       :: checkslpgrad = .false.
      real          :: slpgradthresh = 0.0015 !mb/km, as per Marchok tracker

      logical       :: warm_core_only = .false.   !Backwards compatability, if true warm_core_check is enabled
      

!!!!! flist options: Names of variables in input file(s)
      character(len=40) :: name_lon = 'lon' ! Longitude dimension axis
      character(len=40) :: name_lat = 'lat' ! Latitude dimension axis
      character(len=40) :: name_lon_2D = 'grid_lont'
      character(len=40) :: name_lat_2D = 'grid_latt'
      character(len=40) :: name_area_2D = 'area' 
      character(len=40) :: name_time = 'time'
      character(len=40) :: name_land_mask = 'land_mask'
      character(len=40) :: name_zsurf = 'zsurf'
      character(len=40) :: name_SLP = 'slp'
      character(len=40) :: name_TM  = 'tm'
      character(len=40) :: name_u_ref = 'u_ref' 
      character(len=40) :: name_v_ref = 'v_ref' 
      character(len=40) :: name_u_steer = 'u500' 
      character(len=40) :: name_v_steer = 'v500' 
      character(len=40) :: name_vort850 = 'vort850'
      character(len=40) :: name_precip  = 'precip'

      character(len=40) :: new_name_lon, new_name_lat

      character*100 :: timestring
      character*100 :: errstr
      character*400 :: current_filename, current_avg_filename = '', current_idfilename
      integer :: ios = 0, length = 0

      integer :: outtimes = 0
      integer :: centerid, rgnid, londimid, idfileid, latdimid, timedimid, timeid, lonid, latid

      logical :: errflag
      integer :: numfilled

      logical :: notmasked

      logical :: use_avg_file

      integer :: nargc
      character*400 :: nmlfile = 'namelist.dat'
      integer, external ::  iargc

!$    integer :: omp_get_thread_num
      integer :: num_threads = 1


      integer*8 timecenters_unit ! Pointer to udunits "unit" type
      integer*8 timecenters_unit_avg
      !The size of the pointer variable appears to be system-dependent. Check your
      !udunits.inc for the exact definition (look for 'UD_POINTER')

      !!! UDUNITS ROUTINES
      !!! To get around incombatability with f77 udunits interface
      integer, external :: utopen, utmake, utdec, utcaltime

#ifdef CHECK_MISSING
      print*, '**** Checking for missing values ****'
      print*, '   Tracker will take longer to run   '
#endif

      !!! Constants
      pi = 4.*atan(1.)
      zvir = rvgas/rdgas - 1.
      ginv = 1./grav

      namelist /nlist/infile, infile_avg, outfile, nfiles, &
           dist_threshold, extrap_traj, periodicx, extrap_adv_traj, avg_adv_winds, ncontours, &
           warm_core_check, vort_thresh, lon_s, lon_e, lat_s, lat_e, idfile, &
           topo_file, use_wind_landmask, use_wind_oceanmask, useregionmask, regionmask, cint_slp, extrap_rk, &
           maxSLPthreshold, &
           checkslpgrad, slpgradthresh, do_sort, latlon_grid, debug, grid_file, save_untracked_storms, &
           warm_core_only, max_cyrad, max_warm_offset, max_warm_cyrad, &
           julian_calendar_fix, one_variable_per_file, precip_max_radius

      namelist /flist/ name_lon, name_lat, name_lon_2D, name_lat_2D, name_area_2D, name_time, &
           name_land_mask, name_zsurf, name_SLP, name_TM, &
           name_u_ref, name_v_ref, name_u_steer, name_v_steer, name_vort850, name_precip

      namelist /tlist/ num_threads

!**************
! Namelist var:
!**************
      outfile       = 'track.txt'
      infile(1)     = 'infile.nc'
      infile_avg(1) = ''
      do n=2,MAXFILES
         infile(n)     = ''
         infile_avg(n) = ''
      enddo

!! Process command line options; just namelist file name for now
      !nargc = iargc()
      if (iargc() > 0) then
         call getarg(1,nmlfile)
      endif

      !Initialize udunits; note that Fortran is limited to using udunits v1
      if (utopen('/usr/local/udunits-1.12.11/etc/udunits.dat') /= 0) then
         print*, 'COULD NOT OPEN UDUNITS FILE. Stop.'
         stop 150
      endif

!****************************************************************************
!
! Get control parameters from namelist
!
      open(10,file=trim(nmlfile),status='old')
      read(10,nml=nlist)
      if (nfiles /= -1) print*, 'Namelist option nfiles no longer required or used.'
      if (warm_core_only) then
         warm_core_check = .TRUE.
         print*, 'Namelist option warm_core_only has been renamed warm_core_check.'
         print*, 'Setting warm_core_check to .TRUE.'
      endif
      !Optional namelists
      read(10,nml=flist,iostat=ios,iomsg=errstr)
      if (ios > 0) then
         print*, errstr
         stop 200
      endif
      read(10,nml=nlist,iostat=ios,iomsg=errstr)
      if (ios > 0) then
         print*, errstr
         stop 200
      endif
      close(10)

!$omp parallel
!$    call omp_set_num_threads(num_threads)
!$omp end parallel

      !! Convert km into meters
      max_cyrad = max_cyrad*1000.
      max_warm_offset = max_warm_offset*1000.
      max_warm_cyrad = max_warm_cyrad*1000.
      precip_max_radius = precip_max_radius*1000.

      ! Output text file
      open(13,file=outfile,form='formatted',status='unknown')

      if (useregionmask) then
         open(14,file=trim(outfile)//'.wind',form='formatted',status='unknown')
         open(15,file=trim(outfile)//'.snow',form='formatted',status='unknown')
      endif

      !Open up the first file so we can get ntimes. We will assume all files have roughly the same number of times in them.
      if (one_variable_per_file) then
         call variable_file(current_filename, infile(1), name_SLP)
      else
         current_filename = infile(1)
      endif
      call open_ncfile( current_filename, fid )
      call get_ncdim1( fid, name_time, ntimes )


      allocate(allstorms(MAXSTORM,4)) !i,j,time,track ID
      allocate(allstormval(MAXSTORM,NUMFIELDS)) !lon, lat, minslp, maxwind, maxprecip, maxrain,
      !maxsnow, maxvort, blizzard-cells, mean 500mb u, mean 500 mb v, area
      allocate(allstorms_m1(MAXSTORM,4)) !i,j,time,track ID
      allocate(allstormval_m1(MAXSTORM,NUMFIELDS)) !lon, lat, minslp, maxwind, maxprecip, maxrain,
      !maxsnow, maxvort, blizzard-cells, mean 500mb u, mean 500 mb v, area
      allocate(allstorms_m2(MAXSTORM,4)) !i,j,time,track ID
      allocate(allstormval_m2(MAXSTORM,NUMFIELDS)) !lon, lat, minslp, maxwind, maxprecip, maxrain,
      !maxsnow, maxvort, blizzard-cells, mean 500mb u, mean 500 mb v, area
      allocate(claimed(MAXSTORM))

      !Note that this does not mask lakes
      if (use_wind_landmask) then
         write(*,*) 'Topography file: ', trim(topo_file)
         call open_ncfile( topo_file, fidtopo )
         call get_ncdim1( fidtopo, name_lon, im_t )
         call get_ncdim1( fidtopo, name_lat, jm_t )

         allocate ( land(im_t,jm_t) )
         call get_real2( fidtopo, name_land_mask, im_t, jm_t, land )

         if (use_wind_landmask .and. use_wind_oceanmask) then
            print*, 'ERROR both use_wind_landmask and use_wind_oceanmask cannot be .true. Stop.'
            stop 201
         endif

         if (use_wind_landmask .or. use_wind_oceanmask) then
            allocate ( landtopomask_wind(im_t,jm_t) )
            if (use_wind_landmask) then
               do j=1,jm_t
                  do i=1,im_t
                     landtopomask_wind(i,j) = land(i,j) >= 1.
                  enddo
               enddo
            else if (use_wind_oceanmask) then
               do j=1,jm_t
                  do i=1,im_t
                     landtopomask_wind(i,j) = land(i,j) <= 0.
                  enddo
               enddo
            endif
         endif

         deallocate(land)


      endif

      if (useregionmask) then
         write(*,*) 'Region mask file: ', trim(regionmask)
         call open_ncfile( regionmask, fidrgn )
         call get_ncdim1( fidrgn, 'lon', im_r )
         call get_ncdim1( fidrgn, 'lat', jm_r )

         allocate ( regions(im_r,jm_r) )
         call get_int2( fidrgn, 'regions', im_r, jm_r, regions )

         allocate(regions_wind(MAXSTORM, maxval(regions)))
         allocate(regions_snow(MAXSTORM, maxval(regions)))

         call close_ncfile( fidrgn )

         write(region_format,'(A, I4, A)') '(I6, ', maxval(regions), 'F8.2 )'
         print*, 'region_format = ', region_format
      endif

      if (len(trim(grid_file)) > 0) then
         use_grid_file = .true.
         write(*,*) 'Using grid file ', trim(grid_file)
         call open_ncfile( grid_file, fidgrid )
         call get_ncdim1( fidgrid, name_lon, im_g )
         call get_ncdim1( fidgrid, name_lat, jm_g )

         allocate(lon(im_g,jm_g))
         allocate(lat(im_g,jm_g))
         allocate(area(im_g,jm_g))

         call get_var2_real ( fidgrid, name_lon_2d,  im_g, jm_g, lon )
         call get_var2_real ( fidgrid, name_lat_2d,  im_g, jm_g, lat )
         call get_var2_real ( fidgrid, name_area_2d, im_g, jm_g, area )
      endif

      nstormtot = 0
      nstorm_m1 = -1
      nstorm_m2 = -1

      nfile = 0
      FILELOOP: do

         nfile = nfile + 1
         if (one_variable_per_file) then
            call variable_file(current_filename, infile(nfile), name_SLP)
         else
            current_filename = infile(nfile)
         endif

         if (trim(current_filename) == '') exit
         write(*,'(A12, A)') 'Input file: ', trim(current_filename) 

         if (nfile > 1) then
            call open_ncfile( current_filename, fid )
            call get_ncdim1( fid, name_time, ntimes )
         end if
         current_filename = infile(nfile)

         call get_ncdim1( fid, name_lon, im, new_name_lon )
         call get_ncdim1( fid, name_lat, jm, new_name_lat )
         name_lon = new_name_lon
         name_lat = new_name_lat
         allocate(timeval(ntimes))
         call get_var1_real ( fid, name_time, ntimes, timeval )

         !Get time info
         !Code from http://www.unidata.ucar.edu/cgi-bin/man-cgi?udunits+-s3f
         ! and
         ! http://www.esrl.noaa.gov/psd/data/gridded/readgeneral.f
         ! see also http://www.unidata.ucar.edu/software/netcdf/time/recs.html
         TIMECENTERS_UNIT = UTMAKE()
         call get_var_att_str(fid, name_time, 'units', timestring, nt == 1 .and. nfile == 1)
         STATUS = UTDEC(trim(timestring), TIMECENTERS_UNIT)
         if (status /= 0) then
            print*, 'UTDEC FAILED ', status
            stop 151
         endif

         !!! Average file (just precip for now)
         if (trim(infile_avg(nfile)) /= '') then
            if (one_variable_per_file) then
               call variable_file(current_avg_filename, infile_avg(nfile), name_precip)
            else
               current_avg_filename = infile_avg(nfile)
            endif            
            use_avg_file = .true.
            print*, ' Reading averages file ' , current_avg_filename
            call open_ncfile(current_avg_filename, fidavg )
            call get_ncdim1( fidavg, name_lon, im_avg )
            call get_ncdim1( fidavg, name_lat, jm_avg )
            call get_ncdim1( fidavg, name_time, ntimes_avg )
            allocate(timeval_avg(2,ntimes_avg))
            call get_var2_real ( fidavg, 'time_bounds', 2, ntimes_avg, timeval_avg )
            if (im /= im_avg .or. jm /= jm_avg .or. ntimes /= ntimes_avg) then
               print*, 'Dimensions of average file does not match input file (check time dimension too)'
               stop 2017
            endif
            if (any(timeval /= timeval_avg(2,:))) then
               print*, 'Times in average file do not match those of the input file'
               print*, timeval(1), timeval_avg(2,1)
               stop 2018
            endif
         else
            current_avg_filename = ''
            use_avg_file = .false.
         endif

         if ((use_wind_landmask .or. use_wind_oceanmask )) then
            if (im /= im_t .or. jm /= jm_t) then
               print*, 'Dimensions of topography file does not match input file'
               stop 2019
            endif
         endif
         if (useregionmask) then
            if (im /= im_r .or. jm /= jm_r) then
               print*, 'Dimensions of region mask file does not match input file'
               stop 2020
            endif
         endif
         if (use_grid_file) then
            if (im /= im_g .or. jm /= jm_g) then
               print*, 'Dimensions of grid file do not match input file'
               print*, im, im_g, jm, jm_g
               stop 2021
            endif
         endif

! If we are creating a netcdf output file for the identified locations then set it up here
         if (len(trim(idfile)) > 0) then

            !Get hour
            time8 = real(timeval(1),8)
            if (julian_calendar_fix) time8 = time8 - 13.
            status = UTCALTIME(time8,TIMECENTERS_UNIT,year,month,day,hour,tmin,tsec)
            if (status /= 0) then
               print*, 'UTCALTIME FAILED ', status
               stop 153
            endif

            length = len(trim(idfile))
            !Assumes name ends in .nc
            write(current_idfilename,'(2A, I4.4, 2I2.2, A)') idfile(1:length-3), '.', year, month, day, '.nc'
            status = nf_create(trim(current_idfilename), NF_NOCLOBBER, idfileid)
            if (status .ne. NF_NOERR) then
               print*, 'Could not create idfile ', trim(current_idfilename)
               call handle_err(status)
            endif

            print*, 'Creating idfile ', trim(current_idfilename)
            outtimes = 0

            status = nf_def_dim(idfileid, 'lon', im, londimid) 
            if (status .ne. NF_NOERR) call handle_err(status)
            status = nf_def_dim(idfileid, 'lat', jm, latdimid)
            if (status .ne. NF_NOERR) call handle_err(status)
            status = nf_def_dim(idfileid, 'time', NF_UNLIMITED, timedimid)
            if (status .ne. NF_NOERR) call handle_err(status)

            status = nf_def_var(idfileid, 'time', NF_DOUBLE, 1, (/ timedimid /), timeid)
            if (status .ne. NF_NOERR) call handle_err(status)
            call get_var_att_str(fid, name_time, 'units', units, nt == 1 .and. nfile == 1)
            status =  NF_PUT_ATT_TEXT  (idfileid, timeid,'units', len(trim(units)), trim(units))
            if (status .ne. NF_NOERR) call handle_err(status)

            !Define all our variables
            status = nf_def_var(idfileid, 'ID_regions', NF_INT, 3, (/ londimid, latdimid, timedimid /), rgnid)
            if (status .ne. NF_NOERR) call handle_err(status)
            call set_var_att_int(idfileid, 'ID_regions', '_FillValue', 0)
            !status = nf_def_var(idfileid, 'ID_centers', NF_INT, 3, (/ londimid, latdimid, timedimid /), centerid)
            !if (status .ne. NF_NOERR) call handle_err(status)
            status = nf_def_var(idfileid, 'lon', NF_FLOAT, 1, (/ londimid /), lonid)
            if (status .ne. NF_NOERR) call handle_err(status)
            status = nf_def_var(idfileid, 'lat', NF_FLOAT, 1, (/ latdimid /), latid)
            if (status .ne. NF_NOERR) call handle_err(status)

            status = nf_enddef(idfileid)

         endif

         status = nf_inq_varid (fid, name_slp, varid)
         
         if (.not. use_grid_file ) then
            allocate ( lon(im,jm) )
            allocate ( lat(im,jm) )
            allocate ( area(im,jm) )
         endif
         allocate ( slp(im,jm) )
         if (warm_core_check) allocate (  tm(im,jm) )
         allocate ( vort(im,jm) )
         allocate ( u_ref(im,jm) )
         allocate ( v_ref(im,jm) )
         allocate ( u_steer(im,jm) )
         allocate ( v_steer(im,jm) )
         allocate ( wk(im,jm) )
         allocate ( slp_s(im,jm) )
         allocate ( id(im,jm) )
         if (len(trim(idfile)) > 0) allocate ( id_out(im,jm) )
         allocate ( centers(im,jm) )

         if (warm_core_check) then
            allocate(warm_core_center(MAXSTORM,2))
            allocate ( id_warm(im,jm) )
         endif

         if (use_avg_file) then
            allocate(precip(im,jm))
         endif

         allocate ( lon_in(im) ) 
         allocate ( lat_in(jm) ) 
         if (latlon_grid) then
            call get_var1_real ( fid, name_lon, im, lon_in )
            call get_var1_real ( fid, name_lat, jm, lat_in )
            do j=1,jm
               do i=1,im
                  lon(i,j) = lon_in(i)
                  lat(i,j) = lat_in(j)
                  !uniform lat-lon grid
                  area(i,j) = RADIUS**2 * (2.*pi/real(im)) * 2. * cos(pi/180.*lat(i,j)) * sin( pi/real(jm+1))
               enddo
            enddo
         else
            !Assume a 2d array
            if (.not. use_grid_file) then
               call get_var2_real ( fid, name_lon_2d, im, jm, lon )
               call get_var2_real ( fid, name_lat_2d, jm, jm, lat )
               call get_var2_real ( fid, name_area_2d, jm, jm, area )
            endif
         endif

         if (nfile == 1 .and. len(trim(idfile)) > 0) then
            if (.not. latlon_grid) then
               do j=1,jm
                  do i=1,im
                     lon_in(i) = real(i)
                     lat_in(j) = real(j)
                  enddo
               enddo
            endif
            status = nf_put_vara_real(idfileid, lonid, (/1/), (/im/), lon_in)
            if (status .ne. NF_NOERR) call handle_err(status)
            status = nf_put_vara_real(idfileid, latid, (/1/), (/jm/), lat_in)
            if (status .ne. NF_NOERR) call handle_err(status)
         endif
         deallocate(lon_in)
         deallocate(lat_in)

         if (nfile == 1) then
            call pmaxmin4( 'LON', lon, im, jm, 1. )            ! lon [0,360]
            call pmaxmin4( 'LAT', lat, im, jm, 1. )            ! lat [-90,90]
            call pmaxmin4( 'AREA (km**2)', area, im, jm, 1.e-6 )            
         endif

         if (latlon_grid) then !assuming roughly uniform grid spacing in both directions
            if (periodicx) then
               dphi_avg = 360./real(im) ! degrees
               dx_avg = dphi_avg*RADIUS*pi/180.
            else
               !This estimate is slightly incorrect, since the distance
               !around the latitudinal direction is one point short.
               dphi_avg = (lon(im,1)-lon(1,1))/real(im-1)
               dx_avg = dphi_avg*RADIUS*pi/180.
            endif
         else
            !If not a regular lat-lon grid, take the average of the lengths of the domain edges
            ! and divide that by the number of cells along each edge
            !north and south
            dx_avg = 0.
            do i=1,im-1
               dx_avg = dx_avg + great_circle_dist((/lon(i,1),lat(i,1)/),  (/lon(i+1,1),lat(i+1,1)/),  RADIUS)
               dx_avg = dx_avg + great_circle_dist((/lon(i,jm),lat(i,jm)/),(/lon(i+1,jm),lat(i+1,jm)/),RADIUS)
            enddo
            !east and west
            do j=1,jm-1
               dx_avg = dx_avg + great_circle_dist((/lon(1,j),lat(1,j)/),  (/lon(1,j+1),lat(1,j+1)/),  RADIUS)
               dx_avg = dx_avg + great_circle_dist((/lon(im,j),lat(im,j)/),(/lon(im,j+1),lat(im,j+1)/),RADIUS)
            enddo
            dx_avg = dx_avg/(2.*(im-2)+2.*(jm-2))
         endif

         if (warm_core_check) warm_core_center = 0

         !TIMELOOP: do nt = 91, ntimes ! for testing purposes
         TIMELOOP: do nt = 1, ntimes

            claimed = .false.
            allstorms = 0
            allstormval = 0.
            allstormval(:,3) = maxSLPthreshold

            if (useregionmask) then
               regions_wind = 0.
               regions_snow = 0.
            endif

            if (nfile > 1 .or. nt > 1 ) then
               year_prev = year
               month_prev = month
               day_prev = day
               hour_prev = hour
               tmin_prev = tmin
               tsec_prev = tsec
            endif

            !Get hour
            time8 = real(timeval(nt),8)
            if (julian_calendar_fix) time8 = time8 - 13.
            status = UTCALTIME(time8,TIMECENTERS_UNIT,year,month,day,hour,tmin,tsec)
            if (status /= 0) then
               print*, 'UTCALTIME FAILED ', status
               stop 152
            endif

            !If after the first timestep calculate dt_hr. Since dt_hr
            !is an integer number of hours, we only need to drill down to hours.
            if (nfile > 1 .or. nt > 1 ) then
               if (year > year_prev) then
                  month_prev = month_prev - 12
               endif
               if (month > month_prev) then
                  !We include negative values to take care of negative months due to the year changing
                  select case (month_prev)
                  case (1, 3, 5, 7, 8, 10, 12, -11, -9, -7, -5, -4, -2, 0)
                     day_prev = day_prev - 31
                  case (4, 6, 9, 11, -8, -6, -3, -1)
                     day_prev = day_prev - 30
                  case (2, -10)
                     if (mod(year,4)) then
                        if (mod(year,100) .and. .not. mod(year,400)) then
                           day_prev = day_prev - 28
                        else
                           day_prev = day_prev - 29
                        end if
                     else
                        day_prev = day_prev - 29
                     end if
                  end select
               endif
               if (day > day_prev) then
                  hour_prev = hour_prev - 24
               end if
               dt_hr = hour - hour_prev
               !print*, 'dt_hr = ', dt_hr
               !!! Sanity check. dt_hr should probably be no greater than
               !!! 48, and should not be negative.
               if (dt_hr > 48) then
                  print*, 'ERROR COMPUTING DT_HR: too large', dt_hr
                  print*, 'CURRENT TIME:', year, month, day, hour
                  print*, 'PREV    TIME:', year_prev, month_prev, day_prev, hour_prev, ' (PROCESSED)'
                  stop 155
               endif
               if (dt_hr < 0) then
                  print*, 'ERROR COMPUTING DT_HR: negative value', dt_hr
                  print*, 'CURRENT TIME:', year, month, day, hour
                  print*, 'PREV    TIME:', year_prev, month_prev, day_prev, hour_prev, ' (PROCESSED)'
                  stop 156
               endif
            endif

            call get_real3( current_filename, fid, name_slp, im, jm, nt, slp, one_variable_per_file )
            call get_var_att_str(fid, name_slp, 'units', units, nt == 1 .and. nfile == 1)
#ifdef CHECK_MISSING
            call get_var_att_real(fid, name_slp, '_FillValue', slp_fillvalue, nt == 1 .and. nfile == 1)
#endif
            units = units(1:2) !necessary because apparently get_var_att_str returns
                               !INVISIBLE GARBAGE preventing a proper comparison
            dumunit = 'Pa'
            if (trim(units) .eq. 'Pa' .or. trim(units) .eq. 'pa' .or. trim(units) .eq. 'PA' .or. trim(units) .eq. '  ') then
               slp = slp/100.
#ifdef CHECK_MISSING
               slp_fillvalue = slp_fillvalue/100.
#endif
            end if
            if (warm_core_check) then
               call get_real3( current_filename, fid, name_tm, im, jm, nt, tm, one_variable_per_file )
#ifdef CHECK_MISSING
               call get_var_att_real(fid, name_tm, '_FillValue', tm_fillvalue, nt == 1 .and. nfile == 1)
#endif
            endif

            if (use_avg_file ) then
               call get_real3( current_avg_filename, fidavg, name_precip, im, jm, nt, precip)
            endif

            call get_real3( current_filename, fid, name_u_ref, im, jm, nt, u_ref, one_variable_per_file )
            call get_real3( current_filename, fid, name_v_ref, im, jm, nt, v_ref, one_variable_per_file )
#ifdef CHECK_MISSING
            call get_var_att_real(fid, name_u_ref, '_FillValue', u_fillvalue, nt == 1 .and. nfile == 1)              
            call get_var_att_real(fid, name_v_ref, '_FillValue', v_fillvalue, nt == 1 .and. nfile == 1)              
#endif               
            !Check if input has 500 mb winds available. If not, then we will use
            !surface winds times a fudge factor to start the cyclone extrapolation process
            if (one_variable_per_file) then
               if (nfile == 1 .and. nt == 1) then
                  print*, 'Steering winds not yet compatible with one_variable_per_file.'
                  print*, 'using surface winds*1.5 for first step of cyclone track extrapolation'
               endif
               u_steer = u_ref * 1.5
               v_steer = v_ref * 1.5
            else
               status = nf_inq_varid (fid, name_u_steer, varid)

               if (status == NF_NOERR) then
                  call get_real3( current_filename, fid, name_u_steer, im, jm, nt, u_steer, one_variable_per_file )
                  call get_real3( current_filename, fid, name_v_steer, im, jm, nt, v_steer, one_variable_per_file )
               else
                  if (nt == 1) print*, 'steering winds not found; using surface winds*1.5 for first step of cyclone track extrapolation'
                  u_steer = u_ref * 1.5
                  v_steer = v_ref * 1.5
               endif
            endif

            if (vort_thresh > 0.) then
               call get_real3( current_filename, fid, name_vort850, im, jm, nt, vort, one_variable_per_file )
#ifdef CHECK_MISSING
               call get_var_att_real(fid, 'vort850', '_FillValue', vort_fillvalue, nt == 1 .and. nfile == 1)
#endif
            end if

            if ( nt==1 .and. nfile==1) then
               call pmaxmin4( 'SLP', slp, im, jm, 1. )
               call pmaxmin4( 'u10', u_ref, im, jm, 1. )
               call pmaxmin4( 'v10', v_ref, im, jm, 1. )
               if (vort_thresh > 0.) then
                  call pmaxmin4( 'vort850', vort, im, jm, 1. )
               else
                  print*, 'NO VORTICITY THRESHOLD'
               endif
               if (use_avg_file) call pmaxmin4( 'PRECIP (mm/day)', precip, im, jm, 86400. )
            endif

            !The '.2' and '.4' enables zero-padding of integers
            write(*,'(A, I4.4"-"I2.2"-"I2.2" "I2.2"Z")') 'ANALYZING TIME ', year, month, day, hour 

#ifdef CHECK_MISSING
            call split_smoother(im, jm, slp, slp_s, periodicx, slp_fillvalue)
#else
            call split_smoother(im, jm, slp, slp_s, periodicx)
#endif
            id = 0
#ifdef CHECK_MISSING
            call find_extrema(slp_s, id, im, jm, allstorms(:,1:2), allstormval(:,3), nstorm, .true., maxSLPthreshold, ncontours, &
                 lon, lat, lon_s, lon_e, lat_s, lat_e, max_cyrad, cint_slp, checkslpgrad, slpgradthresh, periodicx, slp_fillvalue)
#else
            call find_extrema(slp_s, id, im, jm, allstorms(:,1:2), allstormval(:,3), nstorm, .true., maxSLPthreshold, ncontours, &
                 lon, lat, lon_s, lon_e, lat_s, lat_e, max_cyrad, cint_slp, checkslpgrad, slpgradthresh, periodicx)
#endif

            print*,' FOUND ', nstorm, ' STORMS'
            if (nstorm > MAXSTORM) then
               print*, ' nstorm = ', nstorm, ' exceeds MAXSTORM = ', MAXSTORM
               print*, ' exiting; increase MAXSTORM, ncontours, or cint_slp'
               stop 104
            endif

            allstorms(:,3) = nt

            !Collect storms
!$omp parallel do default(shared) private(i,j,n)
            do n=1,nstorm

               !! Chern's Method:
               ! fine tuning the center location by assuming the sub-grid distribution around 
               ! the center (i,j) is second order polynomial
               !Be sure to use the SMOOTHED vorticity since that is what is used to find extrema in the first place
               i = allstorms(n,1)
               j = allstorms(n,2)

               if (j == 1 .or. j == jm) then
                  !At the domain boundary (or poles), cannot fine-tune
                     allstormval(n,1) = lon(i,j)
                     allstormval(n,2) = lat(i,j)
               else if (i == 1) then
                  if (periodicx) then
                     if (max(abs(slp_s(im,j)),abs(slp_s(i+1,j)),abs(slp_s(i,j+1)),abs(slp_s(i,j-1))) < 1.e8) then
                        !Assuming domain is 360 degrees around
                        allstormval(n,1) = lon(i,j) + & 
                             0.25*(lon(i+1,j)-lon(im,j)+360.)*(slp_s(im,j)-slp_s(i+1,j))/(slp_s(im,j) - 2.*slp_s(i,j) + slp_s(i+1,j))
                        allstormval(n,2) = lat(i,j) + & 
                             0.25*(lat(i,j+1)-lat(i,j-1))*(slp_s(i,j-1)-slp_s(i,j+1))/(slp_s(i,j-1) - 2.*slp_s(i,j) + slp_s(i,j+1))
                     endif
                  else
                     !At the domain boundary, cannot fine-tune
                     allstormval(n,1) = lon(i,j)
                     allstormval(n,2) = lat(i,j)
                  endif
               else if (i == im) then
                  if (periodicx) then
                     if (max(abs(slp_s(1,j)),abs(slp_s(i-1,j)),abs(slp_s(i,j+1)),abs(slp_s(i,j-1))) < 1.e8) then
                        allstormval(n,1) = lon(i,j) + & 
                             0.25*(lon(1,j)-lon(i-1,j)+360.)*(slp_s(i-1,j)-slp_s(1,j))/(slp_s(i-1,j) - 2.*slp_s(i,j) + slp_s(1,j))
                        allstormval(n,2) = lat(i,j) + & 
                             0.25*(lat(i,j+1)-lat(i,j-1))*(slp_s(i,j-1)-slp_s(i,j+1))/(slp_s(i,j-1) - 2.*slp_s(i,j) + slp_s(i,j+1))
                     endif
                  else
                     !At the domain boundary, cannot fine-tune
                     allstormval(n,1) = lon(i,j)
                     allstormval(n,2) = lat(i,j)
                  end if
               else
                  if (ALL(slp_s(i-1:i+1,j-1:j+1) < 1.e8)) then
                     allstormval(n,1) = lon(i,j) + & 
                          0.25*(lon(i+1,j)-lon(i-1,j))*(slp_s(i-1,j)-slp_s(i+1,j))/(slp_s(i-1,j) - 2.*slp_s(i,j) + slp_s(i+1,j))
                     allstormval(n,2) = lat(i,j) + & 
                          0.25*(lat(i,j+1)-lat(i,j-1))*(slp_s(i,j-1)-slp_s(i,j+1))/(slp_s(i,j-1) - 2.*slp_s(i,j) + slp_s(i,j+1))
                  endif
               endif

               if (debug) print*, 'Storm ', n, ' at ', i, j, ' tuned to ', allstormval(n,1:2)
               if (debug .and.i > 1 .and. i < im .and. j > 1 .and. j < jm) then
                  do jj=j-1,j+1
                     print*, slp(i-1:i+1,jj)
                  enddo
               endif

               !Check for bad lat/lon values
               !Making sure the absolute value is greater than 0.01
               ! avoids issues when doing the formatted write; if the
               ! value is too close to zero Fortran formats write out
               ! asterisks.
               !Using 370 instead of 360 for max lat to preserve values
               !which may have been moved across the prime meridian
               if (allstormval(n,1) > 370.) allstormval(n,1) = allstormval(n,1) - 360
               if (abs(allstormval(n,1)) < 0.01 ) allstormval(n,1) = sign(0.01,allstormval(n,1))
               if (abs(allstormval(n,2)) > 90.  ) allstormval(n,2) = sign(90. ,allstormval(n,2))
               if (abs(allstormval(n,2)) < 0.01 ) allstormval(n,2) = sign(0.01,allstormval(n,2))

               !allstormval(n,1) = lon(storms_1t(1,n))
               !allstormval(n,2) = lat(storms_1t(2,n))

            enddo

            if (vort_thresh > 0.) then
#ifdef CHECK_MISSING
               call split_smoother(im, jm, vort, wk, periodicx, vort_fillvalue)
#else
               call split_smoother(im, jm, vort, wk, periodicx)
#endif
               vort = wk
            end if

            !Get diagnostics
!$omp parallel do default(shared) private(i,j,n,wind,rad)
            do j=1,jm

               do i=1,im

                  if (lat(i,j) > lat_e .or. lat(i,j) < lat_s) cycle

                  if (lon(i,j) > lon_e .or. lon(i,j) < lon_s) cycle

                  if (id(i,j) > 0)  then
                     if (id(i,j) > nstorm) print*, 'ID =', id(i,j), i, j, im, jm
                     n = id(i,j)

                     !Minslp
                     if (slp(i,j) < allstormval(n,3)) then
                        allstormval(n,3) = slp(i,j)
                     endif

                     !Maxwind
                     !Calculate if not doing masking (either land or ocean) or if the
                     !mask variable is true here (both branches have identical code)
                     if ( .not.( use_wind_landmask .or. use_wind_oceanmask) )  then
#ifdef CHECK_MISSING
                        if ( .not. ( u_ref(i,j) == u_fillvalue .or. v_ref(i,j) == v_fillvalue) ) then
#endif
                        wind = sqrt( u_ref(i,j)**2 + v_ref(i,j)**2 )
                        if (wind > allstormval(n,4)) allstormval(n,4) = wind

                        if (useregionmask) then
                           if ( regions(i,j) > 0) regions_wind(n,regions(i,j)) = max(wind,regions_wind(n,regions(i,j)))
                        endif
#ifdef CHECK_MISSING
                        endif
#endif
                     else if (landtopomask_wind(i,j) ) then 
#ifdef CHECK_MISSING
                        if ( .not.( u_ref(i,j) == u_fillvalue .or. v_ref(i,j) == v_fillvalue) ) then
#endif
                        wind = sqrt( u_ref(i,j)**2 + v_ref(i,j)**2 )
                        if (wind > allstormval(n,4)) allstormval(n,4) = wind

                        if (useregionmask) then
                           if (regions(i,j) > 0) regions_wind(n,regions(i,j)) = max(wind,regions_wind(n,regions(i,j)))
                        endif
#ifdef CHECK_MISSING
                        endif
#endif
                     endif

                     !Max smoothed vorticity (note hemispheric dependence)
                     if (vort_thresh > 0.) then
#ifdef CHECK_MISSING
                        if (wk(i,j) /= vort_fillvalue) then
#endif
                        if (lat(i,j) >= 0) then
                           if (wk(i,j) >= allstormval(n,8)) allstormval(n,8) = wk(i,j)
                        else
                           if (wk(i,j) <= allstormval(n,8)) allstormval(n,8) = wk(i,j)
                        endif
#ifdef CHECK_MISSING
                        endif
#endif
                     endif

                     !Area
                     allstormval(n,12) = allstormval(n,12) + area(i,j)
                     !if ( nt==1 .and. nfile==1) then
                     !   print*, 'AREA: ', lon(i), lat(j), area
                     !endif

                     !Mean 500mb u and v; add up, then divide by area measurement
                     allstormval(n,10) = allstormval(n,10) + u_steer(i,j)*area(i,j)
                     allstormval(n,11) = allstormval(n,11) + v_steer(i,j)*area(i,j)

                     rad = great_circle_dist((/lon(i,j), lat(i,j)/), allstormval(n,1:2), RADIUS)
                     !Precipitation, if using. Units are mm*m**2/s Go no more than 500 km from the storm center.
                     if (use_avg_file) then
                        if ( rad <= precip_max_radius) then
                           allstormval(n, 5) = allstormval(n,5) + precip(i,j)*area(i,j)
                        endif
                     endif

                     !Warm core test
                     !TO DO: Need to change from cells to actual distance
                     if (warm_core_check ) then
                        if (tm(i,j) > allstormval(n,13)) then
#ifdef CHECK_MISSING
                           if (tm(i,j) /= tm_fillvalue) then
#endif
                           !Want to make sure we are not too far from storm center
                           if (rad <= max_warm_offset) then
                              allstormval(n,13) = tm(i,j)
                              warm_core_center(n,:) = (/i, j/)
                           endif
#ifdef CHECK_MISSING
                           endif
#endif

                        endif
                     endif

                  endif

               enddo
            enddo

!$omp parallel do default(shared)
            do n=1,nstorm
               !Divide out area to get average 500 mb wind
               allstormval(n,10) = allstormval(n,10)/allstormval(n,12)
               allstormval(n,11) = allstormval(n,11)/allstormval(n,12)
               !do the same for precip, converting to mm/d
               allstormval(n,5)  = allstormval(n,5) /allstormval(n,12) * 86400.
            end do

            !Zero out disturbances that do not meet various criteria.
            if (warm_core_check) then
               wk = tm
#ifdef CHECK_MISSING
               call split_smoother(im, jm, tm, wk, periodicx, tm_fillvalue)
#else
               call split_smoother(im, jm, tm, wk, periodicx)
#endif
            end if

!This part notyet parallelizable due to the flood fill (?)
            do n=1,nstorm
               if (vort_thresh > 0. .and. abs(allstormval(n,8)) < vort_thresh) then
                  allstorms(n,4) = -1
                  if (debug) print*, 'STORM ', n, 'DID NOT MEET VORTICITY THRESHOLD', abs(allstormval(n,8))
               endif
               if (warm_core_check) then
                  !Perform warm-core test: a closed 2C contour in tm
                  id_warm = 0.
                  numfilled = 0

                  if (warm_core_center(n,1) > 0 .and. warm_core_center(n,2) > 0) then
                   if (debug) print*, 'Searching for closed tm contour: ', &
                         lon(warm_core_center(n,1),warm_core_center(n,2)),  &
                          lat(warm_core_center(n,1),warm_core_center(n,2)), allstormval(n,13)
                     call flood_fill(wk,id_warm,warm_core_center(n,1),warm_core_center(n,2), &
                                                warm_core_center(n,1),warm_core_center(n,2), &
                                                im,jm,allstormval(n,13),allstormval(n,13)-2.,&
                                                n,.false., errflag, numfilled, lon, lat, max_warm_cyrad, 1, .false.)
                     if (errflag) then
                        if (debug) print*, 'STORM ', n, 'DID NOT MEET WARM-CORE TEMPERATURE THRESHOLD', allstormval(n,13)
                        allstormval(n,13) = -1.
                     endif
                  else
                     allstormval(n,13) = -1.
                  endif
               else
                  allstormval(n,13) = 0.
               endif
            enddo

            !Sort in order of ascending pressure. Later when we connect tracks the strongest cyclones
            ! will start by looking through the previous time's strongest cyclones
            if (do_sort) then
               if (debug) print*, 'DOING SORTING'
!$omp parallel do default(shared) private(n)
               do n=1,nstorm
                  if (allstorms(n,4) >= 0) then
                     allstorms(n,4) = n
                  else
                     allstorms(n,4) = 0
                  endif
               enddo
!!$                  !!! DEBUG CODE
!!$                  print*, n, allstorms(n,4), allstormval(n,3), allstormval(n,13)
!!$                  !!! END DEBUG CODE
               call heapsort(nstorm,NUMFIELDS,4,3,allstormval(1:nstorm,:),allstorms(1:nstorm,:))
               stormids = -1
!$omp parallel do default(shared) private(n)
               do n=1,nstorm
                  if (allstorms(n,4) > 0) then
                     stormids(allstorms(n,4)) = n
                  endif
               enddo
!!$               !!! DEBUG CODE
!!$               do n=1,nstorm
!!$                  print*, n, stormids(n), allstorms(n,4), allstormval(n,3), allstormval(n,13)
!!$               enddo
!!$               !!! END DEBUG CODE
            endif

            !Perform storm tracking. 
            if (nt == 1 .and. nfile == 1 .or. nstorm_m1 <= 0) then
!$omp parallel do default(shared)
               do n=1,nstorm
                  if (.not. allstorms(n,4) < 0) then
                     allstorms(n,4) = trackcount
                     trackcount = trackcount+1
                  endif
               enddo
            else

            if (debug) print*, 'TIME: ', nfile, nt, 'NSTORM: ', nstorm, nstorm_m1, nstorm_m2

            !! New process: two-step track connection
            !! 1. For all previous time's storms:
            !!    a. compute extrapolated track
            !!    b. sort by minimum pressure (?) 
            !! 2. For all current time's storms:
            !!    a. again sort by minimum pressure
            !!    b. Go down the sorted list, from deepest to least-deep, of
            !!       the current time's storms. Compute the distance to each un-claimed storm


            !WARNING does this mean we cannot parallelize? Might need to re-think this
            do n=1,nstorm
               
               if (allstorms(n,4) < 0) cycle

               distmin = 1.1*dist_threshold
               minstorm = 0
               !search for prior storm at earlier time. That storm should clear the vorticity threshold.
               !do nprev=1,nstorm_m1

               !NOTE: This appears to be extremely inefficient. We should compute the
               !extrapolated trajectory for each of the prior time's storms first, save that list, then perform the search.
               !do nprev=nstorm_m1,1,-1
               do nprev=1,nstorm_m1

                  !If a previous-time storm has not yet been claimed by
                  !another, compute the distance.

                  !We will now compute distances from where we expect the storm to have moved in the past six hours.
                  !1. If the previous time is the first time for this storm, track using h500 winds.
                  !2. If two or more previous times exist, then perform a linear extrapolation.
                  if (.not. claimed(nprev) .and. allstorms_m1(nprev,4) > 0) then

                     if (extrap_adv_traj) then

                        prevstorm = -1

                        if (extrap_traj .and. nstorm_m2 > 0) then
                           !do np2 = 1,nstorm_m2
                           do np2 = nstorm_m2,1,-1
                              if (allstorms_m2(np2,4) == allstorms_m1(nprev,4)) then
                                 prevstorm = np2
                                 !Perform extrapolation
                                 if (extrap_rk) then
                                    print*, 'Runge-Kutta extrapolation not implemented. Stop.'
                                    stop 157
                                 else
                                    !Linear extrapolation
                                    call intp_great_circle(2., allstormval_m2(np2,1:2), allstormval_m1(nprev,1:2), extp(1), extp(2))
                                 endif
                                 dist = great_circle_dist(extp, allstormval(n,1:2), RADIUS)
                                 if (debug) print*, 'EXTRAP DIST = ', dist, nstormtot+n, nprev, np2, &
                                      extp, allstormval_m2(np2,1:2), allstormval_m1(nprev,1:2)
                                 exit
                              endif
                           enddo
                        endif

                        if (prevstorm == -1) then
                           !Perform advection using h500 winds
                           if (avg_adv_winds) then
                              extp(1) = allstormval_m1(nprev,1) + &
                                   allstormval(n,10)/RADIUS/cos(pi/180.*allstormval_m1(nprev,2) )*dt_hr*3600.
                              extp(2) = allstormval_m1(nprev,2) + &
                                   allstormval(n,11)/RADIUS*dt_hr*3600.
                           else
                              !NOTE: Assumes regular lat-lon grid
                              dum = minloc(abs(lon - allstormval_m1(nprev,1)),1)
                              i = dum(1)
                              dum = minloc(abs(lat - allstormval_m1(nprev,2)),2)
                              j = dum(2)
                              extp(1) = allstormval_m1(nprev,1) + &
                                   u_steer(i,j)/RADIUS/cos(pi/180.*allstormval_m1(nprev,2) )*dt_hr*3600.
                              extp(2) = allstormval_m1(nprev,2) + &
                                   v_steer(i,j)/RADIUS*dt_hr*3600.                                 
                           endif
                           dist = great_circle_dist(extp, allstormval(n,1:2), RADIUS)
                           if (debug) print*, 'ADV DIST = ', dist, nstormtot+n, nprev, &
                                extp, allstormval_m1(nprev,1:2), &
                                allstormval(n,10), &
                                allstormval(n,11), &
                                allstormval(n,10)/RADIUS/cos(pi/180.*allstormval_m1(nprev,2) )*dt_hr*3600., &
                                allstormval(n,11)/RADIUS*dt_hr*3600.
                        endif

                     else
                        !No advection or extrapolation
                        dist = great_circle_dist( allstormval_m1(nprev,1:2), allstormval(n,1:2), RADIUS)
                     endif
                     if (debug) print*, 'DIST = ', dist, nstormtot+n, nprev, dist_threshold, distmin
                     if (dist > dist_threshold) cycle
                     if (dist < distmin) then
                        distmin = dist
                        minstorm = nprev
                     endif

                  endif


               enddo

               if (minstorm == 0) then
                  !Found no associated storms at earlier time; start a new track
                  if (debug) print*, 'STARTING NEW TRACK', trackcount, n+nstormtot
                  allstorms(n,4) = trackcount
                  trackcount = trackcount + 1
               else
                  allstorms(n,4) = allstorms_m1(minstorm,4)
                  claimed(minstorm) = .true.
                  if (debug) print*, 'TRACK: ', allstorms(n,4), n, minstorm
               endif

            enddo
         endif


         !FORMAT
         ! 1) record ID
         ! 2--5) year, month, day, hour
         ! 6) Storm ID
         ! 7) center longitude
         ! 8) center latitutde
         ! 9) minimum slp
         ! 10) maximum wind
         ! 11--13) maximum precip/rain/snow rate
         ! 14) greatest smoothed cyclonic vorticity 
         ! 15) Number of cells with blizzard conditions
         ! 16) Warm core flag
         ! 17) Storm area

         !Storm file I/O
         do n=1,nstorm
            if (allstorms(n,4) == -1 .and. .not. save_untracked_storms) cycle
            write(13,'(I8, I6, 3I4, 1x, I7, 7(3x, F10.2), 3x, E8.2, 3x, I5, 3x, F8.2, F16.2)') nstormtot+n, year, month, day, hour, allstorms(n,4), allstormval(n,1:8), max(0,int(allstormval(n,9))), allstormval(n,13), allstormval(n,12)/1.e6
            if (useregionmask) then
               write(14,trim(region_format)) nstormtot+n, regions_wind(n,:)
               write(15,trim(region_format)) nstormtot+n, regions_snow(n,:)
            endif
!!$            !!! DEBUG CODE
!!$            print*, n, stormids(n), allstorms(stormids(n),4)
!!$            !!! END DEBUG CODE
         enddo

!!! Identified region output; replace ID with track numbers
         if (len(trim(idfile)) > 0) then
            id_out = 0
            outtimes = outtimes + 1
            do j=1,jm
               do i=1,im
                  if (id(i,j) > 0) then
                     n = stormids(id(i,j))
                     if (n > 0) id_out(i,j) = allstorms(n,4)
                  endif
               enddo
            enddo
            status = nf_put_vara_int(idfileid, rgnid, (/1, 1, outtimes/), (/im, jm, 1/), id_out)
            if (status .ne. NF_NOERR) call handle_err(status)
            !Also write out storm centers
            !centers = 0
            !do n=1,nstorm
            !   i = allstorms(n,1)
            !   j = allstorms(n,2)
            !   centers(i,j) = allstorms(n,4)
            !end do
            !status = nf_put_vara_int(idfileid, centerid, (/1, 1, outtimes/), (/im, jm, 1/), centers)
            !if (status .ne. NF_NOERR) call handle_err(status)
            time8 = real(timeval(nt),8)
            status = nf_put_vara_double(idfileid, timeid, (/outtimes/), (/1/), time8)
            if (status .ne. NF_NOERR) call handle_err(status)

            !Skipping sync to improve performance
            !status = nf_sync(idfileid)
            !if (status .ne. NF_NOERR) call handle_err(status)
         endif

         nstormtot = nstormtot+nstorm
         nstorm_m2 = nstorm_m1
         nstorm_m1 = nstorm

         allstorms_m2 = allstorms_m1
         allstorms_m1    = allstorms

         allstormval_m2 = allstormval_m1
         allstormval_m1    = allstormval

      end do TIMELOOP


9020  format(i4,1x,i2,1x,i2,1x,i2,1x,f6.1,1x,f5.1,1x,f6.1)

      call close_ncfile( fid )

      if (.not. use_grid_file) then
         deallocate ( lon ) 
         deallocate ( lat ) 
         deallocate ( area )
      endif
      deallocate ( timeval )
      if (use_avg_file) deallocate(timeval_avg)
      deallocate ( slp )
      deallocate ( slp_s )
      if (warm_core_check) deallocate (  tm )
      deallocate ( vort  )
      deallocate ( u_ref )
      deallocate ( v_ref )
      deallocate ( u_steer )
      deallocate ( v_steer )
      deallocate ( wk )
      deallocate ( id )
      if (len(trim(idfile)) > 0) then
         deallocate ( id_out )
         call close_ncfile(idfileid)
      endif
      deallocate ( centers )

      if (warm_core_check) then
         deallocate(id_warm)
         deallocate(warm_core_center)
      endif

      if (use_avg_file) then
         deallocate(precip)
      endif

   end do FILELOOP

   if (use_wind_landmask .or. use_wind_oceanmask) deallocate(landtopomask_wind)

   if (useregionmask) then
      deallocate( regions )
      deallocate( regions_snow )
      deallocate( regions_wind )
   endif

   close(13)
   if ( useregionmask) then
      close(14)
      close(15)
   endif

 contains

      !This subroutine is not used and may be out of date
 subroutine compute_vort(im, jm, us, vs,  vort, lon, lat)

    implicit none
    integer, intent(in):: im, jm
    real*4,    intent(in):: lon(im,jm)
    real*4,    intent(in):: lat(im,jm)
    real*4,    intent(in):: us(im,jm), vs(im,jm)
    real*4,    intent(out):: vort(im,jm)
! Local:
    real:: vdy(im+1), udx(im,jm)
    real:: d_area, deg2rad, pi,dlamda
    integer i,j


     pi = 4.*atan(1.)
     deg2rad = pi/180.

     !NOTE: Assumes regular lat-lon grid
     dlamda = 2.*pi/real(im)

      do j=2,jm
         do i=1,im
            udx(i,j) = 0.5*(us(i,j-1)+us(i,j))*cos( deg2rad*0.5*(lat(i,j-1)+lat(i,j)) )*dlamda
         enddo
      enddo  

      vort = 0.    ! ignore poles (fro TC detection)
      do j=2,jm-1
         do i=2,im
            vdy(i) = 0.5*(vs(i-1,j)+vs(i,j))*0.5*(lat(i,j+1)-lat(i,j-1))*deg2rad
         enddo
         i=1
            vdy(i) = 0.5*(vs(im,j)+vs(i,j))*0.5*(lat(i,j+1)-lat(i,j-1))*deg2rad
         i=im+1
            vdy(i) = 0.5*(vs(i-1,j)+vs(1,j))*0.5*(lat(i,j+1)-lat(i,j-1))*deg2rad
         do i=1,im
            d_area = dlamda*cos( deg2rad*lat(i,j) )*0.5*(lat(i,j+1)-lat(i,j-1))*deg2rad
            vort(i,j) = (udx(i,j)-udx(i,j+1)-vdy(i)+vdy(i+1))/ d_area
         enddo
      enddo  

  end subroutine compute_vort

!This subroutine is not used
 subroutine locate_center_slp(found, im, jm,  slp, tm, lon, lat,  &
                              lon_b, lon_e, lat_b, lat_e, tm_c,   &
                              i0, j0, xx, yy, maxSLPthreshold) 

    implicit none
    integer, intent(in):: im, jm
    real,    intent(in):: lon(im,jm)
    real,    intent(in):: lat(im,jm)
    real*4,    intent(in):: slp(im,jm), tm(im,jm)
    real,    intent(in):: tm_c, lon_b, lon_e, lat_b, lat_e
    real,    intent(inout):: maxSLPthreshold
    real,    intent(out):: xx, yy
    integer, intent(out):: i0, j0
    logical, intent(out):: found
! Local:
    integer i,j

      found = .false.
      i0 = 1;    j0 = 1
      xx = -100.
      yy = -100.

      do j=2,jm-1
         do i=2,im-1
            if ( lat(i,j)>lat_b .and. lat(i,j)<lat_e .and. lon(i,j)> lon_b .and. lon(i,j)<lon_e ) then
!------------------
! Locate TC center:
!------------------
! Using SLP and warm core temp:
! Check if a local min.
               if ( slp(i,j)<min(maxSLPthreshold, slp(i-1,j),slp(i+1,j),slp(i,j-1),slp(i,j+1)) .and. ALL(abs(slp(i-1:i+1,j-1:j+1))<1.e8)) then
                  if ( tm(i,j)>tm_c ) then
                    maxSLPthreshold = slp(i,j)
                      xx = lon(i,j)
                      yy = lat(i,j)
                      i0 = i
                      j0 = j
                      found = .true.
                  endif
               endif
            endif
         enddo     ! i-loop
      enddo        ! j-loop

  end subroutine locate_center_slp

  subroutine find_extrema(field, id, im, jm, storms, stormval, nstorm, minflag, threshold, mincontours, &
#ifdef CHECK_MISSING
       lon, lat, lon_s, lon_e, lat_s, lat_e, max_rad, cint, checkgrad, gradthresh, periodicx, fillvalue)
#else
       lon, lat, lon_s, lon_e, lat_s, lat_e, max_rad, cint, checkgrad, gradthresh, periodicx)
#endif

    !Looking for closed contours: ie contiguous regions 
    implicit none
    integer, intent(in) :: im, jm
    real*4,    intent(in):: field(im,jm)
    integer, intent(out) :: id(im,jm)
    integer, intent(OUT) :: storms(MAXSTORM,2)
    real, intent(OUT)    :: stormval(MAXSTORM)
    real*4,  intent(in) :: threshold
    logical, intent(in) :: minflag !.true. if looking for minimum
    integer, intent(out) :: nstorm
    integer, intent(in) :: mincontours
    real, intent(in) :: lon(im,jm), lat(im,jm), lon_s, lon_e, lat_s, lat_e, max_rad, cint
    logical, intent(in) :: checkgrad, periodicx
    real, intent(IN) :: gradthresh
#ifdef CHECK_MISSING
    real, intent(IN) :: fillvalue
#endif

    integer i,j, n, numfilled, ns
    real*4 :: minf, minf_round
    integer id_tmp(im,jm), id_conf(im,jm)
    logical :: errflag, stillworking
    logical :: finished(MAXSTORM)

    logical, parameter :: DEBUG = .false.
    real:: deg2rad, pi, rdx
    real, dimension(3,3) :: grad
    integer :: istart, iend

     pi = 4.*atan(1.)
     deg2rad = pi/180.

    !Assuming minflag == .true. for now
    if (.not. minflag) then
       print*, 'FIND_EXTREMA: minflag == .false. not implemented'
       stop 800
    endif

    if (periodicx) then
       istart = 1
       iend = im
    else
       istart = 2
       iend = im-1
    endif

    nstorm = 0
    id = 0
    finished = .false.
    stormval = 0.
    storms = 0

    !First pass: identify candidates
    do j=2,jm-1
    do i=istart,iend !periodicity

       if (lat(i,j) > lat_e .or. lat(i,j) < lat_s) cycle
       if (lon(i,j) > lon_e .or. lon(i,j) < lon_s) cycle

       !check to see if center is a local minimum and below our threshold, and optionally that the gradient is sufficiently large
       if (i == 1) then
#ifdef CHECK_MISSING
          if (ANY(field(1:2,j-1:j+1) == fillvalue) .or. ANY(field(im,j-1:j+1) == fillvalue)) cycle
#endif
          if (field(i,j) > min(minval(field(1:2,j-1:j+1)),minval(field(im,j-1:j+1))) .or. field(i,j) > threshold) cycle
       else if (i == im) then
#ifdef CHECK_MISSING
          if (ANY(field(1,j-1:j+1) == fillvalue) .or. ANY(field(im-1:im,j-1:j+1) == fillvalue)) cycle
#endif
          if (field(i,j) > min(minval(field(im-1:im,j-1:j+1)),minval(field(1,j-1:j+1))) .or. field(i,j) > threshold) cycle
       else
#ifdef CHECK_MISSING
          if (ANY(field(i-1:i+1,j-1:j+1) == fillvalue)) cycle
#endif
          if (field(i,j) > minval(field(i-1:i+1,j-1:j+1)) .or. field(i,j) > threshold) cycle
          if (checkgrad) then
             grad = field(i-1:i+1,j-1:j+1) - field(i,j)
             !NOTE: assume dx == dy
             rdx = 1./(RADIUS/1000.*cos(lat(i,j)*deg2rad)*deg2rad*abs(lon(i+1,j)-lon(i-1,j))*0.5)
             !print*, 1./rdx, maxval(grad), gradthresh
             if (maxval(grad)*rdx < gradthresh) cycle
          endif
       endif

       minf = field(i,j)

       !First flood_fill
       id_tmp = id
       numfilled = 0
       call flood_fill(field,id_tmp,i,j,i,j,im,jm,minf+cint,minf,nstorm+1,.true.,errflag, numfilled, lon, lat, max_rad, 1, minflag)
       if (errflag) cycle

       if (debug) print*, 'POSSIBLE STORM FOUND ', i, j, nstorm + 1, minf, minf+cint

       id_conf = id
       !Now look for further filled contours. When these fail that only means
       !that we can find no further closed contours.
       do n=1,mincontours
          !If the next contour is above the threshold, then quit
          if (minf+cint*real(n) > threshold) then
             errflag = .true.
             exit
          endif

          numfilled = 0
          id_tmp = id
          call flood_fill(field,id_tmp,i,j,i,j,im,jm,minf+cint*real(n),minf,nstorm+1,.true.,errflag, numfilled, lon, lat, max_rad, 1, minflag)
          !call flood_fill(field,id_tmp,i,j,i,j,im,jm,minf+cint*real(n),-1.e25,nstorm+1,.true.,errflag, numfilled)
          if (errflag) then
             if (debug) print*, 'CONTOUR NOT FOUND', i, j, n, minf+cint*n
             exit
          endif
          id_conf = id_tmp
          id_tmp = id
          if (debug) print*, 'STORM CONTOUR FOUND ', i, j, n, minf+cint*n

       enddo
!!$       id = id_conf
       if (.not. errflag) then
          id = id_conf
          

          !If the first set of flood_fills succeed, ie we found n closed contours, then we have found a storm.
          nstorm = nstorm + 1
          
          storms(nstorm,1) = i
          storms(nstorm,2) = j
          stormval(nstorm) = minf


          if (debug) print*, 'STORM VERIFIED ', nstorm, i, j, minf, n
          
       else

          if (debug) print*, 'CONTOURS ABOVE THRESHOLD, STORM NOT VERIFIED'
          if (debug) print*, '     ', nstorm, i, j, minf+cint*real(n), threshold

       end if


    enddo
    enddo
    
    !second pass: fill out contours as much as possible
    !How to get flood_fill to handle already-filled areas?
    n = mincontours+1
    do 
       stillworking = .false.
       do ns=1,nstorm


          minf = stormval(ns)

          if (finished(ns)) cycle
          if (minf+cint*real(n) > threshold) then
             finished(ns) = .true.
             cycle
          endif

          id_tmp = id
          do j=1,jm
          do i=1,im
             if (id_tmp(i,j) == ns) id_tmp(i,j) = 0
          enddo
          enddo

          i = storms(ns,1)
          j = storms(ns,2)

          !print*, 'DEBUG: ', n, ns, i, j, minf, minf+cint*real(n), finished(ns), id(i,j)
          numfilled = 0
          errflag = .false.
          call flood_fill(field,id_tmp,i,j,i,j,im,jm, &
               minf+cint*real(n),minf,ns,.false., errflag, numfilled, lon, lat, max_rad, 1, minflag)
          if (errflag .or. numfilled == 0) then
             finished(ns) = .true.
             cycle
          endif
          !print*, 'numfilled = ', numfilled
          stillworking = .true.
          if (debug) print*, 'ADDITIONAL CONTOUR FOUND ', ns, storms(ns,:),  minf+cint*real(n)
          id = id_tmp
       enddo
       if (.not. stillworking) exit
       n = n+1
    enddo

  end subroutine find_extrema

  recursive subroutine flood_fill(field,id,i,j,io,jo,im,jm,upvalue,lowvalue,fillvalue,norefill,errflag,numfilled, &
       lon, lat, max_rad, depth, increasing)
    
    !see http://en.wikipedia.org/wiki/Flood_fill

    implicit none
    integer, intent(in) :: im, jm, i, j, io, jo !io, jo == 'origin'
    real*4,    intent(in):: field(im,jm)
    integer,    intent(inout):: id(im,jm)
    real*4, intent(in) :: upvalue, lowvalue
    real, intent(IN) :: lon(im, jm), lat(im,jm)
    real, intent(IN) :: max_rad
    integer, intent(in) :: fillvalue
    logical, intent(in) :: norefill, increasing
    logical, intent(out) :: errflag
    integer, intent(inout) :: numfilled
    integer, intent(in) :: depth

    real :: rad
    integer :: xdist

    errflag = .false.

    rad = great_circle_dist((/lon(i,j),lat(i,j)/),(/lon(io,jo),lat(io,jo)/),RADIUS)
    if (rad > max_rad) then
       errflag = .true.
       return
    endif

    if (id(i,j) == -fillvalue) then
       id(i,j) = fillvalue
    endif

    if (id(i,j) /= fillvalue) then
       
       if (increasing) then
          if (field(i,j) >= upvalue) return
       else
          if (field(i,j) <= lowvalue) return
       endif



       if (id(i,j) /= 0) then
          !Pixel is already assigned
          !print*, 'REGIONS INTERSECTING: ', i, j, id(i,j), fillvalue, upvalue
          errflag = .true.
          return
       endif


       if (increasing) then
          !If lowvalue = -1.e20. then we turn this check off
          if (field(i,j) < lowvalue .and. lowvalue > -1.e20) then
             !print*, 'BELOW LOWVALUE', i,j, field(i,j), lowvalue
             errflag = .true.
             return
          endif
       else
          if (field(i,j) > upvalue .and. upvalue > 1.e20) then
             !print*, 'ABOVE UPVALUE', i,j, field(i,j), upvalue
             errflag = .true.
             return
          endif
       endif

       !Return an error if at the edge of the domain without finding a contour
       if (((i == 1 .or. i == im) .and. .not. periodicx) .or. j == 1 .or. j == jm .or. abs(field(i,j) > 1.e8)) then
          !print*, 'DOMAIN EDGE ', i, j, fillvalue
          errflag = .true.
          return
       endif

       id(i,j) = fillvalue
       numfilled = numfilled + 1

    endif
    
    !For each adjacent cell: if fill value it is already done; if another
    !nonzero value we have hit another region, an error; else call flood_fill for it\
    !Alternately, if 'refilling' (attempting to expand the fill region) always proceed
    !AWAY from the origin instead of looking for whether the next point has been filled.

    if (i == im .and. periodicx) then

       if (id(1,j) /= fillvalue) then
          call flood_fill(field,id,1,j,io,jo,im,jm,upvalue,lowvalue,fillvalue,norefill,errflag,numfilled, lon, lat, max_rad, depth+1, increasing)
          if (errflag) return
       endif       

    else

       if (id(i+1,j) /= fillvalue) then
          call flood_fill(field,id,i+1,j,io,jo,im,jm,upvalue,lowvalue,fillvalue,norefill,errflag,numfilled, lon, lat, max_rad, depth+1, increasing)
          if (errflag) return
       endif

    endif
       
    if (i == 1 .and. periodicx) then
       if (id(im,j) /= fillvalue) then
          call flood_fill(field,id,im,j,io,jo,im,jm,upvalue,lowvalue,fillvalue,norefill,errflag,numfilled, lon, lat, max_rad, depth+1, increasing)
          if (errflag) return
       endif
    else
       if (id(i-1,j) /= fillvalue) then
          call flood_fill(field,id,i-1,j,io,jo,im,jm,upvalue,lowvalue,fillvalue,norefill,errflag,numfilled, lon, lat, max_rad, depth+1, increasing)
          if (errflag) return
       endif
    endif
       
       if (id(i,j+1) /= fillvalue) then
          call flood_fill(field,id,i,j+1,io,jo,im,jm,upvalue,lowvalue,fillvalue,norefill,errflag,numfilled, lon, lat, max_rad, depth+1, increasing)
          if (errflag) return
       endif
       
       if (id(i,j-1) /= fillvalue) then
          call flood_fill(field,id,i,j-1,io,jo,im,jm,upvalue,lowvalue,fillvalue,norefill,errflag,numfilled, lon, lat, max_rad, depth+1, increasing)
          if (errflag) return
       endif

  end subroutine flood_fill

#ifdef CHECK_MISSING
 subroutine split_smoother(im, jm, ai, ao, periodicx, fillvalue)
#else
 subroutine split_smoother(im, jm, ai, ao, periodicx)
#endif

    implicit none
    integer, intent(in):: im, jm
    real*4,    intent(in):: ai(im,jm)
#ifdef CHECK_MISSING
    real*4,    intent(in):: fillvalue
#endif
    real*4,   intent(out):: ao(im,jm)
    logical, intent(in) :: periodicx
! Local:
    real*4:: wk(im,jm)
    integer i,j

    ao = 0.
    wk = 0.

! X-
      do j=1,jm
         if (periodicx) then
#ifdef CHECK_MISSING
            !! check for fillvalue
            if (ai(1,j) == fillvalue) then
               wk(1,j) = fillvalue
            elseif (ai(im,j) == fillvalue .or. ai(2,j) == fillvalue) then
               wk(1,j) = ai(1,j)
            else
               wk(1,j) = 0.25*(ai(im,j) + 2.*ai(1,j) + ai(2,j))
            end if
#else
            wk(1,j) = 0.25*(ai(im,j) + 2.*ai(1,j) + ai(2,j))
#endif
         endif
         do i=2,im-1
#ifdef CHECK_MISSING
            !! check for fillvalue
            if (ai(i,j) == fillvalue) then
               wk(i,j) = fillvalue
            elseif (ai(i-1,j) == fillvalue .or. ai(i+1,j) == fillvalue) then
               wk(i,j) = ai(i,j)
            else
               wk(i,j) = 0.25*(ai(i-1,j) + 2.*ai(i,j) + ai(i+1,j))
            end if
#else
            wk(i,j) = 0.25*(ai(i-1,j) + 2.*ai(i,j) + ai(i+1,j))
#endif
         enddo     ! i-loop
         if (periodicx) then
#ifdef CHECK_MISSING
            !! check for fillvalue
            if (ai(im,j) == fillvalue) then
               wk(im,j) = fillvalue
            elseif (ai(im-1,j) == fillvalue .or. ai(1,j) == fillvalue) then
               wk(im,j) = ai(im,j)
            else
               wk(im,j) = 0.25*(ai(im-1,j) + 2.*ai(im,j) + ai(1,j))
            end if
#else
            wk(im,j) = 0.25*(ai(im-1,j) + 2.*ai(im,j) + ai(1,j))
#endif
         endif
      enddo        ! j-loop
! Y-
      do j=2,jm-1
         do i=1,im
#ifdef CHECK_MISSING
            if (wk(i,j) == fillvalue) then
               ao(i,j) = fillvalue
            elseif (wk(i,j-1) == fillvalue .or. wk(i,j+1) == fillvalue) then
               ao(i,j) = wk(i,j)
            else
               ao(i,j) = 0.25*(wk(i,j-1) + 2.*wk(i,j) + wk(i,j+1))
            end if
#else
            ao(i,j) = 0.25*(wk(i,j-1) + 2.*wk(i,j) + wk(i,j+1))
#endif
         enddo     ! i-loop
      enddo        ! j-loop


  end subroutine split_smoother

 real function great_circle_dist( q1, q2, radius )
      real, intent(IN)           :: q1(2), q2(2)
      real, intent(IN), optional :: radius

      real*4:: p1(2), p2(2)
      real*4:: beta, pi
      integer n

      pi = 4.*atan(1.)

      !Note that the data used in this program uses DEGREES which must be converted to radians

      do n=1,2
         p1(n) = q1(n)*pi/180.
         p2(n) = q2(n)*pi/180.
      enddo

      beta = asin( sqrt( sin((p1(2)-p2(2))/2.)**2 + cos(p1(2))*cos(p2(2))*   &
                         sin((p1(1)-p2(1))/2.)**2 ) ) * 2.

      if ( present(radius) ) then
           great_circle_dist = radius * beta
      else
           great_circle_dist = beta   ! Returns the angle
      endif

  end function great_circle_dist


 subroutine intp_great_circle(beta, p1, p2, x_o, y_o)
 real, intent(in)::  beta    ! [0,1]
 real, intent(in)::  p1(2), p2(2)
 real, intent(out):: x_o, y_o     ! between p1 and p2 along GC; output to degrees
!------------------------------------------
    real:: pm(2)
    real:: prad1(2), prad2(2) !radians; input is in degrees
    real:: e1(3), e2(3), e3(3)
    real:: s1, s2, s3, dd, alpha
    real*4::pi

      pi = 4.*atan(1.)

      prad1 = p1*pi/180.
      prad2 = p2*pi/180.

      call latlon2xyz(prad1, e1)
      call latlon2xyz(prad2, e2)

       alpha = 1. - beta

       s1 = alpha*e1(1) + beta*e2(1)
       s2 = alpha*e1(2) + beta*e2(2)
       s3 = alpha*e1(3) + beta*e2(3)

       dd = sqrt( s1**2 + s2**2 + s3**2 )

       e3(1) = s1 / dd
       e3(2) = s2 / dd
       e3(3) = s3 / dd

      call cart_to_latlon(1, e3, pm(1), pm(2))

      x_o = pm(1)*180./pi
      y_o = pm(2)*180./pi

 end subroutine intp_great_circle


 subroutine latlon2xyz(p, e)
!
! Routine to map (lon, lat) to (x,y,z)
!
 real, intent(in) :: p(2)
 real, intent(out):: e(3)

 integer n
 real *4:: q(2)
 real *4:: e1, e2, e3

    do n=1,2
       q(n) = p(n)
    enddo

    e1 = cos(q(2)) * cos(q(1))
    e2 = cos(q(2)) * sin(q(1))
    e3 = sin(q(2))
!-----------------------------------
! Truncate to the desired precision:
!-----------------------------------
    e(1) = e1
    e(2) = e2
    e(3) = e3

 end subroutine latlon2xyz




 subroutine cart_to_latlon(np, q, xs, ys)
! vector version of cart_to_latlon1
  integer, intent(in):: np
  real, intent(inout):: q(3,np)
  real, intent(inout):: xs(np), ys(np)
! local
  real, parameter:: esl=1.e-10
  real *4:: p(3)
  real *4:: dist, lat, lon
  integer i,k
    real*4::pi

      pi = 4.*atan(1.)

  do i=1,np
     do k=1,3
        p(k) = q(k,i)
     enddo
     dist = sqrt(p(1)**2 + p(2)**2 + p(3)**2)
     do k=1,3
        p(k) = p(k) / dist
     enddo

     if ( (abs(p(1))+abs(p(2)))  < esl ) then
          lon = 0.
     else
          lon = atan2( p(2), p(1) )   ! range [-pi,pi]
     endif

     if ( lon < 0.) lon = 2.*pi + lon
     lat = asin(p(3))

     xs(i) = lon
     ys(i) = lat
! q Normalized:
     do k=1,3
        q(k,i) = p(k)
     enddo
  enddo

 end  subroutine cart_to_latlon




      subroutine open_ncfile( iflnm, ncid )
      implicit         none
#include <netcdf.inc>
      character*(*), intent(in)::  iflnm
      integer, intent(out)::      ncid
      integer::  status

      status = nf_open (iflnm, NF_NOWRITE, ncid)
      if (status .ne. NF_NOERR) then
         print*, iflnm
         call handle_err(status)
      endif


      end subroutine open_ncfile


      subroutine close_ncfile( ncid )
      implicit         none
#include <netcdf.inc>
      integer, intent(in)::      ncid
      integer::  status

      status = nf_close (ncid)
      if (status .ne. NF_NOERR) call handle_err(status)


      end subroutine close_ncfile

      subroutine get_ncdim1( ncid, var1_name_in, im, var1_name_out )
      implicit         none
#include <netcdf.inc>
      integer, intent(in):: ncid
      character*(*), intent(in)::  var1_name_in
      character*(*), intent(out), optional::  var1_name_out
      integer, intent(out):: im
      character(len=100) :: var1_name
      integer::  status, var1id

      !Here we allow alternate names to be used for some
      !well-known dimensions, if the original isn't found.

      var1_name = var1_name_in

      status = nf_inq_dimid (ncid, var1_name, var1id)
      if (status .eq. NF_EBADDIM) then
         !Find the RIGHT name
         do
            select case(trim(var1_name))
            case('lon')
               var1_name = 'grid_xt'
            case('grid_xt')
               var1_name = 'lon'
            case('lat')
               var1_name = 'grid_yt'
            case('grid_yt')
               var1_name = 'lat'
            end select
            

            if (trim(var1_name) .eq. trim(var1_name_in)) exit !If we come back around to the original name
                                                        !or can't find an altername name, give up;
                                                        !the next handle_err call will then stop execution
            status = nf_inq_dimid (ncid, trim(var1_name), var1id)
            if (status .eq. NF_NOERR) print*, 'Could not find dimension ', trim(var1_name_in), &
                 ', using ', trim(var1_name), ' instead'
            if (status .ne. NF_EBADDIM) exit
         enddo
      endif


      if (status .ne. NF_NOERR) call handle_err(status)

      status = nf_inq_dimlen (ncid, var1id, im)
      if (status .ne. NF_NOERR) call handle_err(status)

      if (present(var1_name_out)) var1_name_out = trim(var1_name)

      end subroutine get_ncdim1

      subroutine get_var1_real( ncid, var1_name, im, var1 )
      implicit         none
#include <netcdf.inc>
      integer, intent(in):: ncid
      character*(*), intent(in)::  var1_name
      integer, intent(in):: im
      real*4, intent(out):: var1(im)

      integer::  status, var1id

      status = nf_inq_varid (ncid, var1_name, var1id)
      if (status .ne. NF_NOERR) call handle_err(status)
!     write(*,*) 'Got var1id', var1id, var1_name

      status = nf_get_var_real (ncid, var1id, var1)
      if (status .ne. NF_NOERR) call handle_err(status)


      end subroutine get_var1_real



      subroutine get_var1_double( ncid, var1_name, im, var1 )
      implicit         none
#include <netcdf.inc>
      integer, intent(in):: ncid
      character*(*), intent(in)::  var1_name
      integer, intent(in):: im
      real*8, intent(out):: var1(im)

      integer::  status, var1id

      status = nf_inq_varid (ncid, var1_name, var1id)
      if (status .ne. NF_NOERR) call handle_err(status)
!     write(*,*) 'Got var1id', var1id

      status = nf_get_var_double (ncid, var1id, var1)
      if (status .ne. NF_NOERR) call handle_err(status)


      end subroutine get_var1_double



      subroutine get_var2_real( ncid, var2_name, im, jm, var2)
      implicit         none
#include <netcdf.inc>
      integer, intent(in):: ncid
      character*(*), intent(in)::  var2_name
      integer, intent(in):: im, jm
      real*4, intent(out):: var2(im,jm)

      integer::  status, var2id

      status = nf_inq_varid (ncid, var2_name, var2id)
      if (status .ne. NF_NOERR) call handle_err(status)
!     write(*,*) 'Got var2id', var2id

      status = nf_get_var_real (ncid, var2id, var2)
      if (status .ne. NF_NOERR) call handle_err(status)


      end subroutine get_var2_real


!!$      !Not used. Commented out to avoid problems
!!$      subroutine get_var2_double( ncid, var2_name, im, jm, var2 )
!!$      implicit         none
!!$#include <netcdf.inc>
!!$      integer, intent(in):: ncid
!!$      character*(*), intent(in)::  var2_name
!!$      integer, intent(in):: im, jm
!!$      real*8, intent(out):: var2(im,jm)
!!$
!!$      integer::  status, var2id
!!$
!!$      status = nf_inq_varid (ncid, var2_name, var2id)
!!$      if (status .ne. NF_NOERR) call handle_err(status)
!!$!     write(*,*) 'Got var2id', var2id
!!$
!!$      status = nf_get_var_double (ncid, var2id, var2)
!!$      if (status .ne. NF_NOERR) call handle_err(status)
!!$
!!$
!!$      end subroutine get_var2_double


      subroutine get_var3_double( ncid, var3_name, im, jm, km, var3 )
      implicit         none
#include <netcdf.inc>
      integer, intent(in):: ncid
      character*(*), intent(in)::  var3_name
      integer, intent(in):: im, jm, km
      real*8, intent(out):: var3(im,jm,km)

      integer::  status, var3id

      status = nf_inq_varid (ncid, var3_name, var3id)

      if (status .ne. NF_NOERR) call handle_err(status)
!     write(*,*) var3_name, ' Got var3id=', var3id

      status = nf_get_var_double (ncid, var3id, var3)
      if (status .ne. NF_NOERR) call handle_err(status)

      end subroutine get_var3_double

!------------------------------------------------------------------------
      subroutine get_var4_double( ncid, var4_name, im, jm, km, nt, var4 )
      implicit         none
#include <netcdf.inc>
      integer, intent(in):: ncid
      character*(*), intent(in)::  var4_name
      integer, intent(in):: im, jm, km, nt
      real*8, intent(out):: var4(im,jm,km,1)
      integer::  status, var4id
!
      integer:: start(4), icount(4) 

      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = nt

      icount(1) = im    ! all range
      icount(2) = jm    ! all range
      icount(3) = km    ! all range
      icount(4) = 1     ! one time level at a time

!     write(*,*) 'Within get_var4_double', var4_name

      status = nf_inq_varid (ncid, var4_name, var4id)
!     write(*,*) var4_name, ' Got var4id=', var4id

      status = nf_get_vara_double(ncid, var4id, start, icount, var4)

      if (status .ne. NF_NOERR) call handle_err(status)

      end subroutine get_var4_double
!------------------------------------------------------------------------

      subroutine get_real2( ncid, var2_name, im, jm, var2 )
! This is for single-time-level 2D var
      implicit         none
#include <netcdf.inc>
      integer, intent(in):: ncid
      character*(*), intent(in)::  var2_name
      integer, intent(in):: im, jm
      real*4, intent(out):: var2(im,jm)
      integer::  status, var2id
      integer:: start(2), icount(2)
      integer:: i,j

      start(1) = 1
      start(2) = 1

      icount(1) = im
      icount(2) = jm

      status = nf_inq_varid (ncid, var2_name, var2id)
      status = nf_get_vara_real(ncid, var2id, start, icount, var2)

      if (status .ne. NF_NOERR) call handle_err(status)

      end subroutine get_real2

      subroutine get_int2( ncid, var2_name, im, jm, var2 )
! This is for single-time-level 2D var
      implicit         none
#include <netcdf.inc>
      integer, intent(in):: ncid
      character*(*), intent(in)::  var2_name
      integer, intent(in):: im, jm
      integer*4, intent(out):: var2(im,jm)
      integer::  status, var2id
      integer:: start(2), icount(2)
      integer:: i,j

      start(1) = 1
      start(2) = 1

      icount(1) = im
      icount(2) = jm

      status = nf_inq_varid (ncid, var2_name, var2id)
      status = nf_get_vara_int(ncid, var2id, start, icount, var2)

      if (status .ne. NF_NOERR) call handle_err(status)

    end subroutine get_int2
!------------------------------------------------------------------------

      subroutine get_real3( fname, ncid, var4_name, im, jm, nt, var4, one_variable_per_file )
! This is for multi-time-level 2D var
      implicit         none
#include <netcdf.inc>
      integer, intent(inout):: ncid
      character*(*), intent(in)::  fname, var4_name
      integer, intent(in):: im, jm, nt
      real*4, intent(out):: var4(im,jm)
      logical, intent(IN), optional :: one_variable_per_file
      character*190 :: fname_var
      integer::  status, var4id
      integer:: start(3), icount(3)
      integer:: i,j

      start(1) = 1
      start(2) = 1
      start(3) = nt

      icount(1) = im
      icount(2) = jm
      icount(3) = 1

      if (present(one_variable_per_file)) then
         if (one_variable_per_file) then
            call close_ncfile(ncid)
            call variable_file(fname_var, fname, var4_name)
            !print*, 'OPENING FILE ', trim(fname_var)
            call open_ncfile(fname_var, ncid)            
         endif
      endif

      status = nf_inq_varid (ncid, var4_name, var4id)
      status = nf_get_vara_real(ncid, var4id, start, icount, var4)

      if (status .ne. NF_NOERR) call handle_err(status)

      end subroutine get_real3
!------------------------------------------------------------------------

!------------------------------------------------------------------------
      subroutine get_var4_real( ncid, var4_name, im, jm, km, nt, var4 )
      implicit         none
#include <netcdf.inc>
      integer, intent(in):: ncid
      character*(*), intent(in)::  var4_name
      integer, intent(in):: im, jm, km, nt
      real*4:: wk4(im,jm,km,4)
      real*4, intent(out):: var4(im,jm)
      integer::  status, var4id
      integer:: start(4), icount(4) 
      integer:: i,j

      start(1) = 1
      start(2) = 1
      start(3) = 1
      start(4) = nt

      icount(1) = im    ! all range
      icount(2) = jm    ! all range
      icount(3) = km    ! all range
      icount(4) = 1     ! one time level at a time

!     write(*,*) nt, 'Within get_var4_double: ', var4_name

      status = nf_inq_varid (ncid, var4_name, var4id)
!     write(*,*) '#1', status, ncid, var4id

      status = nf_get_vara_real(ncid, var4id, start, icount, var4)
!     status = nf_get_vara_real(ncid, var4id, start, icount, wk4)
!     write(*,*) '#2', status, ncid, var4id

      do j=1,jm
      do i=1,im
!        var4(i,j) = wk4(i,j,1,nt)
      enddo
      enddo

      if (status .ne. NF_NOERR) call handle_err(status)

      end subroutine get_var4_real
!------------------------------------------------------------------------

      subroutine get_var_att_str(ncid, var_name, att_name, att, msg)
      implicit none
#include <netcdf.inc>
      integer, intent(in):: ncid
      character*(*), intent(in)::  var_name, att_name
      character*(*), intent(out)::  att
      logical, intent(IN) :: msg

      integer::  status, varid

      status = nf_inq_varid (ncid, var_name, varid)
      status = nf_get_att_text(ncid, varid, att_name, att)

      if (status .ne. NF_NOERR) then
         if (msg) print*, 'Attribute ', trim(att_name), ' not found for variable ', trim(var_name)
         if (msg) print*, ' setting to default value ', att
         return
      end if

      if (status .ne. NF_NOERR) call handle_err(status)

      end subroutine get_var_att_str

      subroutine get_var_att_real(ncid, var_name, att_name, att, msg)
      implicit none
#include <netcdf.inc>
      integer, intent(in):: ncid
      character*(*), intent(in)::  var_name, att_name
      real, intent(inout)::  att
      logical, intent(IN) :: msg

      integer::  status, varid, xtype, len

      status = nf_inq_varid (ncid, var_name, varid)
      status = nf_inq_att   (ncid, varid, att_name, xtype, len) !last 2 args are dummies

      if (status .ne. NF_NOERR) then
         if (msg) print*, 'Attribute ', trim(att_name), ' not found for variable ', trim(var_name)
         if (msg) print*, ' setting to default value ', att
         return
      end if

      status = nf_get_att_real(ncid, varid, att_name, att)

      if (status .ne. NF_NOERR) call handle_err(status)

    end subroutine get_var_att_real

    subroutine set_var_att_int(ncid, var_name, att_name, att)
      implicit none
#include <netcdf.inc>
      integer, intent(in):: ncid
      character*(*), intent(in)::  var_name, att_name
      integer, intent(in)::  att

      integer::  status, varid, xtype, len

      status = nf_inq_varid (ncid, var_name, varid)
      status = nf_put_att_real(ncid, varid, att_name, NF_INT, 1, att)

      if (status .ne. NF_NOERR) call handle_err(status)

    end subroutine set_var_att_int

    subroutine handle_err(status)
      implicit         none
#     include          <netcdf.inc>
      integer          status

      if (status .ne. nf_noerr) then
         print*, status
        print *, nf_strerror(status)
        call tracebackqq
        stop 'Stopped'
      endif

    end subroutine handle_err

 subroutine pmaxmin4( qname, a, imax, jmax, fac )

      character*(*)  qname
      integer imax, jmax
      integer i, j
      real*4 a(imax,jmax)

      real qmin(jmax), qmax(jmax)
      real pmax, pmin
      real fac                     ! multiplication factor

      do j=1,jmax
         pmax = a(1,j)
         pmin = a(1,j)
         do i=2,imax
            pmax = max(pmax, a(i,j))
            pmin = min(pmin, a(i,j))
         enddo
         qmax(j) = pmax
         qmin(j) = pmin
      enddo
!
! Now find max/min of amax/amin
!
            pmax = qmax(1)
            pmin = qmin(1)
         do j=2,jmax
            pmax = max(pmax, qmax(j))
            pmin = min(pmin, qmin(j))
         enddo

      write(*,*) qname, ' max = ', pmax*fac, ' min = ', pmin*fac

 end subroutine pmaxmin4

 subroutine pmaxmin( qname, a, imax, jmax, fac )

      character*(*)  qname
      integer imax, jmax
      integer i, j
      real a(imax,jmax)

      real qmin(jmax), qmax(jmax)
      real pmax, pmin
      real fac                     ! multiplication factor

      do j=1,jmax
         pmax = a(1,j)
         pmin = a(1,j)
         do i=2,imax
            pmax = max(pmax, a(i,j))
            pmin = min(pmin, a(i,j))
         enddo
         qmax(j) = pmax
         qmin(j) = pmin
      enddo
!
! Now find max/min of amax/amin
!
            pmax = qmax(1)
            pmin = qmin(1)
         do j=2,jmax
            pmax = max(pmax, qmax(j))
            pmin = min(pmin, qmin(j))
         enddo

      write(*,*) qname, ' max = ', pmax*fac, ' min = ', pmin*fac

 end subroutine pmaxmin


 subroutine calendar(year, month, day, hour)
      implicit none
      integer, intent(inout) :: year              ! year
      integer, intent(inout) :: month             ! month
      integer, intent(inout) :: day               ! day
      integer, intent(inout) :: hour
!
! Local variables
!
      integer irem4,irem100
      integer mdays(12)                           ! number day of month 
      data mdays /31,28,31,30,31,30,31,31,30,31,30,31/
!
!***********************************************************************
!******         compute current GMT                               ******
!***********************************************************************
!
!**** consider leap year
!
      irem4    = mod( year, 4 )
      irem100  = mod( year, 100 )
      if( irem4 == 0 .and. irem100 /= 0) mdays(2) = 29
!
!!$      do 
!!$         if( month > 12 ) then
!!$            year   = year + 1
!!$            month  = 1
!!$         else
!!$            exit
!!$         end if
!!$      enddo
!!$
!!$      do 
!!$         if( day > mdays(month) ) then
!!$            day    = day - mdays(month)
!!$            month  = month + 1
!!$         else
!!$            exit
!!$         end if
!!$      enddo
!!$
!!$
!!$      do
!!$         if( hour >= 24 ) then
!!$            day    = day + 1
!!$            hour   = hour - 24
!!$         else
!!$            exit
!!$         end if
!!$      end do
         

      if( hour >= 24 ) then
        day    = day + 1
        hour   = hour - 24
      end if

      if( day > mdays(month) ) then
        day    = day - mdays(month)
        month  = month + 1
      end if
      if( month > 12 ) then
        year   = year + 1
        month  = 1
      end if


      return
  end subroutine calendar

 subroutine get_height_field(is, ie, js, je, km, wz, pt, q, peln, zvir, zsurf)
  integer, intent(in):: is, ie, js, je, km
  real, intent(in):: peln(km+1)
  real, intent(in):: pt(is:ie,js:je,km)
  real, intent(in)::  q(is:ie,js:je,km) ! water vapor
  real, intent(in) :: zsurf(is:ie,js:je)
  real, intent(in):: zvir
  real, intent(out):: wz(is:ie,js:je,km+1)
!
  integer i,j,k
  real gg

      gg  = rdgas * ginv

      do j=js,je
         do i=is,ie
            wz(i,j,km+1) = zsurf(i,j)
         enddo
         do k=km,1,-1
            do i=is,ie
               wz(i,j,k) = wz(i,j,k+1) + gg*pt(i,j,k)*(1.+zvir*q(i,j,k))  &
                          *(peln(k+1)-peln(k))
            enddo
         enddo
      enddo

 end subroutine get_height_field

 subroutine get_pressure_given_height(is, ie, js, je, km, wz, kd, height,   &
                                      ts, peln, a2, fac)

 integer,  intent(in):: is, ie, js, je, km
 integer,  intent(in):: kd           ! vertical dimension of the ouput height
 real, intent(in):: wz(is:ie,js:je,km+1)
 real, intent(in):: ts(is:ie,js:je)
 real, intent(in):: peln(km+1)
 real, intent(in):: height(kd)   ! must be monotonically decreasing with increasing k
 real, intent(out):: a2(is:ie,js:je,kd)      ! pressure (pa)
 real, optional, intent(in):: fac

! local:
 integer n, i,j,k
 real ptmp, tm


 do n=1,kd

!$omp parallel do default(shared) private(ptmp, tm)
    do j=js,je

       do 1000 i=is,ie

         if ( height(n) >= wz(i,j,km+1) ) then
!---------------------
! Search from top down
!---------------------
          do k=1,km
             if( height(n) < wz(i,j,k) .and. height(n) >= wz(i,j,k+1) ) then
! Found it!
                 ptmp = peln(k) + (peln(k+1)-peln(k)) *   &
                       (wz(i,j,k)-height(n)) / (wz(i,j,k)-wz(i,j,k+1))
                 a2(i,j,n) = exp(ptmp)
                 go to 500
             endif
          enddo

         else
!-----------------------------------------
! xtrapolation: mean laspe rate 6.5 deg/km
!-----------------------------------------
                tm = rdgas*ginv*(ts(i,j) + 3.25E-3*(wz(i,j,km)-height(n)))
          a2(i,j,n) = exp( peln(km+1) + (wz(i,j,km+1) - height(n))/tm )
         endif
500      if ( present(fac) ) a2(i,j,n) = fac * a2(i,j,n)
1000   continue
    enddo
 enddo

 end subroutine get_pressure_given_height


!             call interpolate_vertical(isc, iec, jsc, jec, npz,   &
!                                       850.e2, Atm(n)%peln, wk, a2)

 subroutine interpolate_vertical(is, ie, js, je, km, plev, peln, a3, a2)

 integer,  intent(in):: is, ie, js, je, km
 real, intent(in):: peln(is:ie,km+1,js:je)
 real, intent(in):: a3(is:ie,js:je,km)
 real, intent(in):: plev
 real, intent(out):: a2(is:ie,js:je)
! local:
 real pm(km)
 real logp
 integer i,j,k

 logp = log(plev)

!$omp parallel do default(shared) private(pm)
 do j=js,je
    do 1000 i=is,ie

       do k=1,km
          pm(k) = 0.5*(peln(i,k,j)+peln(i,k+1,j))
       enddo

       if( logp <= pm(1) ) then
           a2(i,j) = a3(i,j,1)
       elseif ( logp >= pm(km) ) then
           a2(i,j) = a3(i,j,km)
       else
           do k=1,km-1
              if( logp <= pm(k+1) .and. logp >= pm(k) ) then
                  a2(i,j) = a3(i,j,k) + (a3(i,j,k+1)-a3(i,j,k))*(logp-pm(k))/(pm(k+1)-pm(k))
                  go to 1000
              endif
           enddo
       endif
1000   continue
 enddo

 end subroutine interpolate_vertical

subroutine heapsort(n,nfields,nfields_i,sortfield,ra,ia)

  !From numerical recipes in Fortran, Press et al 1992

  integer, intent(IN) :: n, nfields, nfields_i, sortfield
  real, intent(INOUT) :: ra(n,nfields)
  integer, intent(INOUT) :: ia(n,nfields_i)

  integer :: i,ir,j,l
  real rra, rraa(nfields)
  integer ira(nfields_i)

  if (n .lt. 2) return
  l=n/2+1
  ir=n
10 continue
  if (l .gt. 1) then
     l = l-1
     rra=ra(l,sortfield)
     rraa=ra(l,:)
     ira=ia(l,:)
  else
     rra=ra(ir,sortfield)
     rraa=ra(ir,:)
     ira=ia(ir,:)
     ra(ir,:) = ra(1,:)
     ia(ir,:) = ia(1,:)
     ir=ir-1
     if (ir .eq. 1) then
        ra(1,:) = rraa
        ia(1,:) = ira
        return
     endif
  endif
  i=l
  j=l+l
20 if (j .le. ir) then
     if (j .lt. ir) then
        if (ra(j,sortfield) .lt. ra(j+1,sortfield)) j=j+1
     endif
     if (rra .lt. ra(j,sortfield)) then
        ra(i,:) = ra(j,:)
        ia(i,:) = ia(j,:)
        i = j
        j = j+j
     else
        j=ir+1
     endif
     goto 20
  endif

  ra(i,:) = rraa
  ia(i,:) = ira
  goto 10


end subroutine heapsort

subroutine variable_file(filename, filepattern, variable)

  character*(*), intent(IN) :: filepattern, variable
  character*(*), intent(OUT) :: filename

  integer :: stringstart

  stringstart = index(filepattern, '#')
  if (stringstart > 0) then
     filename = filepattern(1:stringstart-1) // trim(variable) // trim(filepattern(stringstart+1:))
  else
     filename = filepattern
  endif

end subroutine variable_file

end program

     
