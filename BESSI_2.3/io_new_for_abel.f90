! Module which is responsible for reading the values from the disk.
!
! I/O for the surface and energy balance model BESSI
! This file contains netcdf writing and reading routines for the climate input and mass and energy balance output
! based on netcdf 4, modeloutput uses standard names following the CF-conventions http://cfconventions.org/
!   Author: Tobias Zolles
!   Developer: Tobias Zolles
!   Mail: tobias.zolles@uib.no 
!   Last Update: 22.09.2020


module io

include 'netcdf.inc'


CONTAINS
    ! =======================================================================================================
    ! READ
    ! =======================================================================================================


    !--------------------------------------------------------------
    ! NetCDF Stuff
    !--------------------------------------------------------------


    function read_variable(filename, NLONS, NLATS, variable)
        ! Reads the variable (8bit real) from a 2D (long x lat, without time) netcdf file.
        !
        ! With inspirations from here:
        ! http://www.unidata.ucar.edu/software/netcdf/examples/programs/simple_xy_rd.f90
        ! http://www.unidata.ucar.edu/software/netcdf/examples/programs/sfc_pres_temp_rd.f

        ! Reads
        use netcdf
        use variables
        implicit none

        character(len=*), intent(in) :: filename
        ! input
        integer, intent(in) :: NLATS
        integer, intent(in) :: NLONS
        character(*), intent(in) :: variable

        ! output
        real(kind=8) :: read_variable(NLONS, NLATS)

        ! Return value, to check if everything was ok
        integer :: retval

        ! This will be the netCDF ID for the file and data variable.
        integer :: ncid, varid

        ! Open the file, Read only
        ! Trim file path: http://stackoverflow.com/questions/15093712/trimming-string-for-directory-path
        retval = nf_open(TRIM(adjustl(filename)), nf_NOWRITE, ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        if(debug > 4 ) then
            print *, "io_read.f90: NetCDF File opened: ", filename
        end if

        !get_the_test
        ! Get the varid of the data variable, based on its name.
        !retval = nf_inq_varid(ncid, 'LANDMASK', varid)
        retval = nf_inq_varid(ncid, variable, varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        if(debug > 4 ) then
            print *, "io_read.f90: NetCDF varid for variable received: ", variable
        end if

        ! Read the data as Array, should also work
        !retval = nf_GET_VARA_DOUBLE(NCID, varid, START, END, landmask_netcdf)
        !if (retval .ne. nf_noerr) call handle_err(retval)
        !read_landmask = real(landmask_netcdf, kind=4)
        ! read the data
        retval = nf_get_var_double(ncid, varid, read_variable)
        if (retval .ne. nf_noerr) call handle_err(retval)
        if(debug > 4 ) then
            print *, "io_read.f90: NetCDF read the data from the file"
        end if

        ! Close the file, freeing all resources.
        retval = nf_close(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        if(debug > 4 ) then
            print *, "io_read.f90: NetCDF file closed"
        end if

        ! If we got this far, everything worked as expected. Yipee!
        if(debug > 0) then
            print *,"*** Successful reading the ",TRIM(adjustl(variable))," from NetCDF file: ", TRIM(adjustl(filename)), ""
        end if
        return
    end function read_variable

    ! Calls in subroutine READ_WATERMASK: 
    ! => clean_watermask (on line <137>)
    function read_watermask(elevation_landmask, NLONS, NLATS, sea_level_local, elevation_bedrock, ice_thickness_local)
        ! Returns a watermask (array, with integers).
        ! elevation_landmask is the elevation of the bedrock. If this elevation is below the sea level, the grid point is marked as water.
        ! But if there is ice on the grid point which is heavier than the water column
        ! (calculated from the difference of the sea level and elevation_bedrock), the grid point will be marked as ice (0)
        ! 
        ! 0 = ice, 1 = water, 2 = no ice, 3 = unstable integration
        ! Variables with the "_local" at the end, are only used in the function.

        use variables
	!use variables_snow
        implicit none
        !character(len=*), intent(in) :: filename
        ! input
        integer, intent(in) :: NLONS
        integer, intent(in) :: NLATS
        real(kind=8), intent(in) :: sea_level_local
        real(kind=8) :: ice_thickness_local(NLONS, NLATS)
        ! output
        integer :: read_watermask(NLONS, NLATS)

        ! Landmask, to get the position of the water
        real(kind=8), intent(in) :: elevation_landmask(NLONS, NLATS)
        real(kind=8), intent(in) :: elevation_bedrock(NLONS, NLATS)
        ! Calculate the ice column in water equivalent
        !real(kind=8), parameter :: rho_ice = 0.910 ! Eismint, Table 1

        ! Loop indexes
        integer :: lat, lon
        ! 0 = ice  1 = water  2 = no ice
        do lat=1,nlats,1
            ! then loop over y-axes
            do lon=1,nlons,1
		! version Michael
                ! first set to water
	        ! 0 = ice  1 = water  2 = no ice		
                read_watermask(lon, lat) = 1
                !if there is ice, check whether it can maintain. if there is no ice, check whether the bedrock is above sealevel
                if(ice_thickness_local(lon, lat).gt. 0.) then
                    if ((ice_thickness_local(lon, lat)*rho_ice/rho_w .gt. sea_level_local - elevation_bedrock(lon, lat)) & 
					.or. (elevation_bedrock(lon, lat) .gt. sea_level_local + shelv_depth )  )then

!.or. (elevation_bedrock(lon, lat)-sea_level_local .gt. shelv_depth )

                        read_watermask(lon, lat) = 0
                    elseif(elevation_bedrock(lon, lat) .gt. sea_level_local+shelv_depth)then ! landmask
                        read_watermask(lon, lat) = 2
                    end if
                elseif(elevation_bedrock(lon, lat) .gt. sea_level_local+shelv_depth)then ! actual bedrock height (results in inland seas)
                !elseif(elevation_landmask(lon, lat) .gt. sea_level_local)then ! height of relaxed initial bedrock
                    read_watermask(lon, lat) = 2
                end if
            end do
        end do

        ! At the end, remove the lakes from the watermask
        call clean_watermask(read_watermask, NLONS, NLATS)

	! read lgm boundary mask (only used to initialize the lgm ice sheets)
	!read_watermask = 0
	!read_watermask = read_variable('/alphadata04/imhof/model/input/topographie/lgmicemask.cdf', NLONS, NLATS, 'LGM_MASK')
	!where(read_watermask(:,:)==2)
	! 	read_watermask(:,:)=1
	!end where



    end function read_watermask

    subroutine clean_watermask(watermask, NLONS, NLATS)
        ! Removes the great lakes from the watermask
        ! Removes isolated water points sourounded by the land.
        implicit none
        integer, intent(in) :: NLONS
        integer, intent(in) :: NLATS
        integer, intent(inout) :: watermask(NLONS, NLATS)
        integer :: ix,iy !local running variables

        ! Remove the lakes in North America
        !if (NLONS .eq. 625) then
        !    watermask(60:300, 120:176) = 2
        !else
        !    watermask(40:150, 60:88) = 2
        !    ! Over whole america
        !    watermask(50:120, 1:80) = 2
        !    watermask(120:130, 42:50) = 2
        !end if
        ! Remove lakes with the size of one grid box
        ! 0 = ice
        ! 1 = water
        ! 2 = no ice
        do ix=2,NLONS-1,1
            do iy=2,NLATS-1,1
                ! In X Direction: 2 water boxes, 1 ice boxes
                if ( (watermask(ix, iy) .eq. 1) .and. &
                    (watermask(ix-1, iy) .ne. 1) .and. &
                    (watermask(ix+1, iy) .ne. 1) .and. &
                    (watermask(ix, iy-1) .ne. 1) .and. &
                    (watermask(ix, iy+1) .ne. 1) ) then

                    watermask(ix:ix, iy:iy) = 2
                endif
            end do
        end do


    end subroutine clean_watermask



    function read_climate(filename, NLONS, NLATS, variable)
        ! Reads the temperature or precipitation data (8bit real) for 365 days out of the given NetCDF File.
        ! variable: prec/temp, units: m/yr and kelvin
        ! Returns 3D Array: Lon, Lat, Day

        ! Reads
        use netcdf
        use variables ! Own Module
        implicit none
        ! input
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: variable
        integer, intent(in) :: NLATS
        integer, intent(in) :: NLONS

        ! Zeitschritte
        ! ndays

        ! output
        ! TODO: Zeitschritt einbinden
        real(kind=8) :: read_climate(NLONS, NLATS, ndays)

        ! Return value, to check if everything was ok
        integer :: retval

        ! This will be the netCDF ID for the file and data variable.
        integer :: ncid, varid

        ! Open the file, Read only
        ! Trim file path: http://stackoverflow.com/questions/15093712/trimming-string-for-directory-path
        retval = nf_open(TRIM(adjustl(filename)), nf_NOWRITE, ncid)

        ! Read Data
        ! Get the varid of the data variable, based on its name.
        retval = nf_inq_varid(ncid, TRIM(adjustl(variable)), varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! read the data
        retval = nf_get_var_double(ncid, varid, read_climate)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Close the file, freeing all resources.
        retval = nf_close(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        if(debug > 0) then
            ! If we got this far, everything worked as expected. Yipee!
            print *,"*** Read ", TRIM(adjustl(variable))," from NetCDF file: ", TRIM(adjustl(filename)), ""
        end if

        return

    end function read_climate



    function read_swrad_func(filename, NLONS, NLATS, variable)
        ! Reads the temperature or precipitation data (8bit real) for 365 days out of the given NetCDF File.
        ! variable: prec/temp, units: m/yr and kelvin
        ! Returns 3D Array: Lon, Lat, Day

        ! Reads
        use netcdf
        use variables ! Own Module
        implicit none
        ! input
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: variable
        integer, intent(in) :: NLATS
        integer, intent(in) :: NLONS

        ! Zeitschritte
        ! ndays

        ! output
        ! TODO: Zeitschritt einbinden
        real(kind=8) :: read_swrad_func(NLONS, NLATS, 12)

        ! Return value, to check if everything was ok
        integer :: retval

        ! This will be the netCDF ID for the file and data variable.
        integer :: ncid, varid

        ! Open the file, Read only
        ! Trim file path: http://stackoverflow.com/questions/15093712/trimming-string-for-directory-path
        retval = nf_open(TRIM(adjustl(filename)), nf_NOWRITE, ncid)

        ! Read Data
        ! Get the varid of the data variable, based on its name.
        retval = nf_inq_varid(ncid, TRIM(adjustl(variable)), varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! read the data
        retval = nf_get_var_double(ncid, varid, read_swrad_func)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Close the file, freeing all resources.
        retval = nf_close(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        if(debug > 0) then
            ! If we got this far, everything worked as expected. Yipee!
            print *,"*** Read ", TRIM(adjustl(variable))," from NetCDF file: ", TRIM(adjustl(filename)), ""
        end if

        return

    end function read_swrad_func





    !--------------------------------------------------------------
    ! Output Directory Stuff
    !--------------------------------------------------------------

    subroutine init_output_directory(directory, identifier)
        ! Creates a supdirectory with a timestamp and the identifier of the experiment
        ! Copies the variable and the initial NetCDF files into it.
        ! Sets the new output directory to the variable directory.

        ! Use the variables
        use variables ! Own Module
        implicit none
        ! input
        character(len=*), intent(inout) :: directory
        character(len=*), intent(in) :: identifier
        integer(kind=4)  stime, tarray(9), time
	integer(kind=4)   VALUES(8)
	character*6 DATE
	character*5 ZONE
	character*10 TTIME


        integer(kind=4)  :: result
	character(len=128) :: constants_string
	character(len=30) :: transient_string

        character(len=256) :: timestamp
        character(len=256) :: command
        character(len=256) :: current_directory
        character(len=256) :: timestart

	stime = time()
        call ltime(stime, tarray)
        ! YYYYMMDDHHMM : http://docs.oracle.com/cd/E19957-01/805-4942/6j4m3r90k/index.html
        write (timestamp, "(I0.4,'_',I0.2,'_',I0.2,'_',I0.2,I0.2,'_',I0.2)") (tarray(6) + 1900), tarray(5)+1, tarray(4),& 
												(tarray(3)), tarray(2),tarray(1)
	!timestamp='pgf00000_'
	!print*,timestamp
	! build directory name
	!CALL DATE_AND_TIME(DATE, TTIME, ZONE, VALUES) 
	!values=1
	!print*, 'val 1-6',values(1) , values(2), values(3), values(5), values(6)
        !write (timestamp, "(I0.4,'_',I0.2,'_',I0.2,'_',I0.2,I0.2,'_',I0.2)") values(1), values(2), values(3), values(5), &
	!												values(6), values(7)

	write (constants_string,"('_ea',F0.2,'_ad',F0.2,'_aw' ,F0.2,'_ai' ,F0.2)") eps_air, albedo_snow_new, albedo_snow_wet, albedo_ice
	!write (eps_string,"F0.2") eps_air


        directory = TRIM(adjustl(directory)) // TRIM(adjustl(identifier)) // '_' // run_foo // '_' // &
        parameters_foo // TRIM(adjustl(timestamp)) //  '/'!& 
			!TRIM(adjustl(timestart)) // TRIM(adjustl(transient_string)) // TRIM(adjustl(constants_string)) //'/'

    print *, TRIM(adjustl(directory))
        ! create directory
        command = 'mkdir ' // directory
        result = system(command)

        if(debug > 0 ) then
            print *, "*** Output directory created: ", TRIM(adjustl(directory))
        end if
	
	! save parameters to txt file
	call save_parameters(directory,timestamp,identifier)


        ! copy variables, depending from which directory the model gets executed.
        if(store_input_variables) then
                CALL getcwd(current_directory)
                if (index(current_directory, 'Debug') .gt. 0) then
                    command = 'cp -p ' // TRIM(adjustl(current_directory)) // '/variables.f90 ' // TRIM(adjustl(directory))
                    result = system(command)
                    command = 'cp -p ' // TRIM(adjustl(current_directory)) // '/io.f90 ' // TRIM(adjustl(directory))
                    result = system(command)
                    command = 'cp -p ' // TRIM(adjustl(current_directory)) // '/smb_emb.f90 ' // TRIM(adjustl(directory))
                    result = system(command)
                    command = 'cp -p ' // TRIM(adjustl(current_directory)) // '/IceBern2D.f90 ' // TRIM(adjustl(directory))
                    result = system(command)
!                     command = 'cp -p ' // TRIM(adjustl(current_directory)) // '/variables_snow.f90 ' // TRIM(adjustl(directory))
!                     result = system(command)
                else if (index(current_directory, 'home') .gt. 0) then
                    command = 'cp -p ' // TRIM(adjustl(current_directory)) // '/variables.f90 ' // TRIM(adjustl(directory))
                    result = system(command)
                    command = 'cp -p ' // TRIM(adjustl(current_directory)) // '/io.f90 ' // TRIM(adjustl(directory))
                    result = system(command)
                    command = 'cp -p ' // TRIM(adjustl(current_directory)) // '/smb_emb.f90 ' // TRIM(adjustl(directory))
                    result = system(command)
                    command = 'cp -p ' // TRIM(adjustl(current_directory)) // '/IceBern2D.f90 ' // TRIM(adjustl(directory))
                    result = system(command)
!                     command = 'cp -p ' // TRIM(adjustl(current_directory)) // '/variables_snow.f90 ' // TRIM(adjustl(directory))
!                     result = system(command)
                else
                    command = 'cp -p ' // TRIM(adjustl(current_directory)) // '/variables.f90 ' // TRIM(adjustl(directory))
                    result = system(command)
                    command = 'cp -p ' // TRIM(adjustl(current_directory)) // '/io.f90 ' // TRIM(adjustl(directory))
                    result = system(command)
                    command = 'cp -p ' // TRIM(adjustl(current_directory)) // '/smb_emb.f90 ' // TRIM(adjustl(directory))
                    result = system(command)
                    command = 'cp -p ' // TRIM(adjustl(current_directory)) // '/IceBern2D.f90 ' // TRIM(adjustl(directory))
                    result = system(command)
!                     command = 'cp -p ' // TRIM(adjustl(current_directory)) // '/variables_snow.f90 ' // TRIM(adjustl(directory))
!                     result = system(command)
                end if
               if(debug > 0 ) then
                    print *, "*** Variables copied: ", TRIM(adjustl(command))
                end if

        end if
    end subroutine init_output_directory

    subroutine copy_to_output_directory(filepath, directory)

        ! Copies the given file to the output directory. This is done via commandline (cp -p)
        implicit none
        character(len=*), intent(in) :: filepath
        character(len=*), intent(in) :: directory

        integer(kind=4)  :: result
        character(len=256) :: command
        
        command = 'cp -p ' // TRIM(adjustl(filepath)) // ' ' // TRIM(adjustl(directory))
        result = system(command)

    end subroutine copy_to_output_directory


    !--------------------------------------------------------------
    ! NetCDF Stuff
    !--------------------------------------------------------------

    subroutine init_netcdf_file(filename, ncid, nlats, nlons, lat_distance, long_distance, &
        height_varid, bed_varid, acc_varid,  diffusivity_varid, &
        discharge_x_varid, discharge_y_varid, assignent_mask_varid)!abl_varid,
        ! Init ncdf file to write into it.
        ! param filename: path and filename to the file
        ! param ncid: filehandle for the ncid file, will be overwritten
        ! param nlats: Number of latititudes (integer)
        ! param nlons: Number of longitued (integer)
        ! param lat_distance: Length in m between the points of latitude
        ! param long_distance: Length in m between the points of longitude
        !
        ! Output
        ! param height_varid: Variable ID if the height (needed to write valued to the netCDF)
        ! param bed_varid: Variable ID of the bedrock values
        ! param acc_varid: Variable ID of the accumulation
        ! param abl_varid: Variable ID of the ablation
        ! param diffusivity_varid: Variable ID for the diffusivity
        ! param discharge_x_varid: Variable ID of the discharge in x direction
        ! param discharge_y_varid: Variable ID of the discharge in y direction
        ! param assignent_mask_varid: Mask for calculations (0 = ice, 1 = water, 2 = land without ice)
        !
        ! Initialise the variables to write in 2D data with time dependency
        ! Prepare the following variables to write into the file:
        ! - Thickness (2D Real Array): varid =
        ! - Bedrock (2D Real Array)
        ! - Elevation Line Position (2D Real Array)
        ! - Mass Balance (in preparation)


        use netcdf
        implicit none
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ncid
        ! input
        integer, intent(in) :: NLATS
        integer, intent(in) :: NLONS
        integer, intent(in) :: lat_distance
        integer, intent(in) :: long_distance
        ! output
        integer, intent(out) :: height_varid        ! Variable for the netCDF Height, parameter
        integer, intent(out) :: bed_varid           ! Variable for the netCDF Bedrock, parameter
        integer, intent(out) :: acc_varid           ! Variable for the netCDF accumulation, parameter
        !integer, intent(out) :: abl_varid           ! Variable for the netCDF ablation, parameter
        integer, intent(out) :: diffusivity_varid   ! Variable for the netCDF Diffusivity, parameter
        integer, intent(out) :: discharge_x_varid   ! Variable for the netCDF discharge in x direction, parameter
        integer, intent(out) :: discharge_y_varid   ! Variable for the netCDF discharge in y direction, parameter
        integer, intent(out) :: assignent_mask_varid! Variable for the netCDF assignment_mask, parameter


        ! return value
        integer :: retval

        ! Copied from Example file!
        !--------------------------
        ! We are writing 2D data with time, We will need 3 netCDF dimensions. (Lats, Long, Time)
        integer, parameter :: NDIMS = 3

        character*(*) :: LAT_NAME, LON_NAME, REC_NAME

        parameter (LAT_NAME='y', LON_NAME='x')
        parameter (REC_NAME = 'time')
        integer lon_dimid, lat_dimid, rec_dimid

        ! The start and count arrays will tell the netCDF library where to
        ! write our data.
        !integer start(NDIMS), count(NDIMS)

        ! In addition to the latitude and longitude dimensions, we will also
        ! create latitude and longitude netCDF variables which will hold the
        ! actual latitudes and longitudes. Since they hold data about the
        ! coordinate system, the netCDF term for these is: "coordinate
        ! variables."
        real :: lats(NLATS)
        real :: lons(NLONS)
        integer lat_varid, lon_varid, rec_varid
        real START_LAT, START_LON
        parameter (START_LAT = 0, START_LON = 0)
        !
        !    ! Variables
        character*(*) HEIGHT_NAME, BED_NAME, ACC_NAME, DIFFUSIVITY_NAME!, ABL_NAME
        character*(*) DISCHARGE_X_NAME, DISCHARGE_Y_NAME, ASSIGNMENT_MASK_NAME
        parameter (HEIGHT_NAME='height')
        parameter (BED_NAME='bedrock')
        parameter (ACC_NAME='accumulation')
        !parameter (ABL_NAME='ablation')
        parameter (DIFFUSIVITY_NAME='diffusivity')
        parameter (DISCHARGE_X_NAME = 'discharge_x')
        parameter (DISCHARGE_Y_NAME = 'discharge_y')
        parameter (ASSIGNMENT_MASK_NAME = 'assignment_mask')
        integer dimids(NDIMS)
        !
        !    ! It's good practice for each variable to carry a "units" attribute.
        character*(*) :: UNITS
        parameter (UNITS = 'units')
        character*(*) HEIGHT_UNITS, BED_UNITS, ACC_UNITS, LAT_UNITS, LON_UNITS, &
             DIFFUSIVITY_UNITS, DISCHARGE_UNITS , REC_UNITS !,ABL_UNITS
	
        parameter (HEIGHT_UNITS = 'm', BED_UNITS = 'm', ACC_UNITS = 'm/yr', &
             DIFFUSIVITY_UNITS = 'm2/yr', &
            DISCHARGE_UNITS = 'm2/yr')!ABL_UNITS = 'm/yr'
        ! Units of Dimensions
        parameter (LAT_UNITS = 'm')
        parameter (LON_UNITS = 'm')
        parameter (REC_UNITS = '1000years')
	!character*(*) :: dummy_string
        !write (dummy_string, "(I4,A5)") netcdf_timesteps,"years"
	!parameter (REC_UNITS = trim(adjustl(dummy_string)))
         integer :: lat,lon !local loop variables
        ! Distance of the lat/long points
        do lat = 1, NLATS
            lats(lat) = START_LAT + (lat - 1) * lat_distance
        end do
        do lon = 1, NLONS
            lons(lon) = START_LON + (lon - 1) * long_distance
        end do

        !print *, '*** Init NetCDF file: ', TRIM(adjustl(filename))

        ! Create the file.
        retval = nf_create(filename, nf_clobber, ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !
        !
        !    ! Define dimensions
        retval = nf_def_dim(ncid, LAT_NAME, NLATS, lat_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, LON_NAME, NLONS, lon_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        ! TODO Time - this should be OK
        retval = nf_def_dim(ncid, REC_NAME, NF_UNLIMITED, rec_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Define the coordinate variables. They will hold the coordinate
        ! information, that is, the latitudes and longitudes. A varid is
        ! returned for each.
        !retval = nf_def_var(ncid, LAT_NAME, NF_REAL, 1, lat_dimid, lat_varid)
        retval = nf_def_var(ncid, LAT_NAME, NF_DOUBLE, 1, lat_dimid, lat_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, LON_NAME, NF_REAL, 1, lon_dimid, lon_varid)
        retval = nf_def_var(ncid, LON_NAME, NF_DOUBLE, 1, lon_dimid, lon_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        ! TODO Time
        !retval = nf_def_var(ncid, REC_NAME, NF_REAL, 1, rec_dimid, rec_varid)
        retval = nf_def_var(ncid, REC_NAME, NF_DOUBLE, 1, rec_dimid, rec_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Assign units attributes to coordinate var data. This attaches a
        ! text attribute to each of the coordinate variables, containing the
        ! units.
        retval = nf_put_att_text(ncid, lat_varid, UNITS, len(LAT_UNITS), LAT_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, lon_varid, UNITS, len(LON_UNITS), LON_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        ! TODO Time
        retval = nf_put_att_text(ncid, rec_varid, UNITS, len(REC_UNITS), REC_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Define the netCDF variables. The dimids array is used to pass the
        ! dimids of the dimensions of the netCDF variables.
        dimids(1) = lon_dimid
        dimids(2) = lat_dimid
        dimids(3) = rec_dimid

        ! Define the netCDF variables for the "real" values
        !retval = nf_def_var(ncid, HEIGHT_NAME, NF_REAL, NDIMS, dimids, height_varid)
        retval = nf_def_var(ncid, HEIGHT_NAME, NF_DOUBLE, NDIMS, dimids, height_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, BED_NAME, NF_REAL, NDIMS, dimids, bed_varid)
        retval = nf_def_var(ncid, BED_NAME, NF_DOUBLE, NDIMS, dimids, bed_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, ACC_NAME, NF_REAL, NDIMS, dimids, ACC_varid)
        retval = nf_def_var(ncid, ACC_NAME, NF_DOUBLE, NDIMS, dimids, acc_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, ABL_NAME, NF_REAL, NDIMS, dimids, ACC_varid)
        !retval = nf_def_var(ncid, ABL_NAME, NF_DOUBLE, NDIMS, dimids, abl_varid)
        !if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, DIFFUSIVITY_NAME, NF_REAL, NDIMS, dimids, diffusivity_varid)
        retval = nf_def_var(ncid, DIFFUSIVITY_NAME, NF_DOUBLE, NDIMS, dimids, diffusivity_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, DISCHARGE_X_NAME, NF_REAL, NDIMS, dimids, discharge_x_varid)
        retval = nf_def_var(ncid, DISCHARGE_X_NAME, NF_DOUBLE, NDIMS, dimids, discharge_x_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, DISCHARGE_Y_NAME, NF_REAL, NDIMS, dimids, discharge_y_varid)
        retval = nf_def_var(ncid, DISCHARGE_Y_NAME, NF_DOUBLE, NDIMS, dimids, discharge_y_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, ASSIGNMENT_MASK_NAME, NF_INT, NDIMS, dimids, assignent_mask_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)


        ! Assign units attributes to the pressure and temperature netCDF
        ! variables.
        retval = nf_put_att_text(ncid, height_varid, UNITS, len(HEIGHT_UNITS), HEIGHT_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, bed_varid, UNITS, len(BED_UNITS), BED_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, acc_varid, UNITS, len(ACC_UNITS), acc_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_put_att_text(ncid, abl_varid, UNITS, len(ABL_UNITS), abl_UNITS)
        !if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, diffusivity_varid, UNITS, len(DIFFUSIVITY_UNITS), DIFFUSIVITY_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, discharge_x_varid, UNITS, len(DISCHARGE_UNITS), DISCHARGE_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, discharge_y_varid, UNITS, len(DISCHARGE_UNITS), DISCHARGE_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! End define mode.
        retval = nf_enddef(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Write the coordinate variable data. This will put the latitudes
        ! and longitudes of our data grid into the netCDF file.
        retval = nf_put_var_real(ncid, lat_varid, lats)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_real(ncid, lon_varid, lons)
        if (retval .ne. nf_noerr) call handle_err(retval)
        print *,'*** Output netCDF file defined: ', TRIM(adjustl(filename))

    end subroutine init_netcdf_file

    ! special function to initialize netcdf file for climate data
    subroutine init_netcdf_climate_file(filename, ncid, nlats, nlons, lat_distance, long_distance, &
        precip_varid, temp_varid, swrad_varid, dev_precip_varid, dev_temp_varid, dev_swrad_varid)

        ! Init ncdf file to write into it.
        ! param filename: path and filename to the file
        ! param ncid: filehandle for the ncid file, will be overwritten
        ! param nlats: Number of latititudes (integer)
        ! param nlons: Number of longitued (integer)
        ! param lat_distance: Length in m between the points of latitude
        ! param long_distance: Length in m between the points of longitude
        !
        ! Output
        ! param height_varid: Variable ID if the height (needed to write valued to the netCDF)
        ! param bed_varid: Variable ID of the bedrock values
        ! param acc_varid: Variable ID of the accumulation
        ! param abl_varid: Variable ID of the ablation
        ! param diffusivity_varid: Variable ID for the diffusivity
        ! param discharge_x_varid: Variable ID of the discharge in x direction
        ! param discharge_y_varid: Variable ID of the discharge in y direction
        ! param assignent_mask_varid: Mask for calculations (0 = ice, 1 = water, 2 = land without ice)
        !
        ! Initialise the variables to write in 2D data with time dependency
        ! Prepare the following variables to write into the file:
        ! - Thickness (2D Real Array): varid =
        ! - Bedrock (2D Real Array)
        ! - Elevation Line Position (2D Real Array)
        ! - Mass Balance (in preparation)


        use netcdf
        implicit none
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ncid
        ! input
        integer, intent(in) :: NLATS
        integer, intent(in) :: NLONS
        integer, intent(in) :: lat_distance
        integer, intent(in) :: long_distance
        ! output
        integer, intent(out) :: precip_varid        ! Variable for the netCDF Height, parameter
        integer, intent(out) :: temp_varid           ! Variable for the netCDF Bedrock, parameter
        integer, intent(out) :: dev_precip_varid           ! Variable for the netCDF accumulation, parameter
        integer, intent(out) :: dev_temp_varid   ! Variable for the netCDF Diffusivity, parameter
        integer, intent(out) :: swrad_varid   ! Variable for the netCDF Diffusivity, parameter
        integer, intent(out) :: dev_swrad_varid   ! Variable for the netCDF Diffusivity, parameter

        ! return value
        integer :: retval

        ! Copied from Example file!
        !--------------------------
        ! We are writing 2D data with time, We will need 3 netCDF dimensions. (Lats, Long, Time)
        integer, parameter :: NDIMS = 3

        character*(*) :: LAT_NAME, LON_NAME, REC_NAME

        parameter (LAT_NAME='y', LON_NAME='x')
        parameter (REC_NAME = 'time')
        integer lon_dimid, lat_dimid, rec_dimid

        ! The start and count arrays will tell the netCDF library where to
        ! write our data.
        !integer start(NDIMS), count(NDIMS)

        ! In addition to the latitude and longitude dimensions, we will also
        ! create latitude and longitude netCDF variables which will hold the
        ! actual latitudes and longitudes. Since they hold data about the
        ! coordinate system, the netCDF term for these is: "coordinate
        ! variables."
        real :: lats(NLATS)
        real :: lons(NLONS)
        integer lat_varid, lon_varid, rec_varid
        real START_LAT, START_LON
        parameter (START_LAT = 0, START_LON = 0)
        !
        !    ! Variables
        character*(*) PRECIP_NAME, TEMP_NAME, DEV_PRECIP_NAME, DEV_TEMP_NAME, SWRAD_NAME, DEV_SWRAD_NAME
        parameter (PRECIP_NAME='precip')
        parameter (TEMP_NAME='temp')
        parameter (DEV_PRECIP_NAME='dev_precip')
        parameter (DEV_TEMP_NAME='dev_temp')
        parameter (SWRAD_NAME='swrad')
        parameter (DEV_SWRAD_NAME='dev_swrad')

        integer dimids(NDIMS)
        !
        !    ! It's good practice for each variable to carry a "units" attribute.
        character*(*) :: UNITS
        parameter (UNITS = 'units')
        character*(*) PRECIP_UNITS, TEMP_UNITS, DEV_PRECIP_UNITS, LAT_UNITS, LON_UNITS, &
             DEV_TEMP_UNITS, REC_UNITS, SWRAD_UNITS, DEV_SWRAD_UNITS

        parameter (PRECIP_UNITS = 'm/s', TEMP_UNITS = 'K', DEV_PRECIP_UNITS = 'm/s', &
             DEV_TEMP_UNITS = 'K', SWRAD_UNITS = 'W/m2', DEV_SWRAD_UNITS = 'W/m2' )
        ! Units of Dimensions
        parameter (LAT_UNITS = 'm')
        parameter (LON_UNITS = 'm')
        parameter (REC_UNITS = '1000years')
	!character*(*) :: dummy_string
        integer :: lat,lon !local loop variables
        !write (dummy_string, "(I4,A5)") netcdf_timesteps,"years"
	!parameter (REC_UNITS = trim(adjustl(dummy_string)))

        ! Distance of the lat/long points
        do lat = 1, NLATS
            lats(lat) = START_LAT + (lat - 1) * lat_distance
        end do
        do lon = 1, NLONS
            lons(lon) = START_LON + (lon - 1) * long_distance
        end do

        !print *, '*** Init NetCDF file: ', TRIM(adjustl(filename))

        ! Create the file.
        retval = nf_create(filename, nf_clobber, ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !
        !
        !    ! Define dimensions
        retval = nf_def_dim(ncid, LAT_NAME, NLATS, lat_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, LON_NAME, NLONS, lon_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        ! TODO Time - this should be OK
        retval = nf_def_dim(ncid, REC_NAME, NF_UNLIMITED, rec_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Define the coordinate variables. They will hold the coordinate
        ! information, that is, the latitudes and longitudes. A varid is
        ! returned for each.
        !retval = nf_def_var(ncid, LAT_NAME, NF_REAL, 1, lat_dimid, lat_varid)
        retval = nf_def_var(ncid, LAT_NAME, NF_DOUBLE, 1, lat_dimid, lat_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, LON_NAME, NF_REAL, 1, lon_dimid, lon_varid)
        retval = nf_def_var(ncid, LON_NAME, NF_DOUBLE, 1, lon_dimid, lon_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        ! TODO Time
        !retval = nf_def_var(ncid, REC_NAME, NF_REAL, 1, rec_dimid, rec_varid)
        retval = nf_def_var(ncid, REC_NAME, NF_DOUBLE, 1, rec_dimid, rec_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Assign units attributes to coordinate var data. This attaches a
        ! text attribute to each of the coordinate variables, containing the
        ! units.
        retval = nf_put_att_text(ncid, lat_varid, UNITS, len(LAT_UNITS), LAT_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, lon_varid, UNITS, len(LON_UNITS), LON_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        ! TODO Time
        retval = nf_put_att_text(ncid, rec_varid, UNITS, len(REC_UNITS), REC_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Define the netCDF variables. The dimids array is used to pass the
        ! dimids of the dimensions of the netCDF variables.
        dimids(1) = lon_dimid
        dimids(2) = lat_dimid
        dimids(3) = rec_dimid

        ! Define the netCDF variables for the "real" values
        !retval = nf_def_var(ncid, PRECIP_NAME, NF_REAL, NDIMS, dimids, precip_varid)
        retval = nf_def_var(ncid, PRECIP_NAME, NF_DOUBLE, NDIMS, dimids, precip_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, TEMP_NAME, NF_REAL, NDIMS, dimids, temp_varid)
        retval = nf_def_var(ncid, TEMP_NAME, NF_DOUBLE, NDIMS, dimids, temp_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, DEV_PRECIP_NAME, NF_REAL, NDIMS, dimids, dev_precip_varid)
        retval = nf_def_var(ncid, DEV_PRECIP_NAME, NF_DOUBLE, NDIMS, dimids, dev_precip_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, DEV_TEMP_NAME, NF_REAL, NDIMS, dimids, dev_temp_varid)
        retval = nf_def_var(ncid, DEV_TEMP_NAME, NF_DOUBLE, NDIMS, dimids, dev_temp_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, SWRAD_NAME, NF_REAL, NDIMS, dimids, swrad_varid)
        retval = nf_def_var(ncid, SWRAD_NAME, NF_DOUBLE, NDIMS, dimids, swrad_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, DEV_SWRAD_NAME, NF_REAL, NDIMS, dimids, dev_swrad_varid)
        retval = nf_def_var(ncid, DEV_SWRAD_NAME, NF_DOUBLE, NDIMS, dimids, dev_swrad_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)



        ! Assign units attributes to the pressure and temperature netCDF
        ! variables.
        retval = nf_put_att_text(ncid, precip_varid, UNITS, len(PRECIP_UNITS), PRECIP_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, temp_varid, UNITS, len(TEMP_UNITS), TEMP_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, dev_precip_varid, UNITS, len(DEV_PRECIP_UNITS), DEV_PRECIP_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, dev_temp_varid, UNITS, len(DEV_TEMP_UNITS), DEV_TEMP_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
	retval = nf_put_att_text(ncid, swrad_varid, UNITS, len(SWRAD_UNITS), SWRAD_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
	retval = nf_put_att_text(ncid, dev_swrad_varid, UNITS, len(DEV_SWRAD_UNITS), DEV_SWRAD_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! End define mode.
        retval = nf_enddef(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Write the coordinate variable data. This will put the latitudes
        ! and longitudes of our data grid into the netCDF file.
        retval = nf_put_var_real(ncid, lat_varid, lats)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_real(ncid, lon_varid, lons)
        if (retval .ne. nf_noerr) call handle_err(retval)
        print *,'*** Output netCDF file defined: ', TRIM(adjustl(filename))

    end subroutine init_netcdf_climate_file
    ! subroutine to open a netcdf file for accumulation data


    subroutine init_netcdf_accum_file(filename, ncid, nlats, nlons, lat_distance, long_distance, &
        accum_varid, rain_varid, melt_varid, refreezing_varid, snow_varid, runoff_varid )

        ! Init ncdf file to write into it.
        ! param filename: path and filename to the file
        ! param ncid: filehandle for the ncid file, will be overwritten
        ! param nlats: Number of latititudes (integer)
        ! param nlons: Number of longitued (integer)
        ! param lat_distance: Length in m between the points of latitude
        ! param long_distance: Length in m between the points of longitude
        !
        ! Output
        ! param height_varid: Variable ID if the height (needed to write valued to the netCDF)
        ! param bed_varid: Variable ID of the bedrock values
        ! param acc_varid: Variable ID of the accumulation
        ! param abl_varid: Variable ID of the ablation
        ! param diffusivity_varid: Variable ID for the diffusivity
        ! param discharge_x_varid: Variable ID of the discharge in x direction
        ! param discharge_y_varid: Variable ID of the discharge in y direction
        ! param assignent_mask_varid: Mask for calculations (0 = ice, 1 = water, 2 = land without ice)
        !
        ! Initialise the variables to write in 2D data with time dependency
        ! Prepare the following variables to write into the file:
        ! - Thickness (2D Real Array): varid =
        ! - Bedrock (2D Real Array)
        ! - Elevation Line Position (2D Real Array)
        ! - Mass Balance (in preparation)


        use netcdf
        implicit none
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ncid
        ! input
        integer, intent(in) :: NLATS
        integer, intent(in) :: NLONS
        integer, intent(in) :: lat_distance
        integer, intent(in) :: long_distance
        ! output
        integer, intent(out) :: accum_varid        ! Variable for the netCDF Height, parameter
        integer, intent(out) :: rain_varid           ! Variable for the netCDF Bedrock, parameter
        integer, intent(out) :: melt_varid           ! Variable for the netCDF accumulation, parameter
        integer, intent(out) :: refreezing_varid   ! Variable for the netCDF Diffusivity, parameter
        integer, intent(out) :: snow_varid   ! Variable for the netCDF Diffusivity, parameter
        integer, intent(out) :: runoff_varid   ! Variable for the netCDF Diffusivity, parameter

        ! return value
        integer :: retval

        ! Copied from Example file!
        !--------------------------
        ! We are writing 2D data with time, We will need 3 netCDF dimensions. (Lats, Long, Time)
        integer, parameter :: NDIMS = 3

        character*(*) :: LAT_NAME, LON_NAME, REC_NAME

        parameter (LAT_NAME='y', LON_NAME='x')
        parameter (REC_NAME = 'time')
        integer lon_dimid, lat_dimid, rec_dimid

        ! The start and count arrays will tell the netCDF library where to
        ! write our data.
        !integer start(NDIMS), count(NDIMS)

        ! In addition to the latitude and longitude dimensions, we will also
        ! create latitude and longitude netCDF variables which will hold the
        ! actual latitudes and longitudes. Since they hold data about the
        ! coordinate system, the netCDF term for these is: "coordinate
        ! variables."
        real :: lats(NLATS)
        real :: lons(NLONS)
        integer lat_varid, lon_varid, rec_varid
        real START_LAT, START_LON
        parameter (START_LAT = 0, START_LON = 0)
        !
        !    ! Variables
        character*(*) accum_NAME, rain_NAME, melt_NAME, refreezing_NAME, snow_NAME, runoff_NAME
        parameter (accum_NAME='accum')
        parameter (rain_NAME='rain')
        parameter (melt_NAME='melt')
        parameter (refreezing_NAME='refreezing')
        parameter (snow_NAME='snow')
        parameter (runoff_NAME='runoff')


        integer dimids(NDIMS)
        !
        !    ! It's good practice for each variable to carry a "units" attribute.
        character*(*) :: UNITS
        parameter (UNITS = 'units')
        character*(*) accum_UNITS, rain_UNITS, melt_UNITS, LAT_UNITS, LON_UNITS, &
             refreezing_UNITS, REC_UNITS, snow_UNITS, runoff_UNITS

        parameter (accum_UNITS = 'kg/m2/a', rain_UNITS = 'kg/m2/a', melt_UNITS = 'kg/m2/a', &
             refreezing_UNITS = 'kg/m2/a', snow_UNITS = 'kg/m2/a', runoff_UNITS = 'kg/m2/a' )
        ! Units of Dimensions
        parameter (LAT_UNITS = 'm')
        parameter (LON_UNITS = 'm')
        parameter (REC_UNITS = '1000years')
	!character*(*) :: dummy_string
        integer :: lat, lon !local loop variables
        !write (dummy_string, "(I4,A5)") netcdf_timesteps,"years"
	!parameter (REC_UNITS = trim(adjustl(dummy_string)))

        ! Distance of the lat/long points
        do lat = 1, NLATS
            lats(lat) = START_LAT + (lat - 1) * lat_distance
        end do
        do lon = 1, NLONS
            lons(lon) = START_LON + (lon - 1) * long_distance
        end do

        !print *, '*** Init NetCDF file: ', TRIM(adjustl(filename))

        ! Create the file.
        retval = nf_create(filename, nf_clobber, ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !
        !
        !    ! Define dimensions
        retval = nf_def_dim(ncid, LAT_NAME, NLATS, lat_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, LON_NAME, NLONS, lon_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        ! TODO Time - this should be OK
        retval = nf_def_dim(ncid, REC_NAME, NF_UNLIMITED, rec_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Define the coordinate variables. They will hold the coordinate
        ! information, that is, the latitudes and longitudes. A varid is
        ! returned for each.
        !retval = nf_def_var(ncid, LAT_NAME, NF_REAL, 1, lat_dimid, lat_varid)
        retval = nf_def_var(ncid, LAT_NAME, NF_DOUBLE, 1, lat_dimid, lat_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, LON_NAME, NF_REAL, 1, lon_dimid, lon_varid)
        retval = nf_def_var(ncid, LON_NAME, NF_DOUBLE, 1, lon_dimid, lon_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        ! TODO Time
        !retval = nf_def_var(ncid, REC_NAME, NF_REAL, 1, rec_dimid, rec_varid)
        retval = nf_def_var(ncid, REC_NAME, NF_DOUBLE, 1, rec_dimid, rec_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Assign units attributes to coordinate var data. This attaches a
        ! text attribute to each of the coordinate variables, containing the
        ! units.
        retval = nf_put_att_text(ncid, lat_varid, UNITS, len(LAT_UNITS), LAT_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, lon_varid, UNITS, len(LON_UNITS), LON_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        ! TODO Time
        retval = nf_put_att_text(ncid, rec_varid, UNITS, len(REC_UNITS), REC_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Define the netCDF variables. The dimids array is used to pass the
        ! dimids of the dimensions of the netCDF variables.
        dimids(1) = lon_dimid
        dimids(2) = lat_dimid
        dimids(3) = rec_dimid

        ! Define the netCDF variables for the "real" values
        !retval = nf_def_var(ncid, accum_NAME, NF_REAL, NDIMS, dimids, accum_varid)
        retval = nf_def_var(ncid, accum_NAME, NF_REAL, NDIMS, dimids, accum_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, rain_NAME, NF_REAL, NDIMS, dimids, rain_varid)
        retval = nf_def_var(ncid, rain_NAME, NF_REAL, NDIMS, dimids, rain_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, melt_NAME, NF_REAL, NDIMS, dimids, melt_varid)
        retval = nf_def_var(ncid, melt_NAME, NF_REAL, NDIMS, dimids, melt_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, refreezing_NAME, NF_REAL, NDIMS, dimids, refreezing_varid)
        retval = nf_def_var(ncid, refreezing_NAME, NF_REAL, NDIMS, dimids, refreezing_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        !retval = nf_def_var(ncid, snow_NAME, NF_REAL, NDIMS, dimids, snow_varid)
        retval = nf_def_var(ncid, snow_NAME, NF_REAL, NDIMS, dimids, snow_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, runoff_NAME, NF_REAL, NDIMS, dimids, runoff_varid)
        retval = nf_def_var(ncid, runoff_NAME, NF_REAL, NDIMS, dimids, runoff_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)


        ! Assign units attributes to the pressure and rainerature netCDF
        ! variables.
        retval = nf_put_att_text(ncid, accum_varid, UNITS, len(accum_UNITS), accum_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, rain_varid, UNITS, len(rain_UNITS), rain_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, melt_varid, UNITS, len(melt_UNITS), melt_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, refreezing_varid, UNITS, len(refreezing_UNITS), refreezing_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, snow_varid, UNITS, len(snow_UNITS), snow_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, runoff_varid, UNITS, len(runoff_UNITS), runoff_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! End define mode.
        retval = nf_enddef(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Write the coordinate variable data. This will put the latitudes
        ! and longitudes of our data grid into the netCDF file.
        retval = nf_put_var_real(ncid, lat_varid, lats)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_real(ncid, lon_varid, lons)
        if (retval .ne. nf_noerr) call handle_err(retval)
        print *,'*** Output netCDF file defined: ', TRIM(adjustl(filename))

    end subroutine init_netcdf_accum_file



    ! initialize netcdf file for snow mask
    subroutine init_netcdf_snow_mask(filename, ncid, nlats, nlons, ndays, lat_distance, long_distance, snow_mask_varid )

        ! Init ncdf file to write into it.
        ! param filename: path and filename to the file
        ! param ncid: filehandle for the ncid file, will be overwritten
        ! param nlats: Number of latititudes (integer)
        ! param nlons: Number of longitued (integer)
        ! param lat_distance: Length in m between the points of latitude
        ! param long_distance: Length in m between the points of longitude
        !
        ! Output
        ! param snow_mask_varid: Mask for calculations (1 = water, 2 = land without snow, 3 = snow, 4 = wet snow)


        use netcdf
        implicit none
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ncid
        ! input
        integer, intent(in) :: NLATS
        integer, intent(in) :: NLONS
        integer, intent(in) :: ndays
        integer, intent(in) :: lat_distance
        integer, intent(in) :: long_distance
        ! output
        integer, intent(out) :: snow_mask_varid! Variable for the netCDF snow_mask, parameter


        ! return value
        integer :: retval

        ! Copied from Example file!
        !--------------------------
        ! We are writing 2D data with time, We will need 3 netCDF dimensions. (Lats, Long, Time)
        integer, parameter :: NDIMS = 3

        character*(*) :: LAT_NAME, LON_NAME, REC_NAME

        parameter (LAT_NAME='y', LON_NAME='x')
        parameter (REC_NAME = 'time')
        integer lon_dimid, lat_dimid, rec_dimid

        ! The start and count arrays will tell the netCDF library where to
        ! write our data.
        !integer start(NDIMS), count(NDIMS)

        ! In addition to the latitude and longitude dimensions, we will also
        ! create latitude and longitude netCDF variables which will hold the
        ! actual latitudes and longitudes. Since they hold data about the
        ! coordinate system, the netCDF term for these is: "coordinate
        ! variables."
        real :: lats(NLATS)
        real :: lons(NLONS)
        integer lat_varid, lon_varid, rec_varid
        real START_LAT, START_LON
        parameter (START_LAT = 0, START_LON = 0)
        !
        !    ! Variables
        character*(*) SNOW_MASK_NAME
        parameter (SNOW_MASK_NAME = 'snow_mask')
        integer dimids(NDIMS)
        !
        !    ! It's good practice for each variable to carry a "units" attribute.
        character*(*) :: UNITS
        parameter (UNITS = 'units')
        character*(*) HEIGHT_UNITS, BED_UNITS, ACC_UNITS, LAT_UNITS, LON_UNITS, &
             DIFFUSIVITY_UNITS, DISCHARGE_UNITS , REC_UNITS !,ABL_UNITS
	
        parameter (HEIGHT_UNITS = 'm', BED_UNITS = 'm', ACC_UNITS = 'm/yr', &
             DIFFUSIVITY_UNITS = 'm2/yr', &
            DISCHARGE_UNITS = 'm2/yr')!ABL_UNITS = 'm/yr'
        ! Units of Dimensions
        parameter (LAT_UNITS = 'm')
        parameter (LON_UNITS = 'm')
        parameter (REC_UNITS = 'days')
        
        integer :: lat, lon !local running variables
	!character*(*) :: dummy_string
        !write (dummy_string, "(I4,A5)") netcdf_timesteps,"years"
	!parameter (REC_UNITS = trim(adjustl(dummy_string)))

        ! Distance of the lat/long points
        do lat = 1, NLATS
            lats(lat) = START_LAT + (lat - 1) * lat_distance
        end do
        do lon = 1, NLONS
            lons(lon) = START_LON + (lon - 1) * long_distance
        end do

        !print *, '*** Init NetCDF file: ', TRIM(adjustl(filename))

        ! Create the file.
        retval = nf_create(filename, nf_clobber, ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !
        !
        !    ! Define dimensions
        retval = nf_def_dim(ncid, LAT_NAME, NLATS, lat_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, LON_NAME, NLONS, lon_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        ! TODO Time - this should be OK
        retval = nf_def_dim(ncid, REC_NAME, NF_UNLIMITED, rec_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Define the coordinate variables. They will hold the coordinate
        ! information, that is, the latitudes and longitudes. A varid is
        ! returned for each.
        !retval = nf_def_var(ncid, LAT_NAME, NF_REAL, 1, lat_dimid, lat_varid)
        retval = nf_def_var(ncid, LAT_NAME, NF_DOUBLE, 1, lat_dimid, lat_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, LON_NAME, NF_REAL, 1, lon_dimid, lon_varid)
        retval = nf_def_var(ncid, LON_NAME, NF_DOUBLE, 1, lon_dimid, lon_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        ! TODO Time
        !retval = nf_def_var(ncid, REC_NAME, NF_REAL, 1, rec_dimid, rec_varid)
        retval = nf_def_var(ncid, REC_NAME, NF_DOUBLE, 1, rec_dimid, rec_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Assign units attributes to coordinate var data. This attaches a
        ! text attribute to each of the coordinate variables, containing the
        ! units.
        retval = nf_put_att_text(ncid, lat_varid, UNITS, len(LAT_UNITS), LAT_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, lon_varid, UNITS, len(LON_UNITS), LON_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        ! TODO Time
        retval = nf_put_att_text(ncid, rec_varid, UNITS, len(REC_UNITS), REC_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Define the netCDF variables. The dimids array is used to pass the
        ! dimids of the dimensions of the netCDF variables.
        dimids(1) = lon_dimid
        dimids(2) = lat_dimid
        dimids(3) = rec_dimid

        ! Define the netCDF variables for the "real" values
        retval = nf_def_var(ncid, SNOW_MASK_NAME, NF_INT, NDIMS, dimids, snow_mask_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! End define mode.
        retval = nf_enddef(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Write the coordinate variable data. This will put the latitudes
        ! and longitudes of our data grid into the netCDF file.
        retval = nf_put_var_real(ncid, lat_varid, lats)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_real(ncid, lon_varid, lons)
        if (retval .ne. nf_noerr) call handle_err(retval)
        print *,'*** Output netCDF file defined: ', TRIM(adjustl(filename))

    end subroutine init_netcdf_snow_mask




    ! subroutine to close an open netcdf file
    subroutine closeNCDFFile(ncid)
        use netcdf
        implicit none
        integer, intent(in) :: ncid

        integer :: retval
        
        ! Close the file.
        retval = nf_close(ncid)
        if (retval .ne. nf_noerr) then
            print *, 'Could not close NetCDF'
            call handle_err(retval)
        end if

        ! If we got this far, everything worked as expected. Yipee!
        print *,'*** netCDF file written!'

    end subroutine closeNCDFFile



    subroutine writeNCDFGridValues(ncid, year, varid, values, nlats, nlons)
        ! Writes the Values (Real(kind=8)) into the netCDF File
        ! Param ncid: File handle of the netCDF File
        ! Param year: The year of the symmulation (integer)
        ! Param varid: Variable ID (a, bedrock, ...) the values belong to
        ! Param values: Array (lon, lat) with values as REAL(kind=8)
        ! param nlats: Number of latititudes (integer)
        ! param nlons: Number of longitued (integer)
        !
        ! Flushs the data to the disk afterwards
        use netcdf
        implicit none
        ! Parameters
        integer, intent(in) :: ncid
        integer, intent(in) :: year
        integer, intent(in) :: varid
        integer, intent(in) :: nlats
        integer, intent(in) :: nlons
        real(kind=4), dimension(nlons,nlats), intent(in) :: values

        ! Error handling.
        integer :: retval

        ! We are writing 2D data with time, We will need 3 netCDF dimensions. (Lats, Long, Time)
        integer, parameter :: NDIMS = 3

        ! The start and count arrays will tell the netCDF library where to
        ! write our data.
        integer start(NDIMS), count(NDIMS)


        ! These settings tell netcdf to write one timestep of data. (The
        ! setting of start(3) inside the loop below tells netCDF which
        ! timestep to write.)
        count(1) = NLONS
        count(2) = NLATS
        count(3) = 1
        start(1) = 1
        start(2) = 1
        start(3) = year

        ! Write the pretend data. This will write the data.
        ! The arrays only hold one timestep worth of data.
        !retval = nf_put_vara_real(ncid, varid, start, count, values)
        retval = nf_put_vara_double(ncid, varid, start, count, values)
        if (retval .ne. nf_noerr) then
            print *, 'Could not write 2D variable to NetCDF: ', varid
            call handle_err(retval)
        end if

        ! Flush data to the disk
        retval = NF_SYNC(ncid)
        if (retval .ne. nf_noerr) then
            print *, 'Could not sync NetCDF'
            call handle_err(retval)
        end if

    end subroutine writeNCDFGridValues


    subroutine writeNCDFGridIntegerValues(ncid, year, varid, values, nlats, nlons)
        ! Writes the Values (Integer) into the netCDF File
        ! Param ncid: File handle of the netCDF File
        ! Param year: The year of the symmulation (integer)
        ! Param varid: Variable ID (a, bedrock, ...) the values belong to
        ! Param values: Array (lon, lat) with values as Integer
        ! param nlats: Number of latititudes (integer)
        ! param nlons: Number of longitued (integer)
        !
        ! Flushs the data to the disk afterwards
        use netcdf
        implicit none
        ! Parameters
        integer, intent(in) :: ncid
        integer, intent(in) :: year
        integer, intent(in) :: varid
        integer, intent(in) :: nlats
        integer, intent(in) :: nlons
        integer, dimension(nlons,nlats), intent(in) :: values

        ! Error handling.
        integer :: retval

        ! We are writing 2D data with time, We will need 3 netCDF dimensions. (Lats, Long, Time)
        integer, parameter :: NDIMS = 3

        ! The start and count arrays will tell the netCDF library where to
        ! write our data.
        integer start(NDIMS), count(NDIMS)


        ! These settings tell netcdf to write one timestep of data. (The
        ! setting of start(3) inside the loop below tells netCDF which
        ! timestep to write.)
        count(1) = NLONS
        count(2) = NLATS
        count(3) = 1
        start(1) = 1
        start(2) = 1
        start(3) = year

        ! Write the pretend data. This will write the data.
        ! The arrays only hold one timestep worth of data.
        !retval = nf_put_vara_real(ncid, varid, start, count, values)
        retval = nf_put_vara_int(ncid, varid, start, count, values)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Flush data to the disk
        retval = NF_SYNC(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

    end subroutine writeNCDFGridIntegerValues

    !--------------------------------------------------------------
    ! NetCDF Stuff
    !--------------------------------------------------------------
    subroutine handle_err(errcode)
        implicit none
        integer errcode
        
        print *, 'Error: ', nf_strerror(errcode)
        stop 4
    end subroutine handle_err

    !33333333333333333333333333333333333333333333333333333333333333333333333333333333333333333
    ! ----------------------- 3D NetCDF  -----------------------------------------------------
    !33333333333333333333333333333333333333333333333333333333333333333333333333333333333333333

    subroutine init_netcdf_3D_file(filename, ncid, nlats, nlons, N_LAYER, lat_distance, long_distance, &
        z_distance, snowman_varid, lwmass_varid, rho_snow_varid, snow_temp_varid, albedo_dynamic_varid)
        ! Init ncdf file to write into it.
        ! param filename: path and filename to the file
        ! param ncid: filehandle for the ncid file, will be overwritten
        ! param nlats: Number of latititudes (integer)
        ! param nlons: Number of longitued (integer)
        ! param N_LAYER: Number of layers (integer)
        ! param lat_distance: Length in m between the points of latitude
        ! param long_distance: Length in m between the points of longitude
        ! param z_distance: Length in m. dummy does not contain physical values
        !
        ! Output
        ! param height_varid: Variable ID if the height (needed to write valued to the netCDF)
        ! param snowman_varid:      Variable ID of the snow mass values
        ! param lwmass_varid:       Variable ID of the lw mass values
        ! param rho_snow_varid:     Variable ID of the density of snow values
        ! param snow_temp_varid:    Variable ID of the snow temperature values

        !
        ! Initialise the variables to write in 3D data with time dependency
        ! Prepare the following variables to write into the file:
        ! - "quantity" (3D Real Array): varid =
        

        
        use netcdf
        implicit none
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ncid
        ! input
        integer, intent(in) :: NLATS
        integer, intent(in) :: NLONS
        integer, intent(in) :: N_LAYER
        integer, intent(in) :: lat_distance
        integer, intent(in) :: long_distance
        integer, intent(in) :: z_distance

        ! output
        integer, intent(out) :: snowman_varid        ! Variable for the netCDF snow mass, parameter
        integer, intent(out) :: lwmass_varid           ! Variable for the netCDF liquid water, parameter
        integer, intent(out) :: rho_snow_varid           ! Variable for the netCDF density, parameter
        integer, intent(out) :: snow_temp_varid           ! Variable for the netCDF temperature, parameter
        integer, intent(out) :: albedo_dynamic_varid

        ! return value
        integer :: retval

        ! Copied from Example file!
        !--------------------------
        ! We are writing 3D data with time, We will need 4 netCDF dimensions. (Lats, Long, z, Time)
        integer, parameter :: NDIMS = 4

        character*(*) :: LAT_NAME, LON_NAME, REC_NAME, Z_NAME, ANNUAL_NAME
        parameter (LAT_NAME='y', LON_NAME='x',Z_NAME='layer')
        parameter (REC_NAME = 'time')
        parameter (ANNUAL_NAME = 'annual')
        integer lon_dimid, lat_dimid, z_dimid, rec_dimid, annual_dimid


        ! The start and count arrays will tell the netCDF library where to
        ! write our data.
        !integer start(NDIMS), count(NDIMS)

        ! In addition to the latitude and longitude dimensions, we will also
        ! create latitude and longitude netCDF variables which will hold the
        ! actual latitudes and longitudes. Since they hold data about the
        ! coordinate system, the netCDF term for these is: "coordinate
        ! variables."
        real :: lats(NLATS)
        real :: lons(NLONS)
        real :: zs(N_LAYER)
        integer lat_varid, lon_varid, rec_varid, z_varid, annual_varid
        real START_LAT, START_LON
        parameter (START_LAT = 0, START_LON = 0)
        !
        !    ! Variables
        character*(*) snowman_NAME, lwmass_NAME, rho_snow_NAME, snow_temp_NAME, albedo_dynamic_NAME
        parameter (snowman_NAME='snowmass')
        parameter (lwmass_NAME='lwmass')
        parameter (rho_snow_NAME='snowdensity')
        parameter (snow_temp_NAME='snowtemp')
        parameter (albedo_dynamic_NAME='albedo')

        integer dimids(NDIMS)
        integer dimids2(NDIMS-1)
        integer dimids_annual(NDIMS)
        !
        !    ! It's good practice for each variable to carry a "units" attribute.
        character*(*) :: UNITS
        parameter (UNITS = 'units')
        character*(*) snowman_UNITS, lwmass_UNITS, rho_snow_UNITS, snow_temp_UNITS, albedo_dynamic_UNITS, LAT_UNITS, LON_UNITS, &
            REC_UNITS, z_UNITS, ANNUAL_UNITS
        parameter ( snowman_UNITS = 'kg/m2', lwmass_UNITS = 'kg/m2', rho_snow_UNITS = 'kg/m3', snow_temp_UNITS = 'K', &
            albedo_dynamic_UNITS = '')
        ! Units of Dimensions
        parameter (LAT_UNITS = 'm')
        parameter (LON_UNITS = 'm')
        parameter (z_UNITS = 'box')
        parameter (ANNUAL_UNITS = 'y')
        parameter (REC_UNITS = '1000years')
        
        integer :: lat,lon,zz !local run variable

        ! Distance of the lat/long points
        do lat = 1, NLATS
            lats(lat) = START_LAT + (lat - 1) * lat_distance
        end do
        do lon = 1, NLONS
            lons(lon) = START_LON + (lon - 1) * long_distance
        end do
        do zz = 1, N_LAYER
            zs(zz) = 1 + (real(zz) - 1) * z_distance
        end do
        !print *, '*** Init NetCDF file: ', TRIM(adjustl(filename))

        ! Create the file.
        retval = nf_create(filename, nf_clobber, ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !
        !
        !    ! Define dimensions
        retval = nf_def_dim(ncid, LAT_NAME, NLATS, lat_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, LON_NAME, NLONS, lon_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, z_NAME, N_LAYER, z_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        retval = nf_def_dim(ncid, ANNUAL_NAME,1, ANNUAL_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        
        ! TODO Time - this should be OK
        retval = nf_def_dim(ncid, REC_NAME, NF_UNLIMITED, rec_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Define the coordinate variables. They will hold the coordinate
        ! information, that is, the latitudes and longitudes. A varid is
        ! returned for each.
        !retval = nf_def_var(ncid, LAT_NAME, NF_REAL, 1, lat_dimid, lat_varid)
        retval = nf_def_var(ncid, LAT_NAME, NF_DOUBLE, 1, lat_dimid, lat_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, LON_NAME, NF_REAL, 1, lon_dimid, lon_varid)
        retval = nf_def_var(ncid, LON_NAME, NF_DOUBLE, 1, lon_dimid, lon_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
         !retval = nf_def_var(ncid, z_NAME, NF_REAL, 1, z_dimid, z_varid)
        retval = nf_def_var(ncid, Z_NAME, NF_DOUBLE, 1, z_dimid, z_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_def_var(ncid, ANNUAL_NAME, NF_DOUBLE, 1, annual_dimid, annual_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! TODO Time
        !retval = nf_def_var(ncid, REC_NAME, NF_REAL, 1, rec_dimid, rec_varid)
        retval = nf_def_var(ncid, REC_NAME, NF_DOUBLE, 1, rec_dimid, rec_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Assign units attributes to coordinate var data. This attaches a
        ! text attribute to each of the coordinate variables, containing the
        ! units.
        retval = nf_put_att_text(ncid, lat_varid, UNITS, len(LAT_UNITS), LAT_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, lon_varid, UNITS, len(LON_UNITS), LON_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, z_varid, UNITS, len(Z_UNITS), Z_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_put_att_text(ncid, annual_varid, UNITS, len(ANNUAL_UNITS), ANNUAL_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)


        ! TODO Time
        retval = nf_put_att_text(ncid, rec_varid, UNITS, len(REC_UNITS), REC_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Define the netCDF variables. The dimids array is used to pass the
        ! dimids of the dimensions of the netCDF variables.
        dimids(1) = lon_dimid
        dimids(2) = lat_dimid
        dimids(3) = z_dimid
        dimids(4) = rec_dimid
        dimids2(1) = lon_dimid
        dimids2(2) = lat_dimid

        dimids2(3) = rec_dimid
        
        dimids_annual(1)=lon_dimid
        dimids_annual(2)=lat_dimid
        dimids_annual(3)=z_dimid
        dimids_annual(4)=annual_dimid

        
        ! Define the netCDF variables for the "real" values
        retval = nf_def_var(ncid, snowman_NAME, NF_REAL, NDIMS, dimids,snowman_varid )
        !retval = nf_def_var(ncid, snowman_NAME, NF_DOUBLE, NDIMS, dimids, snowman_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, lwmass_NAME, NF_REAL, NDIMS, dimids, lwmass_varid)
        !retval = nf_def_var(ncid, lwmass_NAME, NF_DOUBLE, NDIMS, dimids, lwmass_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, rho_snow_NAME, NF_REAL, NDIMS, dimids, rho_snow_varid)
        !retval = nf_def_var(ncid, rho_snow_NAME, NF_DOUBLE, NDIMS, dimids, rho_snow_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, snow_temp_NAME, NF_REAL, NDIMS, dimids, snow_temp_varid)
        !retval = nf_def_var(ncid, snow_temp_NAME, NF_DOUBLE, NDIMS, dimids, snow_temp_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, albedo_dynamic_NAME, NF_REAL, NDIMS-1, dimids2, albedo_dynamic_varid)
        !retval = nf_def_var(ncid, snow_temp_NAME, NF_DOUBLE, NDIMS, dimids, snow_temp_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)




        ! Assign units attributes to the pressure and temperature netCDF
        ! variables.
        retval = nf_put_att_text(ncid, snowman_varid, UNITS, len(snowman_UNITS), snowman_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, lwmass_varid, UNITS, len(lwmass_UNITS), lwmass_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, rho_snow_varid, UNITS, len(rho_snow_UNITS), rho_snow_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, snow_temp_varid, UNITS, len(snow_temp_UNITS), snow_temp_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, albedo_dynamic_varid, UNITS, len(albedo_dynamic_UNITS), albedo_dynamic_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)


        ! End define mode.
        retval = nf_enddef(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Write the coordinate variable data. This will put the latitudes
        ! and longitudes of our data grid into the netCDF file.
        retval = nf_put_var_real(ncid, lat_varid, lats)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_real(ncid, lon_varid, lons)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_real(ncid, z_varid, zs)
        if (retval .ne. nf_noerr) call handle_err(retval)
        print *,'*** Output netCDF file defined: ', TRIM(adjustl(filename))

    end subroutine init_netcdf_3D_file

    !=====================================================================================================
    
    subroutine init_netcdf_snow_file(filename, ncid, nlats, nlons, N_LAYER, n_time, lat_distance, long_distance, &
        z_distance, snowman_annual_varid, lwmass_annual_varid, rho_snow_annual_varid, snow_temp_annual_varid, &
        albedo_dynamic_annual_varid, latent_heat_annual_varid, accum_annual_varid, rain_annual_varid, &
        melt_annual_varid, refreezing_annual_varid, snow_annual_varid, runoff_annual_varid, &
        melt_ice_annual_varid, snow_mask_annual_varid, snow_temp_ave_surface_annual_varid, regridding_annual_varid,&
        real_mass_balance_annual_varid,new_input_file3)
        ! Init ncdf file to write into it.
        ! param filename: path and filename to the file
        ! param ncid: filehandle for the ncid file, will be overwritten
        ! param nlats: Number of latititudes (integer)
        ! param nlons: Number of longitued (integer)
        ! param N_LAYER: Number of layers (integer)
        ! param lat_distance: Length in m between the points of latitude
        ! param long_distance: Length in m between the points of longitude
        ! param z_distance: Length in m. dummy does not contain physical values
        !
        ! Output
        ! param height_varid: Variable ID if the height (needed to write valued to the netCDF)
        ! param snowman_varid:      Variable ID of the snow mass values
        ! param lwmass_varid:       Variable ID of the lw mass values
        ! param rho_snow_varid:     Variable ID of the density of snow values
        ! param snow_temp_varid:    Variable ID of the snow temperature values

        !
        ! Initialise the variables to write in 3D data with time dependency
        ! Prepare the following variables to write into the file:
        ! - "quantity" (3D Real Array): varid =


        use variables
        use netcdf
        implicit none        
        
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ncid
        ! input
        integer, intent(in) :: NLATS
        integer, intent(in) :: NLONS
        integer, intent(in) :: N_LAYER
        integer, intent(in) :: n_time
        integer, intent(in) :: lat_distance
        integer, intent(in) :: long_distance
        integer, intent(in) :: z_distance
        character(128), intent(in) :: new_input_file3
        ! output
        integer, intent(out) :: snowman_annual_varid        ! Variable for the netCDF snow mass, parameter
        integer, intent(out) :: lwmass_annual_varid           ! Variable for the netCDF liquid water, parameter
        integer, intent(out) :: rho_snow_annual_varid           ! Variable for the netCDF density, parameter
        integer, intent(out) :: snow_temp_annual_varid           ! Variable for the netCDF temperature, parameter
        
        integer, intent(out) :: albedo_dynamic_annual_varid
        integer, intent(out) :: latent_heat_annual_varid

                ! output
        integer, intent(out) :: accum_annual_varid        ! Variable for the netCDF Height, parameter
        integer, intent(out) :: rain_annual_varid           ! Variable for the netCDF Bedrock, parameter
        integer, intent(out) :: melt_annual_varid           ! Variable for the netCDF accumulation, parameter
        integer, intent(out) :: refreezing_annual_varid   ! Variable for the netCDF Diffusivity, parameter
        integer, intent(out) :: snow_annual_varid   ! Variable for the netCDF Diffusivity, parameter
        integer, intent(out) :: runoff_annual_varid   ! Variable for the netCDF Diffusivity, parameter
        
        integer, intent(out) :: snow_mask_annual_varid
        
        integer, intent(out) :: melt_ice_annual_varid
        integer, intent(out) :: snow_temp_ave_surface_annual_varid
        integer, intent(out) :: real_mass_balance_annual_varid
        
        integer, intent(out) :: regridding_annual_varid
        ! return value
        integer :: retval

        ! Copied from Example file!
        !--------------------------
        ! We are writing 3D data with time, We will need 4 netCDF dimensions. (Lats, Long, z, Time)
        integer, parameter :: NDIMS = 4

        character*(*) :: LAT_NAME, LON_NAME, REC_NAME, Z_NAME, ANNUAL_NAME
        parameter (LAT_NAME='y', LON_NAME='x',Z_NAME='layer')
        parameter (REC_NAME = 'time')
        parameter (ANNUAL_NAME = 'annual')
        integer lon_dimid, lat_dimid, z_dimid, rec_dimid, annual_dimid


        ! The start and count arrays will tell the netCDF library where to
        ! write our data.
        !integer start(NDIMS), count(NDIMS)

        ! In addition to the latitude and longitude dimensions, we will also
        ! create latitude and longitude netCDF variables which will hold the
        ! actual latitudes and longitudes. Since they hold data about the
        ! coordinate system, the netCDF term for these is: "coordinate
        ! variables."
        real :: lats(NLATS)
        real :: lons(NLONS)
        real :: zs(N_LAYER)
        real, dimension(n_time) :: timecounter
        integer lat_varid, lon_varid, rec_varid, z_varid, annual_varid
        real START_LAT, START_LON
        parameter (START_LAT = 0, START_LON = 0)
        !
        !    ! Variables
        
        !    ! Variables
        character*(*) snowman_annual_NAME, lwmass_annual_NAME, rho_snow_annual_NAME, snow_temp_annual_NAME
        parameter (snowman_annual_NAME='snowmass')
        parameter (lwmass_annual_NAME='lwmass')
        parameter (rho_snow_annual_NAME='snowdensity')
        parameter (snow_temp_annual_NAME='snowtemp')
        
        character*(*) albedo_dynamic_annual_NAME,latent_heat_annual_NAME, melt_ice_annual_NAME
        character*(*) snow_temp_ave_surface_annual_NAME, regridding_annual_NAME, real_mass_balance_annual_NAME
        parameter (albedo_dynamic_annual_NAME='albedo')
        parameter (latent_heat_annual_NAME='latent_heat_flux')
        parameter (melt_ice_annual_NAME='melt_of_ice')
        parameter (snow_temp_ave_surface_annual_NAME='averaged_surface_temp')
        parameter (regridding_annual_NAME='amount_of_regrids')
        parameter (real_mass_balance_annual_NAME='mass_balance')
        
        character*(*) SNOW_MASK_annual_NAME
        parameter (SNOW_MASK_annual_NAME = 'snow_mask')
        
                !    ! Variables 2
        character*(*) accum_annual_NAME, rain_annual_NAME, melt_annual_NAME, &
                    refreezing_annual_NAME, snow_annual_NAME, runoff_annual_NAME
        parameter (accum_annual_NAME='accum_ice_model')
        parameter (rain_annual_NAME='rain')
        parameter (melt_annual_NAME='melt')
        parameter (refreezing_annual_NAME='refreezing')
        parameter (snow_annual_NAME='snow')
        parameter (runoff_annual_NAME='runoff')
        
        character*(*) snowman_annual_STDNAME, lwmass_annual_STDNAME, rho_snow_annual_STDNAME, &
        snow_temp_annual_STDNAME
        parameter (snowman_annual_STDNAME='surface_snow_amount')
        parameter (lwmass_annual_STDNAME='liquid_water_content_of_snow_layer')
        parameter (rho_snow_annual_STDNAME='snow_density')
        parameter (snow_temp_annual_STDNAME='snow_temperature')
        
        character*(*) albedo_dynamic_annual_STDNAME,latent_heat_annual_STDNAME, melt_ice_annual_STDNAME
        character*(*) snow_temp_ave_surface_annual_STDNAME, regridding_annual_STDNAME, &
        real_mass_balance_annual_STDNAME
        parameter (albedo_dynamic_annual_STDNAME='surface_albedo')
        parameter (latent_heat_annual_STDNAME='turbulent_latent_heat_mass_flux')
        parameter (melt_ice_annual_STDNAME='ice_melt_flux')
        parameter (snow_temp_ave_surface_annual_STDNAME='yearly_average_surface_temperature')
        parameter (regridding_annual_STDNAME='amount_of_regrids')
        parameter (real_mass_balance_annual_STDNAME='land_ice_surface_specific_mass_balance_flux')
        
        character*(*) SNOW_MASK_annual_STDNAME
        parameter (SNOW_MASK_annual_STDNAME = 'surface_snow_mask')
        
                !    ! Variables 2
        character*(*) accum_annual_STDNAME, rain_annual_STDNAME, melt_annual_STDNAME, &
                    refreezing_annual_STDNAME, snow_annual_STDNAME, runoff_annual_STDNAME
        parameter (accum_annual_STDNAME='accumulation_ice_model')
        parameter (rain_annual_STDNAME='rainfall_flux')
        parameter (melt_annual_STDNAME='surface_snow_and_ice_melt_flux')
        parameter (refreezing_annual_STDNAME='surface_snow_and_ice_refreezing_flux')
        parameter (snow_annual_STDNAME='snowfall_flux')
        parameter (runoff_annual_STDNAME='runoff_flux')
        
        integer dimids3D(NDIMS)
        integer dimids2D(NDIMS-1)
        integer dimids_annual(NDIMS)
        !
        
        
        character*(*) :: FillValue
        parameter (FillValue = '_FillValue')
        character*(*) :: STANDARDNAME
        parameter (STANDARDNAME = 'standard_name')
        !    ! It's good practice for each variable to carry a "units" attribute.
        character*(*) :: UNITS
        parameter (UNITS = 'units')
        character*(*) snowman_annual_UNITS, lwmass_annual_UNITS, rho_snow_annual_UNITS, snow_temp_annual_UNITS, &
            albedo_dynamic_annual_UNITS, latent_heat_annual_UNITS, melt_ice_annual_UNITS, &
            snow_temp_ave_surface_annual_UNITS, regridding_annual_UNITS, real_mass_balance_annual_UNITS, &
            LAT_UNITS, LON_UNITS, REC_UNITS, z_UNITS, ANNUAL_UNITS 
        character*(*) accum_annual_UNITS, rain_annual_UNITS, melt_annual_UNITS, &
            refreezing_annual_UNITS, snow_annual_UNITS, runoff_annual_UNITS
        parameter ( snowman_annual_UNITS = 'kg/m2', lwmass_annual_UNITS = 'kg/m2',regridding_annual_UNITS='#', &
            rho_snow_annual_UNITS = 'kg/m3', snow_temp_annual_UNITS = 'K', snow_temp_ave_surface_annual_UNITS='K', &
            albedo_dynamic_annual_UNITS = '', latent_heat_annual_UNITS='kg/m2/a', melt_ice_annual_UNITS='kg/m2/a')
        parameter (accum_annual_UNITS = 'kg/m2/a', rain_annual_UNITS = 'kg/m2/a', melt_annual_UNITS = 'kg/m2/a', &
             refreezing_annual_UNITS = 'kg/m2/a', snow_annual_UNITS = 'kg/m2/a', runoff_annual_UNITS = 'kg/m2/a',&
             real_mass_balance_annual_UNITS='kg/m2/a')
        ! Units of Dimensions
        parameter (LAT_UNITS = 'm')
        parameter (LON_UNITS = 'm')
        parameter (z_UNITS = 'box')
        parameter (ANNUAL_UNITS = 'a')
        parameter (REC_UNITS = '1000years')
        
        integer :: tttt, lat, lon, zz
        character(256)::parameters
        character(512)::settings
        ! Distance of the lat/long points
        print*,filename
        do lat = 1, NLATS
            lats(lat) = START_LAT + (lat - 1) * lat_distance
        end do
        do lon = 1, NLONS
            lons(lon) = START_LON + (lon - 1) * long_distance
        end do
        do zz = 1, N_LAYER
            zs(zz) = 1 + (real(zz) - 1) * z_distance
        end do
        do tttt = 1, n_time
            timecounter(tttt)= real(tttt)
        end do
        !print *, '*** Init NetCDF file: ', TRIM(adjustl(filename))

        ! Create the file.
        retval = nf_create(filename, nf_clobber, ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !
        !
        !    ! Define dimensions
        retval = nf_def_dim(ncid, LAT_NAME, NLATS, lat_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, LON_NAME, NLONS, lon_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, z_NAME, N_LAYER, z_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        retval = nf_def_dim(ncid, ANNUAL_NAME,n_time, ANNUAL_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        
        ! TODO Time - this should be OK
        retval = nf_def_dim(ncid, REC_NAME, NF_UNLIMITED, rec_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Define the coordinate variables. They will hold the coordinate
        ! information, that is, the latitudes and longitudes. A varid is
        ! returned for each.
        !retval = nf_def_var(ncid, LAT_NAME, NF_REAL, 1, lat_dimid, lat_varid)
        retval = nf_def_var(ncid, LAT_NAME, NF_REAL, 1, lat_dimid, lat_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, LON_NAME, NF_REAL, 1, lon_dimid, lon_varid)
        retval = nf_def_var(ncid, LON_NAME, NF_REAL, 1, lon_dimid, lon_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
         !retval = nf_def_var(ncid, z_NAME, NF_REAL, 1, z_dimid, z_varid)
        retval = nf_def_var(ncid, Z_NAME, NF_REAL, 1, z_dimid, z_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_def_var(ncid, ANNUAL_NAME, NF_REAL, 1, annual_dimid, annual_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! TODO Time
        !retval = nf_def_var(ncid, REC_NAME, NF_REAL, 1, rec_dimid, rec_varid)
        retval = nf_def_var(ncid, REC_NAME, NF_REAL, 1, rec_dimid, rec_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Assign units attributes to coordinate var data. This attaches a
        ! text attribute to each of the coordinate variables, containing the
        ! units.
        retval = nf_put_att_text(ncid, lat_varid, UNITS, len(LAT_UNITS), LAT_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, lon_varid, UNITS, len(LON_UNITS), LON_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, z_varid, UNITS, len(Z_UNITS), Z_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_put_att_text(ncid, annual_varid, UNITS, len(ANNUAL_UNITS), ANNUAL_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)

        retval = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'title', len(TRIM(adjustl(title))), &
        TRIM(adjustl(title)))
         if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'history', len(TRIM(adjustl(history))), &
        TRIM(adjustl(history)))
         if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'source', len(TRIM(adjustl(source))), &
        TRIM(adjustl(source)))
         if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'institution',len(TRIM(adjustl(institution))),&
        TRIM(adjustl(institution)))
         if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'contact author', len(TRIM(adjustl(contact))), &
        TRIM(adjustl(contact)))
         if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'references', len(TRIM(adjustl(references))), &
        TRIM(adjustl(references)))
         if (retval .ne. nf_noerr) call handle_err(retval)
         
        retval = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'forcing_data', len(TRIM(adjustl(new_input_file3))),&
        TRIM(adjustl(new_input_file3)))
         if (retval .ne. nf_noerr) call handle_err(retval)
         
        retval = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'comment', len(TRIM(adjustl(comment))), &
        TRIM(adjustl(comment)))
         if (retval .ne. nf_noerr) call handle_err(retval)
        
		write (parameters,'(A5,F8.4,A6,F8.4,A5,F8.4,A6,F8.4,A11,F8.4,A6,F8.4,A9,F8.4)') &
		"a_fs:",albedo_snow_new_spam," a_fi:",albedo_snow_wet_spam," a_i:",albedo_ice_spam,&
		" D_SH:",D_sf_spam," D_LH/D_SH:",ratio,&
		" e_air:",eps_air_spam," max_lwc:",max_lwc_spam
		
        retval=NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'parameters', len(TRIM(adjustl(parameters))), &
        TRIM(adjustl(parameters)))
         if (retval .ne. nf_noerr) call handle_err(retval)
        
        write (settings,'(A15,I1,A21,L1,A24,L1,A22,L1,A15,L1,A23,L1,A20,L1)') &
		" albedo_module:",albedo_module, " latent_heat_flux_on:",latent_heat_flux_on, &
		" longwave_from_air_temp:",longwave_from_air_temp, &
		" longwave_downscaling:", longwave_downscaling, &
		" densification:", densification_model, &
		" temperature_lapserate:",active_temperature_lapsrate, &
		" dewpoint_lapserate:",active_dewpoint_lapserate
		
        retval=NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'settings', len(TRIM(adjustl(settings))), &
        TRIM(adjustl(settings)))
         if (retval .ne. nf_noerr) call handle_err(retval)
         
        !Calculater name
        !contact details 
        if (retval .ne. nf_noerr) call handle_err(retval)
        ! TODO Time
        retval = nf_put_att_text(ncid, rec_varid, UNITS, len(REC_UNITS), REC_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Define the netCDF variables. The dimids array is used to pass the
        ! dimids of the dimensions of the netCDF variables.
        dimids3D(1) = lon_dimid
        dimids3D(2) = lat_dimid
        dimids3D(3) = z_dimid
        dimids3D(4) = rec_dimid
        
        dimids2D(1) = lon_dimid
        dimids2D(2) = lat_dimid
        dimids2D(3) = rec_dimid
        
!         dimids_annual(1)=lon_dimid
!         dimids_annual(2)=lat_dimid
!         dimids_annual(3)=z_dimid
!         dimids_annual(4)=annual_dimid
!         print*,'setup vars'
        
        ! Define the netCDF variables for the "real" values
        retval = nf_def_var(ncid, snowman_annual_NAME, NF_REAL, NDIMS, dimids3D,snowman_annual_varid )
        !retval = nf_def_var(ncid, snowman_NAME, NF_DOUBLE, NDIMS, dimids, snowman_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, lwmass_annual_NAME, NF_REAL, NDIMS, dimids3D, lwmass_annual_varid)
        !retval = nf_def_var(ncid, lwmass_NAME, NF_DOUBLE, NDIMS, dimids, lwmass_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, rho_snow_annual_NAME, NF_REAL, NDIMS, dimids3D, rho_snow_annual_varid)
        !retval = nf_def_var(ncid, rho_snow_NAME, NF_DOUBLE, NDIMS, dimids, rho_snow_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, snow_temp_annual_NAME, NF_REAL, NDIMS, dimids3D, snow_temp_annual_varid)
        !retval = nf_def_var(ncid, snow_temp_NAME, NF_DOUBLE, NDIMS, dimids, snow_temp_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
!         print*,'3Dvars'
        retval = nf_def_var(ncid, albedo_dynamic_annual_NAME, NF_REAL, NDIMS-1, dimids2D, albedo_dynamic_annual_varid)
        !retval = nf_def_var(ncid, snow_temp_NAME, NF_DOUBLE, NDIMS, dimids, snow_temp_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_def_var(ncid, latent_heat_annual_NAME, NF_REAL, NDIMS-1, dimids2D, latent_heat_annual_varid)
        !retval = nf_def_var(ncid, snow_temp_NAME, NF_DOUBLE, NDIMS, dimids, snow_temp_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
!         print*,'Newvars'
        ! Define the netCDF variables for the "real" values
        !retval = nf_def_var(ncid, accum_NAME, NF_REAL, NDIMS, dimids, accum_varid)
        retval = nf_def_var(ncid, accum_annual_NAME, NF_REAL, NDIMS-1, dimids2D, accum_annual_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, rain_NAME, NF_REAL, NDIMS, dimids, rain_varid)
        retval = nf_def_var(ncid, rain_annual_NAME, NF_REAL, NDIMS-1, dimids2D, rain_annual_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, melt_NAME, NF_REAL, NDIMS, dimids, melt_varid)
        retval = nf_def_var(ncid, melt_annual_NAME, NF_REAL, NDIMS-1, dimids2D, melt_annual_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, refreezing_NAME, NF_REAL, NDIMS, dimids, refreezing_varid)
        retval = nf_def_var(ncid, refreezing_annual_NAME, NF_REAL, NDIMS-1, dimids2D, refreezing_annual_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, snow_NAME, NF_REAL, NDIMS, dimids, snow_varid)
        retval = nf_def_var(ncid, snow_annual_NAME, NF_REAL, NDIMS-1, dimids2D, snow_annual_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, runoff_NAME, NF_REAL, NDIMS, dimids, runoff_varid)
        retval = nf_def_var(ncid, runoff_annual_NAME, NF_REAL, NDIMS-1, dimids2D, runoff_annual_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
!         print*,'oldvars'
       
        retval = nf_def_var(ncid, SNOW_MASK_annual_NAME, NF_INT, NDIMS-1, dimids2D, snow_mask_annual_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
!         print*,'snowmask'
        !melt of ice melt_ice_annual_NAME
        retval = nf_def_var(ncid, melt_ice_annual_NAME, NF_REAL, NDIMS-1, dimids2D, melt_ice_annual_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !average_surface_temp
        retval = nf_def_var(ncid, snow_temp_ave_surface_annual_NAME, NF_REAL, NDIMS-1, dimids2D, snow_temp_ave_surface_annual_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, regridding_annual_NAME, NF_INT, NDIMS-1, dimids2D, regridding_annual_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, real_mass_balance_annual_NAME, NF_REAL, NDIMS-1, dimids2D, real_mass_balance_annual_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)        
        
        
!         print*,'3Dvars'
        ! Assign units attributes to the pressure and temperature netCDF
        ! Assign standard_name attributes
        ! variables.
        retval = nf_put_att_text(ncid, snowman_annual_varid, UNITS, len(snowman_annual_UNITS), snowman_annual_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, lwmass_annual_varid, UNITS, len(lwmass_annual_UNITS), lwmass_annual_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, rho_snow_annual_varid, UNITS, len(rho_snow_annual_UNITS), rho_snow_annual_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, snow_temp_annual_varid, UNITS, len(snow_temp_annual_UNITS), snow_temp_annual_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_put_att_text(ncid, snowman_annual_varid, STANDARDNAME, len(snowman_annual_STDNAME),&
        snowman_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, lwmass_annual_varid, STANDARDNAME, len(lwmass_annual_STDNAME), &
        lwmass_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, rho_snow_annual_varid, STANDARDNAME, len(rho_snow_annual_STDNAME), &
        rho_snow_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, snow_temp_annual_varid, STANDARDNAME, len(snow_temp_annual_STDNAME), &
        snow_temp_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
! !            print*,'Newvars'
            
        retval = nf_put_att_text(ncid, albedo_dynamic_annual_varid, UNITS, &
        len(albedo_dynamic_annual_UNITS), albedo_dynamic_annual_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, latent_heat_annual_varid, UNITS, &
        len(latent_heat_annual_UNITS), latent_heat_annual_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, accum_annual_varid, UNITS, len(accum_annual_UNITS), accum_annual_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, rain_annual_varid, UNITS, len(rain_annual_UNITS), rain_annual_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, melt_annual_varid, UNITS, len(melt_annual_UNITS), melt_annual_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, refreezing_annual_varid, UNITS, &
        len(refreezing_annual_UNITS), refreezing_annual_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, snow_annual_varid, UNITS, len(snow_annual_UNITS), snow_annual_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, runoff_annual_varid, UNITS, len(runoff_annual_UNITS), runoff_annual_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, melt_ice_annual_varid, &
        UNITS, len(melt_ice_annual_UNITS), melt_ice_annual_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)    
        retval = nf_put_att_text(ncid, snow_temp_ave_surface_annual_varid, &
        UNITS, len(snow_temp_ave_surface_annual_UNITS), snow_temp_ave_surface_annual_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)  
        retval = nf_put_att_text(ncid, regridding_annual_varid, &
        UNITS, len(regridding_annual_UNITS), regridding_annual_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, real_mass_balance_annual_varid, &
        UNITS, len(real_mass_balance_annual_UNITS), real_mass_balance_annual_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
         retval = nf_put_att_text(ncid, albedo_dynamic_annual_varid, STANDARDNAME, &
        len(albedo_dynamic_annual_STDNAME), albedo_dynamic_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, latent_heat_annual_varid, STANDARDNAME, &
        len(latent_heat_annual_STDNAME), latent_heat_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, accum_annual_varid, STANDARDNAME, &
        len(accum_annual_STDNAME), accum_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, rain_annual_varid, STANDARDNAME, &
        len(rain_annual_STDNAME), rain_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, melt_annual_varid, STANDARDNAME,&
        len(melt_annual_STDNAME), melt_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, refreezing_annual_varid, STANDARDNAME, &
        len(refreezing_annual_STDNAME), refreezing_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, snow_annual_varid, STANDARDNAME, &
        len(snow_annual_STDNAME), snow_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, runoff_annual_varid, STANDARDNAME, &
        len(runoff_annual_STDNAME), runoff_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, melt_ice_annual_varid, &
        STANDARDNAME, len(melt_ice_annual_STDNAME), melt_ice_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)    
        retval = nf_put_att_text(ncid, snow_temp_ave_surface_annual_varid, &
        STANDARDNAME, len(snow_temp_ave_surface_annual_STDNAME), snow_temp_ave_surface_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)  
        retval = nf_put_att_text(ncid, regridding_annual_varid, &
        STANDARDNAME, len(regridding_annual_STDNAME), regridding_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, real_mass_balance_annual_varid, &
        STANDARDNAME, len(real_mass_balance_annual_STDNAME), real_mass_balance_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        
      
        ! End define mode.
        retval = nf_enddef(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Write the coordinate variable data. This will put the latitudes
        ! and longitudes of our data grid into the netCDF file.
        retval = nf_put_var_real(ncid, lat_varid, lats)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_real(ncid, lon_varid, lons)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_real(ncid, z_varid, zs)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_real(ncid, rec_varid, timecounter)
        if (retval .ne. nf_noerr) call handle_err(retval)
        print *,'*** Output netCDF file defined: ', TRIM(adjustl(filename))

    end subroutine init_netcdf_snow_file
    
    subroutine init_netcdf_snow_file_time_first(filename, ncid, nlats, nlons, N_LAYER, n_time, lat_distance, long_distance, &
        z_distance, snowman_daily_varid, lwmass_daily_varid, rho_snow_daily_varid, snow_temp_daily_varid, &
        albedo_dynamic_daily_varid, latent_heat_daily_varid, accum_daily_varid, rain_daily_varid, &
        melt_daily_varid, refreezing_daily_varid, snow_daily_varid, runoff_daily_varid, &
        melt_ice_daily_varid, snow_mask_daily_varid, snow_temp_ave_surface_daily_varid, regridding_daily_varid,&
        real_mass_balance_daily_varid,new_input_file3)
        
        use netcdf
        use variables
        implicit none
        !Data file with more than one timestep
        character, intent(in) :: new_input_file3
        character(len=*), intent(in) :: filename
        integer, intent(out) :: ncid
        ! input
        integer, intent(in) :: NLATS
        integer, intent(in) :: NLONS
        integer, intent(in) :: N_LAYER
        integer, intent(in) :: n_time
        integer, intent(in) :: lat_distance
        integer, intent(in) :: long_distance
        integer, intent(in) :: z_distance

        ! output
        integer, intent(out) :: snowman_daily_varid        ! Variable for the netCDF snow mass, parameter
        integer, intent(out) :: lwmass_daily_varid           ! Variable for the netCDF liquid water, parameter
        integer, intent(out) :: rho_snow_daily_varid           ! Variable for the netCDF density, parameter
        integer, intent(out) :: snow_temp_daily_varid           ! Variable for the netCDF temperature, parameter
        
        integer, intent(out) :: albedo_dynamic_daily_varid
        integer, intent(out) :: latent_heat_daily_varid

                ! output
        integer, intent(out) :: accum_daily_varid        ! Variable for the netCDF Height, parameter
        integer, intent(out) :: rain_daily_varid           ! Variable for the netCDF Bedrock, parameter
        integer, intent(out) :: melt_daily_varid           ! Variable for the netCDF accumulation, parameter
        integer, intent(out) :: refreezing_daily_varid   ! Variable for the netCDF Diffusivity, parameter
        integer, intent(out) :: snow_daily_varid   ! Variable for the netCDF Diffusivity, parameter
        integer, intent(out) :: runoff_daily_varid   ! Variable for the netCDF Diffusivity, parameter
        
        integer, intent(out) :: snow_mask_daily_varid
        
        integer, intent(out) :: melt_ice_daily_varid
        integer, intent(out) :: snow_temp_ave_surface_daily_varid
        integer, intent(out) :: real_mass_balance_daily_varid
        
        integer, intent(out) :: regridding_daily_varid
        ! return value
        integer :: retval

        ! Copied from Example file!
        !--------------------------
        ! We are writing 3D data with time, We will need 4 netCDF dimensions. (Lats, Long, z, Time)
        integer, parameter :: NDIMS = 4

        character*(*) :: LAT_NAME, LON_NAME, REC_NAME, Z_NAME, ANNUAL_NAME
        parameter (LAT_NAME='y', LON_NAME='x',Z_NAME='layer')
        parameter (REC_NAME = 'time')
        parameter (ANNUAL_NAME = 'annual')
        integer lon_dimid, lat_dimid, z_dimid, rec_dimid, annual_dimid


        ! The start and count arrays will tell the netCDF library where to
        ! write our data.
        !integer start(NDIMS), count(NDIMS)

        ! In addition to the latitude and longitude dimensions, we will also
        ! create latitude and longitude netCDF variables which will hold the
        ! actual latitudes and longitudes. Since they hold data about the
        ! coordinate system, the netCDF term for these is: "coordinate
        ! variables."
        real :: lats(NLATS)
        real :: lons(NLONS)
        real :: zs(N_LAYER)
        real, dimension(n_time) :: timecounter
        integer lat_varid, lon_varid, rec_varid, z_varid, annual_varid
        real START_LAT, START_LON
        parameter (START_LAT = 0, START_LON = 0)
        !
        !    ! Variables
        character*(*) snowman_daily_NAME, lwmass_daily_NAME, rho_snow_daily_NAME, snow_temp_daily_NAME
        parameter (snowman_daily_NAME='snowmass')
        parameter (lwmass_daily_NAME='lwmass')
        parameter (rho_snow_daily_NAME='snowdensity')
        parameter (snow_temp_daily_NAME='snowtemp')
        
        character*(*) albedo_dynamic_daily_NAME,latent_heat_daily_NAME, melt_ice_daily_NAME
        character*(*) snow_temp_ave_surface_daily_NAME, regridding_daily_NAME, real_mass_balance_daily_NAME
        parameter (albedo_dynamic_daily_NAME='albedo')
        parameter (latent_heat_daily_NAME='latent_heat_flux')
        parameter (melt_ice_daily_NAME='melt_of_ice')
        parameter (snow_temp_ave_surface_daily_NAME='cummulated_SMB')
        parameter (regridding_daily_NAME='amount_of_regrids')
        parameter (real_mass_balance_daily_NAME='mass_balance')
        
        character*(*) SNOW_MASK_daily_NAME
        parameter (SNOW_MASK_daily_NAME = 'snow_mask')
        
                !    ! Variables 2
        character*(*) accum_daily_NAME, rain_daily_NAME, melt_daily_NAME, &
                    refreezing_daily_NAME, snow_daily_NAME, runoff_daily_NAME
        parameter (accum_daily_NAME='accum_ice_model')
        parameter (rain_daily_NAME='rain')
        parameter (melt_daily_NAME='melt')
        parameter (refreezing_daily_NAME='refreezing')
        parameter (snow_daily_NAME='snow')
        parameter (runoff_daily_NAME='runoff')
        
        integer dimids3D(NDIMS)
        integer dimids2D(NDIMS-1)
        integer dimids_daily(NDIMS)
        !
         
        character*(*) :: FillValue
        parameter (FillValue = '_FillValue')
        character*(*) :: STANDARDNAME
        parameter (STANDARDNAME = 'standard_name')
        !    ! It's good practice for each variable to carry a "units" attribute.
        character*(*) :: UNITS
        parameter (UNITS = 'units')
        character*(*) snowman_daily_UNITS, lwmass_daily_UNITS, rho_snow_daily_UNITS, snow_temp_daily_UNITS, &
            albedo_dynamic_daily_UNITS, latent_heat_daily_UNITS, melt_ice_daily_UNITS, &
            snow_temp_ave_surface_daily_UNITS, regridding_daily_UNITS, real_mass_balance_daily_UNITS, &
            LAT_UNITS, LON_UNITS, REC_UNITS, z_UNITS, ANNUAL_UNITS 
        character*(*) accum_daily_UNITS, rain_daily_UNITS, melt_daily_UNITS, &
            refreezing_daily_UNITS, snow_daily_UNITS, runoff_daily_UNITS
        parameter ( snowman_daily_UNITS = 'kg/m2', lwmass_daily_UNITS = 'kg/m2',regridding_daily_UNITS='#', &
            rho_snow_daily_UNITS = 'kg/m3', snow_temp_daily_UNITS = 'K', &
            snow_temp_ave_surface_daily_UNITS='kg/m2', &
            albedo_dynamic_daily_UNITS = '', latent_heat_daily_UNITS='kg/m2/d', melt_ice_daily_UNITS='kg/m2/d')
        parameter (accum_daily_UNITS = 'kg/m2/d', rain_daily_UNITS = 'kg/m2/d', melt_daily_UNITS = 'kg/m2/d', &
             refreezing_daily_UNITS = 'kg/m2/d', snow_daily_UNITS = 'kg/m2/d', runoff_daily_UNITS = 'kg/m2/d',&
             real_mass_balance_daily_UNITS='kg/m2/d')
        ! Units of Dimensions
        parameter (LAT_UNITS = 'm')
        parameter (LON_UNITS = 'm')
        parameter (z_UNITS = 'box')
        parameter (ANNUAL_UNITS = 'd')
        parameter (REC_UNITS = '1000years')
        
        
        character*(*) snowman_annual_STDNAME, lwmass_annual_STDNAME, rho_snow_annual_STDNAME, &
        snow_temp_annual_STDNAME
        parameter (snowman_annual_STDNAME='surface_snow_amount')
        parameter (lwmass_annual_STDNAME='liquid_water_content_of_snow_layer')
        parameter (rho_snow_annual_STDNAME='snow_density')
        parameter (snow_temp_annual_STDNAME='snow_temperature')
        
        character*(*) albedo_dynamic_annual_STDNAME,latent_heat_annual_STDNAME, melt_ice_annual_STDNAME
        character*(*) snow_temp_ave_surface_annual_STDNAME, regridding_annual_STDNAME, &
        real_mass_balance_annual_STDNAME
        parameter (albedo_dynamic_annual_STDNAME='surface_albedo')
        parameter (latent_heat_annual_STDNAME='turbulent_latent_heat_mass_flux')
        parameter (melt_ice_annual_STDNAME='ice_melt_flux')
        parameter (snow_temp_ave_surface_annual_STDNAME='land_ice_surface_cummulated_mass_balance')
        parameter (regridding_annual_STDNAME='amount_of_regrids')
        parameter (real_mass_balance_annual_STDNAME='land_ice_surface_specific_mass_balance_flux')
        
        character*(*) SNOW_MASK_annual_STDNAME
        parameter (SNOW_MASK_annual_STDNAME = 'surface_snow_mask')
        
                !    ! Variables 2
        character*(*) accum_annual_STDNAME, rain_annual_STDNAME, melt_annual_STDNAME, &
                    refreezing_annual_STDNAME, snow_annual_STDNAME, runoff_annual_STDNAME
        parameter (accum_annual_STDNAME='accumulation_ice_model')
        parameter (rain_annual_STDNAME='rainfall_flux')
        parameter (melt_annual_STDNAME='surface_snow_and_ice_melt_flux')
        parameter (refreezing_annual_STDNAME='surface_snow_and_ice_refreezing_flux')
        parameter (snow_annual_STDNAME='snowfall_flux')
        parameter (runoff_annual_STDNAME='runoff_flux')
        
        
        character(256)::parameters
        character(512)::settings
        
        integer :: tttt, lat, lon, zz
        ! Distance of the lat/long points
        do lat = 1, NLATS
            lats(lat) = START_LAT + (lat - 1) * lat_distance
        end do
        do lon = 1, NLONS
            lons(lon) = START_LON + (lon - 1) * long_distance
        end do
        do zz = 1, N_LAYER
            zs(zz) = 1 + (real(zz) - 1) * z_distance
        end do
        do tttt = 1, n_time
            timecounter(tttt)= real(tttt)
        end do
        !print *, '*** Init NetCDF file: ', TRIM(adjustl(filename))

        ! Create the file.
         retval = nf90_create(path=filename, cmode=or(nf90_clobber,nf90_64bit_offset), ncid=ncid)
!         retval = nf_create(filename, nf_clobber, ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !
        !
        !    ! Define dimensions
        retval = nf_def_dim(ncid, LAT_NAME, NLATS, lat_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, LON_NAME, NLONS, lon_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_dim(ncid, z_NAME, N_LAYER, z_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        retval = nf_def_dim(ncid, ANNUAL_NAME,1, ANNUAL_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        
        ! TODO Time - this should be OK
        retval = nf_def_dim(ncid, REC_NAME, n_time, rec_dimid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Define the coordinate variables. They will hold the coordinate
        ! information, that is, the latitudes and longitudes. A varid is
        ! returned for each.
        !retval = nf_def_var(ncid, LAT_NAME, NF_REAL, 1, lat_dimid, lat_varid)
        retval = nf_def_var(ncid, LAT_NAME, NF_REAL, 1, lat_dimid, lat_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, LON_NAME, NF_REAL, 1, lon_dimid, lon_varid)
        retval = nf_def_var(ncid, LON_NAME, NF_REAL, 1, lon_dimid, lon_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
         !retval = nf_def_var(ncid, z_NAME, NF_REAL, 1, z_dimid, z_varid)
        retval = nf_def_var(ncid, Z_NAME, NF_REAL, 1, z_dimid, z_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_def_var(ncid, ANNUAL_NAME, NF_REAL, 1, annual_dimid, annual_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! TODO Time
        !retval = nf_def_var(ncid, REC_NAME, NF_REAL, 1, rec_dimid, rec_varid)
        retval = nf_def_var(ncid, REC_NAME, NF_REAL, 1, rec_dimid, rec_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Assign units attributes to coordinate var data. This attaches a
        ! text attribute to each of the coordinate variables, containing the
        ! units.
        retval = nf_put_att_text(ncid, lat_varid, UNITS, len(LAT_UNITS), LAT_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, lon_varid, UNITS, len(LON_UNITS), LON_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, z_varid, UNITS, len(Z_UNITS), Z_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_put_att_text(ncid, annual_varid, UNITS, len(ANNUAL_UNITS), ANNUAL_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)


        ! TODO Time
        retval = nf_put_att_text(ncid, rec_varid, UNITS, len(REC_UNITS), REC_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Define the netCDF variables. The dimids array is used to pass the
        ! dimids of the dimensions of the netCDF variables.
        dimids3D(1) = rec_dimid
        dimids3D(2) = z_dimid
        dimids3D(3) = lon_dimid
        dimids3D(4) = lat_dimid
        
        dimids2D(1) = rec_dimid
        dimids2D(2) = lon_dimid
        dimids2D(3) = lat_dimid
        
     !   print*,'setup vars'
        
        ! Define the netCDF variables for the "real" values
        retval = nf_def_var(ncid, snowman_daily_NAME, NF_REAL, NDIMS, dimids3D,snowman_daily_varid )
        !retval = nf_def_var(ncid, snowman_NAME, NF_DOUBLE, NDIMS, dimids, snowman_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, lwmass_daily_NAME, NF_REAL, NDIMS, dimids3D, lwmass_daily_varid)
        !retval = nf_def_var(ncid, lwmass_NAME, NF_DOUBLE, NDIMS, dimids, lwmass_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, rho_snow_daily_NAME, NF_REAL, NDIMS, dimids3D, rho_snow_daily_varid)
        !retval = nf_def_var(ncid, rho_snow_NAME, NF_DOUBLE, NDIMS, dimids, rho_snow_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, snow_temp_daily_NAME, NF_REAL, NDIMS, dimids3D, snow_temp_daily_varid)
        !retval = nf_def_var(ncid, snow_temp_NAME, NF_DOUBLE, NDIMS, dimids, snow_temp_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
! !         
     !   print*,'3Dvars'
        retval = nf_def_var(ncid, albedo_dynamic_daily_NAME, NF_REAL, NDIMS-1, dimids2D, albedo_dynamic_daily_varid)
        !retval = nf_def_var(ncid, snow_temp_NAME, NF_DOUBLE, NDIMS, dimids, snow_temp_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = nf_def_var(ncid, latent_heat_daily_NAME, NF_REAL, NDIMS-1, dimids2D, latent_heat_daily_varid)
        !retval = nf_def_var(ncid, snow_temp_NAME, NF_DOUBLE, NDIMS, dimids, snow_temp_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
    !    print*,'Newvars'
        ! Define the netCDF variables for the "real" values
        !retval = nf_def_var(ncid, accum_NAME, NF_REAL, NDIMS, dimids, accum_varid)
        retval = nf_def_var(ncid, accum_daily_NAME, NF_REAL, NDIMS-1, dimids2D, accum_daily_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, rain_NAME, NF_REAL, NDIMS, dimids, rain_varid)
        retval = nf_def_var(ncid, rain_daily_NAME, NF_REAL, NDIMS-1, dimids2D, rain_daily_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, melt_NAME, NF_REAL, NDIMS, dimids, melt_varid)
        retval = nf_def_var(ncid, melt_daily_NAME, NF_REAL, NDIMS-1, dimids2D, melt_daily_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, refreezing_NAME, NF_REAL, NDIMS, dimids, refreezing_varid)
        retval = nf_def_var(ncid, refreezing_daily_NAME, NF_REAL, NDIMS-1, dimids2D, refreezing_daily_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, snow_NAME, NF_REAL, NDIMS, dimids, snow_varid)
        retval = nf_def_var(ncid, snow_daily_NAME, NF_REAL, NDIMS-1, dimids2D, snow_daily_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !retval = nf_def_var(ncid, runoff_NAME, NF_REAL, NDIMS, dimids, runoff_varid)
        retval = nf_def_var(ncid, runoff_daily_NAME, NF_REAL, NDIMS-1, dimids2D, runoff_daily_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
     !   print*,'oldvars'
       
        retval = nf_def_var(ncid, SNOW_MASK_daily_NAME, NF_INT, NDIMS-1, dimids2D, snow_mask_daily_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
   !     print*,'snowmask'
        !melt of ice melt_ice_daily_NAME
        retval = nf_def_var(ncid, melt_ice_daily_NAME, NF_REAL, NDIMS-1, dimids2D, melt_ice_daily_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        !average_surface_temp
        retval = nf_def_var(ncid, snow_temp_ave_surface_daily_NAME, NF_REAL, NDIMS-1, dimids2D, snow_temp_ave_surface_daily_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, regridding_daily_NAME, NF_INT, NDIMS-1, dimids2D, regridding_daily_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_def_var(ncid, real_mass_balance_daily_NAME, NF_REAL, NDIMS-1, dimids2D, real_mass_balance_daily_varid)
        if (retval .ne. nf_noerr) call handle_err(retval)        
        
        
        !print*,'3Dvars'
        ! Assign units attributes to the pressure and temperature netCDF
        ! variables.
        retval = nf_put_att_text(ncid, snowman_daily_varid, UNITS, len(snowman_daily_UNITS), snowman_daily_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, lwmass_daily_varid, UNITS, len(lwmass_daily_UNITS), lwmass_daily_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, rho_snow_daily_varid, UNITS, len(rho_snow_daily_UNITS), rho_snow_daily_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, snow_temp_daily_varid, UNITS, len(snow_temp_daily_UNITS), snow_temp_daily_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
          !  print*,'Newvars'
            
        retval = nf_put_att_text(ncid, albedo_dynamic_daily_varid, UNITS, &
        len(albedo_dynamic_daily_UNITS), albedo_dynamic_daily_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
  !  print*,'albedo'
        retval = nf_put_att_text(ncid, latent_heat_daily_varid, UNITS, &
        len(latent_heat_daily_UNITS), latent_heat_daily_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)

      !  print*,'oldvars'
        ! Assign units attributes to the pressure and rainerature netCDF
        ! variables.
        retval = nf_put_att_text(ncid, accum_daily_varid, UNITS, len(accum_daily_UNITS), accum_daily_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, rain_daily_varid, UNITS, len(rain_daily_UNITS), rain_daily_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, melt_daily_varid, UNITS, len(melt_daily_UNITS), melt_daily_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, refreezing_daily_varid, UNITS, len(refreezing_daily_UNITS), refreezing_daily_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, snow_daily_varid, UNITS, len(snow_daily_UNITS), snow_daily_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, runoff_daily_varid, UNITS, len(runoff_daily_UNITS), runoff_daily_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, melt_ice_daily_varid, &
        UNITS, len(melt_ice_daily_UNITS), melt_ice_daily_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)    
        retval = nf_put_att_text(ncid, snow_temp_ave_surface_daily_varid, &
        UNITS, len(snow_temp_ave_surface_daily_UNITS), snow_temp_ave_surface_daily_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)  
        retval = nf_put_att_text(ncid, regridding_daily_varid, &
        UNITS, len(regridding_daily_UNITS), regridding_daily_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, real_mass_balance_daily_varid, &
        UNITS, len(real_mass_balance_daily_UNITS), real_mass_balance_daily_UNITS)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        
        retval = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'title', len(TRIM(adjustl(title))), &
        TRIM(adjustl(title)))
         if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'history', len(TRIM(adjustl(history))), &
        TRIM(adjustl(history)))
         if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'source', len(TRIM(adjustl(source))), &
        TRIM(adjustl(source)))
         if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'institution',len(TRIM(adjustl(institution))),&
        TRIM(adjustl(institution)))
         if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'contact author', len(TRIM(adjustl(contact))), &
        TRIM(adjustl(contact)))
         if (retval .ne. nf_noerr) call handle_err(retval)
        
        retval = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'references', len(TRIM(adjustl(references))), &
        TRIM(adjustl(references)))
         if (retval .ne. nf_noerr) call handle_err(retval)
         
        retval = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'forcing_data', len(TRIM(adjustl(new_input_file3))),&
        TRIM(adjustl(new_input_file3)))     
         if (retval .ne. nf_noerr) call handle_err(retval)
         
        retval = NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'comment', len(TRIM(adjustl(comment))), &
        TRIM(adjustl(comment)))
         if (retval .ne. nf_noerr) call handle_err(retval)
        
		write (parameters,'(A5,F8.4,A6,F8.4,A5,F8.4,A6,F8.4,A11,F8.4,A6,F8.4,A9,F8.4)') &
		"a_fs:",albedo_snow_new_spam," a_fi:",albedo_snow_wet_spam," a_i:",albedo_ice_spam,&
		" D_SH:",D_sf_spam," D_LH/D_SH:",ratio,&
		" e_air:",eps_air_spam," max_lwc:",max_lwc_spam
		
        retval=NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'parameters', len(TRIM(adjustl(parameters))), &
        TRIM(adjustl(parameters)))
         if (retval .ne. nf_noerr) call handle_err(retval)
        
       write (settings,'(A15,I1,A21,L1,A24,L1,A22,L1,A15,L1,A23,L1,A20,L1)') &
		" albedo_module:",albedo_module, " latent_heat_flux_on:",latent_heat_flux_on, &
		" longwave_from_air_temp:",longwave_from_air_temp, &
		" longwave_downscaling:", longwave_downscaling, &
		" densification:", densification_model, &
		" temperature_lapserate:",active_temperature_lapsrate, &
		" dewpoint_lapserate:",active_dewpoint_lapserate
		
        retval=NF_PUT_ATT_TEXT (NCID, NF_GLOBAL, 'settings', len(TRIM(adjustl(settings))), &
        TRIM(adjustl(settings)))
         if (retval .ne. nf_noerr) call handle_err(retval)
        
        
         retval = nf_put_att_text(ncid, snowman_daily_varid, STANDARDNAME, len(snowman_annual_STDNAME),&
        snowman_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, lwmass_daily_varid, STANDARDNAME, len(lwmass_annual_STDNAME), &
        lwmass_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, rho_snow_daily_varid, STANDARDNAME, len(rho_snow_annual_STDNAME), &
        rho_snow_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, snow_temp_daily_varid, STANDARDNAME, len(snow_temp_annual_STDNAME), &
        snow_temp_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
         retval = nf_put_att_text(ncid, albedo_dynamic_daily_varid, STANDARDNAME, &
        len(albedo_dynamic_annual_STDNAME), albedo_dynamic_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, latent_heat_daily_varid, STANDARDNAME, &
        len(latent_heat_annual_STDNAME), latent_heat_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, accum_daily_varid, STANDARDNAME, &
        len(accum_annual_STDNAME), accum_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, rain_daily_varid, STANDARDNAME, &
        len(rain_annual_STDNAME), rain_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, melt_daily_varid, STANDARDNAME,&
        len(melt_annual_STDNAME), melt_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, refreezing_daily_varid, STANDARDNAME, &
        len(refreezing_annual_STDNAME), refreezing_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, snow_daily_varid, STANDARDNAME, &
        len(snow_annual_STDNAME), snow_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, runoff_daily_varid, STANDARDNAME, &
        len(runoff_annual_STDNAME), runoff_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, melt_ice_daily_varid, &
        STANDARDNAME, len(melt_ice_annual_STDNAME), melt_ice_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)    
        retval = nf_put_att_text(ncid, snow_temp_ave_surface_daily_varid, &
        STANDARDNAME, len(snow_temp_ave_surface_annual_STDNAME), snow_temp_ave_surface_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)  
        retval = nf_put_att_text(ncid, regridding_daily_varid, &
        STANDARDNAME, len(regridding_annual_STDNAME), regridding_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_att_text(ncid, real_mass_balance_daily_varid, &
        STANDARDNAME, len(real_mass_balance_annual_STDNAME), real_mass_balance_annual_STDNAME)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        ! End define mode.
        retval = nf_enddef(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Write the coordinate variable data. This will put the latitudes
        ! and longitudes of our data grid into the netCDF file.
        retval = nf_put_var_real(ncid, lat_varid, lats)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_real(ncid, lon_varid, lons)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_real(ncid, z_varid, zs)
        if (retval .ne. nf_noerr) call handle_err(retval)
        retval = nf_put_var_real(ncid, rec_varid, timecounter)
        if (retval .ne. nf_noerr) call handle_err(retval)
        print *,'*** Output netCDF file defined: ', TRIM(adjustl(filename))

    end subroutine init_netcdf_snow_file_time_first

    !=====================================================================================================

      subroutine writeNCDFSNOW3DValues(ncid, year, varid, values, nlats, nlons, N_LAYER)
        ! Writes the Values (Real(kind=8)) into the netCDF File
        ! Param ncid: File handle of the netCDF File
        ! Param year: The year of the symmulation (integer)
        ! Param varid: Variable ID (a, bedrock, ...) the values belong to
        ! Param values: Array (lon, lat) with values as REAL(kind=8)
        ! param nlats: Number of latititudes (integer)
        ! param nlons: Number of longitued (integer)
        !
        ! Flushs the data to the disk afterwards
        use netcdf
        implicit none
        ! Parameters
        integer, intent(in) :: ncid
        integer, intent(in) :: year
        integer, intent(in) :: varid
        integer, intent(in) :: nlats
        integer, intent(in) :: nlons
        integer, intent(in) :: N_LAYER
        real(kind=4), dimension(nlons,nlats,N_LAYER), intent(in) :: values

        ! Error handling.
        integer :: retval

        ! We are writing 3D data with time, We will need 4 netCDF dimensions. (Lats, Long, z, Time)
        integer, parameter :: NDIMS = 4

        ! The start and count arrays will tell the netCDF library where to
        ! write our data.
        integer start(NDIMS), count(NDIMS)


        ! These settings tell netcdf to write one timestep of data. (The
        ! setting of start(3) inside the loop below tells netCDF which
        ! timestep to write.)
        count(1) = NLONS
        count(2) = NLATS
        count(3) = N_LAYER
        count(4) = 1
        start(1) = 1
        start(2) = 1
        start(3) = 1
        start(4) = year
        
        
        ! Write the pretend data. This will write the data.
        ! The arrays only hold one timestep worth of data.
        retval = nf_put_vara_real(ncid, varid, start, count, values)
        !retval = nf_put_vara_double(ncid, varid, start, count, values)
        if (retval .ne. nf_noerr) then
            print *, 'Could not write 3D variable to NetCDF: ', varid
            call handle_err(retval)
        end if

        ! Flush data to the disk
        retval = NF_SYNC(ncid)
        if (retval .ne. nf_noerr) then
            print *, 'Could not sync NetCDF'
            call handle_err(retval)
        end if

    end subroutine writeNCDFSNOW3DValues
    
    subroutine writeNCDFSNOW2DValues(ncid, year, varid, values, nlats, nlons,N_LAYER)
        ! Writes the Values (Real(kind=8)) into the netCDF File
        ! Param ncid: File handle of the netCDF File
        ! Param year: The year of the symmulation (integer)
        ! Param varid: Variable ID (a, bedrock, ...) the values belong to
        ! Param values: Array (lon, lat) with values as REAL(kind=8)
        ! param nlats: Number of latititudes (integer)
        ! param nlons: Number of longitued (integer)
        !
        ! Flushs the data to the disk afterwards
        !used for the Annual data
        use netcdf
        implicit none
        ! Parameters
        integer, intent(in) :: ncid
        integer, intent(in) :: year
        integer, intent(in) :: varid
        integer, intent(in) :: nlats
        integer, intent(in) :: nlons
        integer, intent(in) :: N_LAYER
        real(kind=4), dimension(nlons,nlats,N_LAYER), intent(in) :: values

        ! Error handling.
        integer :: retval

        ! We are writing 3D data with time, We will need 4 netCDF dimensions. (Lats, Long, z, Time)
        integer, parameter :: NDIMS = 4

        ! The start and count arrays will tell the netCDF library where to
        ! write our data.
        integer start(NDIMS), count(NDIMS)


        ! These settings tell netcdf to write one timestep of data. (The
        ! setting of start(3) inside the loop below tells netCDF which
        ! timestep to write.)
        count(1) = NLONS
        count(2) = NLATS
        count(3) = 1
        start(1) = 1
        start(2) = 1
        start(3) = year

        ! Write the pretend data. This will write the data.
        ! The arrays only hold one timestep worth of data.
        retval = nf_put_vara_real(ncid, varid, start, count, values)
        !retval = nf_put_vara_double(ncid, varid, start, count, values)
        if (retval .ne. nf_noerr) then
            print *, 'Could not write 3D variable to NetCDF: ', varid
            call handle_err(retval)
        end if

        ! Flush data to the disk
        retval = NF_SYNC(ncid)
        if (retval .ne. nf_noerr) then
            print *, 'Could not sync NetCDF'
            call handle_err(retval)
        end if

    end subroutine writeNCDFSNOW2DValues
    
        subroutine writeNCDFSNOW2D_pointwise(ncid, year, varid, values, iy, ix, N_LAYER)
        ! Writes the Values (Real(kind=8)) into the netCDF File
        ! Param ncid: File handle of the netCDF File
        ! Param year: The year of the symmulation (integer)
        ! Param varid: Variable ID (a, bedrock, ...) the values belong to
        ! Param values: Array (lon, lat) with values as REAL(kind=8)
        ! param nlats: Number of latititudes (integer)
        ! param nlons: Number of longitued (integer)
        !
        ! Flushs the data to the disk afterwards
        ! This routine saves the data pointwise 
        use netcdf
        implicit none
        ! Parameters
        integer, intent(in) :: ncid
        integer, intent(in) :: year
        integer, intent(in) :: varid
        integer, intent(in) :: iy
        integer, intent(in) :: ix
        integer, intent(in) :: N_LAYER
        real(kind=4), intent(in) :: values

        ! Error handling.
        integer :: retval

        ! We are writing 3D data with time, We will need 4 netCDF dimensions. (Lats, Long, z, Time)
        integer, parameter :: NDIMS = 4

        ! The start and count arrays will tell the netCDF library where to
        ! write our data.
        integer start(NDIMS), count(NDIMS)

!         print*,values
        ! These settings tell netcdf to write one timestep of data. (The
        ! setting of start(3) inside the loop below tells netCDF which
        ! timestep to write.)
        count(1) = 1
        count(2) = 1
        count(3) = 1
        start(1) = ix
        start(2) = iy
        start(3) = year

        ! Write the pretend data. This will write the data.
        ! The arrays only hold one timestep worth of data.
        retval = nf_put_vara_real(ncid, varid, start, count, values)
        !retval = nf_put_vara_double(ncid, varid, start, count, values)
        if (retval .ne. nf_noerr) then
            print *, 'Could not write 3D variable to NetCDF: ', varid
            call handle_err(retval)
        end if

        ! Flush data to the disk
        retval = NF_SYNC(ncid)
        if (retval .ne. nf_noerr) then
            print *, 'Could not sync NetCDF'
            call handle_err(retval)
        end if

    end subroutine writeNCDFSNOW2D_pointwise
    
    
subroutine writeNCDFSNOW2D_pointwise_t_xy(ncid, starttime, varid, values, iy, ix, endtime)
        ! Writes the Values (Real(kind=8)) into the netCDF File
        ! Param ncid: File handle of the netCDF File
        ! Param year: The year of the symmulation (integer)
        ! Param varid: Variable ID (a, bedrock, ...) the values belong to
        ! Param values: Array (lon, lat) with values as REAL(kind=8)
        ! param nlats: Number of latititudes (integer)
        ! param nlons: Number of longitued (integer)
        !
        ! Flushs the data to the disk afterwards
        use netcdf
        implicit none
        ! Parameters
        integer, intent(in) :: ncid
        integer, intent(in) :: starttime
        integer, intent(in) :: endtime
        integer, intent(in) :: varid
        integer, intent(in) :: iy
        integer, intent(in) :: ix
        real(kind=4), dimension(endtime-starttime), intent(in) :: values

        ! Error handling.
        integer :: retval

        ! We are writing 3D data with time, We will need 4 netCDF dimensions. (Lats, Long, z, Time)
        integer, parameter :: NDIMS = 4
        
        !local variable timestep
        integer :: timestep

        ! The start and count arrays will tell the netCDF library where to
        ! write our data.
        integer start(NDIMS), count(NDIMS)

        timestep=endtime-starttime+1
!         print*,values
        ! These settings tell netcdf to write one timestep of data. (The
        ! setting of start(3) inside the loop below tells netCDF which
        ! timestep to write.)
        count(1) = timestep
        count(2) = 1
        count(3) = 1
        start(1) = starttime
        start(2) = ix
        start(3) = iy

        ! Write the pretend data. This will write the data.
        ! The arrays only hold one timestep worth of data.
        retval = nf_put_vara_real(ncid, varid, start, count, values)
        !retval = nf_put_vara_double(ncid, varid, start, count, values)
        if (retval .ne. nf_noerr) then
            print *, 'Could not write 3D variable to NetCDF: ', ix, iy
            call handle_err(retval)
        end if

        ! Flush data to the disk
        retval = NF_SYNC(ncid)
        if (retval .ne. nf_noerr) then
            print *, 'Could not sync NetCDF'
            call handle_err(retval)
        end if

    end subroutine writeNCDFSNOW2D_pointwise_t_xy
    
    subroutine writeNCDFSNOW3D_pointwise_t_zxy(ncid, starttime, varid, values, iy, ix, N_LAYER, endtime)
        ! Writes the Values (Real(kind=8)) into the netCDF File
        ! Param ncid: File handle of the netCDF File
        ! Param year: The year of the symmulation (integer)
        ! Param varid: Variable ID (a, bedrock, ...) the values belong to
        ! Param values: Array (lon, lat) with values as REAL(kind=8)
        ! param nlats: Number of latititudes (integer)
        ! param nlons: Number of longitued (integer)
        !
        ! Flushs the data to the disk afterwards
        ! used for the daily and monthly data output, reorder of variables reduces the writing time significantly
        use netcdf
        implicit none
        ! Parameters
        integer, intent(in) :: ncid
        integer, intent(in) :: starttime
        integer, intent(in) :: endtime
        integer, intent(in) :: varid
        integer, intent(in) :: iy
        integer, intent(in) :: ix
        integer, intent(in) :: N_LAYER
        real(kind=4), dimension(365,N_LAYER), intent(in) :: values

        ! Error handling.
        integer :: retval

        ! We are writing 3D data with time, We will need 4 netCDF dimensions. (Lats, Long, z, Time)
        integer, parameter :: NDIMS = 4
        
        !local variable timestep
        integer :: timestep

        ! The start and count arrays will tell the netCDF library where to
        ! write our data.
        integer start(NDIMS), count(NDIMS)

        timestep=endtime-starttime+1
!         print*,values
        ! These settings tell netcdf to write one timestep of data. (The
        ! setting of start(3) inside the loop below tells netCDF which
        ! timestep to write.)
        count(1) = timestep
        count(2) = N_LAYER
        count(3) = 1
        count(4) = 1
        start(1) = starttime
        start(2) = 1
        start(3) = ix
        start(4) = iy

        ! Write the pretend data. This will write the data.
        ! The arrays only hold one timestep worth of data.
        retval = nf_put_vara_real(ncid, varid, start, count, values)
        !retval = nf_put_vara_double(ncid, varid, start, count, values)
        if (retval .ne. nf_noerr) then
            print *, 'Could not write 3D variable to NetCDF: ', ix, iy
            call handle_err(retval)
        end if

        ! Flush data to the disk
        retval = NF_SYNC(ncid)
        if (retval .ne. nf_noerr) then
            print *, 'Could not sync NetCDF'
            call handle_err(retval)
        end if

    end subroutine writeNCDFSNOW3D_pointwise_t_zxy

    subroutine writeNCDFSNOW3B_pointwise_time_interval(ncid, year, varid, values, iy, ix, N_LAYER, endtime)
        ! Writes the Values (Real(kind=8)) into the netCDF File
        ! Param ncid: File handle of the netCDF File
        ! Param year: The year of the symmulation (integer)
        ! Param varid: Variable ID (a, bedrock, ...) the values belong to
        ! Param values: Array (lon, lat) with values as REAL(kind=8)
        ! param nlats: Number of latititudes (integer)
        ! param nlons: Number of longitued (integer)
        !
        ! Flushs the data to the disk afterwards
        use netcdf
        implicit none
        ! Parameters
        integer, intent(in) :: ncid
        integer, intent(in) :: year
        integer, intent(in) :: endtime
        integer, intent(in) :: varid
        integer, intent(in) :: iy
        integer, intent(in) :: ix
        integer, intent(in) :: N_LAYER
        real(kind=4), dimension(N_LAYER,365), intent(in) :: values

        ! Error handling.
        integer :: retval

        ! We are writing 3D data with time, We will need 4 netCDF dimensions. (Lats, Long, z, Time)
        integer, parameter :: NDIMS = 4
        
        !local variable timestep
        integer :: timestep

        ! The start and count arrays will tell the netCDF library where to
        ! write our data.
        integer start(NDIMS), count(NDIMS)

        timestep=endtime-year+1
!         print*,values
        ! These settings tell netcdf to write one timestep of data. (The
        ! setting of start(3) inside the loop below tells netCDF which
        ! timestep to write.)
        count(1) = 1
        count(2) = 1
        count(3) = N_LAYER
        count(4) = endtime
        start(1) = ix
        start(2) = iy
        start(3) = 1
        start(4) = 1

        ! Write the pretend data. This will write the data.
        ! The arrays only hold one timestep worth of data.
        retval = nf_put_vara_real(ncid, varid, start, count, values)
        !retval = nf_put_vara_double(ncid, varid, start, count, values)
        if (retval .ne. nf_noerr) then
            print *, 'Could not write 3D variable to NetCDF: ', varid
            call handle_err(retval)
        end if

        ! Flush data to the disk
        retval = NF_SYNC(ncid)
        if (retval .ne. nf_noerr) then
            print *, 'Could not sync NetCDF'
            call handle_err(retval)
        end if

    end subroutine writeNCDFSNOW3B_pointwise_time_interval
    
    subroutine writeNCDF3DGridValues(ncid, year, varid, values, nlats, nlons, N_LAYER)
        ! Writes the Values (Real(kind=8)) into the netCDF File
        ! Param ncid: File handle of the netCDF File
        ! Param year: The year of the symmulation (integer)
        ! Param varid: Variable ID (a, bedrock, ...) the values belong to
        ! Param values: Array (lon, lat) with values as REAL(kind=8)
        ! param nlats: Number of latititudes (integer)
        ! param nlons: Number of longitued (integer)
        !
        ! Flushs the data to the disk afterwards
        use netcdf
        implicit none
        ! Parameters
        integer, intent(in) :: ncid
        integer, intent(in) :: year
        integer, intent(in) :: varid
        integer, intent(in) :: nlats
        integer, intent(in) :: nlons
        integer, intent(in) :: N_LAYER
        real(kind=4), dimension(nlons,nlats,N_LAYER), intent(in) :: values

        ! Error handling.
        integer :: retval

        ! We are writing 3D data with time, We will need 4 netCDF dimensions. (Lats, Long, z, Time)
        integer, parameter :: NDIMS = 4

        ! The start and count arrays will tell the netCDF library where to
        ! write our data.
        integer start(NDIMS), count(NDIMS)


        ! These settings tell netcdf to write one timestep of data. (The
        ! setting of start(3) inside the loop below tells netCDF which
        ! timestep to write.)
        count(1) = NLONS
        count(2) = NLATS
        count(3) = N_LAYER
        count(4) = 1
        start(1) = 1
        start(2) = 1
        start(3) = 1
        start(4) = year

        ! Write the pretend data. This will write the data.
        ! The arrays only hold one timestep worth of data.
        retval = nf_put_vara_real(ncid, varid, start, count, values)
        !retval = nf_put_vara_double(ncid, varid, start, count, values)
        if (retval .ne. nf_noerr) then
            print *, 'Could not write 3D variable to NetCDF: ', varid
            call handle_err(retval)
        end if

        ! Flush data to the disk
        retval = NF_SYNC(ncid)
        if (retval .ne. nf_noerr) then
            print *, 'Could not sync NetCDF'
            call handle_err(retval)
        end if

    end subroutine writeNCDF3DGridValues
    
    function read_climate_long(filename, NLONS, NLATS, variable)
        ! Reads the temperature or precipitation data (8bit real) for 365 days out of the given NetCDF File.
        ! variable: prec/temp, units: m/yr and kelvin
        ! Returns 3D Array: Lon, Lat, Day

        ! Reads
        use netcdf
        use variables ! Own Module
        implicit none
        ! input
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: variable
        integer, intent(in) :: NLATS
        integer, intent(in) :: NLONS

        ! Zeitschritte
        ! ndays

        ! output
        ! TODO: Zeitschritt einbinden
        real(kind=8) :: read_climate_long(NLONS, NLATS, 365)


        ! Return value, to check if everything was ok
        integer :: retval

        ! This will be the netCDF ID for the file and data variable.
        integer :: ncid, varid
        
        print*,variable
        print*,filename
        ! Open the file, Read only
        ! Trim file path: http://stackoverflow.com/questions/15093712/trimming-string-for-directory-path
        retval = nf_open(TRIM(adjustl(filename)), nf_NOWRITE, ncid)

        ! Read Data
        ! Get the varid of the data variable, based on its name.
        retval = nf_inq_varid(ncid, TRIM(adjustl(variable)), varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! read the data
        retval = nf_get_var_double(ncid, varid, read_climate_long)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! Close the file, freeing all resources.
        retval = nf_close(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        if(debug > 0) then
            ! If we got this far, everything worked as expected. Yipee!
            print *,"*** Read ", TRIM(adjustl(variable))," from NetCDF file: ", TRIM(adjustl(filename)), ""
        end if

        return

    end function read_climate_long
    
    function read_snow_data(filename, NLONS, NLATS, NLAYER, variable)
        ! Reads an firn initialization file from an output file of this model for restarts etc.
        ! used in case of restart = .true.
        ! Reads
        use netcdf
        use variables ! Own Module
        implicit none
        ! input
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: variable
        integer, intent(in) :: NLATS
        integer, intent(in) :: NLONS
        integer, intent(in) :: NLAYER

        ! Zeitschritte
        ! ndays

        ! output
        ! TODO: Zeitschritt einbinden
        real(kind=8) :: read_snow_data(NLONS, NLATS, NLAYER)


        ! Return value, to check if everything was ok
        integer :: retval

        ! This will be the netCDF ID for the file and data variable.
        integer :: ncid, varid
        
        
        !print*,filename
        retval = nf_open(TRIM(adjustl(filename)), nf_NOWRITE, ncid)

        ! Read Data
        ! Get the varid of the data variable, based on its name.
        retval = nf_inq_varid(ncid, TRIM(adjustl(variable)), varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! read the data
        retval = nf_get_var_double(ncid, varid, read_snow_data)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
        ! Close the file, freeing all resources.
        retval = nf_close(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        if(debug > 0) then
            ! If we got this far, everything worked as expected. Yipee!
            print *,"*** Read ", TRIM(adjustl(variable))," as initial profile from NetCDF file: ", TRIM(adjustl(filename)), ""
        end if

        return

    end function read_snow_data

        subroutine read_climate_once(filename, NLONS, NLATS, variable1, variable2, variable3, &
        variable4, variable5, variable6, read_climate_long_1, read_climate_long_2, &
        read_climate_long_3, read_climate_long_4, read_climate_long_5,read_climate_long_6, rate)
        ! Reads the temperature or precipitation data (8bit real) for 365 days out of the given NetCDF File.
        ! variable: prec/temp, units: m/yr and kelvin
        ! Returns 3D Array: Lon, Lat, Day

        ! Reads
        use netcdf
        use variables ! Own Module
        implicit none
        ! input
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: variable1
        character(len=*), intent(in) :: variable2
        character(len=*), intent(in) :: variable3
        character(len=*), intent(in) :: variable4
        character(len=*), intent(in) :: variable5
        character(len=*), intent(in) :: variable6
        integer, intent(in) :: NLATS
        integer, intent(in) :: NLONS

        ! Zeitschritte
        ! ndays

        ! output
        ! TODO: Zeitschritt einbinden
        real(kind=8) :: read_climate_long_1(NLONS, NLATS, 365)
        real(kind=8) :: read_climate_long_2(NLONS, NLATS, 365)
        real(kind=8) :: read_climate_long_3(NLONS, NLATS, 365)
        real(kind=8) :: read_climate_long_4(NLONS, NLATS, 365)
        real(kind=8) :: read_climate_long_5(NLONS, NLATS, 365)
        real(kind=8) :: read_climate_long_6(NLONS, NLATS, 365)
        ! Return value, to check if everything was ok
        integer :: retval
        real :: rate

        ! This will be the netCDF ID for the file and data variable.
        integer :: ncid, varid, c3, c4, csw1,csw2
        
        CALL SYSTEM_CLOCK(c3)
        CALL SYSTEM_CLOCK(csw1)
        ! Open the file, Read only
        ! Trim file path: http://stackoverflow.com/questions/15093712/trimming-string-for-directory-path
        retval = nf_open(TRIM(adjustl(filename)), nf_NOWRITE, ncid)
        CALL SYSTEM_CLOCK(csw2)
        if(debug > 1) then
        print*,'opening', (csw2-csw1)/rate
        end if
        CALL SYSTEM_CLOCK(csw1)
        ! Read Data
        ! Get the varid of the data variable, based on its name.
        !print*,variable1
        !print*,filename
        retval = nf_inq_varid(ncid, TRIM(adjustl(variable1)), varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! read the data
        retval = nf_get_var_double(ncid, varid, read_climate_long_1)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
                if(debug > 0) then
            ! If we got this far, everything worked as expected. Yipee!
            print *,"*** Read ", TRIM(adjustl(variable1))," from NetCDF file: ", TRIM(adjustl(filename)), ""
        end if
                        CALL SYSTEM_CLOCK(csw2)
                        if(debug > 1) then
        print*,TRIM(adjustl(variable1)), (csw2-csw1)/rate
        end if
         CALL SYSTEM_CLOCK(csw1)
        !print*,variable2
        !print*,filename
        retval = nf_inq_varid(ncid, TRIM(adjustl(variable2)), varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! read the data
        retval = nf_get_var_double(ncid, varid, read_climate_long_2)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
                if(debug > 0) then
            ! If we got this far, everything worked as expected. Yipee!
            print *,"*** Read ", TRIM(adjustl(variable2))," from NetCDF file: ", TRIM(adjustl(filename)), ""
        end if
                        CALL SYSTEM_CLOCK(csw2)
                        if(debug > 1) then
        print*,TRIM(adjustl(variable2)), (csw2-csw1)/rate
        end if
         CALL SYSTEM_CLOCK(csw1)
        !print*,variable3
        !print*,filename
        retval = nf_inq_varid(ncid, TRIM(adjustl(variable3)), varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! read the data
        retval = nf_get_var_double(ncid, varid, read_climate_long_3)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
                if(debug > 0) then
            ! If we got this far, everything worked as expected. Yipee!
            print *,"*** Read ", TRIM(adjustl(variable3))," from NetCDF file: ", TRIM(adjustl(filename)), ""
        end if
                        CALL SYSTEM_CLOCK(csw2)
        if(debug > 1) then
        print*,TRIM(adjustl(variable3)), (csw2-csw1)/rate
        end if
         CALL SYSTEM_CLOCK(csw1)
        !print*,variable4
        !print*,filename
        retval = nf_inq_varid(ncid, TRIM(adjustl(variable4)), varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! read the data
        retval = nf_get_var_double(ncid, varid, read_climate_long_4)
        if (retval .ne. nf_noerr) call handle_err(retval)

        
                if(debug > 0) then
            ! If we got this far, everything worked as expected. Yipee!
            print *,"*** Read ", TRIM(adjustl(variable4))," from NetCDF file: ", TRIM(adjustl(filename)), ""
        end if
                CALL SYSTEM_CLOCK(csw2)
                if(debug > 1) then
        print*,TRIM(adjustl(variable4)), (csw2-csw1)/rate
        end if
        CALL SYSTEM_CLOCK(csw1)
!         print*,variable5
!         print*,filename
        retval = nf_inq_varid(ncid, TRIM(adjustl(variable5)), varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! read the data
        retval = nf_get_var_double(ncid, varid, read_climate_long_5)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
                if(debug > 0) then
            ! If we got this far, everything worked as expected. Yipee!
            print *,"*** Read ", TRIM(adjustl(variable5))," from NetCDF file: ", TRIM(adjustl(filename)), ""
        end if
                CALL SYSTEM_CLOCK(csw2)
                if(debug > 1) then
        print*,TRIM(adjustl(variable5)), (csw2-csw1)/rate
        end if
        CALL SYSTEM_CLOCK(csw1)
!         print*,variable6
!         print*,filename
        retval = nf_inq_varid(ncid, TRIM(adjustl(variable6)), varid)
        if (retval .ne. nf_noerr) call handle_err(retval)

        ! read the data
        retval = nf_get_var_double(ncid, varid, read_climate_long_6)
        if (retval .ne. nf_noerr) call handle_err(retval)
        
                if(debug > 0) then
            ! If we got this far, everything worked as expected. Yipee!
            print *,"*** Read ", TRIM(adjustl(variable6))," from NetCDF file: ", TRIM(adjustl(filename)), ""
        end if
        CALL SYSTEM_CLOCK(csw2)
        if(debug > 1) then
        print*,TRIM(adjustl(variable6)), (csw2-csw1)/rate
        end if

        CALL SYSTEM_CLOCK(csw1)
        ! Close the file, freeing all resources.
        retval = nf_close(ncid)
        if (retval .ne. nf_noerr) call handle_err(retval)
        CALL SYSTEM_CLOCK(csw2)
        if(debug > 1) then
        print*,'closing', (csw2-csw1)/rate
        end if
        
        CALL SYSTEM_CLOCK(c4)
        if(debug > 0) then
        print*,'*** Reading routine time duration', (c4-c3)/rate
        end if
        return

    end subroutine read_climate_once
    ! --------------------------------------------------------------
    ! SAVE SOME PARAMETERS TO TXT FILE
    ! --------------------------------------------------------------
    subroutine save_parameters(path,date,description)

	use variables
	!use variables_snow
    implicit none
	character(len=*), intent(in) :: path
	character(len=*), intent(in) :: date
	character(len=*), intent(in) :: description

	open (1337, file = trim(adjustl(path)) // 'values.txt', form='formatted')

	! writing to the file
	write(1337,*) "Date: ",  trim(adjustl(date))
	write(1337,*) "Experiment description: ", trim(adjustl(description))
	write(1337,*)
	write(1337,*) "Logicals used in this simulation:"
	write(1337,*)
	!	write(1337,"(' = ',L1 )")  
	write(1337,"('write_netcdf = ',L1 )") write_netcdf
	write(1337,"('eraiterim_climate averaged = ',L1 )")  eraiterim_climate
	write(1337,"('eraiterim_climate transient = ',L1 )")  erai_backandforth_climate
	write(1337,*)'Input climate transient path ', netcdf_input_eraiterim_directory
	write(1337,*)'Input climate averaged', netcdf_input_eraiterim_climate
	write(1337,*)'ERAi_period = ', erai_year_begin, '-',  erai_year_turn
	write(1337,*)'Resolution =' , dx, '*', dy
	write(1337,*)'Gridsize =', L, '*', W
	write(1337,*)'Gridpoints =', nx, '*', ny
! 	write(1337,"('initial_gis = ',L1 )") initial_gis
	write(1337,*)'Chose massbalance model, 1= Positive degree day, 2 = Energyflux mass balance, 3 = Troll model ', smb_model    
	write(1337,*)'Albedo_module, 0 = albedo input file, 1 = constant, 2 = harmonic, 3 = oerlemans, 4 = Aoki 5 Bougamont', Albedo_module
	if (albedo_module==0) then
	write(1337,*)'Albedo_file_type=',albedo_file_variable_name
	write(1337,*)'Albedo',albedo_file_path 
	end if
	write(1337,*)'latent_heat_flux_on=', latent_heat_flux_on
    write(1337,*)'D_lf analog to D_sf:', latent_heat_flux_analog_to_sensible_heat_flux
	write(1337,*)
	write(1337,*) "Output parameters of simulation:"
	write(1337,*)
	!	write(1337,*)' = ', 
	write(1337,*)'maxyears = ', maxyears
	write(1337,*)'write to netcdf every x years = ', netcdf_timesteps
	write(1337,*)'max number of time entries in one file = ', yearly_netcdf_file_freq
	write(1337,*)'Which data is written = ','annual:',annual_data, 'monthly:',monthly_data, 'daily:', daily_data
	write(1337,*)'write annual data ever x years = ', annual_data_frequency
	write(1337,*)'write daily data ever x years = ', daily_data_frequency
	write(1337,*)'write monthly data ever x years = ', monthly_data_frequency
    write(1337,*)'write detailed data = ', daily_data_start, '-', daily_data_end
    write(1337,*)
	write(1337,*)'n_snowlayer = ', n_snowlayer
	write(1337,*)'soll_mass = ', soll_mass
	write(1337,*)'lower_massbound = ', lower_massbound
	write(1337,*)'upper_massbound = ', upper_massbound
	write(1337,*)
	write(1337,*) "Physical constants used in this simulation:"
	write(1337,*)
	write(1337,*)'albedo_snow_new = ', albedo_snow_new
	write(1337,*)'albedo_snow_wet = ', albedo_snow_wet
	write(1337,*)'albedo_ice = ', albedo_ice
	write(1337,*)'D_sf = ', D_sf
	write(1337,*)'ratio = ', ratio
    write(1337,*)'effective D_lf = ', ratio*D_sf/1003.*0.622*(L_v+L_lh)/10000.
	write(1337,*)'eps_air = ', eps_air
	write(1337,*)'latent_heat_flux = ', latent_heat_flux_on
	write(1337,*)'albedo_module = ', albedo_module
	write(1337,*)'max_lwc = ', max_lwc
	write(1337,*)'A_Eismint = ', A_Eismint
	write(1337,*)'n_Eismint = ', n_Eismint
	write(1337,*)
	write(1337,*) "Input_files"
	write(1337,*)
	write(1337,*)'Grid', netcdf_input_bedrock
	



	! we are done with writing the file
	close(1337)
	print*,"Parameters saved at: ", trim(adjustl(path)) // 'spam.txt'
    end subroutine save_parameters

end module io
