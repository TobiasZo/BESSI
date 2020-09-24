! 2D ice sheet model in fortran. Two-dimensional flow line ice sheet model using Glen's law
!
!
! Initial condition for first EISMINT Tests:
! - fixed margin, fixed accumulation rate (0.3m/yr)
! http://homepages.vub.ac.be/~phuybrec/eismint.html
!
! Some (not all!) references:
!------------------------------------------------------
! Oerlemans, J. (1981), Some basic experiments with a vertically-integrated ice sheet model, Tellus 33, 1-11
! Oerlemans, J. (1982), Glacial Cycles and Ice-Sheet Modelling, Climatic Change 4, 353-374
! Huybrechts et al. (1996), The EISMINT benchmarks for testing ice-sheet models, Annals of Glaciology, 23, 1-12
!
! Author:       neff@climate.unibe.ch, born@climate.unibe.ch
! Date:         2015/05
! Revision:     1.1
! SVN Info:     $Id: $
!

! 
!   Co-author: Imhof Michael
!   Co-Developer: Imhof Michael
!   Mail: imhof@vaw.baug.ethz.ch or imhof@climate.unibe.ch
!   parts changed or added by Michael:
!   - Initial ice for PD and LGM (later one unstable)
!   - Several different climate forcings (probably no more working)
!   - Loading and using Bern3D climate data for -136,050 BP until 1950 BC, Bern3D climate is applied as a deviation relative to 1950. Needs an update


! Possible Improvements
!----------------------
! - NetCDF Output: New Output file every 100'000 years, otherwise the NetCDF File gets huge and is difficult to handle.
! - Dynamic watermask: Do not take the initial bedrock, use the current one but forget the ice above it (and it should still increase, even if it is below the water). These will get major Problems in Northamerica, that there will not grow a sea.
!
! Problems
!----------
! - "Todays" greenland as initial state does not work with the surface mass balance.
! - Some years are not written into the NetCDF file, not yet clear why (fixed with real in "if condition"?).
! - Problems on Ubelix: Last year is not written to the NetCDF: Error: NetCDF: Variable not found
! - Integrated Mass Balance: Isostatic melt is not yet considered correctly
! - Heartbeat calculations with the simulated years per hour is probably wrong.
!

! Version 3.0.0a (2016/09)
! - Ice can be allowed to flow into the ocean until a certain depth under the current sealevel
! - 



! Version 2.0.0a (2015/12) Testing routines from Basil copied but not tested
!-------------------------
! - First Version after submission, Version used for paper review
! - BUG Fixed: Calculations in the Flux are now correct. Before it was possible, that the model gains or loses ice without physical background.
! - Possible to store Integrated Mass Balance (IMB) information in external CSV.
! - Problems with the assignment_mask and the changing sea level fixed
! - ELRA Bedrock Model: First version of the ELRA bedrock model. Very slow and not big differences to the local version.
! - Hysteresis mode with sharp temperature change introduced (hysteresis jump)
!
! Version 1.2 (2015/5) Michael Imhof
!----------------------
! - The surfacemassbalance has been outsourced in modules to allow different models for the surface
! - density of accumulation adjusted in pdd
! - ablation adjusted in pdd with densite (beta uses density of water) eventually can be deleted if beta is retuned for adjusted accumulation
! - sealeveladjustment adjusted too with densities
! - snowmass of the firn included in sealevel adjustment
!
! Version 1.1 (2015/3)
!----------------------
! - Clean source code: Remove unused functions (especialy SMB functions), compact modules.
! - Surface Mass Balance is fixed to a temperature feedback at every 10 years (smb_type = 4).
! - Flow Type fixed to EISMINT (flow_type = 2)
! - Introduced hysteresis mode
! - BUG [removed]: Ice density was not converted to ice, instead the density of water was used!
!
! Version 1.0a (2014/2)
!----------------------
! - First Version with sea level adjustments which runs trough. Version 1.0a is here!
! - BUG [from version before,solved]: If the sea level increases again, grid points with a low bedrock are converted back to water, even if there is ice on it.
!       This leads to unstabilities. The idea is to check if the available ice is able to displace the water column, in this case the ice will stay as it is.
! - Restore initial state prepared to handle the sea level.
!
! Version 0.82 (2014/2)
!----------------------
! - Problem with isolated islands: If an island is surounded by water, the hgradient gets 0 which leads to FN,FS, ... to be zero. It accumulates forever. In the case of a small hgradient, the elevation of the bedrock, instead of the sea level is taken.
! - BUG [solved]: Set the height of the bedrock on water to the sea level instead of 0 and adjust it during runtime.
! - Remove the reset of the ice thickness on the isolated islands. The simulation runs now stable with the two adjustments and now artificial reset is needed.
! - Isolated waterpoints are treated as land, if they are sourounded by land. Cause this could lead to problems.
! - Limit the ice flux to the available ice (introduced in v0.8, now in production, cause it makes a difference and does not break)
! - BUG [not solved]: If the sea level increases again, grid points with a low bedrock are converted back to water, even if there is ice on it. This leads to unstabilities. The idea is to check if the available ice is able to displace the water column.
! - Increased the length of the folder name.
! - There are still some DEBUG lines in the code, which are commented. They will be removed later, when the problem with the increasing sea level is solved.
!
! Version 0.81
!-------------
! - Ignore the accumulation in the himalaya region, since the glacier regime is totaly different in this area.
! - Set the great lakes as land and not as water in the assignment mask.
! - Dynamic landmask, ocean surface depends on the amount of ice on land, calculated every 50 years.
! - Set the ice_thickness on the isolated land points to a maximum value. (removed later)
! - Ready for multiple ensemble runs, the tuning parameters are all in the variables.f90 file (ensemble.sh in separate directory)
! - Introduced the ice_flux_adjustments parameter, to adjust the ice flux for the ensemble runs.
!
! Version 0.8
!------------
! - Elevation of the given climate is independent of the bedrock (etopo) elevation
! - Ready for CCSM input data, 20km and 40km
! - Writes all NetCDF variables out, if an unstable integration occurs.
! - Temperature feedback bug for the first 500 years. Solved.
! - Unstable integration with more detailed information.
! - SOLVED: If timestep is below one year, the netcdf and information is written several times during the netcdf_timesteps.
! - Possible to limit the flux to the maximum amount of available ice (commented out, cause it is not used anymore, also stable without this).
! - Accumulation and Ablation are retreived independent from each other and written to the netcdf.
! - Rename all things belonging to accumluation_type to smb_type.
!
! Version 0.72
!-------------
! - More detailed model description in the folder name
! - Renamed accumulation in the netcdf to smb (surface mass balance)
! - Add flag for the precipitaion units in the input netcdf file (m/s or m/yr)
!
! Version 0.71
!-------------
! - Ready for the northern hemisphere. Memory Model: -mcmodel=medium
! - H_ts() is now in the units of m/m^2 instead of the m/grid. Makes it easier to do integration checks (independent of the grid size).
! - BUG: NetCDF was written in the x/y loop. Moved it outside, now it is much faster.
! - Performace improved: Surface mass balanced calculated with array functions (where function) instead of two loops. Improved performance by 1/4.
! - Elevation feedback only all 10 timestep instead of every timestep. Improved performance.
!
! Version 0.7
!------------
! - Calculate the flux only where the height of the ice is greater than 0: For this, the Watermask is extended to an assignment_mask. mask where the flux should be calculated (if thickness = 0, no calculations need to be done)
!
! Version 0.666
!--------------
! - Changed it for a 64bit system
! - Uses now 8byte real values, instead of 4byte.
! - Problems with long directory names solved, when copying files to the output directory.
! - Removed CSV functions and libraries
! - Included NetCDF Libraries from the System instead of own ones
!
! Version 0.66
!-------------
! - Variables in an own fortran file. The variables are stored with the output, the get a fast overview of the initial conditions.
! - Own outputdirectory for every run with the execution date as directory name. Variables and netcdf input files are stored into this directory.
! - Temperature feedback with elevation as accumulation type (smb_type = 4).
! - Read initial state from an old netcdf file.
! - Debug mode with different levels.
!
! Version 0.6
!-------------
! - BUG fixed: Bedrock sinking is now relative to the ice height and noth to the elevation above sea level.
! - Changed variable names, to get a better understanding what is stored in the variable.
! - Get external climate conditions and computer surface mass balance out of it. In this version, there is no temperature feedback for the increasing elevation. Will be realized in the next version.
! - Debug mode introduced
! - Runtime is printed at the end.
!
!
! Version 0.55
!-------------
! - Load Bedrock and Watermask from NetCDF File (io_read.f90) (Watermask is created, from all land values below 0)
! - Ice can't grow on the water, therefore Height is set to 0 where the watermask has 1
! - Landmask can have negative values, which are written out. But for internal calculations the bedrock on water is set to 0 (for ice, it doesn't matter if there bedrock is deep below the seasurface).
! - Change the dimension, to fullfill the input from the NetCDF file (Greenland is turned in this case)
!

program IceModel

    !=========================
    ! Include Own Modules
    !=========================
    use variables       ! Module with all variables
    !use variables_snow
    use io              ! Own module to read the values
    !use smb_pdd         ! Module for positive degree day surface mass balance (Basil Neff)
    use smb_emb         ! Module for energy mass balance  (Michael Imhof)
    !use smb_troll       ! Module troll model (Michael Imhof)

    !use OMP_LIB
    INTEGER :: c1,c2,cr,cm,c_start,c_end, reorder_year
    REAL :: rate
        print *,"Spam"
    ! output directory
    if (store_input_netcdf .or. write_netcdf .or. annual_data .or. daily_data .or. monthly_data) then
        call init_output_directory(output_directory, TRIM(adjustl(experiment_description)))
    endif



    !=========================
    ! initialze Variable
    !=========================
    if(debug > 0 ) then
        print *, "Initialize Variables from external files"
        print *, "----------------------------------------"
    end if
if (erai_reorder_climate) then
    print*,"reading order vector"
    OPEN(UNIT=11, FILE=year_vector_file)
    do kk=1,maxyears_spam 
	 read(11,*) erai_vector(kk)
        if( iostat < 0 )then
        print*,'Warning: Year order file shorter than simulation period'
            print*,kk
            exit
        else if( iostat > 0 )then
            print*,'Error: error reading file'
            stop
        end if
    end do
!     do kk=1,234,1
!         erai_vector(kk)=1979
!     end do
!     do kk=10,50,1
!         erai_vector(kk)=2011
!     end do
!     erai_vector=(/ /)
    print*,erai_vector
end if

    ! for parallelisazion
    !result = system(TRIM(adjustl(set_core_command)))
    !result = system('echo $OMP_NUM_THREADS')
!$OMP PARALLEL
	print *,"Hello"
!$OMP END PARALLEL

    CALL system_clock(count_rate=cr)
    CALL system_clock(count_max=cm)
    rate = REAL(cr)
    WRITE(*,*) "system_clock rate ",rate


    ! If the sea level gets not adjusted, set the offset to 0
    if (.not.adjust_sea_level) then
        sea_level_offset = 0.
    end if

    ! Read the relaxed Bedrock from a NetCDF file
    !--------------------------------------------
    if(debug > 1 ) then
        print *, "*** Read '", TRIM(adjustl(netcdf_input_bedrock_variable)) ,"' from file: ", &
                    TRIM(adjustl(netcdf_input_bedrock))
    end if
    ! Bedrock_initial is how the bedrock would be without iceload.
    Bedrock_initial = read_variable(netcdf_input_bedrock, nx, ny, TRIM(adjustl(netcdf_input_bedrock_variable)))

    ! set the initial value for ice free hemisphere
    elevation_netcdf = Bedrock_initial
    bedrock_netcdf = Bedrock_initial


    if(store_input_netcdf) then
        if(debug > 0) then
            print *, "*** Copy bedrock input (", TRIM(adjustl(netcdf_input_bedrock)) , &
                     ") file to the output directory: ", TRIM(adjustl(output_directory))
        end if
        call copy_to_output_directory(netcdf_input_bedrock, output_directory)
    end if


    ! 1 = water, 0 = normal case, 3 = unstable grid points
    ! assignment_mask = read_watermask(Bedrock_initial, nx, ny, (sea_level + sea_level_offset), Bedrock_netcdf, ice_thickness)




    ! load a present day climate icesheet
	if(initial_gis) then
		! set bedrock to state with ice		
		bedrock_netcdf = read_variable(netcdf_input_bedrock, nx, ny, netcdf_input_pd_bedrock_variable) ! TODO remove _special
		
		! set surface to state with ice
		elevation_netcdf = read_variable(netcdf_input_bedrock, nx, ny, netcdf_input_icesurf_variable) ! TODO remove _special
		print*,'*** PD ice sheets loaded'
	end if

    ! load a LGM climate icesheet
	if(initial_lgm_ice) then
		!load LGM topography
		bedrock_netcdf = read_variable(netcdf_input_lgm_ice, nx, ny, initial_bedrock_variable_name )
		elevation_netcdf = read_variable(netcdf_input_lgm_ice, nx, ny, initial_height_variable_name )
		print*,'*** LGM ice sheets loaded'
	end if


	! avoid negative ice thicknesses
	Ice_thickness = elevation_netcdf - bedrock_netcdf
	where(Ice_thickness(:,:).lt. 0.)
		elevation_netcdf(:,:) = bedrock_netcdf(:,:)
	end where	

	!elevation_netcdf = elevation
	!bedrock_netcdf = Bedrock

	Ice_thickness = elevation_netcdf - bedrock_netcdf


	! Restore Sea level
	!sea_level = get_sea_level(Ice_thickness, nx, ny, dx, dy, ocean_area)

	! Assign the assignment_mask
	if(debug > 0) then
		print *, "*** Get assignment_mask from the initial bedrock."
	end if
	assignment_mask = read_watermask(Bedrock_initial, nx, ny, (sea_level + sea_level_offset),&
	 bedrock_netcdf, Ice_thickness) ! 1 = water, 0 = normal case, 3 = unstable grid points
	if(debug > 0) then
		print *, "*** Created assignment_mask from the initial bedrock."
	end if

	! Set the height of the bedrock on the water to the sea level, Otherwise there will be an unstable integration,
	! cause the flux at the coast will get very high

	elevation = elevation_netcdf
	bedrock = bedrock_netcdf

	where(assignment_mask(:,:) .eq. 1)
		Bedrock(:,:) = sea_level + sea_level_offset
		elevation(:,:) = sea_level + sea_level_offset
		! integrated mass balance calculated before
		ice_thickness(:,:) = 0d0
		elevation_netcdf(:,:) = bedrock_netcdf(:,:)
	end where

	! Elevation of ice surface above sea level (over the water: Bedrock and ice_thickness are 0 -> elevation = 0)
	!elevation = Bedrock + ice_thickness

    	! Distance of each grid box from bottom domain boundary
    	y = y * dy
    	! Distance of each grid box from left domain boundary
    	X = x * dx



    ! Read the initial state out of an NetCDF File
    ! For this the watermask from the bedrock is used!
    ! THIS HAS NOT BEEN TESTED!!!!	
    if (read_initial_state) then
        ! Store the input files
        if(store_input_netcdf) then
            call copy_to_output_directory(initial_netcdf_file, output_directory)
        end if

        ! Bedrock_initial is how the bedrock would be without iceload.
        ! already done before, so not needed twice
        ! Bedrock_initial = read_variable(netcdf_input_bedrock, nx, ny, TRIM(adjustl(netcdf_input_bedrock_variable)))

        if(debug > 0 ) then
            print *, '*** Read bedrock state from file: ', TRIM(adjustl(initial_netcdf_file))
        end if
        Bedrock = read_variable(initial_netcdf_file, nx, ny, initial_bedrock_variable_name)

        elevation_netcdf = read_variable(initial_netcdf_file, nx, ny, initial_height_variable_name)
        elevation = elevation_netcdf

        ice_thickness = elevation - Bedrock

        ! Restore Sea level
        sea_level = get_sea_level(ice_thickness, nx, ny, dx, dy, ocean_area)

        ! Assign the assignment_mask
        assignment_mask = read_watermask(Bedrock_initial, nx, ny, (sea_level + sea_level_offset), bedrock_netcdf, ice_thickness) ! 1 = water, 0 = normal case, 3 = unstable grid points

        ! Set the height of the bedrock on the water to the sea level, Otherwise there will be an unstable integration,
        ! cause the flux at the coast will get very high
        do ix=1,nx,1
            do iy=1,ny,1
                if (assignment_mask(ix,iy) == 1) then
                    Bedrock(ix,iy) = sea_level + sea_level_offset
                    ice_thickness(ix,iy) = 0d0
                endif
            end do
        end do
    end if  ! ENDIF: Read initial state




    ! ELRA - Kugelfunktion
    if (active_elra) then
        open(5, iostat=ios, file=TRIM(adjustl(elra_kei_file)), status='old')
        if (ios /= 0) stop ' Error when opening the kei file!'
        do elra_ii=1,1059,1
            read(5,*) (kei(elra_ii,jj), jj=1,2)
        end do
        close(5, status='keep')
        print *, "*** Read kei file from ", TRIM(adjustl(elra_kei_file)), " successful."
    end if



	! READING THE CLIMATE INPUT DATA
	!-------------------------------


	! read climate iterim base data. 
	! potential_temperature and initial_climate_precipitation contain the base climate of eraiterim. 
	! temerature and precipitation can be adaptet to topography and time.

	if((eraiterim_climate).or.(erai_backandforth_climate).or.(erai_reorder_climate))then	
		! Load climate reference elevation (includes ice)
	    	initial_climate_elevation = read_variable(netcdf_input_eraiterim_initial_climate_elevation, nx, ny, &
					     TRIM(adjustl(netcdf_input_eraiterim_initial_climate_elevation_variable)))
					     print*,'here'
	
	end if

	if(eraiterim_climate)then
		inp_temp =   read_climate_long(netcdf_input_eraiterim_climate, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_temp_variable)) )
		inp_precip = read_climate_long(netcdf_input_eraiterim_climate, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_precip_variable)) )
        inp_dewpT = read_climate_long(netcdf_input_eraiterim_climate, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_dewpT_variable)) )
        inp_wind = read_climate_long(netcdf_input_eraiterim_climate, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_wind_variable)) )
        inp_lwrd = read_climate_long(netcdf_input_eraiterim_climate, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_lwrd_variable)) )

        
    
    
		! Remove errorous negative precipitation
		where(inp_precip(:,:,:) .lt. 0.)
			inp_precip(:,:,:) = 0.
		end where


		if(ndays==365) then
			do id=1,ndays,1
				temperature(:,:,id)= inp_temp(:,:,id)+kelvin ! TODO TODO remove test cooling
				precipitation(:,:,id)= inp_precip(:,:,id)/3600./24./1000. ! from mm/day to m/sec
                DewpT(:,:,id)= inp_dewpT(:,:,id)
                wind(:,:,id)= inp_wind(:,:,id)
                if (longwave_downscaling) then
                        lwrd(:,:,id)= inp_lwrd(:,:,id)/temperature(:,:,id)**4*(temperature(:,:,id) +&
                        ((initial_climate_elevation-elevation) * temperature_lapse_rate))**4
                        else 
                        lwrd(:,:,id)= inp_lwrd(:,:,id)
                    end if 
                    if (lwrd_unit==2) then
                        lwrd(:,:,id)= lwrd(:,:,id)/3600./24.
                end if
				! Calculate potential temperature (at sea level)
				potential_temperature(:,:,id) = temperature(:,:,id) + (initial_climate_elevation * temperature_lapse_rate)

			end do

			!print*,'test temp',temperature(placex,placey,3)
			initial_climate_precipitation = precipitation



			! load short wave radiation
			inp_temp =   read_climate_long(netcdf_input_eraiterim_climate, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_swradboa_variable)) )

			! reduce to 96 Timesteps
			do id=1,ndays,1
				P_sun0(:,:,id)= inp_temp(:,:,id)
				P_sun(:,:,id) = P_sun0(:,:,id)
			end do
		end if





		if(short_wave_damping) then
			call swrad_damping(P_sun,P_sun0, elevation, initial_climate_elevation)
		end if
		print*,'*** ERA-I climate loaded'

	    	


	end if
	



            ! Calculate Elevation Desertification (Budd and Smith (1979)), logic from Vizcaino et al. (2009)
		precipitation(:,:,:) = initial_climate_precipitation(:,:,:)
		print*,'random'

            if(active_elevation_desertification) then
                do id=1,ndays,1
                    precipitation(:,:,id) = (initial_climate_precipitation(:,:,id) )
                    print*,'random2'
		    ! where the initial elevation is below the threshold
                    where (initial_climate_elevation(:,:) .lt. precipitation_threshold)
                        ! Only where the elevation is above the threshold, otherwise precipitation stays the same
                        where (elevation(:,:) .gt. precipitation_threshold)
                            precipitation(:,:,id) = (initial_climate_precipitation(:,:,id) )&
                                * exp( - precipitation_lapse_rate &
                                * (elevation(:,:) - precipitation_threshold) )
                        end where
                    ! where the initial climate elevation is above the threshold
                    ! (in this case, the precipitation can get amplified in lower areas)
                    elsewhere
                        ! Where the elevation is above the threshold, the precipitation gets lower
                        where (elevation(:,:) .gt. precipitation_threshold)
                            precipitation(:,:,id) = (initial_climate_precipitation(:,:,id)  )&
                                * exp( - precipitation_lapse_rate &
                                * (elevation(:,:) - initial_climate_elevation(:,:)) )
                        ! If the elevation falls below the threshold, amplify the precipitation
                        elsewhere ! turn of precipitation amplification
                            precipitation(:,:,id) = (initial_climate_precipitation(:,:,id) )&
                                * exp( - precipitation_lapse_rate &
                                * (precipitation_threshold - initial_climate_elevation(:,:)) )
                        end where
                    end where

			! add precipitation deviation
			precipitation(:,:,id) = precipitation(:,:,id) + deviation_precip(:,:,id)

		    where(precipitation(:,:,id).lt. 0.)
			precipitation(:,:,id) = 0.
		    end where

                end do
 	        !print*,'minimum precip',minval(precipitation)
			
            end if



	! calculate elevation feedback on short wave radiation
	P_sun = P_sun0
	if(short_wave_damping) then
		call swrad_damping(P_sun, P_sun0, elevation, initial_climate_elevation )
	end if

	! add short wave radiation deviation
    	do id=1,ndays,1
		P_sun(:,:,id) = P_sun(:,:,id) + deviation_P_sun(:,:,id)

		where(P_sun(:,:,id).lt. 0.)
			P_sun(:,:,id) = 0.
		end where
	end do


    if(restart)then

        snowman=read_snow_data(restart_file, nx, ny,n_snowlayer, TRIM(adjustl('snowmass')))
        lwmass(1:nx,1:ny,1:n_snowlayer)=read_snow_data(restart_file, nx, ny, n_snowlayer, TRIM(adjustl('lwmass')) )
        snow_temp(1:nx,1:ny,1:n_snowlayer)=read_snow_data(restart_file, nx, ny, n_snowlayer, TRIM(adjustl('snowtemp')) )
        rho_snow(1:nx,1:ny,1:n_snowlayer)=read_snow_data(restart_file, nx, ny, n_snowlayer, TRIM(adjustl('snowdensity')) )
    end if

	! save the climate if wanted
	if((save_climate_forcing).or.(save_initial_climate)) then
		write (netcdf_output_climate_filename, '( "Climate_", I7.7, ".nc" )') 0
		call init_netcdf_climate_file(TRIM(adjustl(output_directory)) // netcdf_output_climate_filename,&
 				filehandle_netcdf3, ny, nx, int(dy), int(dx), precip_varid, temp_varid, swrad_varid,&
				 dev_precip_varid, dev_temp_varid, dev_swrad_varid )
		do id=1,ndays,1
			call writeNCDFGridValues(filehandle_netcdf3, id, precip_varid, real(precipitation(:,:,id)), ny, nx)
			call writeNCDFGridValues(filehandle_netcdf3, id, temp_varid, real(temperature(:,:,id)), ny, nx)
			call writeNCDFGridValues(filehandle_netcdf3, id, dev_precip_varid, real(deviation_precip(:,:,id)) , ny, nx)
			call writeNCDFGridValues(filehandle_netcdf3, id, dev_temp_varid, real(deviation_temp(:,:,id)), ny, nx)
			call writeNCDFGridValues(filehandle_netcdf3, id, swrad_varid, real(P_sun(:,:,id)), ny, nx)
			call writeNCDFGridValues(filehandle_netcdf3, id, dev_swrad_varid, real(deviation_P_sun(:,:,id)), ny, nx)
		end do
		call closeNCDFFile(filehandle_netcdf3)
	end if


    ! initialize netcdf output files
    !--------------------------------
    if (write_netcdf) then

        ! create names for files
        write (netcdf_output_filename , '( "IceBern2D_", I7.7, ".nc" )') 1
        !write (netcdf_output_filename2, '( "Snowcover3D_", I7.7, ".nc" )') 1

	
        ! creat netcdf file for the ice
        call init_netcdf_file(TRIM(adjustl(output_directory)) // netcdf_output_filename, filehandle_netcdf, ny, nx, &
		int(dy), int(dx), height_varid, bedrock_varid, acc_varid, diffusivity_varid, &
            discharge_x_varid, discharge_y_varid, assignent_mask_varid)!, abl_varid

        ! initialize netcdf file for snowdata
        !call init_netcdf_3D_file(TRIM(adjustl(output_directory)) // netcdf_output_filename2, filehandle_netcdf2, &
        !    ny, nx, n_snowlayer, int(dy), int(dx), 1, snowman_varid, lwmass_varid, rho_snow_varid, snow_temp_varid )

    end if

    ! store the data
    if(store_input_netcdf) then
        call copy_to_output_directory(netcdf_input_eraiterim_climate, output_directory)
        !call copy_to_output_directory(netcdf_input_precip_calib, output_directory)
        !call copy_to_output_directory(netcdf_input_swradboa, output_directory)
        !call copy_to_output_directory(netcdf_input_albedo, output_directory)
    end if

    ! NETCDF
    if (write_netcdf) then
        !        ! creat netcdf file for the ice
        !        call init_netcdf_file(TRIM(adjustl(output_directory)) // netcdf_output_filename, filehandle_netcdf, ny, nx, & 
	!		int(dy), int(dx), height_varid, bedrock_varid, acc_varid, abl_varid, diffusivity_varid, &
        !            discharge_x_varid, discharge_y_varid, assignent_mask_varid)

        ! Write initial values as year 0 to the netcdf
        ! but only height and bedrock, since the other values do not exist
        call writeNCDFGridValues(filehandle_netcdf, 1, &
            bedrock_varid, real(bedrock_netcdf(1:nx,1:ny)), ny, nx)

        call writeNCDFGridValues(filehandle_netcdf, 1, &
            height_varid,real(elevation_netcdf(1:nx,1:ny)), ny, nx)

        call writeNCDFGridIntegerValues(filehandle_netcdf, 1, &
            assignent_mask_varid,assignment_mask(1:nx,1:ny), ny, nx)
        ! Accumulation and Ablation
        call writeNCDFGridValues(filehandle_netcdf,  1, &
            acc_varid, real(accumulation(1:nx, 1:ny)) , ny, nx)
        !call writeNCDFGridValues(filehandle_netcdf,  1, &
        !    abl_varid, ablation(1:nx, 1:ny) , ny, nx)

        !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        ! fill snowdata into snowfile

        !call writeNCDF3DGridValues(filehandle_netcdf2, nc_counter, snowman_varid, real(snowman,4), ny, nx, n_snowlayer)
        !call writeNCDF3DGridValues(filehandle_netcdf2, nc_counter, lwmass_varid, real(lwmass,4), ny, nx, n_snowlayer)
        !call writeNCDF3DGridValues(filehandle_netcdf2, nc_counter, rho_snow_varid, real(rho_snow,4), ny, nx, n_snowlayer)
        !call writeNCDF3DGridValues(filehandle_netcdf2, nc_counter, snow_temp_varid, real(snow_temp,4), ny, nx, n_snowlayer)

        nc_counter = nc_counter + 1

    endif
    
    ! Integrated Mass Balance Output
    if (store_integrated_mass_balance) then
        open (unit=imb_filehandle,file=TRIM(adjustl(output_directory)) // imb_output_filename,action="write",status="replace")
        ! CSV Header
        ! year, accumulation, ablation, calving, isostaticmelt, totalmass, mass_change, net_change
        write (imb_filehandle,*) "year;accumulation;ablation;calving;isostaticmelt;totalmass;mass_change,net_change;sealevel"
    end if


    !=========================
    ! Lets go, do the loop
    !=========================
    if(debug > 0) then
        print *, "Loop over time steps"
        print *, "--------------------"
    endif
	call cpu_time(clock_start)
	CALL SYSTEM_CLOCK(c_start)
    do 
        ! Time series of diagnostics
        myyear(it) = it * int(real(dt)/(3600.*24.*365.)) ! the year that will be calculatet now
        ! Check if the loop conditions are at the end
        if ((myyear(it) > maxyears).or.((erai_backandforth_climate).and.(erai_year == erai_end_year+1)) ) then
            print *, 'The end is near (last year calculated): ', myyear(it)-1
            exit ! Jumps out of the loop. Does not exit the application (call exit(1)), otherwise the script is not terminated correctly
        end if


        !=====================================
        ! NEW OUTPUT FILES
        !=====================================
        ! start new netcdf files if time is up
        if((write_netcdf).and.(nc_counter > yearly_netcdf_file_freq))then
            ! close old netcdf files, generate new names and open new netcdf files
            call closeNCDFFile(filehandle_netcdf )
            !call closeNCDFFile(filehandle_netcdf2)

            write (netcdf_output_filename , '( "IceBern2D_", I7.7, ".nc" )') myyear(it)
            !write (netcdf_output_filename2, '( "Snowcover3D_", I7.7, ".nc" )') myyear(it)

            ! creat new netcdf file for the ice
            call init_netcdf_file(TRIM(adjustl(output_directory)) // netcdf_output_filename, filehandle_netcdf, ny, nx, & 
		int(dy), int(dx), height_varid, bedrock_varid, acc_varid, diffusivity_varid, &
                discharge_x_varid, discharge_y_varid, assignent_mask_varid)!, abl_varid

            ! initialize new netcdf file for snowdata
            !call init_netcdf_3D_file(TRIM(adjustl(output_directory)) // netcdf_output_filename2, filehandle_netcdf2, &
            !    ny, nx, n_snowlayer, int(dy), int(dx), 1, snowman_varid, lwmass_varid, rho_snow_varid, snow_temp_varid )
            ! reset next time slot to be written
            nc_counter=1
        end if


        !=====================================
        ! UPDATING CLIMATE AND MAP
        !=====================================
        ! this includes also elevation feedbacks and loading of new data

        ! Adjust Sea Level every 50 years
        if (store_integrated_mass_balance) then
            imb_isostaticmelt = 0d0
        end if


        ! load new climate data and update climate and swradboa
	!---------------------------------------------------------------

	if(erai_backandforth_climate) then
        
		! creat string for the next year
		CALL SYSTEM_CLOCK(c1)
		
		write (spec_format, '("(A", I0,",I",I0,".", I0,",A", I0,")")') &
		len_trim(netcdf_input_name_leading_string), netcdf_input_digit_specification, &
		netcdf_input_digit_specification , len_trim(netcdf_input_name_end_string) 
		if(debug > 1) then
            print*,'input_file_format',spec_format
        end if
		
		
		write (new_input_file,spec_format) &
		netcdf_input_name_leading_string,erai_year,netcdf_input_name_end_string
		
        print*,"*** Calculating mass balance with backandforth climate"
		!write (new_input_file , '( "ERAinterim_",I4.4,".interp.cdf" )') erai_year
! 		write (new_input_file , '( "LGM_global",I4.4,".interp.cdf" )') erai_year

		new_input_file = TRIM(adjustl(netcdf_input_eraiterim_directory)) // TRIM(adjustl(new_input_file))

		! set number of next input file
		! 1979 1980 ... 2015 2016 2016 2015 ... 1980 1979 1979 1980 ...
		if (erai_climate_backwards)then
			erai_year = erai_year-1
		else
			erai_year = erai_year+1
		end if

		if((erai_year == erai_year_turn+1).and.(erai_current_iteration .le. erai_iterations_max ))then
			erai_climate_backwards = .true.
			erai_year = erai_year_turn
		end if

		if(erai_year == erai_year_begin-1)then
			erai_current_iteration = erai_current_iteration+1
			erai_climate_backwards = .false.
			erai_year = erai_year_begin
		end if

		if (( (      erai_climate_backwards).and.(erai_year .ne. erai_year_turn-1 ) ) .or. &
		    ( (.not. erai_climate_backwards).and.(erai_year .ne. erai_year_begin+1) ) .or. &
		    (myyear(it) == 1)	) then

			! load new year
            call read_climate_once(new_input_file, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_temp_variable)), &
             TRIM(adjustl(netcdf_input_eraiterim_precip_variable)), TRIM(adjustl(netcdf_input_eraiterim_dewpT_variable)), &
            TRIM(adjustl(netcdf_input_eraiterim_wind_variable)), TRIM(adjustl(netcdf_input_eraiterim_lwrd_variable)), &
            TRIM(adjustl(netcdf_input_eraiterim_swradboa_variable)),&
            inp_temp, inp_precip, inp_dewpT, inp_wind, inp_lwrd, inp_swrd, rate)	

		!	! Remove errorous negative precipitation
		!	where(inp_precip(:,:,:) .lt. 0.)
		!		inp_precip(:,:,:) = 0.
		!	end where
		
			if(ndays==365) then
				do id=1,ndays,1
					temperature(:,:,id)= inp_temp(:,:,id)+kelvin ! TODO TODO remove test cooling
					precipitation(:,:,id)= inp_precip(:,:,id)/3600./24./1000. ! from mm/day to m/sec
                    DewpT(:,:,id)= inp_dewpT(:,:,id)
                    wind(:,:,id)= inp_wind(:,:,id)
                    if (longwave_downscaling) then
                        lwrd(:,:,id)= inp_lwrd(:,:,id)/temperature(:,:,id)**4*(temperature(:,:,id) +&
                        ((initial_climate_elevation-elevation) * temperature_lapse_rate))**4
                        else 
                        lwrd(:,:,id)= inp_lwrd(:,:,id)
                    end if 
                    if (lwrd_unit==2) then
                        lwrd(:,:,id)= lwrd(:,:,id)/3600./24.
                    end if
                    P_sun0(:,:,id)= inp_swrd(:,:,id)
					P_sun(:,:,id) = P_sun0(:,:,id)
					! Calculate potential temperature (at sea level)
					potential_temperature(:,:,id) = temperature(:,:,id) + (initial_climate_elevation * temperature_lapse_rate)

				end do

				initial_climate_precipitation = precipitation

			
			end if

			if(short_wave_damping) then
				call swrad_damping(P_sun,P_sun0, elevation, initial_climate_elevation)
			end if
			print*,'*** ERA-I ',trim(adjustl(new_input_file)),' climate loaded'

		end if

        CALL SYSTEM_CLOCK(c2)
        WRITE(*,*) "system_clock : ",(c2 - c1)/rate
	end if
	
	if(erai_reorder_climate) then
        
        
        reorder_year=erai_vector(erai_year)
        print*,"*** Reoder climate forcing using year order vector:", reorder_year
		! creat string for the next year
		CALL SYSTEM_CLOCK(c1)
		
		write (spec_format, '("(A", I0,",I",I0,".", I0,",A", I0,")")') &
		len_trim(netcdf_input_name_leading_string), netcdf_input_digit_specification, &
		netcdf_input_digit_specification , len_trim(netcdf_input_name_end_string) 
		if(debug > 1) then
            print*,'input_file_format',spec_format
        end if
		
		write (new_input_file,spec_format) &
		netcdf_input_name_leading_string,reorder_year,netcdf_input_name_end_string
		!write (new_input_file , '( "ERAinterim_",I4.4,".interp.cdf" )') reorder_year
! 		write (new_input_file , '( "LGM_global",I4.4,".interp.cdf" )') reorder_year

		new_input_file = TRIM(adjustl(netcdf_input_eraiterim_directory)) // TRIM(adjustl(new_input_file))

		! set number of next input file
		! 1979 1980 ... 2015 2016 2016 2015 ... 1980 1979 1979 1980 ...
		erai_year=erai_year+1

! 		if (erai_climate_backwards)then
! 			erai_year = erai_year-1
! 		else
! 			erai_year = erai_year+1
! 		end if
! 
! 		if((erai_year == erai_year_turn+1).and.(erai_current_iteration .le. erai_iterations_max ))then
! 			erai_climate_backwards = .true.
! 			erai_year = erai_year_turn
! 		end if
! 
! 		if(erai_year == erai_year_begin-1)then
! 			erai_current_iteration = erai_current_iteration+1
! 			erai_climate_backwards = .false.
! 			erai_year = erai_year_begin
! 		end if

		if (( (      erai_climate_backwards).and.(erai_year .ne. erai_year_turn-1 ) ) .or. &
		    ( (.not. erai_climate_backwards).and.(erai_year .ne. erai_year_begin+1) ) .or. &
		    (myyear(it) == 1)	) then

			! load new year
            call read_climate_once(new_input_file, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_temp_variable)), &
             TRIM(adjustl(netcdf_input_eraiterim_precip_variable)), TRIM(adjustl(netcdf_input_eraiterim_dewpT_variable)), &
            TRIM(adjustl(netcdf_input_eraiterim_wind_variable)), TRIM(adjustl(netcdf_input_eraiterim_lwrd_variable)), &
            TRIM(adjustl(netcdf_input_eraiterim_swradboa_variable)),&
            inp_temp, inp_precip, inp_dewpT, inp_wind, inp_lwrd, inp_swrd, rate)
            
            !optional reading of each variable individually but from the same file 
! 			inp_temp =   read_climate_long(new_input_file, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_temp_variable)) )
! 			inp_precip = read_climate_long(new_input_file, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_precip_variable)) )
!             inp_dewpT = read_climate_long(new_input_file, nx, ny, &
!             TRIM(adjustl(netcdf_input_eraiterim_dewpT_variable)) )
!             inp_wind = read_climate_long(new_input_file, nx, ny, &
!             TRIM(adjustl(netcdf_input_eraiterim_wind_variable)) )
!             inp_lwrd = read_climate_long(new_input_file, nx, ny, &
!             TRIM(adjustl(netcdf_input_eraiterim_lwrd_variable)) )

		!	! Remove errorous negative precipitation
		!	where(inp_precip(:,:,:) .lt. 0.)
		!		inp_precip(:,:,:) = 0.
		!	end where

			if(ndays==365) then
				do id=1,ndays,1
					temperature(:,:,id)= inp_temp(:,:,id)+kelvin ! TODO TODO remove test cooling
					precipitation(:,:,id)= inp_precip(:,:,id)/3600./24./1000. ! from mm/day to m/sec
                    DewpT(:,:,id)= inp_dewpT(:,:,id)
                    wind(:,:,id)= inp_wind(:,:,id)
                    if (longwave_downscaling) then
                        lwrd(:,:,id)= inp_lwrd(:,:,id)/temperature(:,:,id)**4*(temperature(:,:,id) +&
                        ((initial_climate_elevation-elevation) * temperature_lapse_rate))**4
                        else 
                        lwrd(:,:,id)= inp_lwrd(:,:,id)
                    end if 
                    if  (lwrd_unit==2) then
                        lwrd(:,:,id)= lwrd(:,:,id)/3600./24.
                    end if 
                    P_sun0(:,:,id)= inp_swrd(:,:,id)
					P_sun(:,:,id) = P_sun0(:,:,id)
					! Calculate potential temperature (at sea level)
					potential_temperature(:,:,id) = temperature(:,:,id) + (initial_climate_elevation * temperature_lapse_rate)

				end do

				!print*,'test temp',temperature(placex,placey,3)
				initial_climate_precipitation = precipitation


			end if

			if(short_wave_damping) then
				call swrad_damping(P_sun,P_sun0, elevation, initial_climate_elevation)
			end if
			print*,'*** ERA-I ',trim(adjustl(new_input_file)),' climate loaded'

		end if

        CALL SYSTEM_CLOCK(c2)
        WRITE(*,*) "system_clock : ",(c2 - c1)/rate
	end if
    
        ! Calculate the elevation feedback every 20 years.
	if ( (mod(myyear(it), calc_rate_climate) == 0) .or. (erai_backandforth_climate) .or. &
			(erai_reorder_climate))then  		!TODO set back to 50 years and remove that heating
		! if((myyear(it).ge. 101).and.(myyear(it).le. 200)) then	!TODO remove
		! 	potential_temperature = potential_temperature + 0.01	!TODO remove
		! end if 							!TODO remove

!		if(adjust_sea_level) then
!			sea_level = get_sea_level(ice_thickness, nx, ny, dx, dy, ocean_area)
!			! TODO: Do not take the initial bedrock, use the current one but forget the ice above it (and it should still increase, even if it is below the water
!			! But this could get to complicated, cause we do not want any sea in the middle of america.
!			assignment_mask = read_watermask(Bedrock_initial, nx, ny, (sea_level + sea_level_offset), &
!							! (use bedrock_netcdf, cause the bedrock is not equal to the sea level)
!							 Bedrock_netcdf, ice_thickness)
!		        
!			! integrated mass balance
!			! Die IF-Abfrage kann nicht innerhalb vom where integriert werden, daher einzeln
!			if (store_integrated_mass_balance) then
!				do ix=1,nx,1
!					do iy=1,ny,1
!						if ((assignment_mask(ix,iy) .eq. 1) .and. (ice_thickness(ix,iy) .gt. 0d0)) then
!							imb_isostaticmelt = imb_isostaticmelt - ice_thickness(ix,iy)
!						end if
!					end do
!				end do
!		        end if
!		        
!			! Set the sea level in the netcdf file for every water point
!			! 0 = ice, 1 = water, 2 = no ice, 3 = unstable integration
!			where(assignment_mask(:,:) .eq. 1)
!				Bedrock(:,:) = sea_level + sea_level_offset
!				elevation(:,:) = sea_level + sea_level_offset
!				! integrated mass balance calculated before
!				ice_thickness(:,:) = 0d0
!			end where
!			sunken_snow=0
!			do ix=1,nx,1
!				do iy=1,ny,1
!			    		! reset snow where there is sea
!			    		if((sum(snowman(ix,iy,:)).gt. 0.).and.(assignment_mask(ix,iy) .eq. 1)) then
!						sunken_snow = sunken_snow + sum(snowman(ix,iy,:)) + sum(lwmass(ix,iy,:))
!						do ii=1,n_snowlayer,1
!			    				snowman(ix,iy,ii)=0.
!							lwmass(ix,iy,ii)=0.
!							rho_snow(ix,iy,ii)=rho_s
!							snow_temp(ix,iy,ii)=0.
!						end do
!				    	end if
!				end do
!			end do
!		end if ! end adjust_sea_level


		! add the lapse rate for every day
		if(active_hysteresis) then
			do id = 1,ndays,1
			
				! It's still getting colder
				if ( myyear(it) .lt. hysteresis_return_point_in_time ) then
					! Check if the temperature is hold constant
					if((hysteresis_stop) .and. (hysteresis_stop_year .lt. myyear(it))) then
						! Do nothing, hysteresis_stop_temperature is reached
					else
					! 10 because it is calculated every 10 years
						potential_temperature(:,:,id) = potential_temperature(:,:,id) - (10 * hysteresis_temperature_factor)
					end if
					! Its getting warmer again
				else
					!if((hysteresis_stop_at_specific_temperature) == 2 .and. &
					!    (hysteresis_stop_temperature <= ((hysteresis_inital_temperature_offset - hysteresis_temperature_delta) &
					!    + ((myyear(it) - hysteresis_return_point_in_time) * hysteresis_temperature_factor)) ) ) then
					if((hysteresis_stop) .and. (hysteresis_stop_year .lt. myyear(it))) then
						! Do nothing, hysteresis_stop_temperature is reached
					else
						! 10 because it is calculated every 10 years
						potential_temperature(:,:,id) = potential_temperature(:,:,id) + (10. * hysteresis_temperature_factor)
					end if
				end if
			end do

		end if


		! calculate temperature lapsrate
! 		print*,(temperature(:,:,55)-potential_temperature(:,:,55))/elevation
!         print*,initial_climate_elevation-elevation
        inp_temp(:,:,:)=temperature(:,:,:)
		
		
		if(active_temperature_lapsrate) then
		temperature(:,:,:) = potential_temperature(:,:,:)
			do id = 1,ndays,1
				temperature(:,:,id) = potential_temperature(:,:,id)+ deviation_temp(:,:,id) - &
							(elevation * temperature_lapse_rate)
                if(active_dewpoint_lapserate)then
                    dewpT(:,:,id) = dewpT(:,:,id)+ (initial_climate_elevation-elevation) * dewpoint_lapse_rate
                end if
			end do
			
		end if
        

		! calculate elevation feedback on short wave radiation
		P_sun = P_sun0
		if(short_wave_damping) then
			call swrad_damping(P_sun, P_sun0, elevation, initial_climate_elevation )
		end if
		! add short wave radiation deviation
            	!do id=1,ndays,1
!			P_sun(:,:,id) = P_sun(:,:,id) + deviation_P_sun(:,:,id)
!			where(P_sun(:,:,id).lt. 0.)
!				P_sun(:,:,id) = 0.
!			end where
!		end do
		


            ! Calculate Elevation Desertification (Budd and Smith (1979)), logic from Vizcaino et al. (2009)
	    precipitation(:,:,:) = initial_climate_precipitation(:,:,:)
            if(active_elevation_desertification) then
                do id=1,ndays,1
                    precipitation(:,:,id) = (initial_climate_precipitation(:,:,id) )
		    ! where the initial elevation is below the threshold
                    where (initial_climate_elevation(:,:) .lt. precipitation_threshold)
                        ! Only where the elevation is above the threshold, otherwise precipitation stays the same
                        where (elevation(:,:) .gt. precipitation_threshold)
                            precipitation(:,:,id) = (initial_climate_precipitation(:,:,id) )&
                                * exp( - precipitation_lapse_rate &
                                * (elevation(:,:) - precipitation_threshold) )
                        end where
                    ! where the initial climate elevation is above the threshold
                    ! (in this case, the precipitation can get amplified in lower areas)
                    elsewhere
                        ! Where the elevation is above the threshold, the precipitation gets lower
                        where (elevation(:,:) .gt. precipitation_threshold)
                            precipitation(:,:,id) = (initial_climate_precipitation(:,:,id)  )&
                                * exp( - precipitation_lapse_rate &
                                * (elevation(:,:) - initial_climate_elevation(:,:)) )
                        ! If the elevation falls below the threshold, amplify the precipitation
                        elsewhere ! turn of precipitation amplification
                            precipitation(:,:,id) = (initial_climate_precipitation(:,:,id) )&
                                * exp( - precipitation_lapse_rate &
                                * (precipitation_threshold - initial_climate_elevation(:,:)) )
                        end where
                    end where

			precipitation(:,:,id) = precipitation(:,:,id) + deviation_precip(:,:,id)

		    where(precipitation(:,:,id).lt. 0.)
			precipitation(:,:,id) = 0.
		    end where
                end do
 	        !print*,'minimum precip',minval(precipitation)
            end if

            ! pdd accumulation model
            ! the accumulation/ablation i.e. surfacebalance is calculatet for one year
            ! the pdd model smb is then used for 10 years
            ! pdd model
            if(smb_model==1)then
                accumulation = temperature_dependent_accumulation(nx, ny, temperature, precipitation, &
                    precipitation_unit, accumulation_daily_temperature_threshold, ndays)
                ablation = temperature_dependent_ablation(nx, ny, temperature, beta, ndays)
               

                surface_mass_balance = (accumulation - ablation)/seconds_per_year

		where( ( elevation_netcdf(:,:).lt. sea_level + sea_level_offset ) )
			surface_mass_balance(:,:) = min( surface_mass_balance(:,:),0.)
		end where


            end if

                ! fun model
            if(smb_model==3)then
                ! insert yearly data
                surface_mass_balance  = get_troll_accumulation(nx,ny,ndays,temperature, precipitation, assignment_mask, P_sun)
                !print *, "can you see me? "
                !print *, 'Year: ', myyear(it)

		where( ( elevation_netcdf(:,:).lt. sea_level + sea_level_offset ) )
			surface_mass_balance(:,:) = min( surface_mass_balance(:,:),0.)
		end where

                accumulation = surface_mass_balance
                surface_mass_balance = surface_mass_balance/seconds_per_year
            end if
		
        end if ! mod 20 year if for elevation feedbacks
          


	! write climate forcing data to netcdf
	!-------------------------------------
        if((mod(myyear(it), monthly_data_freq)==0).and.(save_climate_forcing)) then 
		write (netcdf_output_climate_filename, '( "Climate_", I7.7, ".nc" )') myyear(it)
		call init_netcdf_climate_file(TRIM(adjustl(output_directory)) // netcdf_output_climate_filename, &
 				filehandle_netcdf3, ny, nx, int(dy), int(dx), precip_varid, temp_varid, swrad_varid,&
				 dev_precip_varid, dev_temp_varid, dev_swrad_varid )
		do id=1,ndays,1
			call writeNCDFGridValues(filehandle_netcdf3, id, precip_varid, real(precipitation(:,:,id)), ny, nx)
			call writeNCDFGridValues(filehandle_netcdf3, id, temp_varid, real(temperature(:,:,id)), ny, nx)
			call writeNCDFGridValues(filehandle_netcdf3, id, dev_precip_varid, real(deviation_precip(:,:,id)) , ny, nx)
			call writeNCDFGridValues(filehandle_netcdf3, id, dev_temp_varid, real(deviation_temp(:,:,id)), ny, nx)
			call writeNCDFGridValues(filehandle_netcdf3, id, swrad_varid, real(P_sun(:,:,id)), ny, nx)
! 			call writeNCDFGridValues(filehandle_netcdf3, id, dev_swrad_varid, real(deviation_P_sun(:,:,id)), ny, nx)
		end do
		call closeNCDFFile(filehandle_netcdf3)
	end if

	
        !=====================================
        ! SNOW MODEL OR SMB STARTS HERE
        !=====================================

        ! the emb model accumulation must be calculatet for every year seperatly due to theyr descrete accumulation
        ! the models calculate the cumulative smb for one year and then the icesheet is adjusted
        ! emb model
        if(smb_model==2)then
            !print *,'myyear(it)=', myyear(it)

		call get_accumulation_snowman(nx, ny, ndays, n_snowlayer, temperature, precipitation, P_sun, assignment_mask,&
		        	elevation_netcdf, sea_level + sea_level_offset, snowman, lwmass, rho_snow, snow_temp, &
				surface_mass_balance, nc_counter, myyear(it), sunken_snow, adaptive_timestep, fast_calculation,&
				albedo_dynamic, DewpT, lwrd, wind, elevation,new_input_file) ! topography and S_BOA

		sunken_snow=0.
		
		! surface_mass_balance in m_ice/second
		accumulation = surface_mass_balance*seconds_per_year ! in m_ice/ year
		    	!ablation = sum(snowman,3) ! snow in kg/m2 




!		if ((myyear(it).lt. start_speed_up ).or.(myyear(it).gt. end_speed_up)  )then
!			spinup=.false.
!		    	call get_accumulation_snowman(nx, ny, ndays, n_snowlayer, temperature, precipitation, P_sun, assignment_mask,&
!		        	snowman, lwmass, rho_snow, snow_temp, surface_mass_balance, nc_counter, myyear(it),sunken_snow, spinup  ) ! topography and S_BOA
!			sunken_snow=0
!		    	
!		    	accumulation = surface_mass_balance ! in m_ice/ year
!		    	!ablation = sum(snowman,3) ! snow in kg/m2
!		    	surface_mass_balance = surface_mass_balance/seconds_per_year ! in m_ice/second
!		
!
!		elseif ((mod(myyear(it)-50,)==0))then !.or.(mod(myyear(it)-50,100)==1))then
!			spinup=.true.
!		    	call get_accumulation_snowman(nx, ny, ndays, n_snowlayer, temperature, precipitation, P_sun, assignment_mask,&
!		        	snowman, lwmass, rho_snow, snow_temp, surface_mass_balance, nc_counter, myyear(it),sunken_snow, spinup ) ! topography and S_BOA
!			sunken_snow=0
!		    	
!		    	accumulation = surface_mass_balance ! in m_ice/ year
!		    	!ablation = sum(snowman,3) ! snow in kg/m2
!		    	surface_mass_balance = surface_mass_balance/seconds_per_year ! in m_ice/second
!		end if
        end if  

                print*,'*** end of snow part of ice model'

        !------------------------smb part ends here-------------------------------------

        !=====================================
        ! ICE STARTS HERE
        !=====================================

        ! //////////////////////// ICE PHYSICS ENDS HERE //////////////////////


        ! Print Heartbeat
        if ((mod(myyear(it), netcdf_timesteps) == 0) .and. (debug > 0)  ) then !.and. (last_netcdf_year .lt. myyear(it))
            !call ETIME(execution_time, runtime)
	     call cpu_time(clock_end)
	     CALL SYSTEM_CLOCK(c_end)
	    runtime = runtime + clock_end - clock_start
            Write( heartbeat, '(i8)' ) myyear(it)
            print *, 'Year: ', TRIM(adjustl(heartbeat)), ', Runtime [s]: ', int(runtime), &
                !', End in [s]: ', int(( (runtime/myyear(it)) * maxyears) - runtime), &
                !', s/1000yr: ', (runtime/myyear(it) * 1000), &
                ', current yr/hour: ', (netcdf_timesteps/(clock_end - clock_start)*3600.), &		
                ', ave yr/hour: ', (myyear(it)/runtime * 60*60), &
                ', Sea level: ', (sea_level + sea_level_offset), &
                ', NetCDF TS: ', nc_counter-1 !((myyear(it)/netcdf_timesteps) + 1)
                print*,'system_time',(c_end-c_start)/rate
          CALL SYSTEM_CLOCK(c_start)
	     call cpu_time(clock_start)
        end if


	    call cpu_time(clock_end)
	    CALL SYSTEM_CLOCK(c_end)
	    runtime = runtime + clock_end - clock_start
            print *, "Runtime [mm:ss]: ", int(runtime/60), ':', mod(int(runtime),60)
            print *, 'Year: ', myyear(it)
           
            print *, 'TEST'
        
 


    if(debug > 0 ) then
        print *, '========================'
        print *, "Successfully terminated!"
        print *, "We are at the end"
        print *, "-----------------"
        print *, "Some statistics:"
        !call ETIME(execution_time, runtime)
	call cpu_time(clock_end)
	CALL SYSTEM_CLOCK(c_end)
	runtime = runtime + clock_end - clock_start
        print *, "Runtime [mm:ss]: ", int(runtime/60), ':', mod(int(runtime),60)
        print *, "Seconds per 1000 years: ", (runtime/maxyears * 1000)
        print *, 'yr/hour: ', (maxyears/runtime * 60*60)
        print *, "-----------------"
        if (store_input_netcdf .or. write_netcdf .or. annual_data .or. daily_data .or. monthly_data) then
            print *, 'The output is stored in the directory: ', TRIM(adjustl(output_directory))
        else
            print *, 'No output was stored from this run!'
        endif
    end if
    print*,it
     it = int(it) + 1
end do

end program
