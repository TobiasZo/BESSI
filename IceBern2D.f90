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
    INTEGER :: c1,c2,cr,cm,c_start,c_end
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
	sea_level = get_sea_level(Ice_thickness, nx, ny, dx, dy, ocean_area)

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

	if((eraiterim_climate).or.(erai_backandforth_climate))then	
		! Load climate reference elevation (includes ice)
	    	initial_climate_elevation = read_variable(netcdf_input_eraiterim_initial_climate_elevation, nx, ny, &
					     TRIM(adjustl(netcdf_input_eraiterim_initial_climate_elevation_variable)))
	
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

		! reduce to 96 Timesteps
		if(ndays==96) then
			do id=1,ndays,1
				temperature(:,:,id)= inp_temp(:,:,nint(real(id)*3.8))+kelvin ! TODO TODO remove test cooling
				precipitation(:,:,id)= inp_precip(:,:,nint(real(id)*3.8))/3600./24./1000. ! from mm/day to m/sec

				! Calculate potential temperature (at sea level)
				potential_temperature(:,:,id) = temperature(:,:,id) + (initial_climate_elevation * temperature_lapse_rate)

			end do

			!print*,'test temp',temperature(placex,placey,3)
			initial_climate_precipitation = precipitation



			! load short wave radiation
			inp_temp =   read_climate_long(netcdf_input_eraiterim_climate, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_swradboa_variable)) )

			! reduce to 96 Timesteps
			do id=1,ndays,1
				P_sun0(:,:,id)= inp_temp(:,:,nint(real(id)*3.8))
				P_sun(:,:,id) = P_sun0(:,:,id)
			end do
		end if


		if(ndays==365) then
			do id=1,ndays,1
				temperature(:,:,id)= inp_temp(:,:,id)+kelvin ! TODO TODO remove test cooling
				precipitation(:,:,id)= inp_precip(:,:,id)/3600./24./1000. ! from mm/day to m/sec
                DewpT(:,:,id)= inp_dewpT(:,:,id)
                wind(:,:,id)= inp_wind(:,:,id)
                lwrd(:,:,id)= inp_lwrd(:,:,id)
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

	    	!! read PD bern3d reference data and Insolation
	     	!call read_climate_bern3d(b3d_pd_temp, b3d_pd_precip, P_sun, netcdf_input_b3d_pd_climate, thematrix  )
		!if(debug > 0 ) then
		! 	print *, '*** Interpolation done, PD b3d climate loaded '
		!end if



	!else
	    	!initial_climate_precipitation = precipitation
	    	!potential_temperature = temperature
		!initial_climate_elevation = 0.
	end if
	


	! for debuging only
	!do ix=1,nx,1	
	!	do iy=1,ny,1
	!		P_sun(ix,iy,:) = sum(P_sun(ix,iy,:))/96
	!	end do
	!end do


	! calculate bern3d climate deviation for ancient climate
	! this is simply be added to the topographically corrected climate data (temperature and precipitation). 
	! maybe make this to subroutine...
	if(transient_ancient_climate) then


	 	print *, "Build interpolation matrix"
		call load_mat(thematrix,interpol_deg)
		print *, "Interpolation matrix crafted"

		!! read PD bern3d reference data and Insolation
	     	call read_climate_bern3d(b3d_pd_temp, b3d_pd_precip, b3d_pd_P_sun, netcdf_input_b3d_pd_climate, thematrix  )
		print*,'*** Bern3D reference climate (1950 AD)'

		! creat string for b3d initial file path depending 
		if(year1+(input_file_no)*new_climate<0.) then
			write (new_input_file , '( "EEM2PIC",I3.3,".", I7.6, "_full_inst.nc" )') input_file_no, (year1+(input_file_no)*new_climate)
		else
			write (new_input_file , '( "EEM2PIC",I3.3,".", I7.7, "_full_inst.nc" )') input_file_no, (year1+(input_file_no)*new_climate)
		end if
		input_file_no = input_file_no+1
 
		netcdf_input_b3d_timestep_climate_file = TRIM(adjustl(netcdf_input_b3d_timestep_climate_path)) // TRIM(adjustl(new_input_file))
		print *, "load climtology B3D from file ", trim(adjustl(new_input_file))
		call read_climate_bern3d(deviation_temp, deviation_precip, deviation_P_sun, netcdf_input_b3d_timestep_climate_file , thematrix )


		do id=1,ndays,1
			deviation_temp(:,:,id) = (deviation_temp(:,:,id)-b3d_pd_temp(:,:,id))
			deviation_P_sun(:,:,id) = (deviation_P_sun(:,:,id)-b3d_pd_P_sun(:,:,id))
			if(turnoff_dprecip)then
				deviation_precip(:,:,id) = 0.
			else
				deviation_precip(:,:,id) = (deviation_precip(:,:,id)-b3d_pd_precip(:,:,id))*b3d_precip_multiplier
				! calibrate precipitation deviations to altitude
				!where (elevation(:,:) .gt. precipitation_threshold)
                            		!deviation_precip(:,:,id) = deviation_precip(:,:,id) &
                            		!* exp( - precipitation_lapse_rate &
                            		!* (elevation(:,:) - precipitation_threshold) )
                    		!end where
			end if
		end do
		print*,'ave temp rel 1450AD :', sum(deviation_temp(:,:,:))/313./313./96.
	end if
	
	if(ccsm_eem_climate)then
		!! load solar iradiation from bern3d eemian (123000 BC)
    		!call read_climate_bern3d(b3d_pd_temp, b3d_pd_precip, P_sun, netcdf_input_b3d_eem_climate, thematrix  )
    		!if(debug > 0 ) then
        	!	print *, '*** Interpolation done, EEM b3d swradboa loaded '
    		!end if

		! load solar radiation from ccsm3 125000 BP
		call read_swradboa_ccsm(P_sun0, netcdf_input_ccsm_eem_swradboa , netcdf_input_ccsm_lgm_swradboa_variable )
		P_sun = P_sun0	
		if(short_wave_damping) then
			call swrad_damping(P_sun,P_sun0, elevation, initial_climate_elevation )
		end if


		! load and calculate ccsm climate deviation EEM-PI

		! initial_ccsm_elevation
    		initial_climate_elevation = read_variable(netcdf_input_ccsm_initial_climate_elevation, nx, ny, &
				     	  TRIM(adjustl(netcdf_input_ccsm_initial_climate_elevation_variable)))
	
		inp_temp =   read_climate_long(netcdf_input_ccsm_eem_climate, nx, ny, TRIM(adjustl(netcdf_input_ccsm_temp_variable)) ) - &
				read_climate_long(netcdf_input_ccsm_pi_climate, nx, ny, TRIM(adjustl(netcdf_input_ccsm_temp_variable)) )
		inp_precip = read_climate_long(netcdf_input_ccsm_eem_climate, nx, ny, TRIM(adjustl(netcdf_input_ccsm_precip_variable)) ) - &
				read_climate_long(netcdf_input_ccsm_pi_climate, nx, ny, TRIM(adjustl(netcdf_input_ccsm_precip_variable)) )

				!! reduce to 96 Timesteps
				!do id=1,ndays,1
				!	temperature(:,:,id)= inp_temp(:,:,nint(real(id)*3.8))!+kelvin
				!	precipitation(:,:,id)= inp_precip(:,:,nint(real(id)*3.8))!/3600./24./1000. ! from mm/day to m/sec

				!	! Calculate potential temperature (at sea level)
				!	potential_temperature(:,:,id) = temperature(:,:,id) + (initial_climate_elevation * temperature_lapse_rate)

				!	!! Adjust temperature to elevation!
			    	!	!temperature(:,:,id) = potential_temperature(:,:,id) - (elevation * temperature_lapse_rate)
				!end do

				!print*,'test temp',temperature(placex,placey,3)prect_125ka_c1_daily-clim_3-33.nc
				!initial_climate_precipitation = precipitation


		do id=1,ndays,1
			! reduce to 96 Timesteps
			deviation_temp(:,:,id)= inp_temp(:,:,nint(real(id)*3.8))
			deviation_precip(:,:,id)= inp_precip(:,:,nint(real(id)*3.8)) ! from m/sec to m/sec

			! scale devtiation_precipitation from initial_ccsm_elevation to initial_climate_elevation
			!(scale 	   precipitation from initial_climate_elevation to elevation)
			
                	! where the initial ccsm elevation is below the threshold
                	where (initial_ccsm_elevation(:,:) .lt. precipitation_threshold)
                    	! Only where the initial_climate_elevation is above the threshold, otherwise precipitation stays the same
                    		where (initial_climate_elevation(:,:) .gt. precipitation_threshold)
                        		precipitation(:,:,id) = deviation_precip(:,:,id) &
                            		* exp( - precipitation_lapse_rate &
                            		* (initial_climate_elevation(:,:) - precipitation_threshold) )
                    		end where
                		! where the initial ccsm elevation is above the threshold
                		! (in this case, the precipitation can get amplified in lower areas)
                	elsewhere
            			! Where the initial_climate_elevation is above the threshold, the precipitation gets lower
            			where (initial_climate_elevation(:,:) .gt. precipitation_threshold)
                			precipitation(:,:,id) = deviation_precip(:,:,id)  &
                    			* exp( - precipitation_lapse_rate &
                    			* (initial_climate_elevation(:,:) - initial_ccsm_elevation(:,:)) )
            			! If the initial_climate_elevation falls below the threshold, amplify the precipitation
            			elsewhere
                			precipitation(:,:,id) = deviation_precip(:,:,id) &
                    			* exp( - precipitation_lapse_rate &
                    			* (precipitation_threshold - initial_ccsm_elevation(:,:)) )
            			end where
                	end where

		end do
		print*,'ave temp rel 1950AD :', sum(deviation_temp(:,:,:))/313./313./96.
		! add to potential_temperature and initial_climate_precipitation

		potential_temperature = potential_temperature + deviation_temp
		initial_climate_precipitation = initial_climate_precipitation + deviation_precip !precipitation
		where(initial_climate_precipitation(:,:,:)<0.)
			initial_climate_precipitation(:,:,:)=0.
		end where
		precipitation = initial_climate_precipitation 
		!write (netcdf_output_climate_filename, '( "Climate_eem123000BC.nc" )') 
		!call init_netcdf_climate_file(TRIM(adjustl(output_directory)) // netcdf_output_climate_filename,&
 		!		filehandle_netcdf3, ny, nx, int(dy), int(dx), precip_varid, temp_varid, dev_precip_varid, dev_temp_varid )
		!do id=1,ndays,1
		!	call writeNCDFGridValues(filehandle_netcdf3, id, precip_varid, precipitation(:,:,id), ny, nx)
		!	call writeNCDFGridValues(filehandle_netcdf3, id, temp_varid, temperature(:,:,id), ny, nx)
		!	call writeNCDFGridValues(filehandle_netcdf3, id, dev_precip_varid, deviation_precip(:,:,id) , ny, nx)
		!	call writeNCDFGridValues(filehandle_netcdf3, id, dev_temp_varid, deviation_temp(:,:,id), ny, nx)
		!end do
		!call closeNCDFFile(filehandle_netcdf3)
		deviation_temp = 0.
		deviation_precip = 0.
	end if
	
	if(dev_b3d_eem_climate)then
		! load solar iradiation and climate from bern3d eemian (123000 BC)	
		call read_climate_bern3d(deviation_temp, deviation_precip, P_sun, netcdf_input_b3d_eem_climate , thematrix )
		print *, "load climtology B3D from file ", trim(adjustl(netcdf_input_b3d_eem_climate))		
		do id=1,ndays,1
			deviation_temp(:,:,id) = (deviation_temp(:,:,id)-b3d_pd_temp(:,:,id))
			if(turnoff_dprecip)then
				deviation_precip(:,:,id) = 0.
			else
				deviation_precip(:,:,id) = (deviation_precip(:,:,id)-b3d_pd_precip(:,:,id))
				! calibrate precipitation deviations to altitude
				!where (elevation(:,:) .gt. precipitation_threshold)
                            		!deviation_precip(:,:,id) = deviation_precip(:,:,id) &
                            		!* exp( - precipitation_lapse_rate &
                            		!* (elevation(:,:) - precipitation_threshold) )
                    		!end where
			end if
		end do
		print*,'ave temp rel 1950AD :', sum(deviation_temp(:,:,:))/313./313./96.
	
		! add to potential_temperature and initial_climate_precipitation

		!potential_temperature = potential_temperature + deviation_temp
		!initial_climate_precipitation = initial_climate_precipitation + deviation_precip !precipitation
		!where(initial_climate_precipitation(:,:,:)<0)
		!	initial_climate_precipitation(:,:,:)=0
		!end where
		!precipitation = initial_climate_precipitation 
		!write (netcdf_output_climate_filename, '( "Climate_eem123000BC.nc" )') 
		!call init_netcdf_climate_file(TRIM(adjustl(output_directory)) // netcdf_output_climate_filename,&
 		!		filehandle_netcdf3, ny, nx, int(dy), int(dx), precip_varid, temp_varid, dev_precip_varid, dev_temp_varid )
		!do id=1,ndays,1
		!	call writeNCDFGridValues(filehandle_netcdf3, id, precip_varid, precipitation(:,:,id), ny, nx)
		!	call writeNCDFGridValues(filehandle_netcdf3, id, temp_varid, temperature(:,:,id), ny, nx)
		!	call writeNCDFGridValues(filehandle_netcdf3, id, dev_precip_varid, deviation_precip(:,:,id) , ny, nx)
		!	call writeNCDFGridValues(filehandle_netcdf3, id, dev_temp_varid, deviation_temp(:,:,id), ny, nx)
		!end do
		!call closeNCDFFile(filehandle_netcdf3)
		!deviation_temp = 0
		!deviation_precip = 0
	end if



	if (ccsm_lgmpt_climate) then

		! load solar iradiation from bern3d eemian (19000 BC)
    		call read_climate_bern3d(b3d_pd_temp, b3d_pd_precip, P_sun, netcdf_input_b3d_lgm_climate, thematrix  )
    		if(debug > 0 ) then
        		print *, '*** Interpolation done, LGM b3d swradboa loaded '
    		end if
		! load lgm ccsm climate 

		! 
    		initial_climate_elevation = read_variable(netcdf_input_ccsm_initial_climate_elevation, nx, ny, &
				     	  TRIM(adjustl(netcdf_input_ccsm_initial_climate_elevation_variable)))
	
		inp_temp =   read_climate_long(netcdf_input_ccsm_lgmpt_climate, nx, ny, TRIM(adjustl(netcdf_input_ccsm_lgm_temp_variable)) ) 
		inp_precip = read_climate_long(netcdf_input_ccsm_lgmpt_climate, nx, ny, TRIM(adjustl(netcdf_input_ccsm_lgm_precip_variable)) ) 

		! reduce to 96 Timesteps
		do id=1,ndays,1
			temperature(:,:,id)= inp_temp(:,:,nint(real(id)*3.8))!+kelvin
			precipitation(:,:,id)= inp_precip(:,:,nint(real(id)*3.8))!/3600./24./1000. ! from mm/day to m/sec

			! Calculate potential temperature (at sea level)
			potential_temperature(:,:,id) = temperature(:,:,id) + (initial_climate_elevation * temperature_lapse_rate)

			!! Adjust temperature to elevation!
            		!temperature(:,:,id) = potential_temperature(:,:,id) - (elevation * temperature_lapse_rate)
		end do

		!print*,'test temp',temperature(placex,placey,3)
        	initial_climate_precipitation = precipitation
		
	end if

	if (ccsm_lgmlt_climate) then

		!! load solar iradiation from bern3d eemian (19000 BC)
    		!call read_climate_bern3d(b3d_pd_temp, b3d_pd_precip, P_sun, netcdf_input_b3d_lgm_climate, thematrix  )
    		!if(debug > 0 ) then
        	!	print *, '*** Interpolation done, LGM b3d swradboa loaded '
    		!end if


		! Load elevation with ice
    		initial_climate_elevation = read_variable(netcdf_input_ccsm_lgm_initial_climate_elevation, nx, ny, &
				     	  TRIM(adjustl(netcdf_input_ccsm_lgm_initial_climate_elevation_variable)))
	
		inp_temp =   read_climate_long(netcdf_input_ccsm_lgmlt_climate, nx, ny, TRIM(adjustl(netcdf_input_ccsm_lgm_temp_variable)) )
		inp_precip = read_climate_long(netcdf_input_ccsm_lgmlt_climate, nx, ny, TRIM(adjustl(netcdf_input_ccsm_lgm_precip_variable)) )


		! reduce to 96 Timesteps
		do id=1,ndays,1
			temperature(:,:,id)= inp_temp(:,:,nint(real(id)*3.8))! -10.  !+kelvin
			precipitation(:,:,id)= inp_precip(:,:,nint(real(id)*3.8))!/3600./24./1000. ! from mm/day to m/sec

			! Calculate potential temperature (at sea level)
			potential_temperature(:,:,id) = temperature(:,:,id) + (initial_climate_elevation * temperature_lapse_rate)


		end do

		!print*,'test temp',temperature(placex,placey,3)
        	initial_climate_precipitation = precipitation
	

		call read_swradboa_ccsm(P_sun0, netcdf_input_ccsm_lgm_swradboa , netcdf_input_ccsm_lgm_swradboa_variable )
		P_sun = P_sun0
		if(short_wave_damping) then
			call swrad_damping(P_sun,P_sun0, elevation, initial_climate_elevation )
		end if


		! load and calculate ccsm climate deviation EEM-PI

		!initial_ccsm_elevation = read_variable(netcdf_input_ccsm_lgm_initial_climate_elevation, nx, ny, &
		!		     	  TRIM(adjustl(netcdf_input_ccsm_lgm_initial_climate_elevation_variable)))
	
		!inp_temp =   read_climate_long(netcdf_input_ccsm_lgmlt_climate, nx, ny, TRIM(adjustl(netcdf_input_ccsm_lgm_temp_variable)) ) - &
		!		read_climate_long(netcdf_input_ccsm_pi_climate, nx, ny, TRIM(adjustl(netcdf_input_ccsm_temp_variable)) )

		!inp_precip = read_climate_long(netcdf_input_ccsm_lgmlt_climate, nx, ny, TRIM(adjustl(netcdf_input_ccsm_lgm_precip_variable)) ) 
				
		!do id=1,365,1
		!
		!	! scale devtiation_precipitation from initial_ccsm_elevation to initial_climate_elevation
		!	!(scale 	   precipitation from initial_climate_elevation to elevation)
		!	
                !	! where the initial ccsm elevation is below the threshold
                !	where (initial_ccsm_elevation(:,:) .lt. precipitation_threshold)
                !    	! Only where the initial_climate_elevation is above the threshold, otherwise precipitation stays the same
                !    		where (initial_climate_elevation(:,:) .gt. precipitation_threshold)
                !        		inp_precip(:,:,id) = inp_precip(:,:,id) &
                !            		* exp( - precipitation_lapse_rate &
                !            		* (initial_climate_elevation(:,:) - precipitation_threshold) )
                !    		end where
                ! 		! where the initial ccsm elevation is above the threshold
                !		! (in this case, the precipitation can get amplified in lower areas)
                !		elsewhere
                !    		! Where the initial_climate_elevation is above the threshold, the precipitation gets lower
                !    			where (initial_climate_elevation(:,:) .gt. precipitation_threshold)
                !        			inp_precip(:,:,id) = inp_precip(:,:,id)  &
                !           			* exp( - precipitation_lapse_rate &
                !            			* (initial_climate_elevation(:,:) - initial_ccsm_elevation(:,:)) )
                !    			! If the initial_climate_elevation falls below the threshold, amplify the precipitation
                !    			elsewhere
                !        			inp_precip(:,:,id) = inp_precip(:,:,id) &
                !            			* exp( - precipitation_lapse_rate &
                !            			* (precipitation_threshold - initial_ccsm_elevation(:,:)) )
                !    			end where
                !		end where
		!end do
 

		!inp_precip = inp_precip - &
		!	read_climate_long(netcdf_input_ccsm_pi_climate, nx, ny, TRIM(adjustl(netcdf_input_ccsm_precip_variable)) )


		!do id=1,ndays,1
		!	! reduce to 96 Timesteps and adjust ccsm lgm temperature to ccsm pi temperature altitude, ie. remove altitude difference
		!	deviation_temp(:,:,id)= inp_temp(:,:,nint(real(id)*3.8)) + &
		!		(initial_ccsm_elevation(:,:)-initial_climate_elevation(:,:))* temperature_lapse_rate
		!	deviation_precip(:,:,id)= inp_precip(:,:,nint(real(id)*3.8)) ! from m/sec to m/sec
		!end do
		!print*,'ave temp rel 1750AD :', sum(deviation_temp(:,:,:))/313./313./96.
		
	end if




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

			! add precipitation deviation
			precipitation(:,:,id) = precipitation(:,:,id) + deviation_precip(:,:,id)

		    where(precipitation(:,:,id).lt. 0.)
			precipitation(:,:,id) = 0.
		    end where

                end do
 	        !print*,'minimum precip',minval(precipitation)
			
            end if

	! calculate temperature lapsrate
	temperature(:,:,:) = potential_temperature(:,:,:)
	if(active_temperature_lapsrate) then
		do id = 1,ndays,1
			temperature(:,:,id) = potential_temperature(:,:,id)+ deviation_temp(:,:,id) - &
						(elevation * temperature_lapse_rate)
		end do	
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


        ! Calculate Elevation Desertification (Budd and Smith (1979)), logic from Vizcaino et al. (2009)
!        if(active_elevation_desertification) then
!            do id=1,ndays,1
!		precipitation(:,:,id) = (initial_climate_precipitation(:,:,id) + deviation_precip(:,:,id))
!                ! where the initial elevation is below the threshold
!                where (initial_climate_elevation(:,:) .lt. precipitation_threshold)
!                    ! Only where the elevation is above the threshold, otherwise precipitation stays the same
!                    where (elevation(:,:) .gt. precipitation_threshold)
!                        precipitation(:,:,id) = (initial_climate_precipitation(:,:,id) + deviation_precip(:,:,id))&
!                            * exp( - precipitation_lapse_rate &
!                            * (elevation(:,:) - precipitation_threshold) )
!                    end where
!                ! where the initial climate elevation is above the threshold
!                ! (in this case, the precipitation can get amplified in lower areas)
!                elsewhere
!                    ! Where the elevation is above the threshold, the precipitation gets lower
!                    where (elevation(:,:) .gt. precipitation_threshold)
!                        precipitation(:,:,id) = (initial_climate_precipitation(:,:,id) + deviation_precip(:,:,id)) &
!                            * exp( - precipitation_lapse_rate &
!                            * (elevation(:,:) - initial_climate_elevation(:,:)) )
!                    ! If the elevation falls below the threshold, amplify the precipitation
!                    elsewhere ! turn of precipitation amplification
!                        precipitation(:,:,id) = (initial_climate_precipitation(:,:,id)+ deviation_precip(:,:,id) )&
!                            * exp( - precipitation_lapse_rate &
!                            * (precipitation_threshold - initial_climate_elevation(:,:)) )
!                    end where
!                end where
!                temperature(:,:,id) = potential_temperature(:,:,id) + deviation_temp(:,:,id) -&
!								 (elevation * temperature_lapse_rate)
!		!precipitation(:,:,id) = precipitation(:,:,id) + deviation_precip(:,:,id)
!
!		where(precipitation(:,:,id).lt.0)
!			precipitation(:,:,id) = 0.
!		end where
 !           end do
  !      end if

	!print*,'minimum precip',minval(precipitation)


	! initialize snow where there is ice
	if(initialize_snow)then
		do ix=1,nx,1
			do iy=1,ny,1
				if(ice_thickness(ix,iy) .gt. 200.) then
					snowman(ix,iy,1:n_snowlayer-1) = soll_mass
					lwmass(ix,iy,1:n_snowlayer-1) = 0.
					snow_temp(ix,iy,1:n_snowlayer-1) = min(temperature(ix,iy,1),kelvin)
					rho_snow(ix,iy,1:n_snowlayer-1) = 550.
				end if
			end do
		end do
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
        call copy_to_output_directory(netcdf_input_b3d_pd_climate, output_directory)
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
	! calculate bern3d climate deviation for ancient climate
	! this is added to the topographically corrected climate data (temperature and precipitation). 
	! maybe make this to subroutine...
	if((transient_ancient_climate).and.(mod(myyear(it), new_climate)==0).and.(input_file_no<=260)) then
		! creat string for b3d initial file path depending 
		
		if(year1+(input_file_no)*new_climate<0.) then
			write (new_input_file , '( "EEM2PIC",I3.3,".", I7.6, "_full_inst.nc" )') input_file_no, (year1+(input_file_no)*new_climate)
		else
			write (new_input_file , '( "EEM2PIC",I3.3,".", I7.7, "_full_inst.nc" )') input_file_no, (year1+(input_file_no)*new_climate)
		end if

		input_file_no = input_file_no+1

		netcdf_input_b3d_timestep_climate_file = TRIM(adjustl(netcdf_input_b3d_timestep_climate_path)) // TRIM(adjustl(new_input_file))
		print *, "load climtology B3D from file ", trim(adjustl(new_input_file))

		call read_climate_bern3d(deviation_temp, deviation_precip, deviation_P_sun, netcdf_input_b3d_timestep_climate_file , thematrix )
		

		do id=1,ndays,1
			deviation_temp(:,:,id) = (deviation_temp(:,:,id)-b3d_pd_temp(:,:,id))
			deviation_P_sun(:,:,id) = (deviation_P_sun(:,:,id)-b3d_pd_P_sun(:,:,id))

			if(turnoff_dprecip)then
				deviation_precip(:,:,id) = 0.
			else
				deviation_precip(:,:,id) = (deviation_precip(:,:,id)-b3d_pd_precip(:,:,id))*b3d_precip_multiplier
				! calibrate precipitation deviations to elevation
				!where (elevation(:,:) .gt. precipitation_threshold)
                            		!deviation_precip(:,:,id) = deviation_precip(:,:,id) &
                            		!* exp( - precipitation_lapse_rate &
                            		!* (elevation(:,:) - precipitation_threshold) )
                    		!end where
			end if
		end do
		print*,'ave temp rel 1950AD :', sum(deviation_temp(:,:,:))/313./313./96.
	end if




	if(erai_backandforth_climate) then
		! creat string for the next year
		CALL SYSTEM_CLOCK(c1)
		write (new_input_file , '( "ERAinterim_",I4.4,".interp.cdf" )') erai_year

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
	
			inp_temp =   read_climate_long(new_input_file, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_temp_variable)) )
			inp_precip = read_climate_long(new_input_file, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_precip_variable)) )
            inp_dewpT = read_climate_long(new_input_file, nx, ny, &
            TRIM(adjustl(netcdf_input_eraiterim_dewpT_variable)) )
            inp_wind = read_climate_long(new_input_file, nx, ny, &
            TRIM(adjustl(netcdf_input_eraiterim_wind_variable)) )
            inp_lwrd = read_climate_long(new_input_file, nx, ny, &
            TRIM(adjustl(netcdf_input_eraiterim_lwrd_variable)) )

		!	! Remove errorous negative precipitation
		!	where(inp_precip(:,:,:) .lt. 0.)
		!		inp_precip(:,:,:) = 0.
		!	end where


			! reduce to 96 Timesteps
			if(ndays==96) then
				do id=1,ndays,1
					temperature(:,:,id)= inp_temp(:,:,nint(real(id)*3.8))+kelvin ! TODO TODO remove test cooling
					precipitation(:,:,id)= inp_precip(:,:,nint(real(id)*3.8))/3600./24./1000. ! from mm/day to m/sec
					DewpT(:,:,id)= inp_dewpT(:,:,id)
                    wind(:,:,id)= inp_wind(:,:,id)
                    lwrd(:,:,id)= inp_lwrd(:,:,id)

					! Calculate potential temperature (at sea level)
					potential_temperature(:,:,id) = temperature(:,:,id) + (initial_climate_elevation * temperature_lapse_rate)

				end do

				!print*,'test temp',temperature(placex,placey,3)
				initial_climate_precipitation = precipitation

				! load short wave radiation
				inp_temp =   read_climate_long(new_input_file, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_swradboa_variable)) )

				! reduce to 96 Timesteps
				do id=1,ndays,1
					P_sun0(:,:,id)= inp_temp(:,:,nint(real(id)*3.8))
					P_sun(:,:,id) = P_sun0(:,:,id)
				end do
			end if


			if(ndays==365) then
				do id=1,ndays,1
					temperature(:,:,id)= inp_temp(:,:,id)+kelvin ! TODO TODO remove test cooling
					precipitation(:,:,id)= inp_precip(:,:,id)/3600./24./1000. ! from mm/day to m/sec
                    DewpT(:,:,id)= inp_dewpT(:,:,id)
                    wind(:,:,id)= inp_wind(:,:,id)
                    lwrd(:,:,id)= inp_lwrd(:,:,id)
					! Calculate potential temperature (at sea level)
					potential_temperature(:,:,id) = temperature(:,:,id) + (initial_climate_elevation * temperature_lapse_rate)

				end do

				!print*,'test temp',temperature(placex,placey,3)
				initial_climate_precipitation = precipitation

				! load short wave radiation
				inp_temp =   read_climate_long(new_input_file, nx, ny, TRIM(adjustl(netcdf_input_eraiterim_swradboa_variable)) )

				! reduce to 96 Timesteps
				do id=1,ndays,1
					P_sun0(:,:,id)= inp_temp(:,:,id)
					P_sun(:,:,id) = P_sun0(:,:,id)
				end do
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
	if ( (mod(myyear(it), calc_rate_climate) == 0) .or. (erai_backandforth_climate) )then  		!TODO set back to 50 years and remove that heating
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
		temperature(:,:,:) = potential_temperature(:,:,:)
		if(active_temperature_lapsrate) then
			do id = 1,ndays,1
				temperature(:,:,id) = potential_temperature(:,:,id)+ deviation_temp(:,:,id) - &
							(elevation * temperature_lapse_rate)
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
                if(ignore_himalaya) then
                    call ignore_himalaya_accumulation(accumulation, nx, ny)
                end if


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
                if(ignore_himalaya) then
                    call ignore_himalaya_accumulation(surface_mass_balance, nx, ny)
                end if

		where( ( elevation_netcdf(:,:).lt. sea_level + sea_level_offset ) )
			surface_mass_balance(:,:) = min( surface_mass_balance(:,:),0.)
		end where

                accumulation = surface_mass_balance
                surface_mass_balance = surface_mass_balance/seconds_per_year
            end if
		
        end if ! mod 20 year if for elevation feedbacks
          


	! write climate forcing data to netcdf
	!-------------------------------------
        if((mod(myyear(it), monthly_data_freq)==0).and.(save_climate_forcing)) then !(mod(myyear(it), new_climate)==0).and.(input_file_no.lt.1e6).and.
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
			call writeNCDFGridValues(filehandle_netcdf3, id, dev_swrad_varid, real(deviation_P_sun(:,:,id)), ny, nx)
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
				surface_mass_balance, nc_counter, myyear(it),& 
				sunken_snow, adaptive_timestep, fast_calculation, albedo_dynamic, DewpT, lwrd, wind, elevation) ! topography and S_BOA

		sunken_snow=0.
		if(ignore_himalaya) then
			call ignore_himalaya_accumulation(surface_mass_balance, nx, ny)
		end if
		! surface_mass_balance in m_ice/second
		accumulation = surface_mass_balance*seconds_per_year ! in m_ice/ year
		    	!ablation = sum(snowman,3) ! snow in kg/m2 




!		if ((myyear(it).lt. start_speed_up ).or.(myyear(it).gt. end_speed_up)  )then
!			spinup=.false.
!		    	call get_accumulation_snowman(nx, ny, ndays, n_snowlayer, temperature, precipitation, P_sun, assignment_mask,&
!		        	snowman, lwmass, rho_snow, snow_temp, surface_mass_balance, nc_counter, myyear(it),sunken_snow, spinup  ) ! topography and S_BOA
!			sunken_snow=0
!		    	if(ignore_himalaya) then
!		        	call ignore_himalaya_accumulation(surface_mass_balance, nx, ny)
!		    	end if
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
!		    	if(ignore_himalaya) then
!		        	call ignore_himalaya_accumulation(surface_mass_balance, nx, ny)
!		    	end if
!		    	accumulation = surface_mass_balance ! in m_ice/ year
!		    	!ablation = sum(snowman,3) ! snow in kg/m2
!		    	surface_mass_balance = surface_mass_balance/seconds_per_year ! in m_ice/second
!		end if
        end if
	
        ! pdd accumulation model
        ! the accumulation/ablation i.e. surfacebalance is calculatet for one year
        ! the pdd model smb is then used for 10 years
        ! pdd model
!        if(smb_model==1)then
!            accumulation = temperature_dependent_accumulation(nx, ny, temperature, precipitation, &
!                precipitation_unit, accumulation_daily_temperature_threshold, ndays)
!            ablation = temperature_dependent_ablation(nx, ny, temperature, beta, ndays)
!
!            if(ignore_himalaya) then
!                call ignore_himalaya_accumulation(accumulation, nx, ny)
!            end if
!            surface_mass_balance = (accumulation - ablation)/seconds_per_year
!        end if

         ! fun model
!        if(smb_model==3)then
!            ! insert yearly data
!            surface_mass_balance  = get_troll_accumulation(nx,ny,temperature, precipitation, assignment_mask, P_sun, ndays,&
!                c_w, rho_w, rho_i, L_lh, D_sf, albedo_ice, sigma, eps_air, eps_snow, snow_fall_temperature)
!            !print *, "can you see me? "
!            !print *, 'Year: ', myyear(it)
!            if(ignore_himalaya) then
!                call ignore_himalaya_accumulation(surface_mass_balance, nx, ny)
!            end if
!            accumulation = surface_mass_balance
!            surface_mass_balance = surface_mass_balance/seconds_per_year
!        end if


                print*,'end snow part'

        !------------------------smb part ends here-------------------------------------

        !=====================================
        ! ICE STARTS HERE
        !=====================================

        ! compute diffusivity
        ! first loop over x-axes
        do ix=1,nx,1
            ! then loop over y-axes
            do iy=1,ny,1
                ! South border
                if (iy == 1) then
                    ! Special Case, the corner in the east and west
                    ! north west corner
                    if (ix == 1) then
                        hgradient = abs( ( ((elevation(ix,iy) - elevation(ix+1,iy)) / real(dx,8))**2. ) &
                                    + ( ((elevation(ix,iy+1) - elevation(ix,iy))/real(dy,8))**2. ) )
                    ! north east corner
                    else if(ix == nx) then
                        hgradient = abs( ( ((elevation(ix-1,iy) - elevation(ix,iy)) / real(dx,8))**2. ) &
                                    + ( ((elevation(ix,iy+1) - elevation(ix,iy))/real(dy,8))**2 ) )
                    ! South border
                    else
                        hgradient = abs( ( ((elevation(ix-1,iy) - elevation(ix+1,iy)) / (2. * dx))**2. ) &
                                    + ( ((elevation(ix,iy+1) - elevation(ix,iy))/dy)**2. ) )
                    endif
                ! North border
                else if (iy == ny) then
                    ! north west corner
                    if (ix == 1) then
                        hgradient = abs( ( ((elevation(ix,iy) - elevation(ix+1,iy)) / dx)**2 ) &
                                    + ( ((elevation(ix,iy) - elevation(ix,iy-1))/dy)**2 ) )
                    ! north east corner
                    else if(ix == nx) then
                        hgradient = abs( ( ((elevation(ix-1,iy) - elevation(ix,iy)) / dx)**2 ) &
                                    + ( ((elevation(ix,iy) - elevation(ix,iy-1))/dy)**2 ) )
                    ! North Border
                    else
                        hgradient = abs( ( ((elevation(ix-1,iy) - elevation(ix+1,iy)) / (2*dx))**2 ) &
                                    + ( ((elevation(ix,iy) - elevation(ix,iy-1))/dy)**2 ) )
                    endif
                ! East Border
                else if (ix == 1) then
                    hgradient = abs( ( ((elevation(ix,iy) - elevation(ix+1,iy)) / dx)**2 ) &
                                + ( ((elevation(ix,iy+1) - elevation(ix,iy-1))/(2*dy))**2 ) )
                ! West Border
                else if (ix == nx) then
                    ! wrong
                    hgradient = abs( ( ((elevation(ix-1,iy) - elevation(ix,iy)) / dx)**2 ) &
                                + ( ((elevation(ix,iy+1) - elevation(ix,iy-1))/(2*dy))**2 ) )
                ! the normal case in the field
                else
                    ! Discretization with full steps (x-1 and x+1)
                    ! Do the check with the assignment_mask only in the normal case, otherwise it will get to confising in the code. And most of the cells are anyway not at the border.
                    ! Calculate the gradient only, if there is at least one grid point with a zero in the 3x3 grid around the center.
                    if(minval(assignment_mask((ix-1):(ix+1),(iy-1):(iy+1))) == 0) then
                        hgradient = abs( ( ((elevation(ix-1,iy) - elevation(ix+1,iy)) / (2*dx))**2 ) &
                                    + ( ((elevation(ix,iy+1) - elevation(ix,iy-1))/(2*dy))**2 ) )

                        ! If the hgradient is zero or very close to it (one point island, cause the ocean has elevation 0), the ice gets accumulated to infinity until the model gives an unstable integration.
                        ! In this case, the gradient is calculated with the elevation of the sea floor.
                        ! This may not be very realistic (since the sea floor does not interact with the ice), but it only affects small islands.
                        if (abs(hgradient) .lt. 1E-014) then
                            hgradient = abs( ( ((elevation_netcdf(ix-1,iy) - elevation_netcdf(ix+1,iy)) / (2*dx))**2 ) &
                                        + ( ((elevation_netcdf(ix,iy+1) - elevation_netcdf(ix,iy-1))/(2*dy))**2 ) )
                        end if
                   else
                       hgradient = 0
                       ! DEBUG
                       if(debug > 6 ) then
                           print *,it,' - hgradient not calculated at grid point: ', ix, iy, &
                                      ', assignment_mask: ', assignment_mask(ix,iy)
                       endif
                   end if
                end if

                ! Diffusivity calculation at point ix, iy
                ! Eismint (Huybrechts, 1996): equation 3 (first part, without H gradient: This part is multiplied later)
                if(hgradient .ne. 0d0) then ! Check if hgradient is calculated before
                    D(ix,iy) = ((2 * (ice_flux_adjustments * A_Eismint) * (rho_ice * g) ** n_Eismint)/(n_Eismint + 2)) &
                                * (ice_thickness(ix,iy)**(n_Eismint + 2)) &
                                * (hgradient**((n_Eismint - 1)/2d0))
                else
                    D(ix, iy) = 0d0
                endif
            end do ! loop over y-axes
        end do ! loop over x-axes

!        ! Write NetCDF with Diffusivity
!        if (write_netcdf) then
!            ! only write in specific timesteps (netcdf_timesteps) to the netCDF File
!            if (mod(myyear(it), netcdf_timesteps) == 0 .and. (last_netcdf_year .lt. myyear(it))) then
!                call writeNCDFGridValues(filehandle_netcdf, (myyear(it)/netcdf_timesteps) + 1, &
!                        diffusivity_varid, D(1:nx, 1:ny) * seconds_per_year, ny, nx)
!            end if
!        end if
        
        ! Compute Ice flux
        ! first loop over x-axes
        do ix=1,nx,1
            ! then loop over y-axes
            do iy=1,ny,1
                ! Default, all are 0, which is the value at the border and corners.
                ! Set them later, if the position is away from the border
                FS(ix,iy) = 0d0
                FE(ix,iy) = 0d0
                FN(ix,iy) = 0d0
                FW(ix,iy) = 0d0

                ! TODO: Flow into north can be calculated, if iy = 1! The same holds true for FS, FE, FW with other boundary conditions.
                ! North/South Direction
                if (iy > 1 .and. iy < ny) then
                  FN(ix,iy) = 0.5 * ((D(ix,iy+1) + D(ix,iy)) * ((elevation(ix,iy+1) - elevation(ix,iy))/real(dy,8) ))
                  ! FN must be less than a quarter of the local ice:
                  if (FN(ix,iy) .gt. 0d0) then ! Positiv = Into the grid point
                    FN(ix,iy) = min((ice_thickness(ix,iy+1)*dy/dt)*.25,FN(ix,iy))
                  else ! Negativ = Away from the gridpoint
                  FN(ix,iy) = max(-(ice_thickness(ix,iy)*dy/dt)*.25,FN(ix,iy))
                  end if
                end if
                if (ix > 1 .and. ix < nx) then
                  FE(ix,iy) = 0.5 * ((D(ix + 1,iy) + D(ix,iy)) * ((elevation(ix + 1,iy) - elevation(ix,iy))/real(dx,8) ))
                  ! FE must be less than a quarter of the local ice:
                  if(FE(ix,iy) .gt. 0d0) then
                    FE(ix,iy) = min((ice_thickness(ix + 1,iy)*dx/dt)*.25,FE(ix,iy))
                  else ! Negativ = Away from the grid point
                    FE(ix,iy) = max(-(ice_thickness(ix,iy)*dx/dt)*.25,FE(ix,iy))
                  end if
                end if
            
            end do ! loop over y-axes
        end do ! loop over x-axes

        ! ABo: remove fluxes where assignment_mask != 0:
        do ix=2,nx-1,1
          do iy=2,ny-1,1
            FS(ix,iy) = FN(ix,iy-1)
            FW(ix,iy) = FE(ix-1,iy)
            if(minval(assignment_mask(ix,iy:(iy+1))) .ne. 0) then
              FN(ix,iy) = 0d0;
            end if
            if(minval(assignment_mask(ix,(iy-1):iy)) .ne. 0) then
              FS(ix,iy) = 0d0;
            end if
            if(minval(assignment_mask(ix:(ix+1),iy)) .ne. 0) then
              FE(ix,iy) = 0d0;
            end if
            if(minval(assignment_mask((ix-1):ix,iy)) .ne. 0) then
              FW(ix,iy) = 0d0;
            end if
          end do
        end do

        ! compute mass flow (discharge)
        ! This is only done for analytical reasons (during the eismint test), to check the values in the netcdf file.
        ! The values are not needed for later calculatieons.
        !if (write_netcdf) then
        !    ! only write in specific timesteps (netcdf_timesteps) to the netCDF File
        !    if (mod(myyear(it), netcdf_timesteps) == 0 .and. (last_netcdf_year .lt. myyear(it))) then
        !        discharge_x(:, :) = FE(:,:) * seconds_per_year !(((real(FE(ix,iy)) - FW(ix,iy))/real(dx)) * seconds_per_year)
        !        discharge_y(:, :) = FN(:,:) * seconds_per_year ! (((real(FN(ix,iy)) - FS(ix,iy))/real(dy)) * seconds_per_year)
        !        ! Write it to the NetCDF file
        !        call writeNCDFGridValues(filehandle_netcdf, (myyear(it)/netcdf_timesteps) + 1, &
        !                                 discharge_x_varid, discharge_x(1:nx, 1:ny) , ny, nx)
        !        call writeNCDFGridValues(filehandle_netcdf, (myyear(it)/netcdf_timesteps) + 1, &
        !                                 discharge_y_varid, discharge_y(1:nx, 1:ny) , ny, nx)
        !    end if
        !end if
        
        ! Hysteresis temperature transition
        if(active_hysteresis .and. hysteresis_temperature_jump .and. &
                                   (myyear(it) .ge. hysteresis_temperature_jump_year) .and. &
                                   (myyear(it) .le. hysteresis_temperature_jump_year + &
                                   (abs(hysteresis_temperature_jump_temperature) * 1000) )) then
            ! add the temperature every day
            do id = 1,365,1
                potential_temperature(:,:,id) = potential_temperature(:,:,id) + &
                                                ! Damit das Vorzeichen stimmt
                                                ( (abs(hysteresis_temperature_jump_temperature)/&
                                                hysteresis_temperature_jump_temperature)/1000)
            end do
            print *,'**** ',it,' *** Hysteresis temperature jump: ', ( (abs(hysteresis_temperature_jump_temperature)/&
                                                                        hysteresis_temperature_jump_temperature)/1000)
        end if


        ! here used to be the smb calculus
        

        imb_tmp = 0d0
        imb_calving = 0d0
        imb_accumulation = 0d0
        imb_ablation = 0d0
        imb_cellcount = 0d0
        ! Compute new ice thickness
        ! first loop over x-axes
        do ix=1,nx,1
            ! then loop over y-axes
            do iy=1,ny,1
                
                ! IMB: Get flux into water.
                ! 0 = ice, 1 = water, 2 = no ice, 3 = unstable integration
                if ((assignment_mask(ix, iy) .eq. 1) .and. store_integrated_mass_balance) then
                    imb_calving = imb_calving - &
                                ( ( (((FN(ix,iy) - FS(ix,iy))/dy) * dt) + (((FE(ix,iy) - FW(ix,iy))/dx) * dt) &
                                  )  * dx * dy * (rho_ice/rho_water) &
                                )
                endif
                
                ! Height
                ! Only do it where there is no water (0 = ice, 1 = water, 2 = no ice, 3 = unstable integration)
                if ((assignment_mask(ix, iy) .ne. 1)  ) then	! TODO: also claculate there where ocean shallower than shelv_depth
                    
!.or.(bedrock_netcdf(ix, iy)+sea_level.gt. shelv_depth)

                    ! Integrated Mass Balance
                    if (store_integrated_mass_balance) then
                        imb_accumulation = imb_accumulation + max(0d0,(surface_mass_balance(ix,iy) &
                                                                   * dx * dy * dt * (rho_ice/rho_water)))
                        ! Ablation, not more ice than available
                        if ( (surface_mass_balance(ix,iy) .le. 0d0)) then
                            ! Ablation, not more ice than available
                            !imb_ablation = imb_ablation +   (max(-(ice_thickness(ix,iy) &
                            !                                    ! Flux in North/South direction
                            !                                    + (((FN(ix,iy) - FS(ix,iy))/dy) * dt) &
                            !                                    ! Flux in East/West direction
                            !                                    + (((FE(ix,iy) - FW(ix,iy))/dx) * dt) ), &
                            !                                    (surface_mass_balance(ix,iy) * dt) ) &
                            !                                 * dx * dy * (rho_ice/rho_water) )
                            ! Available ice
                            imb_tmp = max(0d0, ice_thickness(ix,iy) + (((FN(ix,iy) - FS(ix,iy))/dy) * dt) &
                                                           + (((FE(ix,iy) - FW(ix,iy))/dx) * dt) )
                            ! ABo: limit negative smb to what is locally available:
                            surface_mass_balance(ix,iy) = max(surface_mass_balance(ix,iy),-imb_tmp/dt)
                            imb_ablation = imb_ablation +   (surface_mass_balance(ix,iy) * dt &
                                                                * dx * dy * (rho_ice/rho_water) )
                            
                            !imb_cellcount = imb_cellcount + (imb_tmp  * dx * dy * (rho_ice/rho_water))
                            !imb_cellcount = imb_cellcount + (((FN(ix,iy) - FS(ix,iy))/dy) * dt) &
                            !                              + (((FE(ix,iy) - FW(ix,iy))/dx) * dt)
                            
                            ! Debug
                            !if (imb_tmp .gt. 0) then
                            !    print *,it,' - Max ablation at', ix, iy, ': ', ice_thickness(ix,iy) , '(thick) vs ', &
                            !                                    -imb_tmp, '(av) vs ', (surface_mass_balance(ix,iy) * dt), &
                            !                                    ' (smb) = ', (max( (surface_mass_balance(ix,iy) * dt), -imb_tmp )) &
                            !                                    , assignment_mask(ix,iy)
                            !endif
                        end if
                    end if
                    
                    ! DEBUG
                    imb_tmp = ice_thickness(ix,iy)
                    
                    ! Real Calculations
		    if(ice_dynamics_on)then
		            ice_thickness(ix,iy) = max(0d0, (ice_thickness(ix,iy) &
		                      ! Flux in North/South direction
		                    + (((FN(ix,iy) - FS(ix,iy))/dy) * dt) &
		                    ! Flux in East/West direction
		                    + (((FE(ix,iy) - FW(ix,iy))/dx) * dt) &
		                    ! Accumulation
		                    + (surface_mass_balance(ix,iy) * dt)) )
                    end if

                    ! DEBUG
                    !imb_cellcount = imb_cellcount + ((ice_thickness(ix,iy) - imb_tmp) * dx * dy * (rho_ice/rho_water))
                    ! DEBUG: Only ablation
                    !imb_cellcount = imb_cellcount - max(real(0,8),(surface_mass_balance(ix,iy) &
                    !                                    * dx * dy * dt * (rho_ice/rho_water)))
                    
                    ! If the new ice_thickness is 0, set the assignment_mask to 2
                    ! 0 = ice, 1 = water, 2 = no ice, 3 = unstable integration
                    if( ice_thickness(ix,iy) .le. 0d0) then
                        ! DEBUG
                        if(debug > 6 ) then
                            print *,it,' - Assigned value 2 at grid point: ', ix, iy, ', height: ', ice_thickness(ix,iy)
                        endif
                        assignment_mask(ix,iy) = 2
			!if (bedrock_netcdf(ix,iy) .lt. sea_level) then
                        !	assignment_mask(ix,iy) = 1
			!end if
                    else
                        !if (assignment_mask(ix,iy) .ne. 0) then
                        !    print *,it,' - Assigned value 0 at grid point: ', ix, iy, ', height: ', ice_thickness(ix,iy), &
                        !               ', elevation: ', Bedrock_netcdf(ix,iy), 'assignement_mask: ', assignment_mask(ix,iy)
                        !endif
                        assignment_mask(ix,iy) = 0
                    endif


                endif
            end do ! loop over y-axes
        end do ! loop over x-axes

        ! Integrated Mass Balance Output
        if (store_integrated_mass_balance) then
            ! Mass Change
            imb_totalmass_change = (sum(ice_thickness) * dx * dy) * (rho_ice/rho_water) - imb_totalmass
            ! total ice mass
            imb_totalmass =        (sum(ice_thickness) * dx * dy) * (rho_ice/rho_water)
            
            ! Open questions:
            !    Himalaya,
            !    Bedrock wird nur alle 10 Jahre berechnet
            !    Sealevel wird nur alle 50yr ausgerechnet

            write (imb_filehandle,*) it, ";",imb_accumulation,";",imb_ablation,";",imb_calving &
                                        ,";",imb_isostaticmelt, ";",imb_totalmass,";",imb_totalmass_change &
                                        ,';', (imb_accumulation + imb_ablation + imb_calving &
                                              + imb_isostaticmelt - imb_totalmass_change) &
                                        , ';', (sea_level + sea_level_offset)
                                        !, ';', imb_cellcount
                                        !, ';', imb_cellcount - imb_ablation
        end if
        
        ! Compute bedrock sinking
        if ((active_bedrock).and.(ice_dynamics_on)) then
            ! Bedrock calculation with the ELRA method
            if (active_elra) then
                ! Octave 1D Code from Andreas
                ! w = zeros(size(B));
                ! for jx=1:nx
                !     for ix=1:nx
                !         delta_x = abs(x(ix)-x(jx));% distance of local ice column from target grid box
                !         Vi = rho_i*g*H(ix)*dx;% mass of local ice column
                !         w(jx) = w(jx) + 1/(8*K1) * Vi * alpha_br^3 * exp(-delta_x/alpha_br).*(cos(delta_x/alpha_br)+sin(delta_x/alpha_br));
                !     end
                !     B(jx) = B(jx) - (B(jx)+w(jx))*(dt/tstar);
                ! end
                
                ! first loop over x-axes
                do ix=1,nx,1
                    ! then loop over y-axes
                    do iy=1,ny,1
                        
                        ! ABo, ELRA: calculate equilibrium sinking depth elra_wss (positive downward):
                        ! every 10 years to speed things up:
                        if (mod(myyear(it), 10) == 0) then
                        ! Grid cell in the center: ix, iy
                        
                            elra_wss(ix,iy) = 0d0
                            ! limit internal loops to +-24 because this includes 7*elra_L_r, more efficient:
                            do jx=max(1,ix-24),min(nx,ix+24),1
                                do jy=max(1,iy-24),min(ny,iy+24),1
                                    elra_kei_value = 0d0
                                    elra_delta_distance = (abs(x(ix) - x(jx)) ** 2 + abs(y(iy) - y(jy)) ** 2)**.5/elra_L_r
                                    if (elra_delta_distance .le. (7d0)) then
                                        elra_f_0 = ice_thickness(jx,jy)*rho_ice_g_dx_dy
                                        ! find value of kei for given distance:
                                        do i=1,1059,1
                                            if (kei(i,1) .gt. elra_delta_distance) then
                                                ! weighted average:
                                                elra_kei_value = (kei(i,2)*(kei(i,1)-elra_delta_distance) + kei(i-1,2) &
                                                                 *(elra_delta_distance-kei(i-1,1)))/(kei(i,1)-kei(i-1,1))
                                                exit
                                            end if
                                        end do
                                        elra_wss(ix,iy) = elra_wss(ix,iy) - elra_factor * elra_f_0 * elra_kei_value
                                    end if
                                end do
                            end do
                        end if ! Mod 10
                        
                        ! ABo, ELRA:
                        Bedrock(ix,iy) = bedrock_netcdf(ix,iy) - (elra_wss(ix,iy) &
                                         + (bedrock_netcdf(ix,iy) - Bedrock_Initial(ix,iy)) ) * (real(dt,8)/elra_tstar)
                        ! bedrock_netcdf = Real bedrock, with negative values in the ocean (Bedrock is set to the sealevel on the ocean).
                        ! The variable "Bedrock" is only used for calculation of hgradient. Otherwise bedrock_netcdf is used.
                        ! For understanding reasons, it would be usefull to set the "Bedrock" at a water grid cell right here to sealevel.
                        bedrock_netcdf(ix, iy) = Bedrock(ix, iy)
                        ! Elevation for netcdf output, but with negative values on the water
                        elevation_netcdf(ix, iy) = Bedrock(ix, iy) + ice_thickness(ix, iy)
                        ! if the bedrock is read from a netcdf file, set the ice height on the water to zero
                        if((assignment_mask(ix, iy) == 1)) then
                            if ((assignment_mask(ix,iy) .eq. 1) .and. (ice_thickness(ix,iy) .gt. 0)) then
                                imb_isostaticmelt = imb_isostaticmelt - ice_thickness(ix,iy)
                            end if
                            ice_thickness(ix,iy) = 0d0
                        endif
                        
                        
                    end do
                end do
                
            ! Or Bedrock calculations with local changes (1/3 of ice with a relaxation time)
            else
                ! first loop over x-axes
                do ix=1,nx,1
                    ! then loop over y-axes
                    do iy=1,ny,1
                        Bedrock(ix,iy) = bedrock_netcdf(ix,iy) - ((real(1,8)/3 * ice_thickness(ix,iy)) &
                                          + (bedrock_netcdf(ix,iy) - Bedrock_Initial(ix,iy)) ) * (real(dt,8)/tstar)
                        ! bedrock_netcdf = Real bedrock, with negative values in the ocean (Bedrock is set to the sealevel on the ocean).
                        ! The variable "Bedrock" is only used for calculation of hgradient. Otherwise bedrock_netcdf is used.
                        ! For understanding reasons, it would be usefull to set the "Bedrock" at a water grid cell right here to sealevel.
                        bedrock_netcdf(ix, iy) = Bedrock(ix, iy)
                        ! Elevation for netcdf output, but with negative values on the water
                        elevation_netcdf(ix, iy) = Bedrock(ix, iy) + ice_thickness(ix, iy)
                        ! if the bedrock is read from a netcdf file, set the ice height on the water to zero
                        if((assignment_mask(ix, iy) == 1)) then
                            if ((assignment_mask(ix,iy) .eq. 1) .and. (ice_thickness(ix,iy) .gt. 0d0)) then
                                imb_isostaticmelt = imb_isostaticmelt - ice_thickness(ix,iy)
                            end if
                            ice_thickness(ix,iy) = 0d0
                        end if
                    end do
                end do
            end if ! Switch: Elra and other Bedrock model
        endif ! Switch active bedrock
        
        ! Write also the bedrock to the NetCDF file, if bedrock sinking is not enabled


	! Update assignment_mask, bedrock and elevation variable


	! calculate new sealevel every 50 years
	if (mod(myyear(it), calc_rate_climate) == 0) then  				!TODO set back to 50 years
		if(adjust_sea_level) then
			sea_level = get_sea_level(ice_thickness, nx, ny, dx, dy, ocean_area)
			! TODO: Do not take the initial bedrock, use the current one but forget the ice above it (and it should still increase, even if it is below the water
			! But this could get to complicated, cause we do not want any sea in the middle of america.
			assignment_mask = read_watermask(Bedrock_initial, nx, ny, (sea_level + sea_level_offset), &
							! (use bedrock_netcdf, cause the bedrock is not equal to the sea level)
							 Bedrock_netcdf, ice_thickness)
			
			! integrated mass balance
			! Die IF-Abfrage kann nicht innerhalb vom where integriert werden, daher einzeln
			if (store_integrated_mass_balance) then
				do ix=1,nx,1
					do iy=1,ny,1
						if ((assignment_mask(ix,iy) .eq. 1) .and. (ice_thickness(ix,iy) .gt. 0d0)) then
							imb_isostaticmelt = imb_isostaticmelt - ice_thickness(ix,iy)
						end if
					end do
				end do
			end if
			
			! Set the sea level in the netcdf file for every water point
			! 0 = ice, 1 = water, 2 = no ice, 3 = unstable integration
			elevation = elevation_netcdf
			bedrock = bedrock_netcdf
			where(assignment_mask(:,:) .eq. 1)
				Bedrock(:,:) = sea_level + sea_level_offset
				elevation(:,:) = sea_level + sea_level_offset
				! integrated mass balance calculated before
				ice_thickness(:,:) = 0d0
				elevation_netcdf(:,:) = bedrock_netcdf(:,:)
			end where
			sunken_snow=0.
			do ix=1,nx,1
				do iy=1,ny,1
			    		! reset snow where there is sea
			    		if((sum(snowman(ix,iy,:)).gt. 0.).and.(assignment_mask(ix,iy) .eq. 1)) then
						sunken_snow = sunken_snow + sum(snowman(ix,iy,:)) + sum(lwmass(ix,iy,:))
						do ii=1,n_snowlayer,1
			    				snowman(ix,iy,ii)=0.
							lwmass(ix,iy,ii)=0.
							rho_snow(ix,iy,ii)=rho_s
							snow_temp(ix,iy,ii)=0.
						end do
				    	end if
				end do
			end do
		end if ! end adjust_sea_level

	else ! update assignment_mask, bedrock and elevation variables
		assignment_mask = read_watermask(Bedrock_initial, nx, ny, (sea_level + sea_level_offset), &
							! (use bedrock_netcdf, cause the bedrock is not equal to the sea level)
							 Bedrock_netcdf, ice_thickness)

		elevation = elevation_netcdf
		bedrock = bedrock_netcdf
		where(assignment_mask(:,:) .eq. 1)
			Bedrock(:,:) = sea_level + sea_level_offset
			elevation(:,:) = sea_level + sea_level_offset
			! integrated mass balance calculated before
			ice_thickness(:,:) = 0d0
			elevation_netcdf(:,:) = bedrock_netcdf(:,:)
		end where

	end if






        ! update prognostic variable
        ! This is done for the output some lines before
        !elevation = Bedrock + ice_thickness


        ! Time series of diagnostics
        ! (max(1,size(ice_thickness)))) = Amount of Grid Boxes
        ! real((real(L,8) * dx * W * dy),8) ) = square meters
        H_ts(it) = ( sum(Ice_thickness))/real((real(L,8) * dx *real( W,8) * dy )) ! Average height at 1m^2: mean(H)


        ! //////////////////////// ICE PHYSICS ENDS HERE //////////////////////


        ! netcdf stuff of icesheets
        !==========================

        ! Write NetCDF with Diffusivity

        if ((write_netcdf).and. (mod(myyear(it), netcdf_timesteps) == 0 )) then
            ! only write in specific timesteps (netcdf_timesteps) to the netCDF File
            call writeNCDFGridValues(filehandle_netcdf, nc_counter, &
                diffusivity_varid, real(D(1:nx, 1:ny) * seconds_per_year), ny, nx)

            call writeNCDFGridValues(filehandle_netcdf, nc_counter, &
                acc_varid, real(surface_mass_balance(1:nx, 1:ny) * seconds_per_year) , ny, nx)

            !call writeNCDFGridValues(filehandle_netcdf, (myyear(it)/netcdf_timesteps) + 1, &
            !    abl_varid, ablation(1:nx, 1:ny) , ny, nx)

            call writeNCDFGridValues(filehandle_netcdf, nc_counter, &
                & bedrock_varid, real(bedrock_netcdf(1:nx,1:ny)), ny, nx)

            call writeNCDFGridValues(filehandle_netcdf, nc_counter, &
                height_varid,real(elevation_netcdf(1:nx,1:ny)), ny, nx)
            ! Write assignment mask to netCDF
            call writeNCDFGridIntegerValues(filehandle_netcdf, nc_counter, &
                assignent_mask_varid,assignment_mask(1:nx,1:ny), ny, nx)
            ! this is the last "official" netCDF access in this loop. Assign the year to the last_netcdf_year variable.


            discharge_x(:, :) = FE(:,:) * seconds_per_year !(((real(FE(ix,iy)) - FW(ix,iy))/real(dx)) * seconds_per_year)
            discharge_y(:, :) = FN(:,:) * seconds_per_year ! (((real(FN(ix,iy)) - FS(ix,iy))/real(dy)) * seconds_per_year)
            ! Write it to the NetCDF file
            call writeNCDFGridValues(filehandle_netcdf, nc_counter, &
                discharge_x_varid, real(discharge_x(1:nx, 1:ny)) , ny, nx)
            call writeNCDFGridValues(filehandle_netcdf, nc_counter, &
                discharge_y_varid, real(discharge_y(1:nx, 1:ny)) , ny, nx)

            ! rise counter for next nc entry
            nc_counter = nc_counter + 1
            last_netcdf_year = myyear(it)
        end if



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



        !        if ((mod(myyear(it), netcdf_timesteps) == 0 ) .and. (write_netcdf)) then
        !            nc_counter = nc_counter + 1
        !        end if

        ! CHECK IF ABORT IS NECESSARY
        !============================
        ! Integration check: H-ts(it) = mean ice thickness at 1m^2
        if (H_ts(it) > 10000 ) then ! 10000m of ice on every 1m^2
            print *, 'unstable integration!'
            print *, '---------------------'
            print *, 'Year: ', myyear(it), it
            print *, 'H_ts(it): ', H_ts(it)
            print *, 'Sea Level: ', sea_level + sea_level_offset

            ! Get unstable points and assign it in the assignment_mask to value 3
            ! x-axes
            do ix=1,nx,1
                ! then loop over y-axes
                do iy=1,ny,1
                    if (ice_thickness(ix,iy) > 6000) then
                        assignment_mask(ix:ix,iy:iy) = 3
                        if(debug > 3) then
                            print *, 'Position [x,y]: ', ix, iy
                            print *, 'Ice Thickness: ', ice_thickness(ix,iy)
                            print *, 'SMB: ', surface_mass_balance(ix,iy)
                            print *, 'FN: ', FN(ix,iy)
                            if (abs(FN(ix,iy)) .gt. 3) then
                                print *, 'D (iy + 1, iy):', D(ix,iy + 1), D(ix, iy)
                                print *, 'Elevation (iy + 1, iy):', elevation(ix,iy + 1), elevation(ix, iy)
                                print *, 'Ice Thickness (iy + 1, iy):', ice_thickness(ix,iy + 1), ice_thickness(ix, iy)
                                print *, 'Bedrock_initial (iy + 1, iy):', Bedrock_initial(ix,iy + 1), Bedrock_initial(ix, iy)
                            end if
                            print *, 'FE: ', FE(ix,iy)
                            if (abs(FE(ix,iy)) .gt. 3) then
                                print *, 'D (ix + 1, ix):', D(ix + 1,iy), D(ix, iy)
                                print *, 'Elevation (ix + 1, ix):', elevation(ix + 1,iy), elevation(ix, iy)
                                print *, 'Ice Thickness (ix + 1, ix):', ice_thickness(ix + 1,iy), ice_thickness(ix, iy)
                                print *, 'Bedrock_initial (ix + 1, ix):', Bedrock_initial(ix + 1,iy), Bedrock_initial(ix, iy)
                            end if
                            print *, 'FS: ', FS(ix,iy)
                            if (abs(FS(ix,iy)) .gt. 3) then
                                print *, 'D (iy - 1, iy):', D(ix,iy - 1), D(ix, iy)
                                print *, 'Elevation (iy - 1, iy):', elevation(ix,iy - 1), elevation(ix, iy)
                                print *, 'Ice Thickness (iy - 1, iy):', ice_thickness(ix,iy - 1), ice_thickness(ix, iy)
                                print *, 'Bedrock_initial (iy - 1, iy):', Bedrock_initial(ix,iy - 1), Bedrock_initial(ix, iy)
                            end if
                            print *, 'FW: ', FW(ix,iy)
                            if (abs(FW(ix,iy)) .gt. 3) then
                                print *, 'D (ix - 1, ix):', D(ix - 1,iy), D(ix, iy)
                                print *, 'Elevation (ix - 1, ix):', elevation(ix - 1,iy), elevation(ix, iy)
                                print *, 'Ice Thickness (ix - 1, ix):', ice_thickness(ix - 1,iy), ice_thickness(ix, iy)
                                print *, 'Bedrock_initial (ix - 1, ix):', Bedrock_initial(ix - 1,iy), Bedrock_initial(ix, iy)
                            end if
                            print *, '--------------------------'
                        endif ! Debug > 3
                    endif
                end do ! y axes
            end do ! x axes

            ! Write last values to the NetCDF file and and close it.
            if (write_netcdf) then
                ! Calculate/Refresh discharge
                discharge_x(:, :) = FE(:,:) * seconds_per_year
                discharge_y(:, :) = FN(:,:) * seconds_per_year

                ! Write all variables for later reanalyse
                call writeNCDFGridValues(filehandle_netcdf, nc_counter, &
                    bedrock_varid, real(bedrock_netcdf(1:nx,1:ny)), ny, nx)
                call writeNCDFGridValues(filehandle_netcdf, nc_counter, &
                    height_varid,real(elevation_netcdf(1:nx,1:ny)), ny, nx)
                call writeNCDFGridIntegerValues(filehandle_netcdf, nc_counter, &
                    assignent_mask_varid,assignment_mask(1:nx,1:ny), ny, nx)
                call writeNCDFGridValues(filehandle_netcdf, nc_counter, &
                    diffusivity_varid, real(D(1:nx, 1:ny) * seconds_per_year), ny, nx)
                call writeNCDFGridValues(filehandle_netcdf, nc_counter, &
                    discharge_x_varid, real(discharge_x(1:nx, 1:ny)) , ny, nx)
                call writeNCDFGridValues(filehandle_netcdf, nc_counter, &
                    discharge_y_varid, real(discharge_y(1:nx, 1:ny)) , ny, nx)
                call writeNCDFGridValues(filehandle_netcdf, nc_counter, &
                    acc_varid, real(surface_mass_balance(1:nx, 1:ny) * seconds_per_year), ny, nx)
                !call writeNCDFGridValues(filehandle_netcdf, (myyear(it)/netcdf_timesteps) + 1, &
                !    abl_varid, ablation(1:nx, 1:ny) * seconds_per_year, ny, nx)
                call closeNCDFFile(filehandle_netcdf)
            endif
            ! Runtime information
            !call ETIME(execution_time, runtime)
	    call cpu_time(clock_end)
	    CALL SYSTEM_CLOCK(c_end)
	    runtime = runtime + clock_end - clock_start
            print *, "Runtime [mm:ss]: ", int(runtime/60), ':', mod(int(runtime),60)
            print *, 'Year: ', myyear(it)
            CALL EXIT(42)
            print *, 'TEST'
        end if ! end of huge bug abort if

        it = int(it) + 1
    end do ! Timestep loop

    if(debug > 0 ) then
        print *, "-----------------"
    endif

    ! NETCDF: Close file with ice data
    if (write_netcdf) then
        call closeNCDFFile(filehandle_netcdf)
    endif
    
    ! Close Integrated Mass Balance Output
    if (store_integrated_mass_balance) then
        close(imb_filehandle)
    end if

    ! NETCDF: Close file with snow data
    if (write_netcdf) then
        !call closeNCDFFile(filehandle_netcdf2)
    endif



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
    endif

! =======================================================================================================
! FUNCTIONS
! =======================================================================================================
CONTAINS


    ! --------------------------------------------------------------
    ! SEA LEVEL
    ! --------------------------------------------------------------
    function get_sea_level(ice_thickness, nx, ny, dx, dy, ocean_area)
        ! calculates the sea level,
        ! with the assumption that no ice is equal to a sea level of 0, ice leads to a negative sea level
        ! use the parameter sea_level_offset in variables.f90 to adjust to sea level without ice

        ! input
        integer, intent(in) :: NX
        integer, intent(in) :: NY
        real(kind=8), intent(in) :: DX
        real(kind=8), intent(in) :: DY
        real(kind=8), intent(in) :: ice_thickness(nx, ny)       ! Ice thickness
        real(kind=8), intent(in) :: ocean_area                  ! Area of the ocean
        ! Output
        !real(kind=8), intent(out) :: get_sea_level
        real(kind=8) :: get_sea_level

        ! Sea level
        get_sea_level = -(sum(Ice_thickness) * dx * dy)/ocean_area*rho_ice/rho_w !- (sum(snowman)*dx*dy)/rho_w/ocean_area
                                                                                 ! Volume/Area = height adjusted density

        return
    end function get_sea_level

    ! --------------------------------------------------------------
    ! SURFACE MASS BALANCE FUNCTIONS
    ! --------------------------------------------------------------

    function temperature_dependent_accumulation(nx, ny, temperature_array, precipitation_array, precipitation_unit,&
                                                    daily_temperature_threshold, ndays)
        ! Used to calculate the accumulation over a specific domain given by the precipitation and temperature over the year.

        ! Static parameters
        real(kind=8), parameter :: kelvin = 273.15

        ! input
        real(kind=8), intent(in) :: daily_temperature_threshold
        integer, intent(in) :: ndays
        integer, intent(in) :: NX
        integer, intent(in) :: NY
        real(kind=8), intent(in) :: temperature_array(nx, ny, ndays)     ! Temperature at the level of the ice elevation
        real(kind=8), intent(in) :: precipitation_array(nx, ny, ndays)   ! Precipitation
        integer, intent(in) :: precipitation_unit
        ! Output
        real(kind=8) :: temperature_dependent_accumulation(nx, ny)


        ! accumulation out of the precipitation
        ! Precipication in m
        temperature_dependent_accumulation = 0

        do id=1,ndays,1
            ! http://www.stanford.edu/class/me200c/tutorial_90/07_arrays.html
            where(temperature_array(:,:,id) .lt. (daily_temperature_threshold + kelvin))
                ! Accumulation: All elements in the array, where the temperature is below 0 degree celsius
                temperature_dependent_accumulation = temperature_dependent_accumulation + precipitation_array(:,:, id)*rho_w/rho_ice
            end where
        end do

        ! Accumulation: m/yr to m/s
        if(precipitation_unit == 1) then
            temperature_dependent_accumulation = temperature_dependent_accumulation / (ndays * 24 * 60 * 60)
        endif
        ! m/s divide through the days of the integrated year
        if(precipitation_unit == 2) then
            !temperature_dependent_accumulation = temperature_dependent_accumulation / ndays
            temperature_dependent_accumulation = temperature_dependent_accumulation *dt_firn
        endif

    end function temperature_dependent_accumulation

    function temperature_dependent_ablation(nx, ny, temperature_array, beta, ndays)
        ! Ablation calculated on the princible of the positive degree days (PDD).
        !
        ! nx and ny: Size of the 2d array
        ! temperature_array: (nx, ny, 365): Temperature of the grid point all over the year at the level of the ice elevation in kelvin.

        ! Static parameters
        !integer, parameter :: ndays = 365
        real(kind=8), parameter :: kelvin = 273.15
        real(kind=8), parameter :: daily_temperature_threshold = 0 ! in Celsius
        ! http://www.igsoc.org/journal/59/218/j13J081.pdf
        ! 3 mm * C^–1 * d^–1 = 3.47e-8 m * C^-1 * s^-1 for snow
        ! 8 mm * C^–1 * d^–1 = 9.26e-8 m * C^-1 * s^-1 for ice
        ! 3 - 8 mm * C^–1 * d^–1 = 3 / (1000 * 24 * 60 * 60) = 3.5e-8 m * C^-1 * s^-1
        real(kind=8), intent(in) :: beta ! = 6e-8/days

        ! input
        integer, intent(in) :: NX
        integer, intent(in) :: NY
        integer, intent(in) :: ndays
        real(kind=8), intent(in) :: temperature_array(nx, ny, ndays)     ! Temperature at the level of the ice elevation
        ! Output
        real(kind=8) :: temperature_dependent_ablation(nx, ny)

        real(kind=8) :: positive_degree_days(nx, ny)    ! Days with a temperature over 0 degrees

        ! accumulation out of the precipitation
        ! Calculate positive degree days
        positive_degree_days = 0

        do id=1,ndays,1
            ! http://www.stanford.edu/class/me200c/tutorial_90/07_arrays.html
            where(temperature_array(:,:,id) .ge. kelvin + daily_temperature_threshold)
                ! Ablation (PDD): All elements, where the temperature is above the temperature threshhold degree Celsius
                positive_degree_days = positive_degree_days + (temperature_array(:,:, id)-kelvin-daily_temperature_threshold)
            end where
        end do

        ! Calculate accumulation/ablation out of it
        !temperature_dependent_ablation = (beta * positive_degree_days / ndays)*rho_w/rho_ice
        temperature_dependent_ablation = (beta * positive_degree_days)*rho_w/rho_ice*dt_firn
        return
    end function temperature_dependent_ablation



    !    function temperature_dependent_accumulation(nx, ny, temperature_array, precipitation_array, precipitation_unit,&
    !                                                    daily_temperature_threshold)
    !        ! Used to calculate the accumulation over a specific domain given by the precipitation and temperature over the year.
    !
    !        ! Static parameters
    !        integer, parameter :: days = 365
    !        real(kind=8), parameter :: kelvin = 273.15
    !
    !        ! input
    !        real(kind=8), intent(in) :: daily_temperature_threshold
    !        integer, intent(in) :: NX
    !        integer, intent(in) :: NY
    !        real(kind=8), intent(in) :: temperature_array(nx, ny, days)     ! Temperature at the level of the ice elevation
    !        real(kind=8), intent(in) :: precipitation_array(nx, ny, days)   ! Precipitation
    !        integer, intent(in) :: precipitation_unit
    !        ! Output
    !        real(kind=8) :: temperature_dependent_accumulation(nx, ny)
    !
    !
    !        ! accumulation out of the precipitation
    !        ! Precipication in m
    !        temperature_dependent_accumulation = 0
    !
    !        do id=1,days,1
    !            ! http://www.stanford.edu/class/me200c/tutorial_90/07_arrays.html
    !            where(temperature_array(:,:,id) .lt. (daily_temperature_threshold + kelvin))
    !                ! Accumulation: All elements in the array, where the temperature is below 0 degree celsius
    !                temperature_dependent_accumulation = temperature_dependent_accumulation + precipitation_array(:,:, id)
    !            end where
    !        end do
    !
    !        ! Accumulation: m/yr to m/s
    !        if(precipitation_unit == 1) then
    !            temperature_dependent_accumulation = temperature_dependent_accumulation / (days * 24 * 60 * 60)
    !        endif
    !        ! m/s divide through the days of the integrated year
    !        if(precipitation_unit == 2) then
    !            temperature_dependent_accumulation = temperature_dependent_accumulation / days
    !        endif
    !
    !    end function temperature_dependent_accumulation
    !
    !    function temperature_dependent_ablation(nx, ny, temperature_array, beta)
    !        ! Ablation calculated on the princible of the positive degree days (PDD).
    !        !
    !        ! nx and ny: Size of the 2d array
    !        ! temperature_array: (nx, ny, 365): Temperature of the grid point all over the year at the level of the ice elevation in kelvin.
    !
    !        ! Static parameters
    !        integer, parameter :: days = 365
    !        real(kind=8), parameter :: kelvin = 273.15
    !        real(kind=8), parameter :: daily_temperature_threshold = 0 ! in Celsius
    !        ! http://www.igsoc.org/journal/59/218/j13J081.pdf
    !        ! 3 mm * C^–1 * d^–1 = 3.47e-8 m * C^-1 * s^-1 for snow
    !        ! 8 mm * C^–1 * d^–1 = 9.26e-8 m * C^-1 * s^-1 for ice
    !        ! 3 - 8 mm * C^–1 * d^–1 = 3 / (1000 * 24 * 60 * 60) = 3.5e-8 m * C^-1 * s^-1
    !        real(kind=8), intent(in) :: beta ! = 6e-8/days
    !
    !        ! input
    !        integer, intent(in) :: NX
    !        integer, intent(in) :: NY
    !        real(kind=8), intent(in) :: temperature_array(nx, ny, days)     ! Temperature at the level of the ice elevation
    !        ! Output
    !        real(kind=8) :: temperature_dependent_ablation(nx, ny)
    !
    !        real(kind=8) :: positive_degree_days(nx, ny)    ! Days with a temperature over 0 degrees
    !
    !        ! accumulation out of the precipitation
    !        ! Calculate positive degree days
    !        positive_degree_days = 0
    !
    !        do id=1,days,1
    !            ! http://www.stanford.edu/class/me200c/tutorial_90/07_arrays.html
    !            where(temperature_array(:,:,id) .ge. kelvin + daily_temperature_threshold)
    !                ! Ablation (PDD): All elements, where the temperature is above the temperature threshhold degree Celsius
    !                positive_degree_days = positive_degree_days + (temperature_array(:,:, id)-kelvin-daily_temperature_threshold)
    !            end where
    !        end do
    !
    !        ! Calculate accumulation/ablation out of it
    !        temperature_dependent_ablation = (beta * positive_degree_days / days)
    !        return
    !    end function temperature_dependent_ablation
    !


    ! --------------------------------------------------------------
    ! IGNOR HIMALAYA
    ! --------------------------------------------------------------
    subroutine ignore_himalaya_accumulation(accumulation, nx, ny)
        ! Sets the accumultion in the himalaya region to 0.
        ! This can also be used for the ablation

        integer, intent(in) :: NX
        integer, intent(in) :: NY
        real(kind=8), intent(inout) :: accumulation(nx, ny)

        if (nx .eq. 625) then
            accumulation(440:610, 480:620) = 0
        else
            accumulation(220:305, 240:310) = 0
        end if

    end subroutine ignore_himalaya_accumulation



    ! --------------------------------------------------------------
    ! GET SOLAR RADIATION
    ! --------------------------------------------------------------
    subroutine get_solarrad(nx,ny,P)
        ! simulate solar irradiance

        integer, intent(in) :: NX
        integer, intent(in) :: NY
        real(kind=8), intent(inout) :: P(nx,ny,365)

        REAL(kind=8), PARAMETER :: pi = 3.1415927
        integer, parameter :: days = 365
        real(kind=8), dimension(nx,ny) :: latitude
        real(kind=8), dimension(days):: elipticpart
        real(kind=8), dimension(days):: obliquity

        ! get latitude of every grid cell
        do iy=1,ny,1
            do ix=1,nx,1
                latitude(ix,iy) = pi/2 - 2*pi/40000000*sqrt(((real(ix)-0.5)*real(dx)-real(L)/2)**real(2) &
                    +  ((real(iy)-0.5)*real(dy)-real(W)/2)**real(2))
            end do
        end do

        do id=1, days,1
            do iy=1,ny,1
                do ix=1,nx,1
                    elipticpart(id)=1+0.033*cos((-3+real(id) )/365*2*pi ) ! -3.3% till 3.3%  -3
                    obliquity(id)= -23.4*cos(pi+(8+real(id) )/365*2*pi )*pi/180 !23.4 8
                    P(ix,iy,id) =  max((sin(latitude(ix,iy)+obliquity(id))*P_sun(ix,iy,id))*elipticpart(id),real(0))
                end do
            end do
        end do
        !P= max(P,real(0))    ! set negative values to 0 (polar night)
    end subroutine get_solarrad


    ! #####################################################################################
    ! LOAD BERN3D DATA AND INTERPOLATE IT ONTO ICEGRID
    ! #####################################################################################

    subroutine read_climate_bern3d(temperature, precipitation, P_sun, path_input, mat )
        real(kind=8), intent(inout) :: P_sun(nx,ny,ndays)
        real(kind=8), intent(inout) :: precipitation(nx,ny,ndays)
        real(kind=8), intent(inout) :: temperature(nx,ny,ndays)
        real(kind=8), intent(inout) :: mat(nx,ny,interpol_deg*3)
	character(len=*), intent(in) :: path_input

        real(kind=8), dimension(bnx,bny,ndays) :: bern_temp = 0
        real(kind=8), dimension(bnx,bny,ndays) :: bern_swradboa = 0
        real(kind=8), dimension(bnx,bny,ndays) :: bern_albedo = 0
        real(kind=8), dimension(bnx,bny,ndays) :: bern_precip = 0
        !real(kind=8), dimension(bnx,bny)       :: bern_precip_mean = 0
        real(kind=8), dimension(bnx,bny)       :: bern_dprecip = 0
        !real(kind=8), dimension(nx,ny,interpol_deg*3) :: mat = 0


        ! make parameter names and load stuff

        !load bern3d data
        !----------------

        !load temperature
        if(debug > 0) then
            !print *, "*** Read 'TEMPERATURE' from file: ", TRIM(adjustl(path_input))
        end if
        bern_temp = read_climate(path_input, nx, ny, TRIM(adjustl(netcdf_input_b3d_temp_variable)) )

        !load precipitation
        if(debug > 0) then
            !print *, "*** Read 'PRECIPITATION' from file: ", TRIM(adjustl(path_input))
        end if
        bern_precip = read_climate(path_input, nx, ny, TRIM(adjustl(netcdf_input_b3d_precip_variable)) )

        ! load swradiation
        if(debug > 0) then
            !print *, "*** Read '",TRIM(adjustl(netcdf_input_b3d_swradboa_variable)),"' from file:", &
            !    TRIM(adjustl(path_input))
        end if
        bern_swradboa = read_climate(path_input, nx, ny, &
            TRIM(adjustl(netcdf_input_b3d_swradboa_variable)))

        ! load albedo
        if(debug > 0) then
            !print *, "*** Read '",TRIM(adjustl(netcdf_input_b3d_albedo_variable)),"' from file:", &
            !    TRIM(adjustl(path_input))
        end if
        bern_albedo = read_climate(path_input, nx, ny, &
            TRIM(adjustl(netcdf_input_b3d_albedo_variable)))

        ! extract swradboa arriving at boa from energy uptaken by soil
        bern_swradboa = bern_swradboa/(1.-bern_albedo)



        ! load l anual mean precipitation
        !if(debug > 0) then
        !    print *, "*** Read '",TRIM(adjustl(netcdf_input_precip_variable)),"' from file:", &
        !        TRIM(adjustl(netcdf_input_precip_calib))
        !end if
        !bern_precip_mean = read_variable(netcdf_input_precip_calib, nx, ny, &
        !    TRIM(adjustl(netcdf_input_precip_variable)))


        ! calibrate precipitation
        !------------------
        ! this part needs to be modified during glacials

        !bern_dprecip = (bern_precip_mean-sum(bern_precip,3)/real(ndays))

        !do ii=1,ndays,1
        !    bern_precip(:,:,ii) = bern_precip(:,:,ii)+bern_dprecip(:,:)
        !end do

        !where (bern_precip(:,:,:) .lt. 0)
        !    bern_precip(:,:,:)=0
        !end where


        ! load interpolation matrix
        !------------------------
        !call load_mat(mat,interpol_deg)


        ! interpolate on new grid
        !------------------------

        temperature = interpolate(bern_temp,mat)
        precipitation = interpolate(bern_precip,mat)
        P_sun = interpolate(bern_swradboa,mat)




    end subroutine read_climate_bern3d


    ! --------------------------------------------------------------
    ! LOAD INTERPOLATION MATIRX
    ! --------------------------------------------------------------
    subroutine load_mat(mat,a)
        ! input
        integer, intent(in) :: a
        real(kind=8), intent(inout) :: mat(nx,ny,a*3)

        ! this function builds the actuial 3D matrix of the unroled
        ! transformation debris txt.

        ! local stuff
        real(kind=8)  :: wgt = 0
        real(kind=8), dimension(nx*ny*a,5) :: debris
        character(128) :: debris_filename
        integer :: max_rows
        integer :: max_cols = 5
        integer :: row
        integer :: col
        max_rows = int(nx*ny*a)
        write (debris_filename, '( "matrix", I1.1, "p_1.txt" )') a

        open(UNIT=11, FILE = trim(adjustl(matrixpath)) // trim(adjustl(debris_filename)) )

        do row = 1,max_rows,1
            read(11,*) (debris(row,col),col=1,max_cols)
        end do

        close(11)


!        print *, 'debris:',debris(1,:)
!        print *, 'debris:',debris(2,:)
!        print *, 'debris:',debris(30000,:)



        do ii=1,nx*ny,1
            do kk=1, a, 1
                ! store bern coords with their weights
                mat( int(debris(kk+a*(ii-1),1)), int(debris(kk+a*(ii-1),2)), 1+(kk-1)*3) = debris(kk+a*(ii-1),3)
                mat( int(debris(kk+a*(ii-1),1)), int(debris(kk+a*(ii-1),2)), 2+(kk-1)*3) = debris(kk+a*(ii-1),4)
                mat( int(debris(kk+a*(ii-1),1)), int(debris(kk+a*(ii-1),2)), 3+(kk-1)*3) = debris(kk+a*(ii-1),5) ! weight
            end do
            ! normalize weights to 1
            wgt=0

            do kk=1,a,1
                wgt=wgt+ mat(int(debris(kk+a*(ii-1),1)), int(debris(kk+a*(ii-1),2)), 3+(kk-1)*3)
            end do

            do kk=1,a,1
                mat(int(debris(kk+a*(ii-1),1)), int(debris(kk+a*(ii-1),2)), 3+(kk-1)*3) = mat(int(debris(kk+a*(ii-1),1)),&
                    int(debris(kk+a*(ii-1),2)), 3+(kk-1)*3)/wgt
            end do

        end do

!        print *, 'mat:',mat(1,1,:)
!        print *, 'mat:',mat(placex,placey,:)


    end subroutine load_mat

    ! --------------------------------------------------------------
    ! INTERPOLATION FUNCTION
    ! --------------------------------------------------------------
    function interpolate(var_in,mat)

        real(kind=8), intent(in) :: var_in(bnx,bny,ndays)
        real(kind=8), intent(in) :: mat(nx,nx,interpol_deg*3)


        real(kind=8) :: interpolate(nx,ny,ndays)
        integer :: tt = 1
        integer :: xx = 1
        integer :: yy = 1
        integer :: ii = 1
        interpolate = 0

        do tt=1,ndays,1
            do xx=1,nx,1
                do yy=1,ny,1
                    do ii=1,3*interpol_deg,3
                        interpolate(xx,yy,tt) = interpolate(xx,yy,tt) + var_in(int(mat(xx,yy,ii)), &
                                        int(mat(xx,yy,ii+1)),tt)*mat(xx,yy,ii+2)
                    end do
                end do
            end do
        end do
!       print *, 'input ber3d:', var_in(35,12,22)
!       print *, 'interpolate:', interpolate(placex,placey,22)
!       print *, 'mat:', mat(placex,placey,:)

    end function interpolate


    ! --------------------------------------------------------------
    ! Troll function
    ! --------------------------------------------------------------

    function get_troll_accumulation(nx,ny, ndays,temperature_array, precipitation_array, assignment_mask,P_sun )
        ! Used to calculate the accumulation over a specific domain given by the precipitation and temperature over the year.


        ! input
        integer, intent(in) :: NX
        integer, intent(in) :: NY
        integer, intent(in) :: assignment_mask(nx, ny)
        integer, intent(in) :: ndays
        real(kind=8), intent(in) :: temperature_array(nx, ny, ndays)     ! Temperature at the level of the ice elevation
        real(kind=8), intent(in) :: precipitation_array(nx, ny, ndays)   ! Precipitation
        real(kind=8), intent(in) :: P_sun(nx, ny, ndays)   ! sw radiation

        ! Output
        real(kind=8) :: get_troll_accumulation(nx, ny)

        ! local variables
        real(kind=8), dimension(nx,ny) :: accum
        real(kind=8) :: dQ_lh
        real(kind=8) :: dQ_sh
        real(kind=8) :: dQ_sw
        real(kind=8) :: dQ_lw
        real(kind=8) :: dQ_tot
        real(kind=8), dimension(nx,ny) :: abla
        ! change units from m/yr to m/s and from K to C
        !precip =  precipitation_array/(365 * 24 * 60 * 60)
        !temp = temperature_array-kelvin
        accum=0
        abla=0
        dQ_tot=0
        dQ_lw=0
        dQ_sw=0
        dQ_sh=0
        dQ_lh=0


        ! ITERATION STARTS HERE
	do ix=1,nx,1
		do iy=1,ny,1
        		do id=1,ndays,1
			    ! calc snow_temp
		!            where(snowman(:,:,1).gt.0)
		!                snow_temp(:,:,1)=snowq(:,:,1)/c_i/snowman(:,:,1)
		!            end where

			    ! first accumulate
			    if( (temperature_array(ix,iy,id) .lt. kelvin ) .and. (assignment_mask(ix,iy) .NE. 1) )then
				accum(ix,iy) = accum(ix,iy)+precipitation_array(ix,iy,id)*dt_firn*rho_w/rho_ice


			    else
				    dQ_lh = ( temperature_array(ix,iy,id) -kelvin)*precipitation_array(ix,iy, id)*c_w*rho_w
				    dQ_sh = D_sf*( temperature_array(ix,iy,id) -kelvin)
				    dQ_sw = P_sun(ix,iy,id)*(1.-albedo_ice) 
				    dQ_lw = sigma*(eps_air*(temperature_array(ix,iy,id)  )**real(4)-eps_snow*(kelvin)**real(4))
				    dQ_tot = dQ_lw+dQ_sh+dQ_lh+dQ_sw
				    abla(ix,iy) = abla(ix,iy) - max(dQ_tot,real(0))/L_lh*dt_firn/rho_ice  ! use ice density of ice model to get the right height
			    end if
        		end do
		end do
	end do
	!print*, abla(111,111)
        get_troll_accumulation=(accum+abla)  ! return change in height in m ice equivalents


    end function get_troll_accumulation

    ! --------------------------------------------------------------
    ! Subroutine to read short wave radiation of CCSM
    ! --------------------------------------------------------------

    subroutine read_swradboa_ccsm(P_sun, path, var_name)
	real(kind=8), intent(inout) :: P_sun(nx, ny, ndays) 
	character(*), intent(in) :: path
	character(*), intent(in) :: var_name

	real(kind=8),dimension(nx,ny,12) :: metadata
	!print*,'before loading swradboa'
	metadata = read_swrad_func(path, nx, ny, var_name) ! if some error occurs, change nx and ny
	!print*,'after loading swradboa'
	!print*, metadata(111,111,:)
	!print*, minval(metadata), maxval(metadata)
	do mm=1,12,1
		do ii=1,8,1
			P_sun(:,:,(mm-1)*8+ii)=metadata(:,:,mm)
		end do
	end do
	!print*,'after expanding data, shoult exit subroutine now'
   end subroutine read_swradboa_ccsm


    ! --------------------------------------------------------------
    ! Subroutine altitude swrad damping
    ! --------------------------------------------------------------

    subroutine swrad_damping(P_sun,P_sun0, height, swrad_topo)
	real(kind=8), intent(inout) :: P_sun(nx, ny, ndays) 
	real(kind=8), intent(in) :: P_sun0(nx, ny, ndays)
	real(kind=8), intent(in) :: height(nx, ny)
	real(kind=8), intent(in) :: swrad_topo(nx, ny)


	do id=1,ndays,1
		do ix=1,nx,1
			do iy=1,ny,1
				P_sun(ix,iy,id) = P_sun0(ix,iy,id)*( exp(k_extinct * &
						( min(height(ix,iy),swrad_topo(ix,iy)) - swrad_topo(ix,iy) ) ) )
			end do
		end do
	end do
	!print*,'after expanding data, shoult exit subroutine now'
   end subroutine swrad_damping


end program
