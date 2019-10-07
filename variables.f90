! Module with all variables
! Author:       neff@climate.unibe.ch
! Date:         2015/05
! Revision:     1.1
! SVN Info:     $Id: $


! variables for new SMB model and different climate data added by Michael Imhof
! Co-author: 	Michael Imhof
! Mail:		imhof@vaw.baug.ethz.ch (or imhof@climate.unibe.ch )


module variables
	!Tuning Parameters

	real(kind=8), parameter	:: A_Eismint_spam =  3.1709792e-24   
	real(kind=8), parameter :: albedo_snow_new_spam = 0.7677    ! albedo of dry snow 0.75 - 0.85 
	real(kind=8), parameter :: albedo_snow_wet_spam = 0.5989!6      ! albedo of wet snow i.e. snow with 0.5-0.7
	real(kind=8), parameter :: albedo_ice_spam = 0.3564!2            ! albedo of ice 0.3 - 0.4
	real(kind=8), parameter :: D_sf_spam = 16.52          ! sensible heat flux snow-air [W/m2/K]
! 	real(kind=8), parameter :: D_lf_spam = D_sf_spam/1003*0.622*(L_v+L_lh)/10000
	real(kind=8), parameter :: ratio = 1.158
	real(kind=8), parameter :: eps_air_spam = 0.6287!85        ! highly dependent on cloubdcover and watervapor preasure [0.75 0.95]
	real(kind=8), parameter :: max_lwc_spam = 0.1373    ! percentage of empty volume that can be filed with water
	integer, parameter :: maxyears_spam =   500      !int(4e5)

	character(256), parameter :: experiment_description_spam = &
!						'snowbern_const_3.8_climgris'
						'ensemble' ! flag for ubelix script

    !=========================
    ! Flags
    !=========================

    ! Flag, if the values should be written into a netCDF-file
    logical, parameter :: write_netcdf = .false.
    
    ! Flag, if averaged ERA-I climate data should be used 
    logical, parameter :: eraiterim_climate = .false.

    ! Flag, if 30 years of ERA-I climate data should be looped forwards and backwards
    logical, parameter :: erai_backandforth_climate = .true.

    ! Flag, if LGM laurentide topo climate shall be simulated with era iterim and ccsm deviations
	! SET THIS TO FALSE DURING TRANSIENT SIMULATIONS!!!!!!!!
    logical, parameter :: ccsm_lgmlt_climate = .false.

    ! Flag, if ancient climate shall be simulated with bern3d data deviation
    logical, parameter :: transient_ancient_climate = .False.

    ! Flag, switch off transient precipitation adjustments
    logical, parameter :: turnoff_dprecip = .false.

    real(kind=8), parameter :: b3d_precip_multiplier= 1.

    ! Flag, if eemian climate shall be simulated with era iterim and ccsm deviations
	! SET THIS TO FALSE DURING TRANSIENT SIMULATIONS!!!!!!!!
    logical, parameter :: ccsm_eem_climate = .FALSE.  	! make sure erai climate is on, not tested yet

    ! Flag, if eemian climate shall be simulated with era iterim and ccsm deviations
	! SET THIS TO FALSE DURING TRANSIENT SIMULATIONS!!!!!!!!
    logical, parameter :: dev_b3d_eem_climate = .FALSE.

    ! Flag, if LGM present topo climate shall be simulated with era iterim and ccsm deviations
	! SET THIS TO FALSE DURING TRANSIENT SIMULATIONS!!!!!!!!
    logical, parameter :: ccsm_lgmpt_climate = .FALSE.

    ! Flag, write climate to netcdf file
    logical, parameter :: save_climate_forcing = .false.

    ! Flag, write only initial climate to netcdf file after loaded
    logical, parameter :: save_initial_climate = .false.

    ! Flag, use Bern3D only
    logical, parameter :: b3d_climatology = .FALSE.

    ! Flag, initial ice sheet in greenland
    logical, parameter :: initial_gis = .true.

    ! Flag, initial ice sheet in greenland
    logical, parameter :: initial_lgm_ice = .false.

    ! Flag, to turn off Ice Dynamics, also shuts down the bedrock model
    logical, parameter :: ice_dynamics_on = .false.

    ! Debug Level. Higher Number = Less Information
    integer, parameter :: debug = 1

    ! Copies the netcdf input files to the output directory.
    logical, parameter :: store_input_netcdf = .FALSE.

    ! Copies the variables.f90 and all the other code files to the output directory.
    logical, parameter :: store_input_variables = .TRUE.

    ! Reads the initial state out of an NetCDF file, which was created with this model.
    logical, parameter :: read_initial_state = .FALSE.

    ! Ignore the SMB in the himalaya region.
    logical, parameter :: ignore_himalaya = .FALSE.

    ! Enable or disable bedrock relaxation
    logical, parameter :: active_bedrock = .true.

    ! Enable or disable ELRA bedrock relaxation
    logical, parameter :: active_elra = .FALSE.

    ! Adjust the sea level to the ice which is frozen on land. If it is not adjusted, the sea level offset will be 0!
    logical, parameter :: adjust_sea_level = .TRUE.

    ! Chose massbalance model, 1= Positive degree day, 2 = Energyflux mass balance, 3= Troll model
    integer, parameter :: smb_model = 2
    
    !Choose albedo module, 1=constant, 2=Bougamont, 3=Oerlemans , 4=Aoki 5=harmonic
    integer,parameter :: albedo_module= 1
    
    character(256), parameter :: albedo_file_variable_name = 'OERLEMANS'
    
    logical, parameter :: latent_heat_flux_on = .true.
    
    logical, parameter :: latent_heat_flux_analog_to_sensible_heat_flux = .true.

    ! Snow densification model
    logical, parameter :: densification_model = .True.

    ! Adjust temperature to different elevations
    logical, parameter :: active_temperature_lapsrate = .TRUE.

    ! Adjust the precipitation with elevation (Budd and Smith (1979))
    ! Above a specific threshold, the precipitation will get reduced by a factor of 0.5 every 1000m.
    logical, parameter :: active_elevation_desertification = .False. 

    ! Adjust shortwave radiation with respect to elevation for changing heights. This only includes damping, no amplification. 10%/1000m 
    ! there is no real physics behind it.
    logical, parameter :: short_wave_damping = .FALSE.


    ! A description of the experiment, which will be used for the foldername.
    character(256), parameter :: experiment_description = experiment_description_spam !&
                                        !'time3e5_LGM_beta7mm_iceflux1_acctemp2_bedrock3000'
                                        !'time5e5_LGM_lgmtopo_tempbias_adjustsealevel_beta7e-8_iceflux0.9_acctemp1'
                                        !'snowman_test_pd'
    !
    logical, parameter :: active_hysteresis = .FALSE.
    
    ! CSV output with the integrated mass balance
    logical, parameter :: store_integrated_mass_balance = .FALSE.

    ! Unit in the climate netcdf file, which is used to calculate the accumulation
    ! (not yet 100% implemented, check the output if you change this value)
    integer, parameter :: precipitation_unit = 2    ! 1 = m/yr, 2 = m/s, 3 = mm/yr

    ! create txt to check mass conservation
    logical, parameter :: mass_checking = .false.

    ! create txt to check energy conservation
    logical, parameter :: energy_checking = .false.

    ! Calculate new climate and sea level every # years
    integer, parameter :: calc_rate_climate = 20 ! 20 or 50 years

    ! Ice model data output
    integer, parameter :: netcdf_timesteps = 5! 500 years
    integer, parameter :: yearly_netcdf_file_freq = 500 ! max number of time entries in one file

! ! ! ! 
! ! ! !     ! Snow model output
! ! ! !     ! Toggle monthly snow data output, the rate in years is further down below
! ! ! ! !     logical, parameter :: monthly_data = .true.
! ! ! !     logical, parameter :: add_accum_output = .false. ! if additional output is desired about accumulation, rain, melt and refreezing of the snowcover
! ! ! !     
! ! ! !     integer, parameter :: monthly_data_freq = 5 ! every 'number' years	! 240
! ! ! !     integer, parameter :: monthly_data_num = 12 ! number of outputs during one year (12,24,48,(96))
! ! ! ! 
! ! ! !     logical, parameter :: daily_snowmask = .true.	! for every year a file is generated containing a mask with 
! ! ! ! 							! 1=Water, 2=Land, 3=Snow, 4=wet Snow, 5=pfa for every timestep
! ! ! ! 
! ! ! !     integer, parameter :: daily_snowmask_start = 1			!240	! year from which on diurnal snow extent output should be produced
! ! ! !     logical, parameter :: write_output_files_imhof = .false.
    
    !Snow model output frequency
    logical, parameter :: annual_data = .true.
    logical, parameter :: monthly_data = .false.
    logical, parameter :: daily_data = .false.
    
    integer, parameter :: annual_data_frequency = 50 !write annual data every X years, extremly high value if no intervalls are wanted
    integer, parameter :: monthly_data_frequency = 10 !write monthly data every X years
    integer, parameter :: daily_data_frequency = 100 !write daily data every X years
    
    integer, parameter :: daily_data_start = 400 !from which day on every data is written
    integer, parameter :: daily_data_end = 500
    
    logical, parameter :: include_values_over_ice = .false. !should melting ice be treated as run off on all cells

    ! Settings for cyclic ERA-I climate
    logical :: erai_climate_backwards = .False. ! do not change this
    integer :: erai_year=1979		! year to start with
    integer :: erai_year_begin=1979	! lower bound for cycled time periode
    integer :: erai_year_turn=2017	! upper bound for cycled time periode
    integer :: erai_iterations_max = 10
    integer :: erai_end_year = 2017
    integer :: erai_current_iteration=1	! do not change this
    logical :: erai_climate_final_iteration = .False. ! do not change this




    ! Initialize snow at the begining of the simulation where ice thickness > 200 m
    logical, parameter :: initialize_snow = .false.


    ! Snow Speed up routine
    logical :: adaptive_timestep = .false.
    integer, parameter :: calc_rate_snow = 20	!calc_rate_climate
    integer, parameter :: start_speed_up = 0
    integer, parameter :: end_speed_up   = 1000000



    !=========================
    ! Constant initialization
    !=========================
    ! maximum simulation duration in years: 2e4 =     20'000 years
    !                                       5e5 =    500'000 years
    !                                       1e7 = 10'000'000 years
    integer, parameter :: maxyears = maxyears_spam       !int(4e5)
    ! distance between two grid points [m] in X-Direction:  20000m = 20km
    !                                                       40000m = 40km
    real(kind=8), parameter :: dx = 10000. ! 40km
    ! distance between two grid points [m] in Y-Direction:  20000m = 20km
    !                                                       40000m = 40km
    real(kind=8), parameter :: dy = 10000. ! 40km
    ! Length of domain: Greenland:   1640km = 1.64e6 = 82 * 20 km
    !                   NH20:       12480km = 1.248e7 = 624 * 20 km
    !                   NH40:       12480km = 1.248e7 = 312 * 40 km
    integer, parameter :: L = int(1.64e6) ! 1.248e7 m
    ! Width of domain:  Greenland:   2800km = 2.8e6 = 140 * 20 km
    !                   NH20:       12480km = 1.248e7 = 624 * 20
    !                   NH40:       12480km = 1.248e7 = 312 * 40 km
    integer, parameter :: W = int(2.80e6) ! 1.248e7 m

    ! number of grid points in X Direction (Longitude). +1 cause 0 is also a grid box!
    integer, parameter :: nx = floor(real(W/dx)) + 1 ! declare and init variable, longitude
    ! number of grid points in Y Direction (Latitude).  +1 cause 0 is also a grid box!
    integer, parameter :: ny = floor(real(L/dy)) + 1 ! declare and init variable, latitude

    !=========================
    ! Tuning parameters
    !
    ! These parameters were used in the MasterThesis of Basil Neff (2014):
    ! IceBern2D: An efficient two-dimensional ice sheet model for paleoclimate studies
    !=========================
    ! Flux adjustments to tune the ice flux. This is a percentual adjustment of A_Eismint: 1 = 100%
    real(kind=8), parameter :: ice_flux_adjustments = 1
    ! http://www.igsoc.org/journal/59/218/j13J081.pdf
    ! Seguinot, Julien (2013): Spatial and seasonal effects of temperature variability in a positive degree-day glacier surface mass-balance model
    ! 3 mm * C^???1 * d^???1 = 3.47e-8 m * C^-1 * s^-1 for snow
    ! 8 mm * C^???1 * d^???1 = 9.26e-8 m * C^-1 * s^-1 for ice
    !
    ! 5 mm * C^???1 * d^???1 = 5.79e-8 m * C^-1 * s^-1
    ! 6 mm * C^???1 * d^???1 = 6.94e-8 m * C^-1 * s^-1 <- suggested by Neff
    ! 7 mm * C^???1 * d^???1 = 8.10e-8 m * C^-1 * s^-1
    ! 8 mm * C^???1 * d^???1 = 9.26e-8 m * C^-1 * s^-1
    ! 9 mm * C^???1 * d^???1 = 1.04e-7 m * C^-1 * s^-1
    !10 mm * C^???1 * d^???1 = 1.16e-7 m * C^-1 * s^-1
    real(kind=8), parameter :: beta = 2.*10.4e-8 ! 3 - 8 mm * C^???1 * d^???1 = 3 / (1000 * 24 * 60 * 60) = 3.5e-8 m * C^-1 * s^-1
    ! accumulation temperature threshhold in Celsius (higher values lead to a lower SMB):
    ! determines at which daily mean temperature the precipitation is treated as accumulation.
    real(kind=8), parameter :: accumulation_daily_temperature_threshold = 2
    ! relaxation time for bedrock sinking [s]:   3'000 yr = 9.4608e10 s
    !                                            6'000 yr = 1.89216e11 s
    real(kind=8), parameter :: tstar = 9.4608e10 ! 3.1536e11

    !=========================
    ! Input from files
    !=========================

    ! Landmask/Assignment mask
    !-------------------------
    character(256), parameter :: netcdf_input_bedrock = &
        ! 40km resolution
        !'/home/neff/MasterThesis/data/northhemisphere/topographie/etopo_nhem40_landmask.cdf'
        !'/local_scratch/imhof/model/input/topographie/etopo_nhem40_landmask.cdf'        
!	'/alphadata04/imhof/model/input/topographie/etopo_nhem40_landmask.cdf'
! 	'/home/tobias/PhD/Programming/NEW/input/etopo_nhem40_landmask.cdf'
	'/work/zolles/input_data/ERA_interim/grl10_interp/etopo_grl10_landmask.cdf'
	

    ! The NetCDF Variable name of the relaxed bedrock
    character(256), parameter :: netcdf_input_bedrock_variable = 'ICE_INTERP' ! 'LANDMASK'

    ! The NetCDF Variable name of the nonrelaxed pd bedrock
    character(256), parameter :: netcdf_input_pd_bedrock_variable = 'BED_INTERP' 

    ! The NetCDF Variable name of the surface with pd ice sheet
    character(256), parameter :: netcdf_input_icesurf_variable = 'SURFACE_INTERP' 

    ! Netcdf file with equilibrium ice for PD and the best parameters for PD
    character(256), parameter :: netcdf_input_bedrock_special = &
! 	'/home/tobias/PhD/Programming/NEW/input/ERAi_topog.interp.cdf'
	'/work/zolles/input_data/ERA_interim/grl10_interp/ERAi_topog.interp.cdf'


    ! The NetCDF path to the LGM ice
    character(256), parameter :: netcdf_input_lgm_ice = &
	'/work/zolles/input_data/ERA_interim/grl10_interp/ERAi_topog.interp.cdf'
	!'/home/ubelix/climate/mimhof/model/input/topographie/ice21kabp/lgm_I6_C.VM5a_sbx5_interp_21kabp.cdf'


	!'/alphadata04/imhof/model/input/topographie/ice26kabp/lgm_I6_C.VM5a_sbx5_interp.cdf'
	!'/home/ubelix/climate/mimhof/model/input/topographie/ice26kabp/lgm_I6_C.VM5a_sbx5_interp.cdf'

	!'/alphadata04/imhof/model/input/topographie/I6_C.VM5a_interp.cdf'
	!'/alphadata04/imhof/model/input/topographie/I6_C.VM5a_cleaned_ice_interp.cdf'
	!'/alphadata04/imhof/model/input/topographie/lgm_init_ice_1.cdf'
	!'/home/ubelix/climate/mimhof/model/input/topographie/lgm_init_ice_1.cdf'

    ! The NetCDF Variable name of the ice during the lgm
    character(256), parameter :: netcdf_input_lgm_ice_variable = 'ICE_LGM' 


    
    character(256), parameter :: albedo_file_path = &
    '/home/tobias/PhD/Programming/NEW/modeloutput/compare/single_albedo.cdf'
    ! Climate Bern3D
    !=================
    ! parameters for loading todays Bern3D climate data
    !file coordinates
    character(128), parameter :: netcdf_input_b3d_pd_climate = &
        !'/local_scratch/imhof/model/input/Bern3D_input/b3d_pd_climate.cdf'
        !'/local_scratch/imhof/model/input/Bern3D_input/climatology1/EEM2PIC260.0001950_full_inst.nc'
	!'/local_scratch/imhof/model/input/Bern3D_input/climatology2/EEM2PIC000.-128050_full_inst.nc'
!         '/home/tobias/PhD/Programming/NEW/input/ERAinterim_1979-2016.ave.interp.cdf'
        '/work/zolles/input_data/ERA_interim/grl10_interp/ERAinterim_2010.interp.cdf'
	!'/home/ubelix/climate/mimhof/model/input/Bern3D_input/climatology1/EEM2PIC260.0001950_full_inst.nc'
        

    ! parameters for loading 123ky BC Bern3D climate data
    !file coordinates
    character(128), parameter :: netcdf_input_b3d_eem_climate = &
        !'/local_scratch/imhof/model/input/Bern3D_input/b3d_pd_climate.cdf'
        !'/local_scratch/imhof/model/input/Bern3D_input/climatology1/EEM2PIC259.0001450_full_inst.nc'
!         '/home/tobias/PhD/Programming/NEW/input/ERAinterim_1979-2016.ave.interp.cdf'	
                '/work/zolles/input_data/ERA_interim/grl10_interp/ERAinterim_2010.interp.cdf'
	!'/home/ubelix/climate/mimhof/model/input/Bern3D_input/climatology1/EEM2PIC010.-123050_full_inst.nc'


    ! parameters for loading 19ky BC Bern3D climate data
    !file coordinates
    character(128), parameter :: netcdf_input_b3d_lgm_climate = &
!         '/home/tobias/PhD/Programming/NEW/input/ERAinterim_1979-2016.ave.interp.cdf'	
                '/work/zolles/input_data/ERA_interim/grl10_interp/ERAinterim_2010.interp.cdf'
	!'/home/ubelix/climate/mimhof/model/input/Bern3D_input/climatology1/EEM2PIC218.-019050_full_inst.nc'


    ! temperature
    character(256), parameter :: netcdf_input_b3d_temp_variable = 'ATMT_CLIM'

    ! precipitation
    character(256), parameter :: netcdf_input_b3d_precip_variable = 'PRECIP_CLIM'

    ! swradboa
    character(256), parameter :: netcdf_input_b3d_swradboa_variable = 'SWRADBOA_CLIM'

    ! albedo
    character(256), parameter :: netcdf_input_b3d_albedo_variable = 'ALBEDO_CLIM'

    ! debris weighting matrix
    character(128), parameter :: matrixpath = &
        !'/local_scratch/imhof/model/input/Bern3D_input/'
        '/alphadata04/imhof/model/input/Bern3D_input/'
	!'/home/ubelix/climate/mimhof/model/input/Bern3D_input/'


    ! Changing climate settings and paths
    !-----------------------------------   

    integer, parameter :: new_climate = 500	! rate for new bern3d climate data
    integer :: input_file_no=0 			! 128000 BC
    !integer :: input_file_no=174 		! 40000 BC
    integer, parameter :: year1 = -128050

    ! path to dynamic climatology direcotry
    character(128), parameter :: netcdf_input_b3d_timestep_climate_path = &
        !'/local_scratch/imhof/model/input/Bern3D_input/climatology1/'
        !'/local_scratch/imhof/model/input/Bern3D_input/climatology2/'
        !'/alphadata04/imhof/model/input/Bern3D_input/climatology2/'
        '/alphadata04/imhof/model/input/Bern3D_input/climatology1/'
	!'/home/ubelix/climate/mimhof/model/input/Bern3D_input/climatology1/'



    ! path to dynamic climatology file
    character(128) :: netcdf_input_b3d_timestep_climate_file = &
        'schneefroueli'
    ! dummy to build new input file name
    character(128) :: new_input_file = &
        'schneemandli'
    ! path to ERA-I timeseries files
    character(128) :: erai_yearly_timeseries_path = &
    '/scratch/zolles/input_data/ERA_interim/grl10_interp/'
!         '/home/tobias/PhD/Programming/NEW/input/ERA_interim/nhem40_interp/'
!not use for transient climate 





    ! ERA iterim Climate average (1979 - 2016)
    character(128), parameter :: netcdf_input_eraiterim_climate = &
!	'/alphadata04/imhof/model/input/erai_climate_input_topo/ERAinterim_1979-2016.ave.interp.cdf'
! 	'/home/tobias/PhD/Programming/NEW/input/ERAinterim_1979-2016.ave.interp.cdf'
	        '/work/zolles/input_data/ERA_interim/grl10_interp/ERAinterim_2010.interp.cdf'

    character(128), parameter :: netcdf_input_eraiterim_swradboa = &
	!'/alphadata04/imhof/model/input/erai_swradboa_daily_ave_interp.cdf'
	'/home/tobias/PhD/Programming/NEW/input/erai_swradboa_daily_ave_interp.cdf'


	! ERA interim climate directory with single years 1979 - 2016
    character(128), parameter :: netcdf_input_eraiterim_directory = &
!		'/alphadata04/imhof/model/input/erai_climate_input_topo/'
!	'/local_scratch/imhof/model/input/erai_climate_input_topo/'
! 	'/home/tobias/PhD/Programming/NEW/input/ERA_interim/nhem40_interp/'
		'/scratch/zolles/input_data/ERA_interim/grl10_interp/'

	! ERA interim input variable names
    character(256), parameter :: netcdf_input_eraiterim_temp_variable = 'T2M_INTERP'
    character(256), parameter :: netcdf_input_eraiterim_precip_variable = 'PRECIP_INTERP'
    character(256), parameter :: netcdf_input_eraiterim_swradboa_variable = 'SWRD_INTERP'
    character(256), parameter :: netcdf_input_eraiterim_dewpT_variable = 'D2M_INTERP'
    character(256), parameter :: netcdf_input_eraiterim_wind_variable = 'WIND_INTERP'
    character(256), parameter :: netcdf_input_eraiterim_lwrd_variable = 'LWRD_INTERP'


    ! ERA iterim climate reference topography
    character(128), parameter :: netcdf_input_eraiterim_initial_climate_elevation = &
!	'/alphadata04/imhof/model/input/erai_climate_input_topo/ERAinterim_topog.interp.cdf'
! 	'/home/tobias/PhD/Programming/NEW/input/ERAi_topog.interp.cdf'
	'/work/zolles/input_data/ERA_interim/grl10_interp/ERAi_topog.interp.cdf'

    character(256), parameter :: netcdf_input_eraiterim_initial_climate_elevation_variable = &
		'TOPOG_INTERP'



    ! CCSM PI Climate
    character(128), parameter :: netcdf_input_ccsm_pi_climate = &
        !'/alphadata04/imhof/model/input/climate/ccsm_pi_eem/CCSM4_PI.interp.cdf'
	'/home/ubelix/climate/mimhof/model/input/climate/ccsm_pi_eem/CCSM4_PI.interp.cdf'
    ! CCSM EEM Climate
    character(128), parameter :: netcdf_input_ccsm_eem_climate = &
        !'/alphadata04/imhof/model/input/climate/ccsm_pi_eem/CCSM4_125ka.interp.cdf'
	'/home/ubelix/climate/mimhof/model/input/climate/ccsm_pi_eem/CCSM4_125ka.interp.cdf'
    ! CCSM PI & EEM Climat variables names
    character(256), parameter :: netcdf_input_ccsm_temp_variable = 'T2M_INTERP'
    character(256), parameter :: netcdf_input_ccsm_precip_variable = 'PRECIP_INTERP'

    ! CCSM PI & EEM topography
    character(128), parameter :: netcdf_input_ccsm_initial_climate_elevation = &
	'/alphadata04/imhof/model/input/topographie/CCSM4_PI_topog.interp.cdf'
	!'/home/ubelix/climate/mimhof/model/input/topographie/CCSM4_PI_topog.interp.cdf'
    character(256), parameter :: netcdf_input_ccsm_initial_climate_elevation_variable = &
		'TOPOG_INTERP'

    ! CCSM EEM short wave radiation (monthly)
    character(128), parameter :: netcdf_input_ccsm_eem_swradboa = &
	!'/alphadata04/imhof/model/input/climate/ccsm_pi_eem/swradboa_eem_interp.cdf'
	'/home/ubelix/climate/mimhof/model/input/climate/ccsm_pi_eem/swradboa_eem_interp.cdf'



    ! CCSM LGM Climate on presentday topography
    character(128), parameter :: netcdf_input_ccsm_lgmpt_climate = &
        !'/alphadata04/imhof/model/input/climate/lgm_present_topo/climate_lgm_present_topo_nh40.nc'
	'/home/ubelix/climate/mimhof/model/input/climate/lgm_present_topo/climate_lgm_present_topo_nh40.nc'


    ! CCSM LGM Climate on lgm topography
    character(128), parameter :: netcdf_input_ccsm_lgmlt_climate = &
        !'/alphadata04/imhof/model/input/climate/lgm/climate_lgm_tempbias_nh40.nc'
	'/home/ubelix/climate/mimhof/model/input/climate/lgm/climate_lgm_tempbias_nh40.nc'

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!'/home/ubelix/climate/mimhof/model/input/climate/lgm/climate_lgm_nh40.nc'
        !'/alphadata04/imhof/model/input/climate/lgm/climate_lgm_nh40.nc'

    ! CCSM LGM short wave radiation (monthly)
    character(128), parameter :: netcdf_input_ccsm_lgm_swradboa = &
	!'/alphadata04/imhof/model/input/climate/lgm/swradboa_lgm_lgmtopo_interp.cdf'
	'/home/ubelix/climate/mimhof/model/input/climate/lgm/swradboa_lgm_lgmtopo_interp.cdf'

    ! CCSM LGM topography
    character(128), parameter :: netcdf_input_ccsm_lgm_initial_climate_elevation = &
	!'/alphadata04/imhof/model/input/topographie/ccsm_nhem40_LGM.nc'
	'/home/ubelix/climate/mimhof/model/input/topographie/ccsm_nhem40_LGM.nc'
    character(256), parameter :: netcdf_input_ccsm_lgm_initial_climate_elevation_variable = &
		'SURFACE'


    ! CCSM LGM Climate variables names
    character(256), parameter :: netcdf_input_ccsm_lgm_temp_variable = 'TEMPERATURE'
    character(256), parameter :: netcdf_input_ccsm_lgm_precip_variable = 'PRECIPITATION'
    character(256), parameter :: netcdf_input_ccsm_lgm_swradboa_variable = 'SWRADBOA_INTERP'



    ! inital climate elevation for downscaling to sealevel character dummies
    character(128) :: netcdf_input_initial_climate_elevation = &
        'leer'

    character(256) :: netcdf_input_initial_climate_elevation_variable = &
		'leer'



!    ! Climate CESM
!    !-------------------------
!    ! Climate data (temperature, precipitation and the elevation according to the initial temperature)
!    character(128), parameter :: netcdf_input_climate = &
!        ! LGM - bias subtracted
!        !'/alphadata04/imhof/input/climate/lgm/climate_lgm_tempbias_nh40.nc'
!        ! LGM
!        !'/alphadata04/imhof/input/climate/lgm/climate_lgm_nh40.nc'
!        ! PI
!        !'/alphadata04/imhof/input/climate/pi/climate_pi_tempbias_nh40.nc'
!        '/local_scratch/imhof/model/input/climate/pi/climate_pi_tempbias_nh40.nc'


!    ! The elevation which is used in the model of the input climate. This is used to calculate the temperature at the sea level.
!    ! For the ReAnalyse data, the elevation of the ice surface is used. Please also adjust the 'netcdf_input_climate_elevation_variable' variable.
!    character(256), parameter :: netcdf_input_climate_elevation = &
!        ! 40km resolution
!        ! LGM
!        !'/alphadata04/imhof/input/topographie/ccsm_nhem40_LGM.nc'
!        ! PI
!        !'/alphadata04/imhof/input/topographie/ccsm_nhem40_present.nc'
!        '/local_scratch/imhof/model/input/topographie/ccsm_nhem40_present.nc'
!
!    ! Variable name in the netcdf_input_climate_elevation file which is used for the elevation.
!    character(256), parameter :: netcdf_input_climate_elevation_variable = 'SURFACE' !'ICE_NHEM20' ! ICE_NHEM20 = Northhemisphere for PD, ICE_grl10 = Greenland


    ! Initial state
    ! Not yet 100% tested
    !=========================
    ! Read initial state (this can be controlled with the flag read_initial_state)
    character(128), parameter :: initial_netcdf_file = '/local_scratch/imhof/model/modelinput/input.nc'  ! The initial netcdf file with the the height and bedrock, no time axes!
    character(128), parameter :: initial_bedrock_variable_name = 'BEDROCK'              ! The name of the bedrock variable
    character(128), parameter :: initial_height_variable_name = 'HEIGHT'                ! The name of the height variable

    !=========================
    ! Output to files
    !=========================

    ! netCDF
    !-------
    ! Root output directory, a subdirectory with a timestamp and experiment_description gets created
    character(512) :: output_directory = & 
!	'/local_scratch/imhof/model/modeloutput/'
!	'/local_scratch/imhof/model/hysterese/'
!	'/alphadata04/imhof/no_backup/'
	'/work/zolles/modeloutput/ERAi_ensemble/'

    ! Filename of the netCDF file with all variables of the ice and masks.
    character(128) :: netcdf_output_filename != 'IceBern2D.nc'

    integer :: nc_counter = 1 !next entry to write in open nc file

    !===================================================================================
    ! Do NOT CHANGE ANYTHING BELOW THIS LINE! ONLY IF YOU REALY KNOW WHAT YOU ARE DOING.
    !===================================================================================

    ! Constants
    !==========
    ! seconds in a year (Eismint Table 1): 31556926
    real(kind=8), parameter :: seconds_per_year = 3600.*24.*365. ! 31556926
    ! time steps [s]: 1 year = 3600*24*365*1
    real(kind=8), parameter :: dt = 3600.*24.*365. ! timesteps of ice model
    ! Diffusivity parameter: Oerlemans from the paper: 3.169 * 10^-8 m^(-3/2) * s^-1 = 1 m^(-3/2) * a^-1 with n = 2.5
    real(kind=8), parameter :: A_Oerlemans = 3.16887646e-8
    ! Diffusivity parameter: Eismint: 10^-16 Pa^-3 a^-1 = 3.1709792e-24 Pa^-3 s^-1 with n = 3
    real(kind=8) :: A_Eismint  =  A_Eismint_spam ! , parameter 
    ! Exponent in Glens's flow law used in Eismint: 3
    real(kind=8), parameter :: n_Eismint = 3
    ! Acceleration of gravity: 9.81 m s^-2
    real(kind=8), parameter :: g = 9.81
    ! Ice density: 910 kg m^-3
    real(kind=8), parameter :: rho_ice = 910 ! Eismint, Table 1
    ! Saturated adiabatic lapse rate
    real(kind=8), parameter :: temperature_lapse_rate = 0.0065 ! degree celsius per m
    ! extinction of short wave radiatin in atmosphere
    real(kind=8), parameter :: k_extinct = 0.0001 	! this value is just a guess, it corresponds to 10% per km
    ! Precipitation lapse rate (in SICOPOLIS gamma_p) above a specific treshhold
    real(kind=8), parameter :: precipitation_lapse_rate =  0.7/1000 ! in m^-1 and positive ! log(real(2))
    ! Elevation threshold, above which elevation should the elevation desertification be active
    integer, parameter :: precipitation_threshold = 2000. ! m
    ! Area of the ocean, needed to calculate the sea level rise
    real(kind=8), parameter :: ocean_area = 3.625e14 ! 3.6e8km^2 = 3.6e14 m^2
    ! At which position would the sea level be without ice [m]
    real(kind=8)  :: sea_level_offset = 7.36
    ! Hysteresis Variables
    !---------------------
    ! Hysteresis Return Point in Time in years
    real(kind=8), parameter :: hysteresis_return_point_in_time = real(5e6)/2 !real(maxyears)/2 ! yr
    ! Hysteresis Temperature Delta
    real(kind=8), parameter :: hysteresis_temperature_delta = 10.
    ! Hysteresis factor, a positive value. First it gets colder and
    real(kind=8), parameter :: hysteresis_temperature_factor = hysteresis_temperature_delta/real(hysteresis_return_point_in_time)  ! = K yr^1
    ! Hysteresis initial offset in Kelvin
    real(kind=8), parameter :: hysteresis_inital_temperature_offset = 5. ! in Kelvin
    !-- Hysteresis Stop
    ! Stop hysteresis after specified years?
    logical, parameter :: hysteresis_stop = .False.
    ! Stop Hysteresis at specific temperature
    integer, parameter :: hysteresis_stop_year = (maxyears - 100000) !5e5
    !-- Hysteresis temperature jump
    ! Flag: temperature jump
    logical, parameter :: hysteresis_temperature_jump = .TRUE.
    ! Temperature jump
    real(kind=8), parameter :: hysteresis_temperature_jump_temperature = -2. ! in Kelvin and delta to the temperature at the specific time
    ! Temperature jump year
    integer, parameter :: hysteresis_temperature_jump_year = hysteresis_stop_year !(maxyears - 100000)
    
    ! ELRA constants
    !---------------
    ! density of asthenosphere (kg m^-3)
    real(kind=8), parameter :: rho_asthenosphere = 3300
    ! flexural rigidity (N m)
    real(kind=8), parameter :: elra_K1 = 1.e25
    !
    ! length scale of lithosphere deformation
    real(kind=8), parameter :: elra_L_r = (elra_K1/rho_asthenosphere/g)**0.25
    ! relaxation time for bedrock sinking [s]:   3'000 yr = 9.4608e10 s
    !                                            6'000 yr = 1.89216e11 s
    real(kind=8), parameter :: elra_tstar = tstar !9.4608e10 ! 3.1536e11
    
    real(kind=8), parameter :: elra_factor = elra_L_r**2/(2d0*3.14159*elra_K1)
    
    real(kind=8), parameter :: rho_ice_g_dx_dy = rho_ice * g * dx * dy
    
    ! KEI values
    character(128), parameter :: elra_kei_file = &
                                    '/alphadata03/neff/input/common/kei.csv'
    real(kind=8), dimension(:) :: kei(1059,2)
    real(kind=8), dimension(nx,ny) :: elra_wss = 0d0
    real(kind=8) :: elra_kei_value, elra_f_0
    integer :: elra_ii
    
    
    ! ELRA variables
    !---------------
    ! Distance from the center of the bedrock calculations
    real(kind=8) :: elra_delta_distance = 0
    ! mass of local ice column
    real(kind=8) :: elra_ice_mass = 0
    ! 
    real(kind=8), dimension(nx,ny) :: elra_w = 0
 

    ! amount of timesteps per year for accumulation
    integer, parameter :: ndays = 365		! 365 only working with ERAi input, use 96 for Bern3D-like time steps
      
    ! arrays for climate data
    real(kind=8), dimension(nx,ny,365) :: inp_temp = 0. ! 
    real(kind=8), dimension(nx,ny,365) :: inp_precip = 0. ! 
    real(kind=8), dimension(nx,ny,ndays) :: P_sun = 0. !1360/4*0.75 ! Shortwaveradiation boa
    real(kind=8), dimension(nx,ny,ndays) :: P_sun0 = 0. ! Shortwaveradiation boa raw
    
    real(kind=8), dimension(nx, ny, ndays) :: inp_dewpT = 0.
    real(kind=8), dimension(nx, ny, ndays) :: inp_wind = 0.
    real(kind=8), dimension(nx, ny, ndays) :: inp_lwrd = 0.
    
    real(kind=8), dimension(nx, ny, ndays) :: DewpT = 0.
    real(kind=8), dimension(nx, ny, ndays) :: wind = 0.
    real(kind=8), dimension(nx, ny, ndays) :: lwrd = 0.

    real(kind=8), dimension(nx,ny,ndays) :: b3d_pd_precip = 0.
    real(kind=8), dimension(nx,ny,ndays) :: b3d_pd_temp = 0.
    real(kind=8), dimension(nx,ny,ndays) :: b3d_pd_P_sun = 0.

    real(kind=8), dimension(nx,ny,ndays) :: deviation_precip = 0.
    real(kind=8), dimension(nx,ny,ndays) :: deviation_temp = 0.
    real(kind=8), dimension(nx,ny,ndays) :: deviation_P_sun = 0.


    ! some places
    integer, parameter :: placex = 133 !80
    integer, parameter :: placey = 166 !155

    integer, parameter :: xmin = 75
    integer, parameter :: xmax = 150
    integer, parameter :: ymin = 100
    integer, parameter :: ymax = 190

    ! Bern3D input data parameters

    integer, parameter :: bnx=41
    integer, parameter :: bny=40
    integer, parameter :: interpol_deg = 5

    ! grid interpolation
    real(kind=8), dimension(nx,ny,interpol_deg*3) :: thematrix = 0. !

    ! Simple shelv depth
    real(kind=8) :: shelv_depth = 0. ! -300.




    ! Runtime Variables
    !==================
    ! Write out status information
    character(128) :: heartbeat = ''
    ! To get the execution time: http://gcc.gnu.org/onlinedocs/gcc-4.0.4/gfortran/ETIME.html
    real, dimension(2) :: execution_time ! user, system
    real :: runtime = 0
    real :: clock_start = 0
    real :: clock_end = 0
    ! Loop counter (primary loop, time)
    integer :: it = 1
    integer, parameter :: timesteps = int(real(maxyears)/real(dt)*(3600.*24.*365.))
    ! Timesteps with the years belonging to each step
    integer, dimension(timesteps) :: myyear = 1
    ! Elevation of the bedrock w/o ice load. This is the relaxed bedrock, the postition where it would be without ice.
    real(kind=8), dimension(nx,ny) :: Bedrock_Initial = 0
    ! Elevation of the bedrock, on water the bedrock is equal to the sea level. This is done for internal calculations.
    real(kind=8), dimension(nx,ny) :: Bedrock = 0 ! B = Bedrock_Initial
    ! Ice thickness
    real(kind=8), dimension(nx,ny) :: ice_thickness = 0
    ! Elevation of ice surface above sea level
    real(kind=8), dimension(nx,ny) :: elevation = 0
    ! assignment_mask: 0 = ice, 1 = water, 2 = no ice, 3 = unstable integration
    integer, dimension(nx, ny) :: assignment_mask = 0   ! "Normal" case (ice) = 0, Water = 1, no ice in this time step = 2
    ! Distance of each grid box from domain boundary at the bottom, will be multiplied with dy at a later point.
    real(kind=8), dimension(ny) :: y = (/ (I, I = 0, ny - 1) /) ! Declare variable
    ! Distance of each grid box from left domain boundary, will be multiplied with dx at a later point.
    real(kind=8), dimension(nx) :: x = (/ (I, I = 0, nx - 1) /) ! Declare variable
    ! Sea leval
    real(kind=8) :: sea_level = 0
    ! Output variables, with "realistic" values, values of the bedrock is not 0 at the sea
    !-------------------------------------------------------------------------------------
    ! The bedrock, with values below 0 on the water
    real(kind=8), dimension(nx, ny) :: bedrock_netcdf
    ! Elevation with ice: variable 'elvevation' on land, bedrock_netcdf in the water
    real(kind=8), dimension(nx, ny) :: elevation_netcdf
    ! (Potential) Temperature and Precipitation during runtime
    !---------------------------------------------------------
    ! Temperature in Kelvin at the ice elevation for 365 days to calculate the accumulation and ablation out of it.
    real(kind=8), dimension(nx, ny, ndays) :: temperature = 0.
    

    ! Temperature in Kelvin at the sea level for 365 days to calculate the accumulation and ablation out of it.
    real(kind=8), dimension(nx, ny, ndays) :: potential_temperature = 0.
    ! Initial elevation of the surface from the climate input, where the temperature is measured and precipiation falls. (h_GCM)
    real(kind=8), dimension(nx,ny) :: initial_climate_elevation = 0.

    ! Initial elevation of the surface from the ccsm climate input, where the temperature is measured and precipiation falls. (h_GCM)
    real(kind=8), dimension(nx,ny) :: initial_ccsm_elevation = 0.

    ! Precipitation (m/yr) for 365 days from the climate model (precipitation at the surface of the initial_climate_elevation)
    real(kind=8), dimension(nx, ny, ndays) :: initial_climate_precipitation = 0
    ! Precipitation (m/yr) for 365 days to calculate the accumulation and ablation out of it.
    ! If active_elevation_desertification is true, the precipitation is adjusted to the elevation above a specific threshold (precipitation_threshold)
    real(kind=8), dimension(nx, ny, ndays) :: precipitation = 0
    ! For diffusivity calculations
    real(kind=8) :: hgradient = 0
    ! For diffusivity calculations
    real(kind=8), dimension(nx,ny) :: D = 0
    ! For ice flux calculations, on the north side of the point.
    real(kind=8), dimension(nx,ny) :: FN = 0
    ! For ice flux calculations, on the east side of the point.
    real(kind=8), dimension(nx,ny) :: FE = 0
    ! For ice flux calculations, on the south side of the point.
    real(kind=8), dimension(nx,ny) :: FS = 0
    ! For ice flux calculations, on the west side of the point.
    real(kind=8), dimension(nx,ny) :: FW = 0
    ! Surface Mass Balance
    real(kind=8), dimension(nx,ny) :: surface_mass_balance = 0
    real(kind=8), dimension(nx,ny) :: accumulation = 0
    real(kind=8), dimension(nx,ny) :: ablation = 0
    ! Discharge
    real(kind=8), dimension(nx, ny) :: discharge_x =0 ! Discharge (mass flow) in x direction at each point in m yr^-1
    real(kind=8), dimension(nx, ny) :: discharge_y = 0      ! Discharge (mass flow) in y direction at each point in m yr^-1
    ! Ice frozen on land (needed for sea level)
    real(kind=8) :: ice_volume = 0
    ! For later reanalysis: Medium ice elevation at 1m^2
    real(kind=8), dimension(timesteps) :: H_ts = 0          ! medium ice height at 1m^2 over time




    ! netCDF IceBern2D
    !=================
    ! cant be a parameter, cause the netcdf functions modifies it.
    integer :: filehandle_netcdf
    ! Last year the netCDF was written. This is used, if the timestep is below one year, that the netCDF is not written several times in this year.
    integer :: last_netcdf_year = 0
    ! Variable IDs: height, bed
    !--------------------------
    ! Variable for the netCDF Height, parameter
    integer :: height_varid
    ! Variable for the netCDF Bedrock, parameter
    integer :: bedrock_varid
    ! Variable for the netCDF accumulation, parameter
    integer :: acc_varid
    ! Variable for the netCDF ablation, parameter
    integer :: abl_varid
    ! Variable for the netCDF diffusivity per year, parameter
    integer :: diffusivity_varid
    ! Variable for the netCDF discharge in x direction per year, parameter
    integer :: discharge_x_varid
    ! Variable for the netCDF discharge in y direction per year, parameter
    integer :: discharge_y_varid
    ! Variable for the netCDF assignent_mask, parameter
    integer :: assignent_mask_varid
    
    ! integrated _mass_ balance (not as ice volume!): imb
    !=======
    ! CSV File with integrated mass balance information: accumulation, ablation, calving, icevolume
    character(128), parameter :: imb_output_filename = 'IntegratedMassBalance.csv'
    ! File handle
    integer, parameter :: imb_filehandle=20
    ! accumulation (positive part of the SMB)
    real(kind=8) :: imb_accumulation = 0
    ! ablation (negative part of SMB)
    real(kind=8) :: imb_ablation = 0
    ! calving
    real(kind=8) :: imb_calving = 0
    ! Increasing sealevel lead to ice melt
    real(kind=8) :: imb_isostaticmelt = 0
    ! total ice mass
    real(kind=8) :: imb_totalmass = 0
    ! ice mass change compared to the year before
    real(kind=8) :: imb_totalmass_change = 0
    ! Temp variable for different calculations
    real(kind=8) :: imb_tmp = 0
    
    real(kind=8) :: imb_cellcount = 0



    ! netCDF Climate
    !===================
    ! cant be a parameter, cause the netcdf functions modifies it.
    integer :: filehandle_netcdf3
    ! Variable IDs: height, bed
    !--------------------------
    ! Variable for the netCDF snow mass, parameter
    integer :: precip_varid
    ! Variable for the netCDF lw mass, parameter
    integer :: temp_varid
    ! Variable for the netCDF density of snow, parameter
    integer :: dev_precip_varid
    ! Variable for the netCDF temperature of snow, parameter
    integer :: dev_temp_varid
    ! Variable for the netCDF temperature of snow, parameter
    integer :: swrad_varid
    ! Variable for the netCDF temperature of snow, parameter
    integer :: dev_swrad_varid

    character(128) :: netcdf_output_climate_filename

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

	! Parallel Computing

	!character(256) :: set_core_command = 'export OMP_NUM_THREADS=4'






	! Filename of the netCDF file with all snow variables
	character(128) :: netcdf_output_filename2 != 'Snowcover3D.nc'


	! logicals for snowmodel setting
	!-----------------------------------
	! model for snow densification for densities > 550 kg/m^3
	integer, parameter :: hl=0		! 0 Barnola Pimienta, 1 Herron Langway

	! model for snow temperature diffusivity
	integer, parameter :: diff_model = 1 ! 1 Yen, 2 Sturm, 3 Dusen



	! new constants for emb models
	! ------------------------

	! determines at which temperature in C snow falls with temperature 0C
	real(kind=8), parameter :: snow_fall_temperature=3.

	! physical conastants
	real(kind=8), parameter :: rho_w = 1000.   ! density of water [kg/m3]
	real(kind=8), parameter :: rho_s = 350.! density of falling snow [kg/m3]
	real(kind=8), parameter :: rho_i = 917.! density of ice [kg/m3] defined twice!!
	real(kind=8), parameter :: rho_e = 830.! density of ice [kg/m3]
	real(kind=8), parameter :: c_w = 4181. ! heatcapacity of water at 25C [J/kg/K]
	real(kind=8), parameter :: c_i = 2110. ! heatcapacity of ice at -10C [J/kg/K]
	real(kind=8), parameter :: L_lh = 3.337e5  ! latent heat of ice [J/kg]
	real(kind=8), parameter :: L_v = 2.501e6 !Latent heat of vaporization
	!real(kind=8), parameter :: L_f = 337000. !Latent heat of fusion 
	real(kind=8), parameter :: K_ice=2.33 ! W/m/K
	real(kind=8), parameter :: P_atm=1e5
    real(kind=8), parameter :: cp_air=1003 
    
	real(kind=8), parameter :: albedo_snow_new = albedo_snow_new_spam  ! albedo of snow in arctic 0.8 0.75
	real(kind=8), parameter :: albedo_snow_wet = albedo_snow_wet_spam  ! albedo of wet snow i.e. snow with q=0 0.6 0.5
	real(kind=8), parameter :: albedo_ice = albedo_ice_spam! albedo of ice 0.4 0.2
	real(kind=8), parameter :: kelvin = 273d0
	real(kind=8), parameter :: D_sf = D_sf_spam   ! sensible heat flux snow-air [W/m2/K]

	real(kind=8), parameter :: sigma = 5.670373e-8 ! stefan-bolzman constant [W/m2/K4]
	real(kind=8), parameter :: eps_air =  eps_air_spam   ! highly dependent on cloubdcover and watervapor preasure [0.75 0.95]
	real(kind=8), parameter :: eps_snow = 0.98   ! 0.98 would be more realistic
	real(kind=8), parameter :: max_lwc = max_lwc_spam  ! percentage of empty volume that can be filed with water

!     real(kind=8), parameter :: D_lf = D_lf_spam
!     real(kind=8), parameter :: D_lf = D_sf/1003*0.622*(L_v+L_lh)/10000
	
	
	real(kind=8), parameter :: dt_firn = 365. * 3600. * 24./real(ndays)   ! time step lenght [s]
	real(kind=8), parameter :: rho_pass_on = rho_e
	integer :: simulstep = 1





	! characterisation of firn
	! -------------------------------------
	integer,  parameter :: n_snowlayer = 15 ! number of snowlayers, min 1
	real(kind=8), parameter ::       soll_mass = 300. ! aimed mass/area per layer
	real(kind=8), parameter :: lower_massbound = 100.
	real(kind=8), parameter :: upper_massbound = 500.


	! new runtime variables for emb models
	! -------------------------------------

	real(kind=8), dimension(nx,ny,n_snowlayer) :: snowman = 0.  ! snow mass kg/m2
	real(kind=8), dimension(nx,ny,n_snowlayer) :: lwmass = 0. ! liquid water in a gridcell in kg/m2
	real(kind=8), dimension(nx,ny,n_snowlayer) :: rho_snow = rho_s  ! density of snow in  kg/m3
	real(kind=8), dimension(nx,ny,n_snowlayer) :: snow_temp = 0. ! temperature of snow in C
    real(kind=8), dimension(nx,ny) :: albedo_dynamic = albedo_snow_new_spam
    
    !real(kind=8), dimension(nx,ny,n_snowlayer) :: albedo_runtime = 0.

	logical :: china_syndrome = .false. ! used for the melting subroutine
	logical, dimension(nx,ny) :: fast_calculation = .false. ! temperature of snow in C


	! mass conservation
	real(kind=8)  :: sunken_snow = 0.
	logical  :: spinup = .false.


	! netCDF Snowcover3D
	!===================
	! cant be a parameter, cause the netcdf functions modifies it.
	integer :: filehandle_netcdf2
	! Variable IDs: height, bed
	!--------------------------
	! Variable for the netCDF snow mass, parameter
	integer :: snowman_varid
	! Variable for the netCDF lw mass, parameter
	integer :: lwmass_varid
	! Variable for the netCDF density of snow, parameter
	integer :: rho_snow_varid
	! Variable for the netCDF temperature of snow, parameter
	integer :: snow_temp_varid
	! Variable for the netCDF temperature of snow, parameter
	integer :: albedo_dynamic_varid







end module variables
