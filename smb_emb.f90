! This module calculates the mass surface balance with the engergyfluxes in a multylayer firn
!
! Created on: May 12, 2015
!   Author: Imhof Michael
!   Developer: Imhof Michael
!   Mail: imhof@vaw.baug.ethz.ch or imhof@climate.unibe.ch
!
! turn on with smb_model = 2

! changed and adapted by Tobias Zolles April, 2018
! changes:
! - latent heat flux included over ice and snow, with the bulk method analog to the sensible heat flux
!       - turn on with latent_heat_flux_on and check latent_heat_flux_analog_to_sensible_heat_flux for same exchange coefficient as the sensible heat flux
! - temporal variation of the surface albedo possible
!       - change with albedo_module 1=constant, 2=harmonic fct, 3=Oerlemans , 4=Aoki 5=Bougamont
! - Output changed to single file for the smb output which contains: 
!   3D variables: 
!   snowman, lwmass, rho_snow, snow_temp
!   2D variables:
!   albedo_dynamic,latent_heat, accum, rain, melt, refreezing, snow, runoff, snow_mask, melt_of_ice, averaged_surface_temperature, amount_of_regrids, mass_balance
! - Output can be written in following time scopes:
!       - annualy
!       - daily (slow currently)


! last modifications August 2016:
! - continous mass conserving routine that passes on snow to the ice model



module smb_emb

	use variables
	!use variables_snow
	use io  ! Own module to read the values
	! use OMP_LIB
contains

subroutine get_accumulation_snowman(nx, ny, ndays, n_snowlayer, air_temp_ice, precip_ice, P_sun, landmask, seafloor, sealevel, &
		snowman, lwmass, rho_snow, snow_temp, smb_ice, nc_entry, year, sunken_snow, spinup, fast_calculation, albedo_dynamic, &
		DewpT, lwrd, wind, elevation)
	! input
	real(kind=4), dimension(nx,ny,n_snowlayer,ndays) :: m_test_array
	integer, intent(in) :: NX
	integer, intent(in) :: NY
	integer, intent(in) :: ndays
	integer, intent(in) :: n_snowlayer
	integer, intent(in) :: year
	integer, intent(in) :: nc_entry
	logical, intent(in) :: spinup

	real(kind=8), intent(in) :: air_temp_ice(nx, ny, ndays) ! Temperature at the level of the ice elevation
	real(kind=8), intent(in) :: precip_ice(nx, ny, ndays)   ! Precipitation meta data
	real(kind=8), intent(in) :: P_sun(nx, ny, ndays)   ! sw radiation
	integer, intent(in) :: landmask(nx, ny) ! ice = 0, water  = 1, land with no ice = 3
	real(kind=8), intent(in) :: seafloor(nx, ny)   ! sw radiation
	real(kind=8), intent(in) :: elevation(nx, ny)
	real(kind=8), intent(in) :: sunken_snow
	real(kind=8), intent(in) :: sealevel
	!real(kind=8), intent(in) :: RHi(nx, ny, ndays) 
    real(kind=8), intent(in) :: DewpT(nx, ny, ndays)
    real(kind=8), intent(in) :: lwrd(nx, ny, ndays)
    real(kind=8), intent(in) :: wind(nx, ny, ndays)
    
   
    
	! output
	real(kind=8), intent(inout) :: snowman(nx, ny, n_snowlayer)
	real(kind=8), intent(inout) :: lwmass(nx, ny, n_snowlayer)
	real(kind=8), intent(inout) :: rho_snow(nx, ny, n_snowlayer)
	real(kind=8), intent(inout) :: snow_temp(nx, ny, n_snowlayer)
	real(kind=8), intent(inout) :: smb_ice(nx, ny) ! height in ice equivalents
	logical, intent(inout) :: fast_calculation(nx, ny)
	real(kind=8), intent(inout) :: albedo_dynamic(nx, ny)
	real(kind=8), dimension(nx,ny) :: albedo_runtime
    real(kind=8), dimension(nx,ny) :: p_air(nx,ny)


! ! ! ! 	! arrays for monthly data output !imhof version
! ! ! ! 	real(kind=4), dimension(:,:,:,:), allocatable :: m_snowman
! ! ! ! 	real(kind=4), dimension(:,:,:,:), allocatable :: m_lwmass
! ! ! ! 	real(kind=4), dimension(:,:,:,:), allocatable :: m_rho_snow
! ! ! ! 	real(kind=4), dimension(:,:,:,:), allocatable :: m_snow_temp
! ! ! !     
! ! ! !     real(kind=4), dimension(:,:,:), allocatable :: m_albedo_dynamic
! ! ! !     
! ! ! !     real(kind=4), dimension(:,:), allocatable :: m_vaporflux
! ! ! ! 	real(kind=4), dimension(:,:), allocatable :: m_accum
! ! ! ! 	real(kind=4), dimension(:,:), allocatable :: m_rain
! ! ! ! 	real(kind=4), dimension(:,:), allocatable :: m_refreezing
! ! ! ! 	real(kind=4), dimension(:,:), allocatable :: m_melt
! ! ! ! 	real(kind=4), dimension(:,:), allocatable :: m_snow
! ! ! ! 	real(kind=4), dimension(:,:), allocatable :: m_runoff
! ! ! ! 
! ! ! ! 	integer,      dimension(:,:,:), allocatable :: snow_mask

!     ! arrays for annual data output
!     !3D snow layers snapshot values
    real(kind=4), dimension(:,:,:), allocatable, save :: m_snowman_annual
	real(kind=4), dimension(:,:,:), allocatable, save :: m_lwmass_annual
	real(kind=4), dimension(:,:,:), allocatable, save :: m_rho_snow_annual
	real(kind=4), dimension(:,:,:), allocatable, save :: m_snow_temp_annual
	
    real(kind=4), dimension(:,:), allocatable, save :: m_albedo_dynamic_annual
    
    real(kind=4), dimension(:,:), allocatable, save :: m_snow_temp_ave_surface_annual
    
    real(kind=4), dimension(:,:), allocatable, save :: m_vaporflux_annual
    real(kind=4), dimension(:,:), allocatable, save :: m_accum_annual
	real(kind=4), dimension(:,:), allocatable, save :: m_rain_annual
	real(kind=4), dimension(:,:), allocatable, save :: m_refreezing_annual
	real(kind=4), dimension(:,:), allocatable, save :: m_melt_annual
	real(kind=4), dimension(:,:), allocatable, save :: m_snow_annual
	real(kind=4), dimension(:,:), allocatable, save :: m_runoff_annual
	real(kind=4), dimension(:,:), allocatable, save :: m_melt_ice_annual

	integer,      dimension(:,:), allocatable, save :: snow_mask_annual
	
	integer,      dimension(:,:), allocatable, save :: m_regridding_annual
	
	real(kind=4), dimension(:,:), allocatable, save :: m_real_mass_balance_annual
	
	!     ! arrays for monthly data output
!     !3D snow layers snapshot values
    real(kind=4), dimension(:,:,:), allocatable, save :: m_snowman_monthly
	real(kind=4), dimension(:,:,:), allocatable, save :: m_lwmass_monthly
	real(kind=4), dimension(:,:,:), allocatable, save :: m_rho_snow_monthly
	real(kind=4), dimension(:,:,:), allocatable, save :: m_snow_temp_monthly
	
    real(kind=4), dimension(:,:), allocatable, save :: m_albedo_dynamic_monthly
    
    real(kind=4), dimension(:,:), allocatable, save :: m_snow_temp_ave_surface_monthly
    
    real(kind=4), dimension(:,:), allocatable, save :: m_vaporflux_monthly
    real(kind=4), dimension(:,:), allocatable, save :: m_accum_monthly
	real(kind=4), dimension(:,:), allocatable, save :: m_rain_monthly
	real(kind=4), dimension(:,:), allocatable, save :: m_refreezing_monthly
	real(kind=4), dimension(:,:), allocatable, save :: m_melt_monthly
	real(kind=4), dimension(:,:), allocatable, save :: m_snow_monthly
	real(kind=4), dimension(:,:), allocatable, save :: m_runoff_monthly
	real(kind=4), dimension(:,:), allocatable, save :: m_melt_ice_monthly

	integer,      dimension(:,:), allocatable, save :: m_snow_mask_monthly
	
	integer,      dimension(:,:), allocatable, save :: m_regridding_monthly
	
	real(kind=4), dimension(:,:), allocatable, save :: m_real_mass_balance_monthly
	!     ! arrays for annual data output
!     !3D snow layers snapshot values
    real(kind=4), dimension(:,:), allocatable, save :: m_snowman_daily
	real(kind=4), dimension(:,:), allocatable, save :: m_lwmass_daily
	real(kind=4), dimension(:,:), allocatable, save :: m_rho_snow_daily
	real(kind=4), dimension(:,:), allocatable, save :: m_snow_temp_daily
	
    real(kind=4), dimension(:), allocatable, save :: m_albedo_dynamic_daily
    
    real(kind=4), dimension(:), allocatable, save :: m_snow_temp_ave_surface_daily
    
    real(kind=4), dimension(:), allocatable, save :: m_vaporflux_daily
    real(kind=4), dimension(:), allocatable, save :: m_accum_daily
	real(kind=4), dimension(:), allocatable, save :: m_rain_daily
	real(kind=4), dimension(:), allocatable, save :: m_refreezing_daily
	real(kind=4), dimension(:), allocatable, save :: m_melt_daily
	real(kind=4), dimension(:), allocatable, save :: m_snow_daily
	real(kind=4), dimension(:), allocatable, save :: m_runoff_daily
	real(kind=4), dimension(:), allocatable, save :: m_melt_ice_daily

	integer,      dimension(:), allocatable, save :: m_snow_mask_daily
	
	integer,      dimension(:), allocatable, save :: m_regridding_daily
	
	real(kind=4), dimension(:), allocatable, save :: m_real_mass_balance_daily

	!local variables
	real(kind=8) :: D_lf
	real(kind=8) :: Q_heat
	real(kind=8) :: dummy
	real(kind=8) :: H_lh
	real(kind=8) :: K_lh
	real(kind=8) :: K_sw
	real(kind=8) :: accum
	real(kind=8) :: rainman
	real(kind=8) :: masssum
	real(kind=8) :: lwc
	real(kind=8) :: vaporflux
	real(kind=8) :: dummy_melt_ice
	real(kind=8) :: dummy_rain_ice
	real(kind=8), dimension(n_snowlayer) :: dz
	integer :: time
	integer:: nday_snowfall
	integer :: dummy_regrid

	integer :: simulstep=1
	integer :: month=1
	integer :: day
	integer :: shift

! ! ! ! 	! file name for monthyly data
! ! ! ! 	integer :: filehandle_netcdf_monthly
! ! ! ! 	character(128) :: netcdf_output_filename_monthly
! ! ! ! 	
! ! ! 
! ! ! 	
! ! ! 
! ! ! ! 	! file name for monthyly data
! ! ! ! 	integer :: filehandle_netcdf_snow_mask
! ! ! ! 	character(128) :: netcdf_output_filename_snow_mask
! ! ! ! 	integer :: snow_mask_varid

	! variables for mass conservation

	real(kind=8), dimension(ndays) :: check_init_mass_snow 		!kg/m2
	real(kind=8), dimension(ndays) :: check_init_mass_water 	!kg/m2
	real(kind=8), dimension(ndays) :: check_accum_snow 		!kg/m2
	real(kind=8), dimension(ndays) :: check_accum_water		!kg/m2
	real(kind=8), dimension(ndays) :: check_melted_snow		!kg/m2
	real(kind=8), dimension(ndays) :: check_runoff_water		!kg/m2
	real(kind=8), dimension(ndays) :: check_refreeze_water		!kg/m2
	real(kind=8), dimension(ndays) :: check_end_mass_snow		!kg/m2
	real(kind=8), dimension(ndays) :: check_end_mass_water		!kg/m2
	real(kind=8), dimension(ndays) :: check_ice2ice			!kg/m2

	real(kind=8):: dummy_melt
	real(kind=8):: dummy_runoff
	real(kind=8):: dummy_refreeze
	real(kind=8):: dummy_ice2ice


	! variables for energy conservation

	real(kind=8), dimension(ndays) :: check_init_energy_snow	!J/m2
	real(kind=8), dimension(ndays) :: check_init_energy_water 	!J/m2
	real(kind=8), dimension(ndays) :: check_end_energy_snow		!J/m2
	real(kind=8), dimension(ndays) :: check_end_energy_water 	!J/m2

	real(kind=8), dimension(ndays) :: check_surface_e_flux_obs 	 	!J/m2
	real(kind=8), dimension(ndays) :: check_surface_e_flux_diag 		!J/m2
	real(kind=8), dimension(ndays) :: check_e_accum_tot 		!J/m2

	real(kind=8), dimension(ndays) :: check_e_melt_tot 		!J/m2
	real(kind=8), dimension(ndays) :: check_e_melt_heat 		!J/m2
	real(kind=8), dimension(ndays) :: check_e_melt_qq 		!J/m2
	real(kind=8), dimension(ndays) :: check_e_melt_runoff		!J/m2

	real(kind=8), dimension(ndays) :: check_e_freeze_tot		!J/m2
	real(kind=8), dimension(ndays) :: check_e_freeze_heat 		!J/m2

	real(kind=8), dimension(ndays) :: check_e_perc_runoff 		!J/m2
	real(kind=8), dimension(ndays) :: check_e_ice2ice 		!J/m2
	real(kind=8) :: dummy_energy					!J/m2
	real(kind=8) :: dummy_heat					!J/m2
	real(kind=8) :: dummy_e_qq					!J/m2

   
	character(1) :: tab

 	real(kind=8) :: mass0
	real(kind=8) :: dmass

! ! ! ! 	! netCDF Accumulation
! ! ! ! 	!===================
! ! ! ! 	! cant be a parameter, cause the netcdf functions modifies it.
! ! ! ! 	integer :: filehandle_netcdf4
! ! ! ! 	! Variable IDs: height, bed
! ! ! ! 	!--------------------------
! ! ! ! 	! Variable for the netCDF snow mass, parameter
! ! ! ! 	integer :: accum_varid
! ! ! ! 		! Variable for the netCDF snow mass, parameter
! ! ! ! 	integer :: rain_varid
! ! ! ! 	! Variable for the netCDF density of snow, parameter
! ! ! ! 	integer :: refreezing_varid
! ! ! ! 	! Variable for the netCDF density of snow, parameter
! ! ! ! 	integer :: melt_varid
! ! ! ! 	! Variable for the netCDF snow, parameter
! ! ! ! 	integer :: snow_varid
! ! ! ! 	! Variable for the netCDF runoff of snow, parameter
! ! ! ! 	integer :: runoff_varid
! ! ! ! 	! Variable for the netCDF runoff of snow, parameter
! ! ! ! 	integer :: latent_heat_varid
! ! ! ! 	! Variable for the netCDF runoff of snow, parameter
! ! ! ! 	integer :: albedo_dynamic_varid

	
    ! file name for annual data
	integer :: filehandle_netcdf_annual
	character(128) :: netcdf_output_filename_annual
	! Variable IDs: for annual data
	!--------------------------
	! Variable for the netCDF snow mass, parameter
	
	integer :: snowman_annual_varid
	integer :: rho_snow_annual_varid
	integer :: snow_temp_annual_varid
	integer :: lwmass_annual_varid
	integer :: accum_annual_varid
		! Variable for the netCDF snow mass, parameter
	integer :: rain_annual_varid
	! Variable for the netCDF density of snow, parameter
	integer :: refreezing_annual_varid
	! Variable for the netCDF density of snow, parameter
	integer :: melt_annual_varid
	! Variable for the netCDF snow, parameter
	integer :: snow_annual_varid
	! Variable for the netCDF runoff of snow, parameter
	integer :: runoff_annual_varid
	! Variable for the netCDF runoff of snow, parameter
	integer :: latent_heat_annual_varid
	! Variable for the netCDF runoff of snow, parameter
	integer :: albedo_dynamic_annual_varid
	! Variable snow_mask
	integer :: snow_mask_annual_varid
	integer :: regridding_annual_varid
	integer :: melt_ice_annual_varid
	integer :: snow_temp_ave_surface_annual_varid
	integer :: real_mass_balance_annual_varid
	
	    ! file name for daily data
	integer :: filehandle_netcdf_daily
	character(128) :: netcdf_output_filename_daily
	! Variable IDs: for annual data
	!--------------------------
	! Variable for the netCDF snow mass, parameter
	
	integer :: snowman_daily_varid
	integer :: rho_snow_daily_varid
	integer :: snow_temp_daily_varid
	integer :: lwmass_daily_varid
	integer :: accum_daily_varid
		! Variable for the netCDF snow mass, parameter
	integer :: rain_daily_varid
	! Variable for the netCDF density of snow, parameter
	integer :: refreezing_daily_varid
	! Variable for the netCDF density of snow, parameter
	integer :: melt_daily_varid
	! Variable for the netCDF snow, parameter
	integer :: snow_daily_varid
	! Variable for the netCDF runoff of snow, parameter
	integer :: runoff_daily_varid
	! Variable for the netCDF runoff of snow, parameter
	integer :: latent_heat_daily_varid
	! Variable for the netCDF runoff of snow, parameter
	integer :: albedo_dynamic_daily_varid
	! Variable snow_mask
	integer :: snow_mask_daily_varid
	integer :: regridding_daily_varid
	integer :: melt_ice_daily_varid
	integer :: snow_temp_ave_surface_daily_varid
	integer :: real_mass_balance_daily_varid

	real(kind=8),dimension(nx,ny)::testi
	integer, dimension(365) :: dummy_vector=(/(i, i=1, 365)/)
	integer, dimension(ndays, n_snowlayer) :: dummy_matrix
! 	integer, dimension(nx,ny,n_snowlayer) :: dummy_test
! 	integer, dimension(n_snowlayer,ny,nx) :: dummy_test2
! 	real(kind=8) :: t1,t2
! 	do i=1,nx,1
!         do k=1,ny,1 
!             do j=1,n_snowlayer,1
!                 dummy_test(i,k,j)=i+k+j
!                 dummy_test2(j,k,i)=i+k+j
!             end do
!             end do
!             end do
    !calculate pressure coordinates of elevation
    p_air(:,:)= 101325*exp(-9.80665*0.0289644*elevation(:,:)/288.15/8.31447)
            
	do i=1,n_snowlayer,1
        dummy_matrix(:,i)=dummy_vector+i
    end do
	tab=char(09)
    albedo_runtime=0.
	! preparations
	mass0 = 0.
	masssum = 0.
	rainman = 0.
	dummy_melt=0.
	dummy_refreeze=0.
	accum = 0.
	K_lh = 0.
	H_lh = 0.
	K_sw = 0.
	dz=0.
	shift=0;
	vaporflux=0.
	!shift = int(floor(real(ndays)/real(monthly_data_num)/2.))
	!print*,shift
	day =   int(floor(real(ndays)/real(monthly_data_num)))
	
! 	albedo_dynamic=0.85
	
	if(albedo_module==0) then
	testi= read_variable(albedo_file_path, nx, ny, albedo_file_variable_name) ! TODO remove _special
	albedo_dynamic(1:nx,1:ny)=testi
	end if
	
    print*,'Albedo_module',albedo_module
    
    if (latent_heat_flux_on) then
        if (latent_heat_flux_analog_to_sensible_heat_flux) then
             D_lf=ratio*D_sf/cp_air*0.622*(L_v+L_lh)
        else
             D_lf=D_sf/cp_air*0.622*(L_v+L_lh)
        end if
    else
        D_lf=0
    end if
       
    print*,'Latent_heat_flux_on',latent_heat_flux_on

	check_init_mass_snow = 0. 		!kg/m2
	check_init_mass_water = 0. 		!kg/m2
	check_accum_snow  = 0.			!kg/m2
	check_accum_water = 0.			!kg/m2
	check_melted_snow = 0.			!kg/m2
	check_runoff_water = 0.			!kg/m2
	check_end_mass_snow = 0.		!kg/m2
	check_end_mass_water = 0.		!kg/m2
	check_refreeze_water = 0.		!kg/m2
	check_ice2ice = 0.			!kg/m2

	dummy_melt = 0.				!kg/m2
	dummy_runoff = 0.			!kg/m2
	dummy_refreeze = 0.			!kg/m2
	dummy_ice2ice = 0.			!kg/m2


	check_init_energy_snow  = 0.		!J/m2
	check_end_energy_snow  = 0.		!J/m2
	check_init_energy_water = 0.		!J/m2
	check_end_energy_water = 0.		!J/m2

	check_surface_e_flux_obs = 0.	 	!J/m2
	check_surface_e_flux_diag = 0.		!J/m2
	check_e_accum_tot = 0.			!J/m2

	check_e_freeze_tot = 0.			!J/m2
	check_e_freeze_heat = 0.		!J/m2

	check_e_melt_tot = 0.			!J/m2
	check_e_melt_heat = 0.			!J/m2
	check_e_melt_qq = 0.			!J/m2
	check_e_melt_runoff = 0.		!J/m2

	check_e_perc_runoff = 0.		!J/m2
	check_e_ice2ice = 0.			!J/m2

	dummy_energy = 0.			!J/m2
	dummy_heat = 0.				!J/m2
	dummy_e_qq = 0.				!J/m2
! ! ! 
! ! ! 	! monthly data? if yes the initialize an outputfile
! ! ! ! 	if ( (((monthly_data).and.(mod(year,monthly_data_freq)==0)) &
! ! ! ! 		.or. ((year .ge. daily_snowmask_start))).and.write_output_files_imhof ) then
! ! ! ! 		 
! ! ! ! 		allocate(m_snowman(nx,ny,n_snowlayer,monthly_data_num))
! ! ! ! 		allocate(m_lwmass(nx,ny,n_snowlayer,monthly_data_num))
! ! ! ! 		allocate(m_rho_snow(nx,ny,n_snowlayer,monthly_data_num))
! ! ! ! 		allocate(m_snow_temp(nx,ny,n_snowlayer,monthly_data_num))
! ! ! ! 		allocate(m_albedo_dynamic(nx,ny,monthly_data_num))
! ! ! ! 
! ! ! ! 		m_snowman = 0.
! ! ! ! 		m_lwmass = 0.
! ! ! ! 		m_rho_snow = rho_s
! ! ! ! 		m_snow_temp = 0.
! ! ! ! 		m_albedo_dynamic = 0.
! ! ! ! 
! ! ! ! 		! creat name for output file
! ! ! ! 		write (netcdf_output_filename_monthly, '( "Snowcover3D_monthly_", I7.7, ".nc" )') year
! ! ! ! 		! initialize output file
! ! ! ! 		call init_netcdf_3D_file(TRIM(adjustl(output_directory)) // netcdf_output_filename_monthly, filehandle_netcdf_monthly, &
! ! ! ! 		ny, nx, n_snowlayer, int(dy), int(dx), 1, snowman_varid, lwmass_varid, rho_snow_varid, snow_temp_varid, &
! ! ! ! 		albedo_dynamic_varid)
! ! ! ! 		
! ! ! ! 
! ! ! ! 		if(add_accum_output)then
! ! ! ! 			allocate(m_accum(nx,ny))
! ! ! ! 			allocate(m_rain(nx,ny))
! ! ! ! 			allocate(m_refreezing(nx,ny))
! ! ! ! 			allocate(m_melt(nx,ny))
! ! ! ! 			allocate(m_snow(nx,ny))
! ! ! ! 			allocate(m_runoff(nx,ny))
! ! ! ! 			allocate(m_vaporflux(nx,ny))
! ! ! ! 
! ! ! ! 			m_accum = 0.
! ! ! ! 			m_rain = 0.
! ! ! ! 			m_refreezing = 0.
! ! ! ! 			m_melt = 0.
! ! ! ! 			m_snow = 0.
! ! ! ! 			m_runoff = 0.
! ! ! ! 
! ! ! ! 
! ! ! ! 			write (netcdf_output_filename_monthly, '( "Internal_snow_balance_", I7.7, ".nc" )') year
! ! ! ! 			! initialize output file
! ! ! ! 			call init_netcdf_accum_file(TRIM(adjustl(output_directory)) // netcdf_output_filename_monthly, filehandle_netcdf4, &
! ! ! ! 			ny, nx, int(dy), int(dx), accum_varid, rain_varid, melt_varid, refreezing_varid, snow_varid, runoff_varid )
! ! ! ! 		end if
! ! ! ! 
! ! ! ! 	end if

! ! ! ! 	if ( ((year .ge. daily_snowmask_start).or.((monthly_data).and.(mod(year,monthly_data_freq)==0))).and.&
! ! ! !         (daily_snowmask).and.write_output_files_imhof  ) then
! ! ! ! 		! allocate variable
! ! ! ! 		allocate(snow_mask(nx,ny,ndays))
! ! ! ! 
! ! ! ! 		! 1 = ocean
! ! ! ! 		! 2 = land
! ! ! ! 		! 3 = snow
! ! ! ! 		! 4 = wet snow
! ! ! ! 		! 5 = perenial firn aquifer
! ! ! ! 
! ! ! ! 		! initialize variable
! ! ! ! 		snow_mask = 1	! 1 = ocean
! ! ! ! 		do id=1,ndays,1
! ! ! ! 			where (landmask(:,:) .ne. 1)
! ! ! ! 				snow_mask(:,:,id) = 2	! 2 = land
! ! ! ! 			end where
! ! ! ! 		end do
! ! ! ! 
! ! ! ! 		! creat name for output file
! ! ! ! 		write (netcdf_output_filename_snow_mask, '( "Mask_snow_", I7.7, ".nc" )') year
! ! ! ! 
! ! ! ! 		! initialize output file
! ! ! ! 		call init_netcdf_snow_mask(TRIM(adjustl(output_directory)) // netcdf_output_filename_snow_mask, filehandle_netcdf_snow_mask, &
! ! ! ! 		ny, nx, ndays, int(dy), int(dx), snow_mask_varid )
! ! ! ! 
! ! ! ! 	end if

	!Allocate daily output variables
	if ( (daily_data) .and.( ((year .ge. daily_data_start).and.(year .le. daily_data_end))&
	.or.(mod(year,daily_data_frequency)==0))  ) then
	print*,'Inititalizing files for daily data',year
	
        write (netcdf_output_filename_daily, '( "DAILY_", I7.7, ".nc" )') year
	
        call init_netcdf_snow_file_time_first(TRIM(adjustl(output_directory)) // netcdf_output_filename_daily, &
        filehandle_netcdf_daily, &
		ny, nx, n_snowlayer, 365, int(dy), int(dx), 1, snowman_daily_varid, lwmass_daily_varid, rho_snow_daily_varid, &
		snow_temp_daily_varid, albedo_dynamic_daily_varid, latent_heat_daily_varid, accum_daily_varid, &
		rain_daily_varid, melt_daily_varid, refreezing_daily_varid, snow_daily_varid, runoff_daily_varid, &
        melt_ice_daily_varid, snow_mask_daily_varid, snow_temp_ave_surface_daily_varid, regridding_daily_varid, &
        real_mass_balance_daily_varid)
		
		
		!allocate memory for annual data
		allocate(m_snowman_daily(ndays,n_snowlayer))
		allocate(m_lwmass_daily(ndays,n_snowlayer))
		allocate(m_rho_snow_daily(ndays,n_snowlayer))
		allocate(m_snow_temp_daily(ndays,n_snowlayer))
	
        allocate(m_albedo_dynamic_daily(ndays))
        allocate(m_snow_temp_ave_surface_daily(ndays))
        
        allocate(m_vaporflux_daily(ndays))
        allocate(m_accum_daily(ndays))
        allocate(m_rain_daily(ndays))
        allocate(m_refreezing_daily(ndays))
        allocate(m_melt_daily(ndays))
        allocate(m_snow_daily(ndays))
        allocate(m_runoff_daily(ndays))
        allocate(m_melt_ice_daily(ndays))
        
        allocate(m_snow_mask_daily(ndays))
        allocate(m_regridding_daily(ndays))
        allocate(m_real_mass_balance_daily(ndays))

        m_snowman_daily=0.
        m_lwmass_daily=0.
        m_rho_snow_daily=0.
        m_snow_temp_daily=0.
        m_albedo_dynamic_daily=0.
        m_snow_temp_ave_surface_daily=0.
        m_vaporflux_daily=0.
        m_accum_daily=0.
        m_rain_daily=0.
        m_refreezing_daily=0.
        m_melt_daily=0.
        m_snow_daily=0.
        m_runoff_daily=0.
        m_snow_mask_daily=0.
        m_melt_ice_daily=0.
        m_real_mass_balance_daily=0.
        m_regridding_daily=0.
        
	
	end if
	
    if ( (monthly_data) .and.( ((year .ge. daily_data_start).and.(year .le. daily_data_end))&
	.or.(mod(year,monthly_data_frequency)==0))  ) then
	print*,'monthly',year
	end if
	
	
	!Initialize annual output file and annual data storage
    if ( (annual_data) .and.( ((year .ge. daily_data_start).and.(year .le. daily_data_end))&
	.or.(mod(year,annual_data_frequency)==0))  ) then
	print*,'Inititalizing files for annual data',year
	
        write (netcdf_output_filename_annual, '( "ANNUAL_", I7.7, ".nc" )') year
        print*,output_directory
        call init_netcdf_snow_file(TRIM(adjustl(output_directory)) // netcdf_output_filename_annual, filehandle_netcdf_annual, &
		ny, nx, n_snowlayer, 1, int(dy), int(dx), 1, snowman_annual_varid, lwmass_annual_varid, rho_snow_annual_varid, &
		snow_temp_annual_varid, albedo_dynamic_annual_varid, latent_heat_annual_varid, accum_annual_varid, &
		rain_annual_varid, melt_annual_varid, refreezing_annual_varid, snow_annual_varid, runoff_annual_varid, &
        melt_ice_annual_varid, snow_mask_annual_varid, snow_temp_ave_surface_annual_varid, regridding_annual_varid, &
        real_mass_balance_annual_varid)
		
		!allocate memory for annual data
		allocate(m_snowman_annual(nx,ny,n_snowlayer))
		allocate(m_lwmass_annual(nx,ny,n_snowlayer))
		allocate(m_rho_snow_annual(nx,ny,n_snowlayer))
		allocate(m_snow_temp_annual(nx,ny,n_snowlayer))
	
        allocate(m_albedo_dynamic_annual(nx,ny))
        allocate(m_snow_temp_ave_surface_annual(nx,ny))
        
        allocate(m_vaporflux_annual(nx,ny))
        allocate(m_accum_annual(nx,ny))
        allocate(m_rain_annual(nx,ny))
        allocate(m_refreezing_annual(nx,ny))
        allocate(m_melt_annual(nx,ny))
        allocate(m_snow_annual(nx,ny))
        allocate(m_runoff_annual(nx,ny))
        allocate(m_melt_ice_annual(nx,ny))
        
        allocate(snow_mask_annual(nx,ny))
        allocate(m_regridding_annual(nx,ny))
        allocate(m_real_mass_balance_annual(nx,ny))

        m_snowman_annual=0.
        m_lwmass_annual=0.
        m_rho_snow_annual=0.
        m_snow_temp_annual=0.
        m_albedo_dynamic_annual=0.
        m_snow_temp_ave_surface_annual=0.
        m_vaporflux_annual=0.
        m_accum_annual=0.
        m_rain_annual=0.
        m_refreezing_annual=0.
        m_melt_annual=0.
        m_snow_annual=0.
        m_runoff_annual=0.
        snow_mask_annual=0.
        m_melt_ice_annual=0.
        m_real_mass_balance_annual=0.
        
	end if
! ))))))))))))))))))) MAIN LOOP (((((((((((((((((((((((((((
	!PRIVATE(ix,iy,H_lh,K_lh,masssum,accum,rainman,china_syndrome,K_sw,dz,month) SHARED(landmask,air_temp_ice,precip_ice,snowman,snow_temp,rho_snow,lwmass,P_sun,smb_ice)
! $OMP PARALLEL 
! $OMP DO
	!if((year.le.20).and.(mod(year,10)==0)) then
	!smb_ice = 0.



	if ( (year .lt. start_speed_up ) .or. (year .ge. end_speed_up) )then 
		fast_calculation = .false.
	end if
	print*, 'starting snowbern...' 
! 	call cpu_time(tn)

	call cpu_time(t1)
	do iy=1,ny,1 		!do ix=placex,placex,1
		do ix=1,nx,1	!do iy=placey,placey,1

			!if((landmask(ix,iy).eq. 0).or. (seafloor(ix,iy) .gt. sealevel) )then
			if((landmask(ix,iy).ne. 1) )then
				if ( (fast_calculation(ix,iy) .eqv. .false.) .or. ( mod(myyear(it),calc_rate_snow)==0) .or. (spinup .eqv. .false.) )then 
								! .or. (year .lt. start_speed_up ) .or. (year .ge. end_speed_up)
					

					smb_ice(ix,iy) = 0.
					mass0 = sum(snowman(ix,iy,:))

					month=1

					! reset fast calculation to false everywhere.
					fast_calculation(ix,iy) = .true.
                    

			! 	LOOP STARTS HERE
			!=============================
					do time=1,ndays,1
! 					if (year>2) then
! 						print *, 'x=',ix,' y=',iy,'t=',time
!                     end if
						!print *,snowman(ix,iy,:)



						check_init_mass_snow(time) = check_init_mass_snow(time) + sum(snowman(ix,iy,:))		
						check_init_mass_water(time) = check_init_mass_water(time) + sum(lwmass(ix,iy,:))
						check_init_energy_snow(time) = check_init_energy_snow(time) + sum(snowman(ix,iy,:)*snow_temp(ix,iy,:))*c_i
 						check_init_energy_water(time) = check_init_energy_water(time) + sum(lwmass(ix,iy,:))*(kelvin*c_i+L_lh)


						!if(landmask(ix,iy).ne.1) then
!                      

						!   if((ix==placex).and.(iy==placey)) then
						!print *, 'mass :', snowman(placex,placey,1:min(int(6),n_snowlayer))
						!print *, 'lwmass :', lwmass(placex,placey,1:min(int(6),n_snowlayer))
						!print *, 'temp :', snow_temp(placex,placey,1:min(int(6),n_snowlayer))
						!print *, 'rho:', rho_snow(placex,placey,1:min(int(6),n_snowlayer))
						!end if



						!-------------------------------------------------------------------------
						!! Accumulation
						!-------------------------------------------------------------------------
						! this part determines the amount [kg/m2] of snow that falls on a gridcell
						! and adjusts the inner energy of the gridcell i.e. T. in case of rain the
						! fraction that can freeze shall be determined and inner energy will be
						! adjusted too.

						H_lh=0.
						K_lh=0.
						accum=0.
						rainman=0.
						dummy_melt = 0.				!kg/m2
                        dummy_runoff = 0.			!kg/m2
                        dummy_refreeze = 0.			!kg/m2
                        !dummy_ice2ice = 0.			!kg/m2
                        dummy_melt_ice=0. !kg/m2
                        dummy_rain_ice=0. !kg/m2
                        dummy_regrid=0.

!                         print*,snowman(ix,iy,:)

						if((air_temp_ice(ix,iy,time)-kelvin<=snow_fall_temperature).and.(precip_ice(ix,iy,time)>0.)) then
							! IN CASE OF SNOW
							! acumulation = precipitation that falls at temperatures smaller than snow_fall_temperature C

							if(snowman(ix,iy,1)==0) then ! seeding temperature for first snow.
								snow_temp(ix,iy,1)=air_temp_ice(ix,iy,time)	! min(0.,airtemp(ix,iy,time))
							end if

							! accumulated snow during one timestep in kg/m2/s
							accum=precip_ice(ix,iy,time)*rho_w
							
! 							if(precip_ice(ix,iy,time)>0.1) then
!                             if (accum*dt_firn>1.) then
!                             print*,accum*dt_firn
							nday_snowfall=time
							if(albedo_module>0) then
							!based on oerlemans depth scale inverse 
							albedo_dynamic(ix,iy)=min(albedo_snow_new,albedo_dynamic(ix,iy)+&
							(albedo_snow_new-albedo_snow_wet)*(1-exp(-accum*dt_firn/3.)))
							!alternative also strange
							!albedo_dynamic(ix,iy,1)=min(albedo_snow_new,albedo_dynamic(ix,iy,1)+(albedo_snow_new-albedo_dynamic(ix,iy,1))*(1-exp(-accum*dt_firn))
							end if
! 							end if
! 							end if
! 							print*,nday_snowfall
							rainman=0.
							masssum=snowman(ix,iy,1)+accum*dt_firn
							rho_snow(ix,iy,1) = masssum / (snowman(ix,iy,1)/rho_snow(ix,iy,1) + accum*dt_firn/rho_s)
                            
							snowman(ix,iy,1) = masssum

							! latent heat flux for Energybalance
							H_lh=precip_ice(ix,iy,time)*rho_w*c_i
							K_lh=precip_ice(ix,iy,time)*rho_w*c_i*air_temp_ice(ix,iy,time)


							! IN CASE OF RAIN
							!
						else if(snowman(ix,iy,1) .gt. 0.) then
							accum=0.
							rainman = precip_ice(ix,iy,time)*rho_w ! kg/m2/s
							lwmass(ix,iy,1) = lwmass(ix,iy,1) + rainman*dt_firn

							! latent heat flux for Energybalance
							H_lh=0.
							K_lh=precip_ice(ix,iy,time)*rho_w*c_w*(air_temp_ice(ix,iy,time)-kelvin)

						end if ! end of accumulation



! ! ! ! 						if (  ((year .ge. daily_snowmask_start).or.((monthly_data).and.(mod(year,monthly_data_freq)==0))).and.(daily_snowmask)   & 
! ! ! ! 								 .and. (add_accum_output).and.write_output_files_imhof)then
! ! ! ! 							m_snow(ix,iy) = m_snow(ix,iy) + accum*dt_firn
! ! ! ! 							m_rain(ix,iy) = m_rain(ix,iy) + rainman*dt_firn
! ! ! ! 						end if
						
						!THIS is a if to write

!                         print*,snowman(ix,iy,:)
						if(snowman(ix,iy,1)>0.) then ! do further calculations only where there is snow
! 							print*,'inside snow loop'
							!lwmass(ix,iy,:)=0
							!rainman=0
							check_accum_snow(time) = check_accum_snow(time) + accum*dt_firn
							check_accum_water(time) = check_accum_water(time) + rainman*dt_firn
							check_e_accum_tot(time) = check_e_accum_tot(time) + dt_firn*accum*c_i*snow_temp(ix,iy,1) &
									+ dt_firn*rainman*(c_i*kelvin+L_lh)

							!----------------------------------------------------------------------
							!! Splitting / fusing of boxes
							!----------------------------------------------------------------------
							if((snowman(ix,iy,1).gt. upper_massbound).or.(snowman(ix,iy,1).lt. lower_massbound)) then
								call go_regrid(ix,iy, snowman(ix,iy,:), snow_temp(ix,iy,:), &
									rho_snow(ix,iy,:), lwmass(ix,iy,:), dummy_regrid)

								do while((snowman(ix,iy,1).gt. upper_massbound))
									call go_regrid(ix,iy, snowman(ix,iy,:), snow_temp(ix,iy,:), &
										rho_snow(ix,iy,:), lwmass(ix,iy,:), dummy_regrid)
									!print*,'heavy snowfall at', ix,iy, 'in Year',year, 'day',time
								end do
							end if


							!----------------------------------------------------------------------
							! Densification
							!----------------------------------------------------------------------
							! calculate new density of every layer
							if ((densification_model) .and. (snowman(ix,iy,3).gt. 0.)) then
								call go_densification(snowman(ix,iy,:), rho_snow(ix,iy,:),&
									snow_temp(ix,iy,:), accum+rainman)
							end if
							!-------------------------------------------------------------------------
							!! Energy Fluxes
							!-------------------------------------------------------------------------
							dz(:)=snowman(ix,iy,:)/rho_snow(ix,iy,:)

!                             if ((ix==39) .and. (iy==74)) then
! 							print*,'albedo', albedo_dynamic(ix,iy,1)
! 							print*,snow_temp(ix,iy,1)
! 							print*,nday_snowfall
! 							end if
                            
                            !calculate lwc
                            lwc=lwmass(ix,iy,1)/snowman(ix,iy,1)/rho_w/(1./rho_snow(ix,iy,1) - 1./rho_i)
                            
!                             print*, 'Starting solar calculation'
							call go_calculate_solar(ix,iy,snow_temp(ix,iy,1), air_temp_ice(nx, ny, time), P_sun(ix,iy,time),&
							albedo_dynamic(ix,iy), nday_snowfall, time, K_sw, albedo_module, lwc)
							!, albedo_runtime(ix,iy,:))
							
! 							if ((ix==39) .and. (iy==74)) then
! 							print*,'albedo', albedo_dynamic(ix,iy)
! 							end if
! 							
							! simulate dry and wet snow classic michael
! 							if(snow_temp(ix,iy,1) .lt. kelvin) then
! 								K_sw=P_sun(ix,iy,time)*(1.-albedo_snow_new)
! 							else
! 								K_sw=P_sun(ix,iy,time)*(1.-albedo_snow_wet)
! 							end if
! 							!print *,'before energy flux'

							dummy_energy= sum(snow_temp(ix,iy,:)*snowman(ix,iy,:))*c_i + sum(lwmass(ix,iy,:))*(kelvin*c_i+L_lh)
                            
                            !Energy flux calculation including or excluding Latent heat flux
!                               print*, 'Starting energy calculation'

                            if (latent_heat_flux_on) then
                                call go_energy_flux_tobias(ix,iy,snow_temp(ix,iy,:), rho_snow(ix,iy,:), china_syndrome, &
								air_temp_ice(ix,iy,time), snowman(ix,iy,1), dz, K_sw, H_lh, K_lh, Q_heat, &
								dummy_heat, DewpT(ix,iy,time), vaporflux, D_lf, p_air(ix,iy))
! 								print*,snowman(ix,iy,1)
                                ! 								print*,vaporflux/(L_v+L_lh)
! 							if(snow_temp(ix,iy,1).lt.kelvin) then	
!                                 print*,'pre',snowman(ix,iy,1)
                                snowman(ix,iy,1)=max(0.,snowman(ix,iy,1)+vaporflux/(L_v+L_lh)*dt_firn)
! !                                 print*,'after',snowman(ix,iy,1)
!                             else 
!                                 lwmass(ix,iy,1)=lwmass(ix,iy,1)+vaporflux/(L_v)*dt_firn
!                             end if
                            else
                                call go_energy_flux_new(ix,iy,snow_temp(ix,iy,:), rho_snow(ix,iy,:), china_syndrome, &
								air_temp_ice(ix,iy,time), snowman(ix,iy,1), dz, K_sw, H_lh, K_lh, Q_heat, dummy_heat)
								vaporflux=0.
                            end if

							check_surface_e_flux_obs(time) = check_surface_e_flux_obs(time) - dummy_energy + &
								sum(snow_temp(ix,iy,:)*snowman(ix,iy,:))*c_i + sum(lwmass(ix,iy,:))*(kelvin*c_i+L_lh)

							check_surface_e_flux_diag(time) = check_surface_e_flux_diag(time) + dummy_heat
							!print*,'after energy flux'

							!----------------------------------------------------------------------
							!! Melting snow
							!----------------------------------------------------------------------
							!if((ix==291).and.(iy==279)) then
							!	print*,'b snowtemp 1', snow_temp(ix,iy,1) , snow_temp(ix,iy,2), snow_temp(ix,iy,3)
							!end if

							if(china_syndrome) then
								fast_calculation(ix,iy) = .false.

								!K_sw=P_sun(ix,iy,time)*(1.-albedo_snow_wet)
 
								dummy_energy = sum(snowman(ix,iy,:)*snow_temp(ix,iy,:))*c_i + sum(lwmass(ix,iy,:))*(kelvin*c_i+L_lh)

					   			call go_melting_snow(ix,iy,time,snow_temp(ix,iy,:), snowman(ix,iy,:), lwmass(ix,iy,:), rho_snow(ix,iy,:),&
					 				 air_temp_ice(ix,iy,time), K_sw, K_lh, H_lh, Q_heat, dummy_melt, dummy_runoff, smb_ice(ix,iy), &
					 				 dummy_e_qq, dummy_heat, vaporflux, D_lf, dummy_melt_ice,  dummy_regrid)

								check_melted_snow(time) = check_melted_snow(time) + dummy_melt
								check_runoff_water(time) = check_runoff_water(time) + dummy_runoff
								check_e_melt_tot(time) = check_e_melt_tot(time) - dummy_energy + ( sum(snowman(ix,iy,:)*snow_temp(ix,iy,:))*c_i + &
									 sum(lwmass(ix,iy,:))*(kelvin*c_i+L_lh)) 
								check_e_melt_heat(time) = check_e_melt_heat(time) + dummy_heat
								check_e_melt_qq(time) = check_e_melt_qq(time) + dummy_e_qq
								check_e_melt_runoff(time) = check_e_melt_runoff(time) + dummy_runoff*(c_i*kelvin+L_lh)

! ! ! ! 								if (  (((year .ge. daily_snowmask_start).or.((monthly_data).and.(mod(year,monthly_data_freq)==0))).and.(daily_snowmask) )  & 
! ! ! ! 									.and. (add_accum_output).and.write_output_files_imhof)then
! ! ! ! 									m_runoff(ix,iy) = m_runoff(ix,iy) + dummy_runoff
! ! ! ! 									m_melt(ix,iy) = m_melt(ix,iy) + dummy_melt
! ! ! ! 								end if


								!lwmass(ix,iy,:)=0
							end if
                            

							!if((ix==291).and.(iy==279)) then
							!	print*,'c snowtemp 1', snow_temp(ix,iy,1) , snow_temp(ix,iy,2), snow_temp(ix,iy,3)
							!end if

							if(maxval(lwmass(ix,iy,:))>0.) then
								!----------------------------------------------------------------------
								! Water Percolation
								!----------------------------------------------------------------------
								fast_calculation(ix,iy) = .false.
								!print*, 'fail'
								call go_percolation(snowman(ix,iy,:), lwmass(ix,iy,:), rho_snow(ix,iy,:), &
									dummy_runoff) !snow_temp(ix,iy,:),
								check_runoff_water(time) = check_runoff_water(time) + dummy_runoff
								check_e_perc_runoff(time) = check_e_perc_runoff(time) + dummy_runoff*(c_i*kelvin+L_lh)
								!----------------------------------------------------------------------
								!! Refreezing
								!----------------------------------------------------------------------

								!if((ix==291).and.(iy==279)) then
								!	print*,'d snowtemp 1', snow_temp(ix,iy,1) , snow_temp(ix,iy,2), snow_temp(ix,iy,3)
								!end if
								dummy_energy = sum(snowman(ix,iy,:)*snow_temp(ix,iy,:))*c_i + sum(lwmass(ix,iy,:))*(kelvin*c_i+L_lh)
								call go_refreezing(lwmass(ix,iy,:), snowman(ix,iy,:), rho_snow(ix,iy,:), &
									snow_temp(ix,iy,:), dummy_refreeze, dummy_heat, ix, iy)

								check_e_freeze_heat(time) = check_e_freeze_heat(time) + dummy_heat
								check_refreeze_water(time) = check_refreeze_water(time) + dummy_refreeze
								check_e_freeze_tot(time) = check_e_freeze_tot(time) - dummy_energy + ( sum(snowman(ix,iy,:)*&
									snow_temp(ix,iy,:))*c_i + sum(lwmass(ix,iy,:))*(kelvin*c_i+L_lh)) 

								!if((ix==291).and.(iy==279)) then
								!	print*,'e snowtemp 1', snow_temp(ix,iy,1) , snow_temp(ix,iy,2), snow_temp(ix,iy,3)
								!end if

! ! ! ! 								if ( (((year .ge. daily_snowmask_start).or.((monthly_data).and.(mod(year,monthly_data_freq)==0))).and.(daily_snowmask) )  & 
! ! ! ! 									.and. (add_accum_output).and.write_output_files_imhof)then
! ! ! ! 									m_runoff(ix,iy) = m_runoff(ix,iy) + dummy_runoff
! ! ! ! 									m_refreezing(ix,iy) = m_refreezing(ix,iy) + dummy_refreeze
! ! ! ! 								end if


							end if



						! $$$$$$$$$$$$$$ FIRN PART ENDS HERE $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
                        ! if landmask is ice melt ice , ignore error of access ice melt in case of snowfall on bareground and entire melting gridcell
						else ! if(landmask(ix,iy) == 0) then 
							!----------------------------------------------------------------------
							!! Melting ice
							!----------------------------------------------------------------------
							! TODO: the following four lines contribute each to 7% of the calculation time 
							! of this SMB subroutine. if I programed everything properly they can be 
							! commented out without consequences. I hope...
							!if ( sum(lwmass(ix,iy,:))**2. .gt. 0.)then
							!	print*, 'lwmass', sum(lwmass(ix,iy,:))
							!end if
							!if ( sum(snowman(ix,iy,:))**2. .gt. 0.)then
							!	print*, 'snowman', sum(snowman(ix,iy,:))
							!end if
							!if ( sum(snow_temp(ix,iy,:))**2. .gt. 0.)then
							!	print*, 'snow_temp', sum(snow_temp(ix,iy,:))
							!end if
							!if ( (sum(rho_snow(ix,iy,:))- real(n_snowlayer,8)*rho_s)**2. .gt. 0.   )then
							!	print*, 'rho_snow', (sum(rho_snow(ix,iy,:))- real(n_snowlayer,8)*rho_s )
							!end if

!                             print*,'ice part'
							!lwmass(ix,iy,:) = 0.
							!snowman(ix,iy,:) = 0.
							!snow_temp(ix,iy,:) = 0.
							!rho_snow(ix,iy,:) = rho_s
							! assume that blank ice has a surface temperature 0C.
							call go_melting_ice(smb_ice(ix,iy), air_temp_ice(ix,iy,time), precip_ice(ix,iy,time),&
								 P_sun(ix,iy,time), vaporflux, D_lf, dummy_melt_ice, dummy_rain_ice, DewpT(ix,iy,time), &
								 p_air(ix,iy))
                            albedo_dynamic(ix,iy)=albedo_ice

						end if ! end of if that processes existant snow OR blank ice in 'else'

						! some safety output in case something doesnt go well
						!if(minval(lwmass(ix,iy,:))<0)then
						!	print*,'minimum lwmass at', ix, iy,' is ', minval(lwmass(ix,iy,:)),'time', time
						!end if
						!if(minval(snowman(ix,iy,:))<0)then
						!	print*,'minimum snowmas  at', ix, iy,  'is ', minval(snowman(ix,iy,:)),'time', time
						!end if

						!end if ! end of the if that looks for lan
! 						print*,'about to calculate saving values'
                        !Creating output datas which are scalars on a point bases
                        !sum up daily values
                        if ( (annual_data) .and.( ((year .ge. daily_data_start).and.(year .le. daily_data_end))&
                            .or.(mod(year,annual_data_frequency)==0))  ) then
!                             print*,'inside daily sums'
                            m_snow_annual(ix,iy)=m_snow_annual(ix,iy) + real(accum*dt_firn)
                            m_rain_annual(ix,iy)=m_rain_annual(ix,iy) + real(rainman*dt_firn)   !rain on ice?
                            m_vaporflux_annual(ix,iy)=m_vaporflux_annual(ix,iy) + real(vaporflux*dt_firn/(L_v+L_lh))
                            m_runoff_annual(ix,iy)=m_runoff_annual(ix,iy) + real(dummy_runoff)
                            m_refreezing_annual(ix,iy)= m_refreezing_annual(ix,iy) + real(dummy_refreeze)
                            m_melt_annual(ix,iy)= m_melt_annual(ix,iy) + real(dummy_melt) !melt on ice
                            m_melt_ice_annual(ix,iy)=m_melt_ice_annual(ix,iy) + real(dummy_melt_ice)
                            albedo_runtime(ix,iy)=albedo_runtime(ix,iy)+albedo_dynamic(ix,iy)
                            m_snow_temp_ave_surface_annual(ix,iy)=m_snow_temp_ave_surface_annual(ix,iy)+snow_temp(ix,iy,1)
                            m_regridding_annual(ix,iy)=m_regridding_annual(ix,iy)+dummy_regrid
                            if (include_values_over_ice) then
                                m_runoff_annual(ix,iy)=m_runoff_annual(ix,iy) + real(dummy_melt_ice)
                                m_rain_annual(ix,iy)=m_rain_annual(ix,iy) + real(dummy_rain_ice*dt_firn)
                            end if
                            if (snowman(ix,iy,1) .gt. 0.) then
                                snow_mask_annual(ix,iy)=snow_mask_annual(ix,iy)+1
                            end if
                            
                            
                        end if
                        
                        if ( (daily_data) .and.( ((year .ge. daily_data_start).and.(year .le. daily_data_end))&
                                .or.(mod(year,daily_data_frequency)==0))  ) then
!                             print*,'inside daily sums'
                            
                            m_snowman_daily(time,1:15)=real(snowman(ix,iy,:))  !dummy_test(ix,iy,:) dummy_test2(:,iy,ix)
                            m_lwmass_daily(time,1:15)=real(lwmass(ix,iy,:))
                            m_rho_snow_daily(time,1:15)=real(rho_snow(ix,iy,:))
                            m_snow_temp_daily(time,1:15)=real(snow_temp(ix,iy,:))
                            
                            
                            m_albedo_dynamic_daily(time)=real(albedo_dynamic(ix,iy))
                            m_vaporflux_daily(time)=real(vaporflux)
                            m_accum_daily(time)=smb_ice(ix,iy)*rho_ice*seconds_per_year
                            m_rain_daily(time)=real(rainman*dt_firn)
                            m_refreezing_daily(time)=real(dummy_refreeze)
                            m_melt_daily(time)=real(dummy_melt)
                            m_snow_daily(time)=real(accum*dt_firn)
                            m_runoff_daily(time)=real(dummy_runoff)
                            m_melt_ice_daily(time)=real(dummy_melt_ice)
                            if (time.eq.1) then
                                m_real_mass_balance_daily(time)=real(sum(snowman(ix,iy,:)))-&
                                mass0-real(dummy_melt_ice)
                            else
                                m_real_mass_balance_daily(time)=real(sum(snowman(ix,iy,:)))-&
                                    m_real_mass_balance_daily(time-1)-real(dummy_melt_ice)
                            end if
                            m_regridding_daily(time)=dummy_regrid
                            m_snow_temp_ave_surface_daily(time)=sum(m_real_mass_balance_daily(1:time))
                            if (snowman(ix,iy,1) .gt. 0.) then !extra for wet snow else useless as snowman, general not very useful
                                m_snow_mask_daily(time)=1
                            end if
                            if (include_values_over_ice) then
                                m_runoff_daily(time)=m_runoff_daily(time) + real(dummy_melt_ice)
                                m_rain_daily(time)=m_rain_daily(time) + real(dummy_rain_ice*dt_firn)
                            end if
!                             print*, snowman(20,78,1), time

!             m_albedo_dynamic_daily=0.
!             m_snow_temp_ave_surface_daily=0.!delete
!             m_vaporflux_daily=0.
!             m_accum_daily=0.
!             m_rain_daily=0.
!             m_refreezing_daily=0.
!             m_melt_daily=0.
!             m_snow_daily=0.
!             m_runoff_daily=0.
!             m_snow_mask_daily=0.
!             m_melt_ice_daily=0.
!             m_real_mass_balance_daily=0.

!                             if (snowman(ix,iy,1) .gt. 0.) then
!                                 snow_mask_annual(ix,iy)=snow_mask_annual(ix,iy)+1
!                             end if
!                             
                            
                        end if
                        
!                         print*,albedo_runtime(ix,iy)
!                         print*,m_vaporflux(ix,iy)
!                         print*,vaporflux
!                         albedo_runtime(ix,iy,1:n_snowlayer)=albedo_runtime(ix,iy,1:n_snowlayer)+albedo_dynamic(ix,iy,1:n_snowlayer)
                        
! ! ! ! 						print*,'about to save snow'
! ! ! ! 						if ( (  ((year .ge. daily_snowmask_start)).or.((monthly_data).and.(mod(year,monthly_data_freq)==0)) ) & 
! ! ! ! 							.and. ( mod((time+shift ),day)==0).and.write_output_files_imhof)then
! ! ! ! 							 
! ! ! ! 							! collect montyly data in cache array
! ! ! !                             
! ! ! ! 							m_snowman(ix,iy,1:n_snowlayer,month) = real(snowman(ix,iy,1:n_snowlayer),4)
! ! ! ! 							m_lwmass(ix,iy,1:n_snowlayer,month) = real(lwmass(ix,iy,1:n_snowlayer),4)
! ! ! ! 							m_rho_snow(ix,iy,1:n_snowlayer,month) = real(rho_snow(ix,iy,1:n_snowlayer),4)
! ! ! ! 							m_snow_temp(ix,iy,1:n_snowlayer,month) = real(snow_temp(ix,iy,1:n_snowlayer),4)
! ! ! ! 							m_albedo_dynamic(ix,iy,month) = real(albedo_runtime(ix,iy),4)
! ! ! ! 		
! ! ! ! 							!call writeNCDF3DGridValues(filehandle_netcdf_monthly, month, snowman_varid, snowman, ny, nx, n_snowlayer)
! ! ! ! 							!call writeNCDF3DGridValues(filehandle_netcdf_monthly, month, lwmass_varid, lwmass, ny, nx, n_snowlayer)
! ! ! ! 							!call writeNCDF3DGridValues(filehandle_netcdf_monthly, month, rho_snow_varid, rho_snow, ny, nx, n_snowlayer)
! ! ! ! 							!call writeNCDF3DGridValues(filehandle_netcdf_monthly, month, snow_temp_varid, snow_temp, ny, nx, n_snowlayer)
! ! ! ! 							month=month+1
! ! ! ! 						end if
! ! ! ! 						print*,'saved snow'

						if(time .lt. ndays) then
							check_end_mass_snow(time) = check_end_mass_snow(time) + sum(snowman(ix,iy,:))
							check_end_mass_water(time) = check_end_mass_water(time) + sum(lwmass(ix,iy,:))

							check_end_energy_snow(time) = check_end_energy_snow(time) + sum(snowman(ix,iy,:)*snow_temp(ix,iy,:))*c_i 
							check_end_energy_water(time) = check_end_energy_water(time) + sum(lwmass(ix,iy,:))*(kelvin*c_i+L_lh)
						end if

! 
! ! ! ! 						if (  ((year .ge. daily_snowmask_start).or.((monthly_data).and.&
! ! ! ! 						(mod(year,monthly_data_freq)==0))).and.(daily_snowmask).and.write_output_files_imhof   ) then
! ! ! ! 
! ! ! ! 							if (snowman(ix,iy,1) .gt. 0.) then
! ! ! ! 								if( maxval(lwmass(ix,iy,:)) .gt. 0.) then
! ! ! ! 									snow_mask(ix,iy,time) = 4	! wet snow
! ! ! ! 								else
! ! ! ! 									snow_mask(ix,iy,time) = 3	! dry snow
! ! ! ! 								end if	
! ! ! ! 							end if
! ! ! 
! ! ! ! 							! if the year is over, check if this snow column contains a perenial firn aquifer
! ! ! ! 							if (time .eq. ndays) then
! ! ! ! 								if (minval(snow_mask(ix,iy,:)) .eq. 4) then
! ! ! ! 									snow_mask(ix,iy,:) = 5		! pfa
! ! ! ! 								end if
! ! ! ! 							end if


! ! ! !                             end if
                        
!                         call writeNCDFSNOW2D_pointwise_time_interval(filehandle_netcdf_annual, time, snow_mask_annual_varid, &
!                                     real(snow_mask(ix,iy,time)), iy, ix, n_snowlayer, time)
                        
!                         print*,'end of timestep'
!                             print*,'vaporflux',vaporflux
					end do ! end of the loop over the timesteps of one year


					!calculation of annual mass balance has to be done prior to the shifing of the mass balance, and it can't nicely be calculated within the timesteps
                    if ( (annual_data) .and.( ((year .ge. daily_data_start).and.(year .le. daily_data_end))&
                            .or.(mod(year,annual_data_frequency)==0))  ) then
                    m_real_mass_balance_annual(ix,iy)=real(sum(snowman(ix,iy,:))-&
                    mass0-m_melt_ice_annual(ix,iy))
                    end if
					!  CALCULATE MASS BALANCE FOR ICE MODEL
					!=======================================

					! at the end of one year, limit the maximum mass in one snowcolumn to upper_massbound*n_snowlayer. 
					! all surplus snow is moved to the ice model
					dmass = sum( snowman(ix,iy,:)) - n_snowlayer * soll_mass *1.5    ! upper_massbound
					if ((dmass .gt. 0.)) then

						smb_ice(ix,iy) = max(dmass,0.)/rho_ice/seconds_per_year !+ smb_ice(ix,iy)

						!snowman(:,:,n_snowlayer) = 0.
						!lwmass(:,:,n_snowlayer) = 0.
						!snow_temp(:,:,n_snowlayer) = 0.
						!rho_snow(:,:,n_snowlayer) = rho_s

						!print*, '1 before continous depleet',ix,iy
						call continous_snowman_depleet(snowman(ix,iy,:), rho_snow(ix,iy,:), &
							snow_temp(ix,iy,:), lwmass(ix,iy,:), dmass, dummy_ice2ice, dummy_runoff, dummy_energy) ! TODO smb_ice

						check_ice2ice(ndays) = check_ice2ice(ndays) + dummy_ice2ice
						check_e_ice2ice(ndays) = check_e_ice2ice(ndays) + dummy_energy

						check_runoff_water(ndays) = check_runoff_water(ndays) + dummy_runoff
						check_e_perc_runoff(ndays) = check_e_perc_runoff(ndays) + dummy_runoff*(c_i*kelvin+L_lh)
						
! ! ! ! 						if ( (   ((year .ge. daily_snowmask_start).or.((monthly_data).and.&
! ! ! !                             (mod(year,monthly_data_freq)==0))).and.(daily_snowmask).and.write_output_files_imhof   ) & 
! ! ! ! 							.and. (add_accum_output))then
! ! ! ! 							m_runoff(ix,iy) = m_runoff(ix,iy) + dummy_runoff
! ! ! ! 							!m_accum(ix,iy) = m_accum(ix,iy) + dummy_ice2ice
! ! ! ! 						end if
						

					else
						fast_calculation(ix,iy) = .false. ! slow calculation if snow grid is not yet filled with snow
						!fast_calculation(ix,iy) = .true.
					end if

					check_end_mass_snow(ndays) = check_end_mass_snow(ndays) + sum(snowman(ix,iy,:))
					check_end_mass_water(ndays) = check_end_mass_water(ndays) + sum(lwmass(ix,iy,:))

					check_end_energy_snow(ndays) = check_end_energy_snow(ndays) + sum(snowman(ix,iy,:)*snow_temp(ix,iy,:))*c_i
					check_end_energy_water(ndays) = check_end_energy_water(ndays) + sum(lwmass(ix,iy,:))*(kelvin*c_i+L_lh)	


					! on land fast coupling, unless there is ice nearby. 
					if ((smb_ice(ix,iy).lt. 0.).and.(landmask(ix,iy).eq. 2) )then  
						fast_calculation(ix,iy) = .true.
						if( (landmask(max(ix-1,1),iy).eq. 0) .or. (landmask(min(ix+1,nx),iy).eq. 0) .or. &
							(landmask(ix,max(iy-1,1)).eq. 0) .or. (landmask(ix,min(iy+1,ny)).eq. 0) ) then
							fast_calculation(ix,iy) = .false.
						end if
					end if
					! seasonal snow on ice slow coupling, 
					if ((smb_ice(ix,iy).lt. 0.).and.(landmask(ix,iy).eq. 0) )then  
						fast_calculation(ix,iy) = .false.
					end if



					! collumns with growing mass but no ice yet, slow
					if ((smb_ice(ix,iy).ge. 0.).and.(landmask(ix,iy).eq. 2) )then  
						fast_calculation(ix,iy) = .false.
					end if

					
					
					
					
					
					
					
					
					
					
					
					
! ! ! ! 					if ( (   ((year .ge. daily_snowmask_start).or.((monthly_data).and.(mod(year,monthly_data_freq)==0))).and.(daily_snowmask)   ) & 
! ! ! ! 						.and. (add_accum_output).and.write_output_files_imhof)then
! ! ! ! ! 						print*,'inside saving daily'
! ! ! ! 						m_accum(ix,iy) = m_accum(ix,iy) + smb_ice(ix,iy)*rho_ice * seconds_per_year
! ! ! ! 						m_vaporflux(ix,iy)=real(vaporflux*dt_firn/(L_v+L_lh))+m_vaporflux(ix,iy)
! ! ! ! 						
! ! ! ! 					end if
				
				end if ! for fast coupling
! 			else 
! 			m_snowman_daily(:,:)=0.
			end if  ! end of the if that looks for land

! 			print*,'x=',ix,'y=',iy

            if ( (daily_data) .and.( ((year .ge. daily_data_start).and.(year .le. daily_data_end))&
                .or.(mod(year,daily_data_frequency)==0))  ) then
!                 print*,'prior writing'
                    !3D values in time
                    call writeNCDFSNOW3D_pointwise_t_zxy(filehandle_netcdf_daily, 1, &
                        snowman_daily_varid, m_snowman_daily, iy, ix, n_snowlayer, ndays)
                    call writeNCDFSNOW3D_pointwise_t_zxy(filehandle_netcdf_daily, 1, &
                        lwmass_daily_varid, m_lwmass_daily, iy, ix, n_snowlayer, ndays)
                    call writeNCDFSNOW3D_pointwise_t_zxy(filehandle_netcdf_daily, 1, &
                        rho_snow_daily_varid, m_rho_snow_daily, iy, ix, n_snowlayer, ndays)
                    call writeNCDFSNOW3D_pointwise_t_zxy(filehandle_netcdf_daily, 1, &
                        snow_temp_daily_varid, m_snow_temp_daily, iy, ix, n_snowlayer, ndays)
                    
                    !2D values in time
                    call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, &
                    albedo_dynamic_daily_varid, m_albedo_dynamic_daily, iy, ix, ndays)
                    call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, &
                    latent_heat_daily_varid, m_vaporflux_daily, iy, ix, ndays)
                    call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, &
                    accum_daily_varid, m_accum_daily, iy, ix, ndays)
                    call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, &
                    rain_daily_varid, m_rain_daily, iy, ix, ndays)
                    call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, &
                    refreezing_daily_varid, m_refreezing_daily, iy, ix, ndays)
                    call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, &
                    melt_daily_varid, m_melt_daily, iy, ix, ndays)
                    call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, &
                    snow_daily_varid, m_snow_daily, iy, ix, ndays)
                    call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, &
                    runoff_daily_varid, m_runoff_daily, iy, ix, ndays)
                    call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, &
                    snow_mask_daily_varid, real(m_snow_mask_daily), iy, ix, ndays)
                    call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, &
                    melt_ice_daily_varid, m_melt_ice_daily, iy, ix, ndays)
                    call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, &
                    real_mass_balance_daily_varid, m_real_mass_balance_daily, iy, ix, ndays)
                    call writeNCDFSNOW2D_pointwise_t_xy(filehandle_netcdf_daily, 1, &
                    snow_temp_ave_surface_daily_varid, m_snow_temp_ave_surface_daily, iy, ix, ndays)
                   
                    
                    
            m_snowman_daily=0.
            m_lwmass_daily=0.
            m_rho_snow_daily=0.
            m_snow_temp_daily=0.
            
            m_albedo_dynamic_daily=0.
            m_snow_temp_ave_surface_daily=0.!delete
            m_vaporflux_daily=0.
            m_accum_daily=0.
            m_rain_daily=0.
            m_refreezing_daily=0.
            m_melt_daily=0.
            m_snow_daily=0.
            m_runoff_daily=0.
            m_snow_mask_daily=0.
            m_melt_ice_daily=0.
            m_real_mass_balance_daily=0.
                    
!                 print*,'post writing'
            end if
           !             print*,'end writing' 
            
		end do  ! end of x loop
	end do  ! end of y loop
!     call cpu_time(tm)
!     print*,'end of xy loop', tm-tn
	                                                call cpu_time(t2)
                            print*,'times', t2-t1

! $OMP END DO
! $OMP END PARALLEL

		print*,'end of year',year
	
		if ((minval(snowman)<0.)) then
			print*,'EOY: lightes snow grid cell', minval(snowman),'heaviest snow grid cell', maxval(snowman)
		end if
		if ((minval(lwmass)<0.)) then
			print*,'EOY: negative lw content', minval(lwmass),'wettest snow grid cell', maxval(lwmass)
		end if


		! save monthly data
! ! ! ! 		if (  (((monthly_data).and.(mod(year,monthly_data_freq)==0)) &
! ! ! !             .or. ((year .ge. daily_snowmask_start))).and.write_output_files_imhof ) then
! ! ! ! 			! creat name for output file
! ! ! ! 			!write (netcdf_output_filename_monthly, '( "Snowcover3D_monthly_", I7.7, ".nc" )') year
! ! ! ! 
! ! ! ! 			! initialize output file
! ! ! ! 			!call init_netcdf_3D_file(TRIM(adjustl(output_directory)) // netcdf_output_filename_monthly, filehandle_netcdf_monthly, &
! ! ! ! 			!ny, nx, n_snowlayer, int(dy), int(dx), 1, snowman_varid, lwmass_varid, rho_snow_varid, snow_temp_varid )
! ! ! !             call writeNCDF3DGridValues(filehandle_netcdf_monthly, 1, snow_temp_varid, m_snow_temp(:,:,:,12), ny, nx,&
! ! ! ! 					 n_snowlayer)
! ! ! ! 			! store data into netcdf file
! ! ! ! 			do id=1,monthly_data_num,1
! ! ! ! 				call writeNCDF3DGridValues(filehandle_netcdf_monthly, id, snowman_varid, m_snowman(:,:,:,id), ny, nx, n_snowlayer)
! ! ! ! 				call writeNCDF3DGridValues(filehandle_netcdf_monthly, id, lwmass_varid, m_lwmass(:,:,:,id), ny, nx, n_snowlayer)
! ! ! ! 				call writeNCDF3DGridValues(filehandle_netcdf_monthly, id, rho_snow_varid, m_rho_snow(:,:,:,id), ny, nx, n_snowlayer)
! ! ! ! 				call writeNCDF3DGridValues(filehandle_netcdf_monthly, id, snow_temp_varid, m_snow_temp(:,:,:,id), ny, nx,&
! ! ! ! 					 n_snowlayer)
! ! ! ! !                 call writeNCDF3DGridValuesB(filehandle_netcdf_monthly, id, albedo_dynamic_annual_varid, real(m_accum(:,:)),&
! ! ! ! !                 ny, nx, monthly_data_num)
! ! ! ! 			end do
! ! ! ! 			! close monthly file
! ! ! ! 			call closeNCDFFile(filehandle_netcdf_monthly)
! ! ! ! 		
! ! ! ! 			if(add_accum_output)then
! ! ! ! 
! ! ! ! 				call writeNCDFGridValues(filehandle_netcdf4, 1, accum_varid, m_accum(:,:), ny, nx)
! ! ! ! 				call writeNCDFGridValues(filehandle_netcdf4, 1, rain_varid, m_rain(:,:) , ny, nx)
! ! ! ! 				call writeNCDFGridValues(filehandle_netcdf4, 1, melt_varid, m_melt(:,:), ny, nx)
! ! ! ! 				call writeNCDFGridValues(filehandle_netcdf4, 1, refreezing_varid, m_refreezing(:,:), ny, nx)
! ! ! ! 				call writeNCDFGridValues(filehandle_netcdf4, 1, snow_varid, m_snow(:,:), ny, nx)
! ! ! ! 				call writeNCDFGridValues(filehandle_netcdf4, 1, runoff_varid, m_runoff(:,:), ny, nx)
! ! ! ! ! 				do id=1,ndays,1
! ! ! !                     call writeNCDFGridValues(filehandle_netcdf4, 1, runoff_varid, m_runoff(:,:), ny, nx)
! ! ! ! ! 				end do
! ! ! ! 				! close monthly file with additional output
! ! ! ! 				call closeNCDFFile(filehandle_netcdf4)
! ! ! ! 			end if
! ! ! ! 			month=1
! ! ! ! 		end if

       


! ! ! ! 		! save snow mask data and close netcdf file
! ! ! ! 		if ( ( ((year .ge. daily_snowmask_start).or.((monthly_data).and.(mod(year,monthly_data_freq)==0))).and.&
! ! ! !             (daily_snowmask).and.write_output_files_imhof  )) then
! ! ! !             call cpu_time(t1)
! ! ! ! 			! wirte snow mask to output file
! ! ! ! 			do id=1,ndays,1
! ! ! !                 call cpu_time(t1)
! ! ! ! 				call writeNCDFGridIntegerValues(filehandle_netcdf_snow_mask, id, snow_mask_varid, snow_mask(:,:,id), ny, nx)
! ! ! ! 				call cpu_time(t2)
! ! ! ! 				            print*,"Time Taken 1-->", real(t2-t1)*ndays
! ! ! ! 			end do
! ! ! !             
! ! ! ! 
! ! ! ! 			! close file
! ! ! ! 			call closeNCDFFile(filehandle_netcdf_snow_mask)
! ! ! ! 		end if
        
        do i=1,ndays,1
        m_test_array(:,:,:,i)=real(snowman(:,:,:))
        end do
        
        !New outputs
        !Annual output
        if ( (annual_data) .and.( ((year .ge. daily_data_start).and.(year .le. daily_data_end))&
            .or.(mod(year,annual_data_frequency)==0))  ) then
            !conversion of annual data which are continued to be used on runtime
            m_snowman_annual(:,:,:)=real(snowman(:,:,:))
            m_lwmass_annual=real(lwmass)
            m_snow_temp_annual=real(snow_temp)
            m_rho_snow_annual=real(rho_snow)
            
            !averageing of annual surface data
            m_albedo_dynamic_annual=real(albedo_runtime/ndays)
            m_snow_temp_ave_surface_annual=m_snow_temp_ave_surface_annual/ndays
            
            m_accum_annual = m_accum_annual + smb_ice*rho_ice * seconds_per_year
            

            call writeNCDFSNOW3DValues(filehandle_netcdf_annual, 1, snowman_annual_varid, m_snowman_annual(:,:,:), ny, nx,&
					 n_snowlayer)


            call writeNCDFSNOW3DValues(filehandle_netcdf_annual, 1, lwmass_annual_varid, m_lwmass_annual(:,:,:), ny, nx,&
					 n_snowlayer)

			call writeNCDFSNOW3DValues(filehandle_netcdf_annual, 1, snow_temp_annual_varid, m_snow_temp_annual(:,:,:), ny, nx,&
					 n_snowlayer)		 

            call writeNCDFSNOW3DValues(filehandle_netcdf_annual, 1, rho_snow_annual_varid, m_rho_snow_annual(:,:,:), ny, nx,&
					 n_snowlayer)
					 
            
            !Write 2D values per grid point		
			
            call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, &
                accum_annual_varid, real(m_accum_annual(:,:)), ny, nx,1)
            print*,'here'
            !sum of the latent heat flux  annually
            call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, snow_annual_varid, real(m_snow_annual(:,:)), ny, nx,1)
            
            !sum of accumulation annually
            call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, rain_annual_varid, real(m_rain_annual(:,:)), ny, nx,1)
			
			!sum of liquid precip  annually
            call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, latent_heat_annual_varid, &
                real(m_vaporflux_annual(:,:)), ny, nx,1)
			
			!sum of melt annually
            call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, albedo_dynamic_annual_varid, &
                real(m_albedo_dynamic_annual(:,:)), ny, nx,1)
            print*,'there'
            call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, runoff_annual_varid, real(m_runoff_annual(:,:)), ny, nx,1)
            
            !sum of the refreeze  annually
            call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, melt_annual_varid, real(m_melt_annual(:,:)), ny, nx,1)
			
			!sum of liquid solid annually
            call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, refreezing_annual_varid, &
!             real(p_air(:,:)), ny, nx,1 )
                real(m_refreezing_annual(:,:)), ny, nx,1)
                
            call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, melt_ice_annual_varid, &
                real(m_melt_ice_annual(:,:)), ny, nx,1)

!             call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, melt_ice_annual_varid, &
!                 real(elevation(:,:)), ny, nx,1)
            !int in file 
            call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, snow_mask_annual_varid, real(snow_mask_annual(:,:)), ny, nx,1)
            
            call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, regridding_annual_varid, &
                real(m_regridding_annual(:,:)), ny, nx,1)
            
            call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, &
                snow_temp_ave_surface_annual_varid, real(m_snow_temp_ave_surface_annual(:,:)), ny, nx,1)
                
            call writeNCDFSNOW2DValues(filehandle_netcdf_annual, 1, &
                real_mass_balance_annual_varid, &
                real(m_real_mass_balance_annual(:,:)), ny, nx,1)
        call closeNCDFFile(filehandle_netcdf_annual)
        
            deallocate(m_snowman_annual,m_lwmass_annual,m_rho_snow_annual,m_snow_temp_annual,m_melt_annual)
            deallocate(m_vaporflux_annual,m_rain_annual,m_snow_annual,m_refreezing_annual,m_albedo_dynamic_annual)
            deallocate(m_runoff_annual,snow_mask_annual,m_snow_temp_ave_surface_annual,m_accum_annual)
            deallocate(m_melt_ice_annual, m_real_mass_balance_annual, m_regridding_annual)
        end if
        
        if ( (daily_data) .and.( ((year .ge. daily_data_start).and.(year .le. daily_data_end))&
            .or.(mod(year,daily_data_frequency)==0))  ) then
            !conversion of annual data which are continued to be used on runtime
!      
            deallocate(m_snowman_daily,m_lwmass_daily,m_rho_snow_daily,m_snow_temp_daily,m_melt_daily)
            deallocate(m_vaporflux_daily,m_rain_daily,m_snow_daily,m_refreezing_daily,m_albedo_dynamic_daily)
            deallocate(m_runoff_daily,m_snow_mask_daily,m_snow_temp_ave_surface_daily,m_accum_daily)
            deallocate(m_melt_ice_daily, m_real_mass_balance_daily, m_regridding_daily)
        end if
         
		!----------------------------------------------------------------------
		!! Depleet lowest gridcells
		!----------------------------------------------------------------------
		! grab the lowest gridboxes with a density larger than 840kg/m3 and put them into smb_ice
		!if(spinup) then ! during spinup delete lowest gridcells
		!	snowman(:,:,n_snowlayer) = 0.
		!	lwmass(:,:,n_snowlayer) = 0.
		!	snow_temp(:,:,n_snowlayer) = 0.
		!	rho_snow(:,:,n_snowlayer) = rho_s
		!else
			!call go_depleet_snowman(snowman, rho_snow, snow_temp, lwmass, smb_ice, dummy1) ! calculate the height of ice accumulated or melted in one year
		!end if


		!smb_ice=0


		! Creating a txt file with all mass fluxes
		!-------------------------------------------

		if(mass_checking) then
			!check_tmass_end(ndays) = sum(snowman(:,:,:)) + sum(lwmass(:,:,:))
		  	if (year==1) then
				open(42, file=trim(adjustl(output_directory)) // 'diagnosed_massflux_Snowman3D.txt',&
						 status="new", action="write")
				write(42,*) 'This file contains the total mass of water and snow as well as massfluxes of snow and water '
				write(42,*) 'inside the model Snoman3D. Total masses are in kg and the fluxes are in kg per timestep.'
				write(42,*) 'initial mass snow', tab, 'end mass snow', tab,'initial mass water', tab, &
						 'end mass water', tab, 'accumulation snow', tab, 'accumulation water', tab, 'snow melt', &
							tab, 'water runoff', tab, 'refreeze', tab, 'snow to ice model'
		
		  	else
				open(42, file=trim(adjustl(output_directory)) // 'diagnosed_massflux_Snowman3D.txt',&
						 status="old", position="append", action="write")
		  	end if


			!print*, maxval(check_ice2ice)
			check_init_mass_snow = check_init_mass_snow * dx * dy
			check_end_mass_snow = check_end_mass_snow * dx * dy

			check_init_mass_water = check_init_mass_water * dx * dy
			check_end_mass_water = check_end_mass_water * dx * dy

			check_accum_snow  = check_accum_snow * dx * dy
			check_accum_water = check_accum_water * dx * dy


			check_melted_snow = check_melted_snow * dx * dy
			check_runoff_water = check_runoff_water * dx * dy
			check_refreeze_water = check_refreeze_water * dx * dy
			check_ice2ice = check_ice2ice * dx * dy
			!print*, maxval(check_ice2ice)

			! write fluxes to txt file
			do ii=1,ndays,1
				write(42,*) check_init_mass_snow(ii), tab, check_end_mass_snow(ii), tab, check_init_mass_water(ii), tab, &
					 check_end_mass_water(ii), tab, check_accum_snow(ii), tab, check_accum_water(ii), tab, &
					check_melted_snow(ii), tab, check_runoff_water(ii), tab, check_refreeze_water(ii), tab, check_ice2ice(ii)
				
			end do
			!print*, check_ice2ice
			close(42)
		end if


		! Creating a txt file with all energy fluxes
		!-------------------------------------------

		if(energy_checking) then
		  	if (year==1) then
				open(420, file=trim(adjustl(output_directory)) // 'diagnosed_energyflux_Snowman3D.txt',&
						 status="new", action="write")
				write(420,*) 'This file contains the total energy stored in water and snow as well as energyfluxes into and out '
				write(420,*) 'of the model Snoman3D. Total energies are in J and the fluxes are in J per timestep.'
				write(420,*) 'initial energy snow', tab, 'end energy snow', tab, 'initial energy water', tab, 'end energy water', tab, & 
						'obs. surface energy flux', tab, 'diag. surface energy flux', tab, 'accum tot energy flux', tab, &
						'energy runoff by percolation', tab, 'energy passed to ice', tab, 'freezing energy diag', tab, &
						'freezing energy heat', tab, 'melting energi diag', tab, 'melting energy heat', tab, &
						'Energy entering snow for melt', tab, 'energy runoff by melt-up'
		  	else
				open(420, file=trim(adjustl(output_directory)) // 'diagnosed_energyflux_Snowman3D.txt',&
						 status="old", position="append", action="write")
		  	end if

			check_init_energy_snow  = check_init_energy_snow * dx * dy 
			check_init_energy_water  = check_init_energy_water * dx * dy 
			check_end_energy_snow  = check_end_energy_snow * dx * dy 
			check_end_energy_water  = check_end_energy_water * dx * dy 
			check_surface_e_flux_obs  = check_surface_e_flux_obs * dx * dy 
			check_surface_e_flux_diag = check_surface_e_flux_diag * dx * dy 
			check_e_accum_tot = check_e_accum_tot * dx * dy 

			check_e_freeze_tot  = check_e_freeze_tot * dx * dy 
			check_e_freeze_heat = check_e_freeze_heat * dx * dy 

			check_e_perc_runoff = check_e_perc_runoff * dx * dy 
			check_e_ice2ice = check_e_ice2ice * dx * dy

			check_e_melt_tot = check_e_melt_tot * dx * dy 
			check_e_melt_heat = check_e_melt_heat * dx * dy 
			check_e_melt_qq = check_e_melt_qq * dx * dy 
			check_e_melt_runoff = check_e_melt_runoff * dx * dy 

			! write fluxes to txt file
			do ii=1,ndays,1
				write(420,*) check_init_energy_snow(ii), tab, check_end_energy_snow(ii), tab, check_init_energy_water(ii), tab, &
					 check_end_energy_water(ii), tab,check_surface_e_flux_obs(ii), tab, &
					 check_surface_e_flux_diag(ii), tab, check_e_accum_tot(ii), tab, check_e_perc_runoff(ii), tab, check_e_ice2ice(ii) &
					, tab, check_e_freeze_tot(ii), tab, check_e_freeze_heat(ii), tab, check_e_melt_tot(ii), tab, check_e_melt_heat(ii),&
					 tab, check_e_melt_qq(ii), tab, check_e_melt_runoff(ii)
				
			end do

			close(420)
		end if




	!else
	

	!end if

	!print*,'------------------**Year**--------------- ', year
	!print*,'Initial snowmass= ', check_tmass_init
	!print*,'Accum Snow  = ', sum(check_accums(:))
	!print*,'Accum Rain  = ', sum(check_accumr(:))
	!print*,'Runoff  = ', sum(check_runoff(:))
	!print*,'To IceBern2D= ', dummy1	
	!print*,'Final snowmass by model = ', check_tmass_end
	!print*,'Final estimated snowmas = ', (check_tmass_init +  sum(check_accums(:)) + sum(check_accumr(:)) - &
	!										  sum(check_runoff(:)) - dummy1)

	
    print*,'end of snowbern'

end subroutine get_accumulation_snowman

!############################################## REGRIDDING SUBROUTINE #########################################
subroutine go_regrid(ix,iy, snowman, snow_temp, rho_snow, lwmass, dummy_regrid)

    implicit none
	! input
	integer, intent(in) :: ix
	integer, intent(in) :: iy
	integer, intent(inout) :: dummy_regrid

	real(kind=8), intent(inout) :: snowman(n_snowlayer)   !
	real(kind=8), intent(inout) :: lwmass(n_snowlayer)   !
	real(kind=8), intent(inout) :: snow_temp(n_snowlayer)   !
	real(kind=8), intent(inout) :: rho_snow(n_snowlayer)   !

	!local variables

	real(kind=8) :: masssum
    integer :: ii
    integer :: kk

	! MAIN ROUTINE
    
	!   written and developed by Michael Imhof 25.05.2015
	!   movebox v1.0a

	! this funcion adjusts boxsizes and moves them arround depending on the
	! soll_mass. its only for one collumn of snow

	!   local variables kk,
    
    !counter for calculation the amount of regridings
     dummy_regrid= dummy_regrid+1
    
	! -------------------------------------------------------------------------------
	! spliting and fusing boxes, build function
	!! SPLITING TOP BOX (Accumulating)
	!if((ix==placex).and.(iy==placey)) then
	!print *, 'in spliting routine mass:', snowman
	!print *, 'in spliting routine lwmass:', lwmass
	!print *, 'in spliting routine rho:', rho_snow
	!print *, 'in spliting routine temp:',snow_temp
	!end if
	if(snowman(1) .gt. upper_massbound) then

		! fuse the two lowest boxes
		!if((ix==placex).and.(iy==placey)) then
		!print *, 'got inside splitting loop'
		!end if
		masssum=snowman(n_snowlayer)+snowman(n_snowlayer-1)

		if(masssum==0) then
			rho_snow(n_snowlayer)=rho_s
			snow_temp(n_snowlayer)=0.
		else if(snowman(n_snowlayer)==0) then
			rho_snow(n_snowlayer) = rho_snow(n_snowlayer-1)
			snow_temp(n_snowlayer) = snow_temp(n_snowlayer-1)
		else
			rho_snow(n_snowlayer) = masssum / (snowman(n_snowlayer)/rho_snow(n_snowlayer) +&
			snowman(n_snowlayer-1)/rho_snow(n_snowlayer-1))
			snow_temp(n_snowlayer) = (snowman(n_snowlayer)*snow_temp(n_snowlayer) + &
			snowman(n_snowlayer-1)*snow_temp(n_snowlayer-1))/masssum
		end if

		snowman(n_snowlayer) = masssum
		lwmass(n_snowlayer) =lwmass(n_snowlayer)+lwmass(n_snowlayer-1)

		! push all layers down, except the toplayer
		kk=n_snowlayer-1;
		do while(kk .gt. 2)
			snowman(kk)=snowman(kk-1)
			rho_snow(kk)=rho_snow(kk-1)
			snow_temp(kk)=snow_temp(kk-1)
			lwmass(kk)=lwmass(kk-1)

			kk=kk-1;
		end do
!         print*,'regridding'

		! split top box

		! move proper amount of stuff down to second box
		snowman(2) = soll_mass
		rho_snow(2) = rho_snow(1)
		snow_temp(2) = snow_temp(1)
		lwmass(2) = lwmass(1)*soll_mass/snowman(1)

		! adjust stuff in topbox
		snowman(1)=snowman(1) - soll_mass
		lwmass(1) = lwmass(1) - lwmass(2)
		! density doesn't change
		! temperature doesn't change



		! end of splitting and moving
		!! FUSING TOP BOXES (Melting)
		! splitting lowest box in a meaningfull way!!!
	else if((snowman(1) .lt. lower_massbound).AND.(snowman(2) .gt. 0.)) then
		! fuse both uppermost boxes if there are two

		masssum=snowman(1)+snowman(2)

		rho_snow(1) = masssum / (snowman(1)/rho_snow(1) + snowman(2)/rho_snow(2))
		snow_temp(1) = (snowman(1)*snow_temp(1) + snowman(2)*snow_temp(2))/masssum
		snowman(1) = masssum
		lwmass(1) = lwmass(1)+lwmass(2)

		! drag all layers up, except the toplayer
		do ii=2,n_snowlayer-1,1
			snowman(ii)=snowman(ii+1)
			snow_temp(ii)=snow_temp(ii+1)
			lwmass(ii)=lwmass(ii+1)
			rho_snow(ii)=rho_snow(ii+1)
		end do

		! split second lowest box in 2. the upper should have the
		! soll_mass, the lowest not less thant the soll_mass

		if(snowman(n_snowlayer-1).ge.2*soll_mass) then
			! move down to lowest box
			snowman(n_snowlayer)=snowman(n_snowlayer-1)-soll_mass
			! densitiy doesn't change
			! temperature doesn't change
			lwmass(n_snowlayer)= (1.- soll_mass/snowman(n_snowlayer-1))  *lwmass(n_snowlayer-1)


			! adjust second lowest box
			snowman(n_snowlayer-1)=soll_mass
			! densitiy doent change
			! temperature doesn't change
			lwmass(n_snowlayer-1)=lwmass(n_snowlayer-1)-lwmass(n_snowlayer)

		else
			! reset lowest box
			rho_snow(n_snowlayer)=rho_s
			snowman(n_snowlayer)=0.
			snow_temp(n_snowlayer)=0.
			lwmass(n_snowlayer)=0.
		end if

	end if ! end of splitting or fusing
	!-----------------------------------------------------------------------------
end subroutine go_regrid



!################################ DENSIFICATION SUBROUTINE ####################################################
subroutine go_densification(snowman, rho_snow, temp, At)
    implicit none
	! input
	real(kind=8), intent(in) :: At !
	real(kind=8), intent(inout) :: snowman(n_snowlayer)  !
	real(kind=8), intent(inout) :: temp(n_snowlayer)   !
	real(kind=8), intent(inout) :: rho_snow(n_snowlayer)   !

	!local variables

	real(kind=8) :: ddens
	!real(kind=8) :: dm
	!real(kind=8) :: drdm
	real(kind=8) :: f
	real(kind=8) :: columnsnow
	real(kind=8) :: P_ice
	real(kind=8) :: P_bubble
	real(kind=8) :: dp
    
    integer :: mm, im !local looping variables
	!! densification starts here

	! Michael Imhof, Mai 2015
	! calc temperature of cells

	! consider liquid water mass for B-P

	! local variables: ddens, mm, drdm, dm, f, columnsnow, P_ice, P_bubble,
	! dp


	! temp(isnan(temp))=0;

	do mm=1,n_snowlayer,1

		if(rho_snow(mm) .lt. 550.) then
			! densification  H-L
			if(snowman(mm)>0.) then

				ddens=0.011*exp(-10160.0/8.13/temp(mm))*(rho_i-rho_snow(mm)) *max(At,real(0)) ! is of order 7e-7

				rho_snow(mm)=max(rho_snow(mm),rho_snow(mm)+ddens*dt_firn)
				rho_snow(mm)=min(rho_snow(mm),rho_i)
			end if



		else if(hl==1) then ! H-L for rho>550
			if(mm>1) then
				if(snowman(mm)>0.) then
					ddens=0.575*exp(-21400.0/8.13/temp(mm))*(rho_i-rho_snow(mm))*(1000./3600./24./365.)**(0.5)*&
						(max(real(0),At))**0.5 ! is of order 2e-7
					rho_snow(mm)=max(rho_snow(mm),rho_snow(mm)+ddens*dt_firn)
					rho_snow(mm)=min(rho_snow(mm),rho_i)
				end if
			end if



		else if((rho_snow(mm) .lt. 800.) .and. (rho_snow(mm) .ge. 550.) ) then
			! densification time for rho 550-800 P-B schwander et all

			f=10.**(-29.166*(rho_snow(mm)/rho_i)**3. +84.422*(rho_snow(mm)/rho_i)**2.-87.425* &
			(rho_snow(mm)/rho_i)+30.673)
			columnsnow=0
			if(mm>1) then
				do im=2,mm,1
					columnsnow=columnsnow+snowman(im-1)
				end do
			end if
			P_ice = 9.81*(columnsnow+snowman(mm)/2.)/1e6
			if(rho_snow(mm)>rho_e) then
				P_bubble = P_atm*((1./rho_e-1./rho_i)/(1./rho_snow(mm)-1./rho_i)-1.)/1e6
			else
				P_bubble = 0
			end if
			dp =P_ice-P_bubble
			ddens=25400.*exp(-60000./8.13/temp(mm))*(rho_snow(mm))*f *(dp)**3.
			rho_snow(mm)=max(rho_snow(mm),rho_snow(mm)+ddens*dt_firn)
			rho_snow(mm)=min(rho_snow(mm),rho_i)

		else if(rho_snow(mm) .gt. 800.) then
			! densification time for rho >800 P-B schwander et all
			f=3./16.*(1.-(rho_snow(mm)/rho_i))/(1.-(1.-(rho_snow(mm)/rho_i))**(1./3.))**3.
			columnsnow=0.
			if(mm>1) then
				do im=2,mm,1
					columnsnow=columnsnow+snowman(im-1)
				end do
			end if
			P_ice = 9.81*(columnsnow+snowman(mm)/2.)/1e6
			if(rho_snow(mm)>rho_e) then
				P_bubble = P_atm*((1./rho_e-1./rho_i)/(1./rho_snow(mm)-1./rho_i)-1.)/1e6
			else
				P_bubble = 0.
			end if

			dp=P_ice-P_bubble

			ddens=25400.*exp(-60000./8.13/temp(mm))*(rho_snow(mm))*f *(dp)**3.
			rho_snow(mm)=max(rho_snow(mm),rho_snow(mm)+ddens*dt_firn)
			rho_snow(mm)=min(rho_snow(mm),rho_i)

		end if
	end do

end subroutine go_densification



!K_sw calcuclation line 443ff 502
subroutine go_calculate_solar(ix,iy,snow_temp, T_air, P_sun, albedo_dynamic, nday_snowfall, time, K_sw, albedo_module, lwc)!, albedo_runtime)
                            !ix,iy, snow_temp(ix,iy,1), T_air, P_sun, albedo_dynamic(ix,iy), nday_snowfall(ix,iy), time, K_sw, variables.f90)
    
    !routine calculation the snow albedo for dry and wet snow and the net solar radiation on a point bases
    !   Author: Tobias Zolles
    !   Developer: Tobias Zolles
                            
    !input
    implicit none 
    
    integer, intent(in) :: ix
	integer, intent(in) :: iy
    integer, intent(in) :: albedo_module

	real(kind=8), intent(in) :: T_air
	real(kind=8), intent(in) :: snow_temp
	real(kind=8), intent(in) :: P_sun
	real(kind=8), intent(in) :: lwc
	integer, intent(in) :: nday_snowfall
	integer, intent(in) :: time
	
	real(kind=8), intent(inout) :: albedo_dynamic
	real(kind=8), intent(inout) :: K_sw
	!real(kind=8), intent(inout) :: albedo_runtime

	
	real(16), parameter :: PI_8 = 4 * atan (1.0_8)
		
	real(8)::tempa
	tempa=albedo_dynamic

	    if(albedo_module==1) then !classic calculation
            if(snow_temp<kelvin) then !what does happen to the albedo if the snow became wet, but it got colder again, but no precipitation
                albedo_dynamic=albedo_snow_new
                K_sw=P_sun*(1-albedo_dynamic)
            else
                albedo_dynamic=albedo_snow_wet
                K_sw=P_sun*(1-albedo_dynamic)
    
            end if    
        elseif(albedo_module==2) then !harmonic, for now simple with a sinus and maxium
            if(snow_temp<kelvin) then !what does happen to the albedo if the snow became wet, but it got colder again, but no precipitation
                albedo_dynamic=min(albedo_dynamic,0.9-sin(time/360.*PI_8)*0.25) !avg of 0.74albedo 
!             print*,'albedo', albedo_dynamic
                K_sw=P_sun*(1-albedo_dynamic)
            else
                albedo_dynamic=albedo_snow_wet
                K_sw=P_sun*(1-albedo_dynamic)
    
            end if   
        
        elseif(albedo_module==3) then !temporal decay temperature INdependent Oerlemans and Kap 1998 neglecting the depth scale, which on longer timescales will be ingnored
            if(snow_temp<kelvin) then !what does happen to the albedo if the snow became wet, but it got colder again, but no precipitation
                albedo_dynamic=min(tempa,albedo_snow_wet + (albedo_snow_new-albedo_snow_wet)*&
                exp((nday_snowfall-time)/30.))!albedo_timescale))
                K_sw=P_sun*(1-albedo_dynamic)
            else
                albedo_dynamic=min(tempa,albedo_snow_wet + (albedo_snow_new-albedo_snow_wet)*&
                exp((nday_snowfall-time)/5.))
                K_sw=P_sun*(1-albedo_dynamic)
    
            end if   
        elseif(albedo_module==5) then !temporal decay Bougamont (Oerlemans and Kap 1998) neglecting the depth scale, LWC based reduced albedo
        ! uses a linear increase in reduction rate depending on the lwc of the snowlayer, exponential possible for 
            if(snow_temp<kelvin-10) then !what does happen to the albedo if the snow became wet, but it got colder again, but no precipitation
                albedo_dynamic=min(tempa,albedo_snow_wet + (albedo_snow_new-albedo_snow_wet)*&
                exp((nday_snowfall-time)/100.))!albedo_timescale))
                K_sw=P_sun*(1-albedo_dynamic)
            elseif(snow_temp<kelvin) then
                albedo_dynamic=min(tempa,albedo_snow_wet + (albedo_snow_new-albedo_snow_wet)*&
                exp((nday_snowfall-time)/(30.+7*abs(snow_temp-kelvin)))) !albedo_timescale
                K_sw=P_sun*(1-albedo_dynamic)
            else
                albedo_dynamic=min(tempa,albedo_snow_wet + (albedo_snow_new-albedo_snow_wet)*&
                exp((nday_snowfall-time)/(15.-14*(lwc/max_lwc)))) !option use even lower value for wet snow with max lwc, problem for lwc/max>1
                K_sw=P_sun*(1-albedo_dynamic)
            end if       
        elseif(albedo_module==4) then !calculate dry snow albedo reduction based on Aoki 2003 suggest weaker effect due to fewer impurities in Antarctica and remote sea ice
                albedo_dynamic=min(tempa,tempa-((snow_temp-kelvin)*1.35e-3+0.0278))
                if(albedo_dynamic<albedo_snow_wet) then
                    albedo_dynamic=albedo_snow_wet
                end if
                if(lwc>0) then
                tempa=albedo_dynamic
                albedo_dynamic=max(albedo_snow_wet,min(tempa,tempa-(tempa-albedo_snow_wet)*(lwc/max_lwc)))
                K_sw=P_sun*(1-albedo_dynamic)
                end if
!               if ((ix==39) .and. (iy==74)) then
!                 print*,'internal', tempalbedo
                K_sw=P_sun*(1-albedo_dynamic)
!               print*,albedo_dynamic
        elseif(albedo_module==0) then
            K_sw=P_sun*(1-albedo_dynamic)
        
        
        end if
   !tempa=albedo_runtime(1)
   !albedo_runtime(1)=tempa+albedo_dynamic


end subroutine go_calculate_solar


!############################## ENERGY FLUX SUBROUTINE ###################################################################
subroutine go_energy_flux_tobias( ix, iy, snow_temp, rho_snow, china_syndrome, T_air, mbox, ddz, &
                                K_sw, H_lh, K_lh, Q_heat, heating, DewpT, vaporflux, D_lf, p_air)

    implicit none
	! input
	integer, intent(in) :: ix
	integer, intent(in) :: iy

	real(kind=8), intent(in) :: T_air

	real(kind=8), intent(in) :: DewpT
	real(kind=8), intent(in) :: mbox
	real(kind=8), intent(in) :: K_sw
	real(kind=8), intent(in) :: H_lh
	real(kind=8), intent(in) :: K_lh
    real(kind=8), intent(in) :: D_lf
    real(kind=8), intent(in) :: p_air
! 	real(kind=8), intent(in) :: albedo_dynamic

	real(kind=8), intent(in) :: ddz(n_snowlayer)
	real(kind=8), intent(in) :: rho_snow(n_snowlayer)   !

	real(kind=8), intent(inout) :: Q_heat
	real(kind=8), intent(inout) :: snow_temp(n_snowlayer)   !
	real(kind=8), intent(inout) :: heating

	real(kind=8), intent(inout) :: vaporflux

	logical, intent(inout) :: china_syndrome   !

	!local variables
	integer :: nn = 0 ! amount of active grid cells in vertical direction
	integer :: inn = 1
	integer :: ii = 0
	integer, dimension(n_snowlayer) ::  istheresnow

	logical :: no_ground = .true.
	real(kind=8) :: H=0.
	real(kind=8) :: dummy
	real(kind=8) :: BB_mid1_backup

	real(kind=8) :: ea
	real(kind=8) :: es
	!write(*,*) 'before allocatable definition'

	! local allocatable variables
	real(kind=8), dimension(:), allocatable :: tempi(:)
	real(kind=8), dimension(:), allocatable :: new_temp(:)
	real(kind=8), dimension(:), allocatable :: backup(:)
	real(kind=8), dimension(:), allocatable :: K(:)
	real(kind=8), dimension(:), allocatable :: K_snow(:)
	real(kind=8), dimension(:), allocatable :: dz1(:)
	real(kind=8), dimension(:), allocatable :: BB_up(:)
	real(kind=8), dimension(:), allocatable :: BB_mid(:)
	real(kind=8), dimension(:), allocatable :: BB_down(:)

	!print *,'after allocatable definition'
	! extract active gridcells and add 2 dummy gridcells at surface and bottom.
	
	heating=0.

	istheresnow=0
	where (snow_temp(:)>0.)
		istheresnow(:)=1
	end where
	nn=int(sum(istheresnow))
	

	if(inn==0)then
		print*,'ERROR inn=',inn,'masse=',mbox,'Temperature=',snow_temp
		read(*,*) dummy
	end if

	!nn=in
	allocate (tempi(nn))
	allocate (new_temp(nn))
	allocate (backup(nn))
	allocate (K(nn))
	allocate (K_snow(nn))
	allocate (dz1(nn))
	allocate (BB_up(nn))
	allocate (BB_mid(nn))
	allocate (BB_down(nn))
	!allocate (BB(nn,nn))

	!print *,'after allocating'

	china_syndrome=.false. ! in case of melting must occure this will become true

	Q_heat=0.
	tempi=0.
	new_temp=0.
	backup=0 ! backup of temperature distribution in case melting occures
	K=0.
	K_snow=0.
	dz1=0.
	BB_up=0.
	BB_mid=0.
	BB_down=0.
	BB_mid1_backup=0.

	! Tobias Zolles, adapted from Michael Imhof, January 2018

	! this subroutine calcualtes the energy fluxes in the snowcap ie the
	! new temperatures. an implicit sceme is used for that.

	! Initialize grid with variables

	! initialize resized dz and temperature

	tempi(1:nn)= snow_temp(1:nn)
	dz1(1:nn) = ddz(1:nn)

	! store a backup of the initial temperature
	backup(1:nn)=tempi(1:nn)


    !Latent heat flux: H=0.622rho_aL_v*C_h_u(ea-es)P???1,Oerlemans Rolstad   C_h ~2e-3
    !rho_a air densitiy, P air pressure, L_v 2.5e6
    !es= 611.20Pa exp(L_s(2834e6)/R(461.5)*(1/kelvin(273.16)-1/T_s))
    !es=6.112exp(22.46*T/(272.62+T)) T in C, es hPa Guide to Meterological Instruments and Methods of Obersvation (WMO, 2008)
    !RH=ea/ews*100
    !ea=0.6108*exp(17.27*DewpT/(DewpT+237.3))
    !sensible: rho_a cp Ch u =A Ch=A/rho_a/cp
    
    !artificial DewpT until forcing data is sorted out
!     DewpT=T_air-kelvin-5
!     if(ix==90.and.iy==168) then
!         print*,DewpT
!     end if
    
    ea=610.8*exp(17.27*DewpT/(DewpT+237.3))
!     ea=611.21*exp(17.502*(DewpT-kelvin)/(DewpT-32.19)
!     print*,Dewpt
!     print*,snow_temp(1)-273
    es=611.2*exp(22.46*(snow_temp(1)-kelvin)/(272.62+(snow_temp(1)-kelvin))) !hpa=6.112
!     print*,ea-es
    !absch??tzungen:
    !es=6.122*22.46*exp(22.46*(tempi(1)-kelvin)/(tempi(1)-kelvin+272.62)+1)/(tempi(1)-kelvin+272.62)/(tempi(1)-kelvin+272.62)
!     K_lf=0.622*rho_a*2.5e6*C_h*u

!     if(new_temp(1).lt.kelvin) then !if at melting point it is considered vaporisation   D_sf/cp_air*0.622*L_v/p_air (elevation 
    !dependent pressure coordinate)
!     D_lf=D_sf/1003*0.622*(L_v+L_lh)/10000
!     else 
!     D_lf=D_sf/1003*0.622*(L_v)/10000
!     end if

!     if((ix==90).and.(iy==180)) then
!       print*,'Dewpt',DewpT
!       print*,'latent',D_lf*(ea-es)
!       print*,'sensible',D_sf*(T_air-snow_temp(1))
!       print*,'solar',K_sw
!       print*,'snowtemp',snow_temp(1)
!       print*,'airtemp',T_air
!       print*,'long',sigma*(eps_air*(T_air)**4.-eps_snow*(tempi(1))**4.)
!       print*,ea
!       print*,es
!       print*,D_lf
!     end if
    
	! set up K vector (Snow Temperature INdependent parts)
	K(1)=dt_firn/c_i/mbox*((T_air)*D_sf+sigma*(eps_air*(T_air)**4.+eps_snow*3.*(tempi(1))**4.)+	&
	K_sw + K_lh+D_lf/p_air*(ea-es*1+es*22.46*272.62*snow_temp(1)/(272.62+(snow_temp(1)-kelvin))**2))

	! K(nn) = dt_firn*Q_geo/c_i/mbox 	! This is aobut how the geothermal heatflux would look like. 
						! A corresponding melting routine would be needed.

	! set up H in BB_mid (Snow Temperature DEpendent parts)
	H = dt_firn/c_i/mbox*(D_sf +sigma*eps_snow*4.*(tempi(1))**3. +H_lh&
	+D_lf/p_air*es*22.46*272.62/(272.62+(snow_temp(1)-kelvin))**2)
!     print*,H
!     print*,snow_temp(1)
!     print*,K(1)

!     print*,

	! set up tridiagonal system

	! there are two cases, 1 and 2+ boxes

	if(nn==1)then
		new_temp(1) = (tempi(1)+K(1))/(1+H)
        
		! if the surface temperature has become larger than 0C, do the
		! calculation again without the surface energyfluxes, but set the
		! surfacetemperature to 273 at the beginning and at the end. this is
		! done this way to avoid unphysical temperature fluxes into the depth of
		! the snowcover. the energy used to heat the snow to 0 degrees is subtractet later in the melting
		if(new_temp(1) .gt. kelvin)  then

			china_syndrome = .true.
			! set temperature in topbox to 0C and ignor incomming energy fluxes. keep topbox temperature at 0C and calculate energy necessary to rise temperature to 0C
			Q_heat= (kelvin-backup(1))*c_i*mbox
			new_temp(1)=kelvin
			heating = Q_heat

		else
			! clean the temperature profile. where snow is at melting point over several layers due to water
			! it can happen that the temperature of one or more boxes can become 273.00000000000006 K due to
			! computational limits.
		!	where(new_temp(:) .gt. kelvin)
		!		new_temp(:) = kelvin
		!	end where

			heating = dt_firn * (new_temp(1)*(-D_sf -sigma*eps_snow*4.*(tempi(1))**3. -H_lh ) + &
					((T_air)*D_sf+sigma*(eps_air*(T_air)**4.+eps_snow*3.*(tempi(1))**4.) + K_sw + K_lh) )

		end if

! 		 vaporflux = D_lf*(6.108*exp(17.27*DewpT/(DewpT+237.3))-&
!         6.112*exp(22.46*(snow_temp(1)-273)/(272.62+(snow_temp(1)-273))))
		! update temperature
		snow_temp(1) = new_temp(1)
		 vaporflux = D_lf/p_air*(6.108*exp(17.27*DewpT/(DewpT+237.3))-&
        6.112*exp(22.46*(snow_temp(1)-kelvin)/(272.62+(snow_temp(1)-kelvin))))
		

	else


		! set up model for the temperature diffusion parameter
		if(diff_model==1) then
			K_snow(1:nn) = K_ice*(rho_snow(1:nn)/1000.)**(1.88) ! yen1981
		elseif(diff_model==2) then
			do ii=1,nn,1
				if (rho_snow(ii).gt. 156.) then
					K_snow(ii) = (0.138-1.01e-3*rho_snow(ii)+ 3.233e-6*rho_snow(ii)**2.) ! Sturm 1997
				else
					K_snow(ii) = 0.023+ 0.234e-3*rho_snow(ii)
				end if
			end do	
		else
			K_snow(1:nn) = (2.1e-2+4.2e-4*rho_snow(1:nn) + 2.2e-9*rho_snow(1:nn)**3.) ! Van Dusen (1929)
		end if


		! set up tridiagonal system

		BB_up(1)   = -2.*dt_firn/rho_snow(1)/c_i/dz1(1)*( K_snow(1)*dz1(1)+K_snow(1+1)*dz1(1+1)) /(dz1(1)+dz1(1+1))**2.
		BB_up(nn)  = -999.

		BB_down(1)  = -999.
		BB_down(nn) = -2.*dt_firn/rho_snow(nn)/c_i/dz1(nn)*( K_snow(nn)*dz1(nn)+K_snow(nn-1)*dz1(nn-1)) /(dz1(nn)+dz1(nn-1))**2.

		BB_mid(1)  = 1. - BB_up(1)
		BB_mid(nn) = 1. - BB_down(nn)


		if(nn .gt. 2)then
			do ii=2,nn-1,1
				BB_down(ii) = -2.*dt_firn/rho_snow(ii)/c_i/dz1(ii)*( K_snow(ii)*dz1(ii)+K_snow(ii-1)*dz1(ii-1)) /(dz1(ii)+dz1(ii-1))**2.
				BB_up(ii)   = -2.*dt_firn/rho_snow(ii)/c_i/dz1(ii)*( K_snow(ii)*dz1(ii)+K_snow(ii+1)*dz1(ii+1)) /(dz1(ii)+dz1(ii+1))**2.
				BB_mid(ii)  = 1. - BB_down(ii) - BB_up(ii)
			end do
		end if



		BB_mid1_backup = BB_mid(1)

		!if ((ix==80).and.(iy==80)) then
		!	print*,BB_down
		!	print*,BB_mid
		!	print*,BB_up
		!end if


		BB_mid(1) = BB_mid(1)+H ! Snow Temperature dependent parts (only outgoing sensible heat and longwave rad and latent heatfux)


		! Tridag solver
		call tridag(BB_down,BB_mid,BB_up,(tempi+K),new_temp,nn)


		! if the surface temperature has become larger than 0C, do the
		! calculation again without the surface energyfluxes, but set the
		! surfacetemperature to 273 at the beginning and at the end. this is
		! done this way to avoid unphysical temperature fluxes into the deep of
		! the snowcover. the energy used to heat the snow to 0 degrees is subtractet later in the melting
		if(new_temp(1) .gt. kelvin)  then

			china_syndrome = .true.
			! set temperature in topbox to 0C and ignor incomming energy fluxes. keep topbox temperature at 0C and calculate energy necessary to rise temperature to 0C
			Q_heat= (kelvin-backup(1))*c_i*mbox
			backup(1)=kelvin

			BB_mid(1) = BB_mid1_backup

		!	if ((ix==97).and.(iy==168)) then
		!		print*,'TRID: backup temp', backup
		!	end if

			! Tridag solver
			call tridag(BB_down(1:nn),BB_mid(1:nn),BB_up(1:nn),backup(1:nn),new_temp(1:nn),nn)

		!	if ((ix==97).and.(iy==168)) then
		!		print*,'TRID: backup temp', new_temp
		!	end if
			! set temperature in topbox to 0C. The energy used for this will be subtracted in the melting routine
			Q_heat = Q_heat + (kelvin-new_temp(1))*c_i*mbox
			new_temp(1)=kelvin


			! due to computational limits it is possible that new_temp(x)= 273.00000000000006. this happens especially when the entire collumn is
			! at the freezing point. in order to avoid that we set temperatures larger than 273.000 K to 273. K

			where(new_temp(:) .gt. kelvin)
				new_temp(:) = kelvin
			end where
			heating = Q_heat

			!new_temp(1)=kelvin

			!if(new_temp(2) .gt. kelvin) then
			!	print*,'new temp 2', new_temp(2)
			!end if
			!QQ=(QQ+(kelvin-tempi(2))*c_i*mbox)/dt_firn ! energy flux due to warming which cant be used for melting
			! subtract energy that has been used for keeping the snow at 0C, note QQ is negative

		else
			! clean the temperature profile. where snow is at melting point over several layers due to water
			! it can happen that the temperature of one or more boxes can become 273.00000000000006 K due to
			! computational limits.
			where(new_temp(:) .gt. kelvin)
				new_temp(:) = kelvin
			end where

			heating = dt_firn * (new_temp(1)*(-D_sf -sigma*eps_snow*4.*(tempi(1))**3. -H_lh ) + &
					((T_air)*D_sf+sigma*(eps_air*(T_air)**4.+eps_snow*3.*(tempi(1))**4.) + K_sw + K_lh) )



		end if
! 		tempi(1)=new_temp(1)-snow_temp(1)
!         print*,'this is',tempi(1)
        
!          vaporflux = D_lf*(610.8*exp(17.27*DewpT/(DewpT+237.3))-&
!         611.2*exp(22.46*(snow_temp(1)-273)/(272.62+(snow_temp(1)-273))))
!         print*,D_lf
!         print*,'es',es
!         print*,'ea',ea
!         
!         print*,'snow',vaporflux*dt_firn
!         print*,'snowtemp2',snow_temp(1)
! 		! update temperature
		snow_temp(1:nn) = new_temp(1:nn)
		  vaporflux = D_lf/p_air*(610.8*exp(17.27*DewpT/(DewpT+237.3))-&
        611.2*exp(22.46*(snow_temp(1)-273)/(272.62+(snow_temp(1)-273))))
! 		print*,vaporflux
	end if

	deallocate (tempi)
	deallocate (new_temp)
	deallocate (backup)
	deallocate (K)
	deallocate (K_snow)
	deallocate (dz1)
	deallocate (BB_up)
	deallocate (BB_mid)
	deallocate (BB_down)


end subroutine go_energy_flux_tobias

!############################## ENERGY FLUX SUBROUTINE ###################################################################
subroutine go_energy_flux_new( ix, iy, snow_temp, rho_snow, china_syndrome, T_air, mbox, ddz, P_sun, H_lh, K_lh, Q_heat, heating)

    implicit none
	! input
	integer, intent(in) :: ix
	integer, intent(in) :: iy

	real(kind=8), intent(in) :: T_air
	real(kind=8), intent(in) :: mbox
	real(kind=8), intent(in) :: P_sun
	real(kind=8), intent(in) :: H_lh
	real(kind=8), intent(in) :: K_lh

	real(kind=8), intent(in) :: ddz(n_snowlayer)
	real(kind=8), intent(in) :: rho_snow(n_snowlayer)   !

	real(kind=8), intent(inout) :: Q_heat
	real(kind=8), intent(inout) :: snow_temp(n_snowlayer)   !
	real(kind=8), intent(inout) :: heating


	logical, intent(inout) :: china_syndrome   !

	!local variables
	integer :: nn = 0 ! amount of active grid cells in vertical direction
	integer :: inn = 1
	integer :: ii = 0
	integer, dimension(n_snowlayer) ::  istheresnow

	logical :: no_ground = .true.
	real(kind=8) :: H=0.
	real(kind=8) :: dummy
	real(kind=8) :: BB_mid1_backup
	!write(*,*) 'before allocatable definition'

	! local allocatable variables
	real(kind=8), dimension(:), allocatable :: tempi(:)
	real(kind=8), dimension(:), allocatable :: new_temp(:)
	real(kind=8), dimension(:), allocatable :: backup(:)
	real(kind=8), dimension(:), allocatable :: K(:)
	real(kind=8), dimension(:), allocatable :: K_snow(:)
	real(kind=8), dimension(:), allocatable :: dz1(:)
	real(kind=8), dimension(:), allocatable :: BB_up(:)
	real(kind=8), dimension(:), allocatable :: BB_mid(:)
	real(kind=8), dimension(:), allocatable :: BB_down(:)

	!print *,'after allocatable definition'
	! extract active gridcells and add 2 dummy gridcells at surface and bottom.
	
	heating=0.

	istheresnow=0
	where (snow_temp(:)>0.)
		istheresnow(:)=1
	end where
	nn=int(sum(istheresnow))
	

	if(inn==0)then
		print*,'ERROR inn=',inn,'masse=',mbox,'Temperature=',snow_temp
		read(*,*) dummy
	end if

	!nn=in
	allocate (tempi(nn))
	allocate (new_temp(nn))
	allocate (backup(nn))
	allocate (K(nn))
	allocate (K_snow(nn))
	allocate (dz1(nn))
	allocate (BB_up(nn))
	allocate (BB_mid(nn))
	allocate (BB_down(nn))
	!allocate (BB(nn,nn))

	!print *,'after allocating'

	china_syndrome=.false. ! in case of melting must occure this will become true

	Q_heat=0.
	tempi=0.
	new_temp=0.
	backup=0 ! backup of temperature distribution in case melting occures
	K=0.
	K_snow=0.
	dz1=0.
	BB_up=0.
	BB_mid=0.
	BB_down=0.
	BB_mid1_backup=0.

	! Michael Imhof, January 2017

	! this subroutine calcualtes the energy fluxes in the snowcap ie the
	! new temperatures. an implicit sceme is used for that.

	! Initialize grid with variables

	! initialize resized dz and temperature

	tempi(1:nn)= snow_temp(1:nn)
	dz1(1:nn) = ddz(1:nn)

	! store a backup of the initial temperature
	backup(1:nn)=tempi(1:nn)





	! set up K vector (Snow Temperature INdependent parts)
	K(1)=dt_firn/c_i/mbox*((T_air)*D_sf+sigma*(eps_air*(T_air)**4.+eps_snow*3.*(tempi(1))**4.) + P_sun + K_lh)

	! K(nn) = dt_firn*Q_geo/c_i/mbox 	! This is aobut how the geothermal heatflux would look like. 
						! A corresponding melting routine would be needed.

	! set up H in BB_mid (Snow Temperature DEpendent parts)
	H = dt_firn/c_i/mbox*(D_sf +sigma*eps_snow*4.*(tempi(1))**3. +H_lh )




	! set up tridiagonal system

	! there are two cases, 1 and 2+ boxes

	if(nn==1)then
		new_temp(1) = (tempi(1)+K(1))/(1+H)

		! if the surface temperature has become larger than 0C, do the
		! calculation again without the surface energyfluxes, but set the
		! surfacetemperature to 273 at the beginning and at the end. this is
		! done this way to avoid unphysical temperature fluxes into the depth of
		! the snowcover. the energy used to heat the snow to 0 degrees is subtractet later in the melting
		if(new_temp(1) .gt. kelvin)  then

			china_syndrome = .true.
			! set temperature in topbox to 0C and ignor incomming energy fluxes. keep topbox temperature at 0C and calculate energy necessary to rise temperature to 0C
			Q_heat= (kelvin-backup(1))*c_i*mbox
			new_temp(1)=kelvin
			heating = Q_heat

		else
			! clean the temperature profile. where snow is at melting point over several layers due to water
			! it can happen that the temperature of one or more boxes can become 273.00000000000006 K due to
			! computational limits.
		!	where(new_temp(:) .gt. kelvin)
		!		new_temp(:) = kelvin
		!	end where

			heating = dt_firn * (new_temp(1)*(-D_sf -sigma*eps_snow*4.*(tempi(1))**3. -H_lh ) + &
					((T_air)*D_sf+sigma*(eps_air*(T_air)**4.+eps_snow*3.*(tempi(1))**4.) + P_sun + K_lh) )

		end if

		! update temperature
		snow_temp(1) = new_temp(1)

	else


		! set up model for the temperature diffusion parameter
		if(diff_model==1) then
			K_snow(1:nn) = K_ice*(rho_snow(1:nn)/1000.)**(1.88) ! yen1981
		elseif(diff_model==2) then
			do ii=1,nn,1
				if (rho_snow(ii).gt. 156.) then
					K_snow(ii) = (0.138-1.01e-3*rho_snow(ii)+ 3.233e-6*rho_snow(ii)**2.) ! Sturm 1997
				else
					K_snow(ii) = 0.023+ 0.234e-3*rho_snow(ii)
				end if
			end do	
		else
			K_snow(1:nn) = (2.1e-2+4.2e-4*rho_snow(1:nn) + 2.2e-9*rho_snow(1:nn)**3.) ! Van Dusen (1929)
		end if


		! set up tridiagonal system

		BB_up(1)   = -2.*dt_firn/rho_snow(1)/c_i/dz1(1)*( K_snow(1)*dz1(1)+K_snow(1+1)*dz1(1+1)) /(dz1(1)+dz1(1+1))**2.
		BB_up(nn)  = -999.

		BB_down(1)  = -999.
		BB_down(nn) = -2.*dt_firn/rho_snow(nn)/c_i/dz1(nn)*( K_snow(nn)*dz1(nn)+K_snow(nn-1)*dz1(nn-1)) /(dz1(nn)+dz1(nn-1))**2.

		BB_mid(1)  = 1. - BB_up(1)
		BB_mid(nn) = 1. - BB_down(nn)


		if(nn .gt. 2)then
			do ii=2,nn-1,1
				BB_down(ii) = -2.*dt_firn/rho_snow(ii)/c_i/dz1(ii)*( K_snow(ii)*dz1(ii)+K_snow(ii-1)*dz1(ii-1)) /(dz1(ii)+dz1(ii-1))**2.
				BB_up(ii)   = -2.*dt_firn/rho_snow(ii)/c_i/dz1(ii)*( K_snow(ii)*dz1(ii)+K_snow(ii+1)*dz1(ii+1)) /(dz1(ii)+dz1(ii+1))**2.
				BB_mid(ii)  = 1. - BB_down(ii) - BB_up(ii)
			end do
		end if



		BB_mid1_backup = BB_mid(1)

		!if ((ix==80).and.(iy==80)) then
		!	print*,BB_down
		!	print*,BB_mid
		!	print*,BB_up
		!end if


		BB_mid(1) = BB_mid(1)+H ! Snow Temperature dependent parts (only outgoing sensible heat and longwave rad and latent heatfux)


		! Tridag solver
		call tridag(BB_down,BB_mid,BB_up,(tempi+K),new_temp,nn)


		! if the surface temperature has become larger than 0C, do the
		! calculation again without the surface energyfluxes, but set the
		! surfacetemperature to 273 at the beginning and at the end. this is
		! done this way to avoid unphysical temperature fluxes into the deep of
		! the snowcover. the energy used to heat the snow to 0 degrees is subtractet later in the melting
		if(new_temp(1) .gt. kelvin)  then

			china_syndrome = .true.
			! set temperature in topbox to 0C and ignor incomming energy fluxes. keep topbox temperature at 0C and calculate energy necessary to rise temperature to 0C
			Q_heat= (kelvin-backup(1))*c_i*mbox
			backup(1)=kelvin

			BB_mid(1) = BB_mid1_backup

		!	if ((ix==97).and.(iy==168)) then
		!		print*,'TRID: backup temp', backup
		!	end if

			! Tridag solver
			call tridag(BB_down(1:nn),BB_mid(1:nn),BB_up(1:nn),backup(1:nn),new_temp(1:nn),nn)

		!	if ((ix==97).and.(iy==168)) then
		!		print*,'TRID: backup temp', new_temp
		!	end if
			! set temperature in topbox to 0C. The energy used for this will be subtracted in the melting routine
			Q_heat = Q_heat + (kelvin-new_temp(1))*c_i*mbox
			new_temp(1)=kelvin


			! due to computational limits it is possible that new_temp(x)= 273.00000000000006. this happens especially when the entire collumn is
			! at the freezing point. in order to avoid that we set temperatures larger than 273.000 K to 273. K

			where(new_temp(:) .gt. kelvin)
				new_temp(:) = kelvin
			end where
			heating = Q_heat

			!new_temp(1)=kelvin

			!if(new_temp(2) .gt. kelvin) then
			!	print*,'new temp 2', new_temp(2)
			!end if
			!QQ=(QQ+(kelvin-tempi(2))*c_i*mbox)/dt_firn ! energy flux due to warming which cant be used for melting
			! subtract energy that has been used for keeping the snow at 0C, note QQ is negative

		else
			! clean the temperature profile. where snow is at melting point over several layers due to water
			! it can happen that the temperature of one or more boxes can become 273.00000000000006 K due to
			! computational limits.
			where(new_temp(:) .gt. kelvin)
				new_temp(:) = kelvin
			end where

			heating = dt_firn * (new_temp(1)*(-D_sf -sigma*eps_snow*4.*(tempi(1))**3. -H_lh ) + &
					((T_air)*D_sf+sigma*(eps_air*(T_air)**4.+eps_snow*3.*(tempi(1))**4.) + P_sun + K_lh) )



		end if

		! update temperature
		snow_temp(1:nn) = new_temp(1:nn)

	end if

	deallocate (tempi)
	deallocate (new_temp)
	deallocate (backup)
	deallocate (K)
	deallocate (K_snow)
	deallocate (dz1)
	deallocate (BB_up)
	deallocate (BB_mid)
	deallocate (BB_down)


end subroutine go_energy_flux_new


!############################## MELTING SNOW SUBROUTINE ###################################################################
subroutine go_melting_snow(ix,iy,time,snow_temp, snowman, lwmass, rho_snow, T_air, K_sw, K_lh, &
									H_lh, Q_heat, melted_snow, runoff_water, ice_melt, used_q, l_heat, vaporflux, D_lf,&
									dummy_melt_ice, dummy_regrid)
    implicit none
	integer, intent(in) :: time
	integer, intent(in) :: iy
	integer, intent(in) :: ix
	real(kind=8), intent(in) :: H_lh
	real(kind=8), intent(in) :: T_air
	!real(kind=8), intent(in) :: DewpT
	real(kind=8), intent(in) :: K_sw
	real(kind=8), intent(in) :: K_lh
	real(kind=8), intent(in) :: D_lf

	real(kind=8), intent(inout) :: snowman(n_snowlayer)
	real(kind=8), intent(inout) :: lwmass(n_snowlayer)
	real(kind=8), intent(inout) :: snow_temp(n_snowlayer)
	real(kind=8), intent(inout) :: rho_snow(n_snowlayer)
	real(kind=8), intent(inout) :: Q_heat
	real(kind=8), intent(inout) :: melted_snow
	real(kind=8), intent(inout) :: runoff_water
	real(kind=8), intent(inout) :: ice_melt
	real(kind=8), intent(inout) :: used_q
	real(kind=8), intent(inout) :: l_heat
	real(kind=8), intent(inout) :: vaporflux
	real(kind=8), intent(inout) :: dummy_melt_ice
	integer, intent(inout) :: dummy_regrid

	!logical, intent(inout) :: china_syndrome

	! local variables

	real(kind=8) :: QQ
	real(kind=8) :: Qp_lw
	real(kind=8) :: Qp_sh
	real(kind=8) :: Qp_lh
	real(kind=8) :: Q_v
	real(kind=8) :: dT
	real(kind=8) :: dm
	real(kind=8) :: DewpT

	melted_snow = 0.
	runoff_water = 0.
	QQ = 0.
	l_heat = 0.
	used_q = 0.

	!ice_melt = 0
	! melt the snow depending on energy upptaken by the snowcover if it
	! reeaches 273K.

	! if ((iy==lat)&&(ix==long)) %%&&(snow_temp(ix,iy,1)<250)
	! fprintf('Pre  top temp %3.3f second temp: %3.3f (year %d, step %d )\n', &
	  !snow_temp(ix,iy,1),snow_temp(ix,iy,2),yy,time )
	! end

	!
	!
	!artificial DewpT
! 	DewpT=temp_of_air-kelvin-10
!   !possibility to calculate vaporflux here, but is done at end of the energy flux calculation instead
!     vaporflux = D_lf*(6.108*exp(17.27*DewpT/(DewpT+237.3))-&
!         6.112*exp(22.46*(snow_temp(1)-273)/(272.62+(snow_temp(1)-273))))
!         print*,vaporflux
! 	      print*,'Dewpt',DewpT
!       print*,'latent',D_lf*(ea-es)
!       print*,'sensible',D_sf*(T_air-snow_temp(1))
!       print*,'solar',K_sw
!       print*,'snowtemp',snow_temp(1)
!       print*,ea
!       print*,es
!       print*,D_lf
! 	
	!if(snow_temp(1)>=kelvin) then
		!snow_temp(1)=kelvin ! set temperature to 0

	Qp_lw = sigma*(eps_air*(T_air)**4.-eps_snow*(kelvin)**4.)
	Qp_sh = D_sf*(T_air-kelvin)
	Qp_lh = K_lh-H_lh*snow_temp(1)
	Q_v =vaporflux
!     Q_v=0.
	
	QQ = max((K_sw+Qp_lw+Qp_sh+Qp_lh+Q_v)*dt_firn-Q_heat,real(0.))
	used_q = QQ
	!print*,'QQ=',QQ,'Q_heat=',Q_heat
	Q_heat = 0.
	do while (QQ .gt. 0.)
		dT=QQ/c_i/snowman(1) ! maximum possible heating
		if(kelvin-snow_temp(1) .gt. dT) then ! QQ too small to melt something
			snow_temp(1) = snow_temp(1)+dT
			QQ = 0.

		else if(snowman(1) .gt. c_i*snowman(1)*(dT-(kelvin-snow_temp(1)))/L_lh) then ! snow is heated to 0C and partially melted
			dm = QQ/L_lh-c_i*snowman(1)*(kelvin-snow_temp(1))/L_lh
			l_heat = l_heat + dm*L_lh
			snow_temp(1) = kelvin
			snowman(1) = snowman(1)-dm
			lwmass(1)  = lwmass(1) + dm
			melted_snow = melted_snow + dm
			QQ = 0.
		
		else if(snowman(2) .gt. 0) then ! the entire grid cell melts 
			QQ = QQ - c_i*snowman(1)*(kelvin-snow_temp(1)) - snowman(1)*L_lh
			lwmass(1)  = lwmass(1) + snowman(1)
			melted_snow = melted_snow + snowman(1)
			l_heat = l_heat + snowman(1)*L_lh
			snowman(1) = 0.
			snow_temp(1) = snow_temp(2)
			rho_snow(1) = rho_snow(2)
			call go_regrid(1,1, snowman(:), snow_temp(:), rho_snow(:), lwmass(:), dummy_regrid)

		else ! entire grid cell melts and no second layer to melt. i.e. reset grid cell

			QQ = QQ - c_i*snowman(1)*(kelvin-snow_temp(1)) - snowman(1)*L_lh
			lwmass(1) = lwmass(1) + snowman(1)
			melted_snow = melted_snow + snowman(1)
			runoff_water = runoff_water + lwmass(1)		! TODO: if balance is not ok, think about this line again
			l_heat = l_heat + snowman(1)*L_lh
			snowman(1) = 0.
			snow_temp(1) = snow_temp(2)
			rho_snow(1) = rho_snow(2)
			
			dummy_melt_ice= max(-QQ/L_lh,0.)
			ice_melt= ice_melt + min(-QQ/L_lh/rho_ice/seconds_per_year,0.) ! this underestimates the melt since the ice is darker than the snow
			used_q =  used_q - QQ
			!if((ice_thickness(ix,iy) .gt. 0.).and.(ix .lt. 160)) then
			!	print*,'Overmelt into ice:', min(-QQ/L_lh/rho_ice,0.), ix, iy, time, 'ice thickness:',ice_thickness(ix,iy)
			!end if
			snowman(:)=0.
			lwmass(:)=0.
	 		snow_temp(:)=0. !Temperature to 273K 
			rho_snow(:)=rho_s
			QQ = 0.
		end if
! 	print*,snowman
	
	end do

	!dmdt= max((K_sw+Qp_lw+Qp_sh+Qp_lh-QQ) /L_lh,real(0))

	!if((ix==placex).and.(iy==placey)) then
	!print*,'QQ=',QQ,'manuel:',K_sw+Qp_lw+Qp_sh+Qp_lh,'effektiv:',K_sw+Qp_lw+Qp_sh+Qp_lh-QQ
	!end if

	!snowman(1)=snowman(1)-dmdt*dt_firn ! adjust mass

	!lwmass(1) = lwmass(1)+ dmdt*dt_firn

	 !   if(snowman(1)<0) then
	!! melt into the next layer and push water
	!! downwards. 

	!if(snowman(2)>0) then
	!if(snowman(2)+snowman(1)<0) then
	!!fprintf('Overovermelting in next layer needed!!! %2.1f at long: %d, lat: %d \n',snowman(ix,iy,2)+ &
	!!snowman(ix,iy,1), ix,iy);
	!print *,'Overovermelting in third layer needed, Melting subroutine ', snowman(2)+snowman(1), 'kg/m2'
	!			print *,'Mass first layer ', snowman(1), 'kg/m2'
	!			print *,'Mass secon layer ', snowman(2), 'kg/m2'
	!			print *,'Mass third layer ', snowman(3), 'kg/m2'
	!			print *,'Mass fourt layer ', snowman(4), 'kg/m2'
	!end if
	!! push meltwater and overmelted snow
	!lwmass(2) = lwmass(2) + lwmass(1) + snowman(1)
	!! adjust lwmass, mass and temperature of
	!! second box
	!! QQ a positive amount of energy that can be used to melt in the second layer
	!QQ= max(real(0),-snowman(1)*L_lh+ c_i*(snow_temp(2)-kelvin)*snowman(2))
	!snow_temp(2) = kelvin
	!snowman(2) = snowman(2) - QQ/L_lh
	!
	!! empty top box
	!snowman(1) =0
	!rho_snow(1) = rho_snow(2)
	!snow_temp(1) = snow_temp(2) ! maybe unnecessary
	!lwmass(1) = 0
	!
	!! Move all layers one level up
	!call go_regrid(1,1, snowman(:), snow_temp(:), &
	!rho_snow(:), lwmass(:))
	!
	!else
	!snowman(1)=0
	!lwmass(1)=0
	!snow_temp(1)=0
	!rho_snow(1)=rho_s
	!end if
	!end if

	!end if !snow_temp(1)>=kelvin

end subroutine go_melting_snow



!############################## PERCOLATION SUBROUTINE ###################################################################

subroutine go_percolation(snowman, lwmass, rho_snow, runoff)
    
    implicit none
	real(kind=8), intent(inout) :: snowman(n_snowlayer)
	real(kind=8), intent(inout) :: lwmass(n_snowlayer)
	real(kind=8), intent(inout) :: rho_snow(n_snowlayer)
	real(kind=8), intent(inout) :: runoff
	!real(kind=8), intent(inout) :: snow_temp(n_snowlayer)

	! local variables

	real(kind=8)  :: lwc
	real(kind=8)  :: percolating
	integer :: ii

	! if the whatercontent of one gridcell has grown over a certain %age
	! of the free bubble space, push everything above downwards

	! maybe its better to use a while loop, one could
	! abbort earlier

	!do ii=1,n_snowlayer,1
	runoff=0.
	ii=1
	do while (ii<=n_snowlayer)
		! calculate water content [fraction of free volume]
		if(snowman(ii)>0.)then

			if(rho_snow(ii)>rho_i-10.)then 
				! in case of very dense snow percolate all water downwards. 
				! we do so in order to avoid division by zero
				percolating = lwmass(ii)
				lwmass(ii) = 0.
				!runoff = percolating
				if(ii<n_snowlayer) then
					if(snowman(ii+1)>0.) then
						lwmass(ii+1) = lwmass(ii+1) + percolating
					else
						runoff = runoff + percolating
					end if
				else
					runoff = runoff + percolating
				end if

			else
				lwc=lwmass(ii)/snowman(ii)/rho_w/(1./rho_snow(ii) - 1./rho_i)
				if(lwc>max_lwc ) then
					! perform percolation
					!fprintf('Liquid water excess!!! %f at long: %d, lat: %d \n',lwc, 1,1);
					percolating = (lwc-max_lwc)*rho_w*snowman(ii)*(1./rho_snow(ii)-1./rho_i)
					lwmass(ii) =lwmass(ii)- percolating
					!runoff = percolating
					if(ii<n_snowlayer) then
						if(snowman(ii+1)>0) then
							lwmass(ii+1) = lwmass(ii+1) + percolating
						else
							runoff = runoff + percolating
						end if
					else
						runoff = runoff + percolating
					end if
			   	end if
			end if

		else
			 !   snowman(ii)=0
			 !   lwmass(ii)=0
			 !   snow_temp(ii)=0
			 !   rho_snow(ii)=rho_s
			ii=n_snowlayer
		end if
		ii=ii+1
	end do
	
end subroutine go_percolation


!############################## REFREEZING SUBROUTINE ###################################################################
subroutine go_refreezing(lwmass, snowman, rho_snow, snow_temp, refreeze, heat_fusion, ix,iy)
    
    implicit none
	real(kind=8), intent(inout) :: snowman(n_snowlayer)
	real(kind=8), intent(inout) :: lwmass(n_snowlayer)
	real(kind=8), intent(inout) :: snow_temp(n_snowlayer)
	real(kind=8), intent(inout) :: rho_snow(n_snowlayer)
	real(kind=8), intent(inout) :: refreeze
	real(kind=8), intent(inout) :: heat_fusion

	integer, intent(in) :: ix
	integer, intent(in) :: iy
	! local variables
    integer :: ii
	real(kind=8)  :: icecube

	refreeze = 0.
	heat_fusion = 0.

	do ii=1,n_snowlayer,1
		if((snowman(ii) .gt. 0.).and.(lwmass(ii) .gt. 0.)) then
			if( (kelvin-snow_temp(ii))*c_i*snowman(ii) .lt. lwmass(ii)*L_lh ) then
				! CASE 1
				! water freezes partly

				icecube = (kelvin-snow_temp(ii))*c_i*snowman(ii)/L_lh

				!if(icecube .lt. 0.)then
				!	print*,'lev ', ii
				!	print*,'x y ', ix,iy
				!	print*,'icecube',icecube
				!	print*,'lwmass',lwmass(ii)
				!	print*,'AB snowtemp',snow_temp(ii)
				!	print*,'snowmass',snowman(ii) 
				!	print*,'-----------'
				!end if



				heat_fusion = heat_fusion + icecube*L_lh
				!if ((ix==291).and.(iy==279)) then
					!print*,'Case 1'
					!print*,'lev ', ii
					!print*,'icecube',icecube
					!print*,'lwmass',lwmass(ii)
					!print*,'AB snowtemp',snow_temp(ii)
					!print*,'snowmass',snowman(ii)

				!end if

				snow_temp(ii) = kelvin
				rho_snow(ii) = rho_snow(ii)*(icecube + snowman(ii)) / snowman(ii)
				snowman(ii) = snowman(ii) + icecube
				lwmass(ii) = lwmass(ii) - icecube
				refreeze = refreeze + icecube
				!if (refreeze .lt. 0.) then

					!print*,refreeze
					!print*,lwmass(ii)
				!end if
			else
				! CASE 2
				! all water freezes
				snow_temp(ii) = (lwmass(ii)*L_lh/c_i + lwmass(ii)*kelvin + snow_temp(ii)*snowman(ii) )&
					/(lwmass(ii)+snowman(ii))
					
				rho_snow(ii) = rho_snow(ii)*(lwmass(ii) + snowman(ii)) / snowman(ii)
				snowman(ii) = lwmass(ii) + snowman(ii)
				refreeze = refreeze + lwmass(ii)
				heat_fusion = heat_fusion + lwmass(ii)*L_lh
				!if(lwmass(ii) .lt. 0.)then
				!	print*,'lev ', ii
				!	print*,'x y ', ix,iy
				!	print*,'icecube',icecube
				!	print*,'lwmass',lwmass(ii)
				!	print*,'AB snowtemp',snow_temp(ii)
				!	print*,'snowmass',snowman(ii) 
				!	print*,'-----------'
				!end if
				!if (refreeze .lt. 0.) then 	! TODO: remove this 
				!	print*,'Case 2'
				!	print*,refreeze
				!	print*,lwmass(ii)
				!end if
				lwmass(ii) = 0.


			end if
		end if
	end do

end subroutine go_refreezing

!############################## MELTING ICE SUBROUTINE ###################################################################
subroutine go_melting_ice( ice_melt, T_air, precipitation, P_sun, vaporflux, D_lf, dummy_melt_ice, dummy_rain_ice, DewpT, p_air)

	! water can melt the ice but doesn't freeze, rain and meltwater run off
    implicit none
	real(kind=8), intent(in) :: T_air   ! Temperature at the level of the ice elevation
	real(kind=8), intent(in) :: precipitation   ! Precipitation
	real(kind=8), intent(in) :: P_sun  ! sw radiation
    real(kind=8), intent(out) :: vaporflux
    real(kind=8), intent(in) :: DewpT!, intent(in) :: DewpT
    real(kind=8), intent(in) :: p_air
    real(kind=8) :: ea!, intent(in) :: ea
    real(kind=8) :: es!, intent(in) :: es
    
	! Output
	real(kind=8), intent(inout) :: ice_melt
	real(kind=8), intent(inout) :: D_lf
	real(kind=8), intent(inout) :: dummy_melt_ice
    real(kind=8), intent(inout) :: dummy_rain_ice
	! local variables

	real(kind=8) :: dQ_lh=0
	real(kind=8) :: dQ_sh=0
	real(kind=8) :: dQ_sw=0
	real(kind=8) :: dQ_lw=0
	real(kind=8) :: dQ_v=0
	real(kind=8) :: dQ_tot=0
	
	!artificial DewpT
! 	DewpT=T_air-kelvin-5
! 	D_lf=D_sf/1003*0.622*(L_v+L_lh)/10000
	ea=610.8*exp(17.27*DewpT/(DewpT+237.3))
    es=611.2*exp(22.46*(kelvin-273)/(272.62+(kelvin-273))) 
!     print*,DewpT
    vaporflux = D_lf/p_air*(ea-es)
!     print*,D_lf
!     print*,'es',es
!     print*,'ea',ea
!     print*,'ice',vaporflux*dt_firn
    

	dQ_lh = (T_air-kelvin)*precipitation*c_w*rho_w !,real(0))
	dQ_sh = D_sf*(T_air-kelvin)
	dQ_sw= P_sun*(1.-albedo_ice)
	dQ_lw=sigma*(eps_air*(T_air)**real(4)-eps_snow*(kelvin)**real(4))
	dQ_v=vaporflux
! 	dQ_v=0.
	vaporflux=dQ_v
	dQ_tot=dQ_lw+dQ_sh+dQ_lh+dQ_sw+dQ_v
    ice_melt= ice_melt - max(dQ_tot,real(0))/L_lh*dt_firn/rho_ice/seconds_per_year+dQ_v/(L_v+L_lh)*dt_firn/rho_ice/seconds_per_year  ! use ice density of ice model to get the right height
!     print*, dQ_tot
    dummy_melt_ice=max(dQ_tot,real(0))/L_lh*dt_firn
    dummy_rain_ice=precipitation*rho_w
!     print*,'end of melting ice'
	! change units in ice equivalents

end subroutine go_melting_ice

!############################## DEPLEETION SUBROUTINE ###################################################################
subroutine go_depleet_snowman(snowman, rho_snow, snow_temp, lwmass, smb_ice, passon_mass)
	! this subroutine removes the lowest gridboxes with densities larger thant rho_e and adds the snow to the smb_ice
    implicit none
	real(kind=8), intent(inout) :: snowman(nx,ny,n_snowlayer)
	real(kind=8), intent(inout) :: lwmass(nx,ny,n_snowlayer)
	real(kind=8), intent(inout) :: snow_temp(nx,ny,n_snowlayer)
	real(kind=8), intent(inout) :: rho_snow(nx,ny,n_snowlayer)
	real(kind=8), intent(inout) :: smb_ice(nx,ny)
	real(kind=8), intent(inout) :: passon_mass

	! local variables
	integer :: mm
	integer :: ix
	integer :: iy
	logical :: move = .true.
	logical :: hit_snow = .false.
		passon_mass = 0

	do ix=1,nx,1
		do iy=1,ny,1
			move = .true.
			hit_snow = .false.
			mm=n_snowlayer+1

			! empty lowest grid box
			smb_ice(ix,iy) = smb_ice(ix,iy) + snowman(ix,iy,n_snowlayer)/rho_ice/seconds_per_year
			!passon_mass = passon_mass + snowman(ix,iy,n_snowlayer) + lwmass(ix,iy,n_snowlayer)
			! reset snowman cell
			snowman(ix,iy,n_snowlayer) = 0.
			snow_temp(ix,iy,n_snowlayer) = 0.
			rho_snow(ix,iy,n_snowlayer) = rho_s
			lwmass(ix,iy,n_snowlayer) = 0.

			!do while (move) ! look for more snow to pass
			do while ((move).and.(mm.gt.1).and.(snowman(ix,iy,1).gt.0.))
				mm=mm-1
				if(snowman(ix,iy,mm) .ne. 0.) then
					hit_snow = .true.
				end if


				if(hit_snow) then
					if(rho_snow(ix,iy,mm) .ge. rho_pass_on) then
						! add ice to icesheet
						smb_ice(ix,iy) = smb_ice(ix,iy) + snowman(ix,iy,mm)/rho_ice/seconds_per_year
									!passon_mass = passon_mass + snowman(ix,iy,mm) + lwmass(ix,iy,mm)
						! reset snowman cell
						snowman(ix,iy,mm) = 0.
						snow_temp(ix,iy,mm) = 0.
						rho_snow(ix,iy,mm) = rho_s
						lwmass(ix,iy,mm) = 0.
					else
						move = .false.
					end if
				end if


			end do ! do in vertical direction

		end do  ! do in y-direction
	end do  ! do x-direction

end subroutine go_depleet_snowman


!############################## DEPLEETION SUBROUTINE ###################################################################
subroutine continous_snowman_depleet(snowman1, rho_snow1, snow_temp1, lwmass1, d_m_in, ice2ice, runoff, energy_ice2ice)
	! this subroutine removes snow at the lower end of the snowmodel and converts it a smb for the ice model
    implicit none
	real(kind=8), intent(inout) :: snowman1(n_snowlayer)
	real(kind=8), intent(inout) :: lwmass1(n_snowlayer)
	real(kind=8), intent(inout) :: snow_temp1(n_snowlayer)
	real(kind=8), intent(inout) :: rho_snow1(n_snowlayer)
	real(kind=8), intent(in) :: d_m_in ! TODO
	real(kind=8), intent(inout) :: ice2ice
	real(kind=8), intent(inout) :: runoff
	real(kind=8), intent(inout) :: energy_ice2ice

	! local variables
	integer :: nn
	real(kind=8) :: d_m
	real(kind=8) :: d_lw
	real(kind=8) :: perc

	!d_m = smb_ice1*seconds_per_year*rho_ice
	nn = n_snowlayer
	perc = 0.
	d_m = d_m_in
	ice2ice = 0.
	runoff = 0.
	energy_ice2ice = 0.

	!print*,'test1'

	do while (d_m .gt. 0.)
		!print*,'test2'
		!print*,'nn = ', nn 
		!print*,'d_m = ', d_m
		!pause
		if (d_m .gt. snowman1(nn) ) then
			!print*,'testa3'
			d_m = d_m - snowman1(nn)
			ice2ice = ice2ice + snowman1(nn)
			runoff = runoff + lwmass1(nn)
			energy_ice2ice =energy_ice2ice+ snowman1(nn)*c_i*snow_temp1(nn)

			rho_snow1(nn) = rho_s
			snow_temp1(nn) = 0.
			lwmass1(nn) = 0.
			snowman1(nn) = 0.
		else
			!print*,'testb3 d_m=',d_m
			!perc = (snowman1(nn) - d_m)/snowman1(nn)
			d_lw = d_m * lwmass1(nn) / snowman1(nn)
			snowman1(nn) = snowman1(nn) - d_m
			ice2ice = ice2ice + d_m
			energy_ice2ice =energy_ice2ice+ d_m*c_i*snow_temp1(nn)
			!print*, ice2ice
			runoff = runoff + d_lw

			!lwmass1(nn) = lwmass1(nn) * perc
			lwmass1(nn) = lwmass1(nn) - d_lw
			d_m = 0.
		end if

		nn = nn - 1
		!print*,'test4'
	end do	

	!print*,'test5'
	!print*, ice2ice
end subroutine continous_snowman_depleet



!############################## SUBROUTINE tridag ###################################################################
SUBROUTINE tridag(a,b,c,r,u,n)
	! from Numerical Recipes in Fortran77 Second Edition 1992
	! input
	implicit none
	integer, intent(in) :: n

	real(kind=8), intent(in) :: a(n)
	real(kind=8), intent(in) :: b(n)
	real(kind=8), intent(in) :: c(n)
	real(kind=8), intent(in) :: r(n)

	! output
	real(kind=8), intent(inout) :: u(n)

	! locals
	integer, parameter :: NMAX = 666

	integer :: j
	real(kind=8) :: bet
	real(kind=8) :: gam(NMAX)
	!Solves for a vector u(1:n) of length n the tridiagonal linear set given by equation (2.4.1).
	!a(1:n), b(1:n), c(1:n), and r(1:n) are input vectors and are not modified.

	!Parameter: NMAX is the maximum expected value of n.

	!One vector of workspace, gam is needed.
	!if(b(1).eq.0.)pause ???tridag: rewrite equations???
	!If this happens then you should rewrite your equations as a set of order N ??? 1, with u2
	!trivially eliminated.


	u=0.
	bet=b(1)
	u(1)=r(1)/bet
	do j=2,n
		!Decomposition and forward substitution.
		gam(j)=c(j-1)/bet
		bet=b(j)-a(j)*gam(j)
		!if(bet.eq.0.)pause ???tridag failed???
		!Algorithm fails; see below.
		u(j)=(r(j)-a(j)*u(j-1))/bet
	end do
	do j=n-1,1,-1
		!Backsubstitution.
		u(j)=u(j)-gam(j+1)*u(j+1)
	enddo
	return
END subroutine tridag

! !timekeeping
! subroutine tic(t1)
!   implicit none
!   real(kind=4), intent(inout)::t1
!   call cpu_time(t1)
! end subroutine tic
! 
! subroutine toc(t1,t2)
!   implicit none
!   real(kind=4), intent(inout)::t1
!   real(kind=4), intent(inout)::t2
!   call cpu_time(t2)
!   ! if (rank==0) print*,"Time Taken -->", real(t2-t1)
!   print*,"Time Taken -->", real(t2-t1)
! end subroutine toc

!############################## SUBROUTINE number of nonzero elements #####################################################
function num_nonz_el(var_in)
    implicit none
	real(kind=8), intent(in) :: var_in(n_snowlayer)
	integer :: istheresnow(n_snowlayer)
	integer :: num_nonz_el

	num_nonz_el = 0
	istheresnow=0
	where (var_in(:) .gt. 0.)
		istheresnow(:)=1
	end where
	num_nonz_el=sum(istheresnow)
end function num_nonz_el

end module smb_emb
