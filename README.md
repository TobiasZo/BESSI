# BESSI
This is the Bergen Snow Simulator version 2 (BESSI 2). It is a firn model designed for long-timescale modelling (Born et al., 2019). The second version of this model has been developed with special attention the Greenland ice sheet, and now includes routines to directly calculate turbulent latent heat flux and the reduction in albedo due to snow aging (Zolles and Born 2019).

The version added here contains an additional subroutine that calculates the albedo and the net solar heat flux with different chosen albedo parametrizations. The energy flux calculation now calculates the turbulent latent heat flux, based on the residual method. The vapor flux has significant impact on the SMB on the local scale, for PD conditions. The manuscript describing the model adaption is currently under preparation.

You can download all necessary files in the tar file: The tar fail contains all the model code, and the scripts to run it. In the version used for the sensitivity analysis ensemble simulations. 
    -   README
    -   IceBern2D.f90
    -   variables.f90
    -   smb_emb.f90
    -   io.f90
    -   basic_compile.sh
    
Furthermore two files to show how the ensemble is created and run:
    -   Create_ensemble_matrix_for_analysis.R
    -   Home_test_jobscript.f90
    -   variables_to_replace_home.f90
    
The plotting routines are in a seperate tar file:
    - plotting_analysis.tar
    
The former are in general all files needed to run the model. smb_emb.f90 is the firn model BESSI, it contains all the necessary subroutines to run the model. It is called from the wrapper IceBern2D. The ice dynamics is switched of as it is not numerically stable for the current grid size. (Hopefully I can supply a simple wrapper soon) In the variables.f90 file all the output and input path and variables are specified. This ranges from the basic parameters subject to optimization like fresh snow albedo to the input and output paths. There are three different switches for specific subroutines concerning the energy balance calculation. The turbulent latent heat flux can be switched on and off, there are 4 different albedo parametrizations to choose from. Incoming longwave radiation can either be calculated from the atmospheric temperature or direct meteorological input can be taken. Based on previous results we recommend including the turbulent latent heat flux, using albedo parametrizations 4 or 5 and taking the incoming long-wave from climate data if it is available. 

You can specifiy at which frequency the output is written and which kind of output. Annual data, daily data and monthly data can be written every X,Y,Z years and a "focus" period were all the data is written can be chosen additionally. 

Parameters for albedo, turbulent exchange coefficients etc.  have been optimized relative to GRACE and 10m firn core temperature data. But we are currently working to include more data and better values can be provided. 


P.s.: The extra scripts mentioned under furthermore show which parameters we have been varying usually and are subject of optimization
