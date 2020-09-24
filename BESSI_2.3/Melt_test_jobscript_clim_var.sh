#for run in {11..15..1}
#do

cd /work/zolles/backup/PhD/Programming/NEW/Ensemble_set_up/scp/

j=214
run=$j


i=$(echo ${run}p)
echo $i
#l=78

albedo_new="`sed -n $i input_samples/albedo_new_short2.sample`"
albedo_wet="`sed -n $i input_samples/albedo_wet_short2.sample`"
#albedo_ice="`sed -n $i input_samples/albedo_ice_short2.sample`"
albedo_ice=0.4
D_sf="`sed -n $i input_samples/D_sf_short2.sample`"
D_lf="`sed -n $i input_samples/D_lf_short2.sample`"
D_lf=1
e_air="`sed -n $i input_samples/e_air_short2.sample`"
e_air=-9
latent_heat_switch="`sed -n $i input_samples/latent_heat_switch_short2.sample`"
latent_heat_switch=1
max_lwc="`sed -n $i input_samples/max_lwc_short2.sample`"
max_lwc=0.1
albedo_module="`sed -n $i input_samples/albedo_module_short2.sample`"
vector=\'/work/zolles/backup/PhD/Programming/NEW/Ensemble_set_up/scp/input_samples/variable_climate_${l}_${k}.txt\'
fileshort=${l}_${k}
vector_file=/work/zolles/backup/PhD/Programming/NEW/Ensemble_set_up/scp/input_samples/variable_climate_${fileshort}.txt

echo run"$i".a_new"$albedo_new".a_wet"$albedo_wet".a_i"$albedo_ice".D_sf"$D_sf".D_lf"$D_lf".eps"$e_air".latent"$latent_heat_switch".a_mod"$albedo_module".lwc"$max_lwc"
parameters1=\'"${fileshort}"a_new"$albedo_new"_a_wet"$albedo_wet"_a_i"$albedo_ice"_D_sf"$D_sf"_D_lf"$D_lf"_eps"$e_air"_latent"$latent_heat_switch"_a_mod"$albedo_module"_lwc"$max_lwc"\'
parameters="${parameters1//./}"
echo $parameters

#convert switch to logic
if [ $latent_heat_switch -eq 1 ]
then
latent_heat_switch=.true.
else
latent_heat_switch=.false.
fi

outputpath=/work/zolles/modeloutput/climate_var3/
cp input_samples/variables_to_replace_melt.f90 insert_variables_file.f90
sed -i 's/a_new_foo/'${albedo_new}'/' insert_variables_file.f90
sed -i 's/a_wet_foo/'${albedo_wet}'/' insert_variables_file.f90
sed -i 's/a_ice_foo/'${albedo_ice}'/' insert_variables_file.f90
sed -i 's/D_sf_foo/'${D_sf}'/' insert_variables_file.f90
sed -i 's/D_lf_foo/'${D_lf}'/' insert_variables_file.f90
sed -i 's/eps_air_foo/'${e_air}'/' insert_variables_file.f90
sed -i 's/latent_heat_switch_foo/'${latent_heat_switch}'/' insert_variables_file.f90
sed -i 's/albedo_module_foo/'${albedo_module}'/' insert_variables_file.f90
sed -i 's/lwc_foo/'${max_lwc}'/' insert_variables_file.f90
sed -i 's|output_directory_foo|'${outputpath}'|' insert_variables_file.f90

sed -i 's|erai_vector_foo|'${vector}'|' insert_variables_file.f90
sed -i 's|sim_length_foo|'${l}'|' insert_variables_file.f90
cp IceBern2D_ERAi.f90 ./IceBern2D.f90
sed -i 's|sim_length_foo|'${l}'|' IceBern2D.f90

run_str=\'$(printf "%05d" $run)\'
cp io_new_for_abel.f90 ./io.f90
#cp IceBern2D_io_new_for_abel.f90 ./IceBern2D.f90

sed -i 's/run_foo/'${run_str}'/' io.f90
sed -i 's/run_foo/'${run_str}'/' smb_emb.f90
sed -i 's/parameters_foo/'${parameters}'/' io.f90

cp insert_variables_file.f90 variables.f90

echo "**** Start compling scripts ****"
exename="SMB_simu_${run}_${k}"
rm *.o *.mod SMB_simulation.x

gfortran -O2 -mcmodel=large -c variables.f90
echo "*** Variables compiled"
#gfortran -O3 -c variables_snow.f90
gfortran -O2 -mcmodel=large -I/work/zolles/libs/system/include -c io.f90 #/work/zolles/backup/PhD/Programming/NEW
echo "*** Io compiled"
gfortran -O2 -mcmodel=large -c smb_emb.f90
echo "*** SMB EMB compiled"
#gfortran -O3 -c smb_pdd.f90
gfortran -O2 -mcmodel=large -L/work/zolles/libs/system/lib -o ${exename}.x variables.o io.o smb_emb.o IceBern2D.f90 -lnetcdff -Wl,-rpath -Wl,/work/zolles/libs/system/lib
# gfortran -O3 -mcmodel=large -L/home/tobias/PhD/Programming/NEW/system/lib -lnetcdff -o ${exename}.x variables.o io.o smb_emb.o IceBern2D.f90
#gfortran -O3 -L/usr/local/netcdf/lib -lnetcdff -o LGM_compare_with_ensemble2_new_swradboa_damped_ampprecip.x variables.o io.o smb_emb.o IceBern2D.f90
#gfortran -O3 -L/usr/local/netcdf/lib -lnetcdff -o foo_gf_pdd_transient.x variables.o io.o smb_pdd.o smb_emb.o IceBern2D.f90
echo "*** Ice Bern compiled"


echo $exename".x"
echo "**** Run model ****"
./${exename}.x > output_${run_str}_$k.txt

run_num=$(printf "%05d" $run)
echo "**** change to output directory ****"
cp /work/zolles/backup/PhD/Programming/NEW/Ensemble_set_up/scp/test_$k.txt ${outputpath}ensemble_00214_${run_num}_$k*/

#for number in {402..500}
#do
#leading_zeros=$(printf "%07d" $number)
#inputname=ANNUAL_$leading_zeros.nc
#ncbo --op_typ=add CENTURY_${run_num}.nc $inputname CENTURY_${run_num}.nc -O
#done
## cp DAILY*.nc /work/projects/nn9591k/modeloutput/ensemble_ERAi/
## rm ANNUAL_??????[!0].nc
#echo "**** Clean up ****"
#cp ANNUAL_0000500.nc ANNUAL_500_${run_num}.nc
# rm ANNUAL_????[!5]??.nc
#rm ANNUAL_0000??[!0].nc
#rm IceBern2D_0000001.nc


# cp variables_to_replace_home.f90 insert_variables_file.f90
# sed -i 's/a_new_foo/'0.75'/' insert_variables_file.f90
# sed -i 's/a_wet_foo/'0.51'/' insert_variables_file.f90
# sed -i 's/a_ice_foo/'0.33'/' insert_variables_file.f90
# sed -i 's/D_sf_foo/'14'/' insert_variables_file.f90
# sed -i 's/D_lf_foo/'0.76'/' insert_variables_file.f90
# sed -i 's/eps_air_foo/'0.78'/' insert_variables_file.f90
# sed -i 's/latent_heat_switch_foo/'.true.'/' insert_variables_file.f90
# sed -i 's/albedo_module_foo/'4'/' insert_variables_file.f90
# sed -i 's/lwc_foo/'0.1039'/' insert_variables_file.f90
# sed -i 's/output_directory_foo/'/home/tobias/PhD/Programming/NEW/modeloutput/'/' insert_variables_file.f90
# 
# 
# cp insert_variables_file.f90 variables.f90

