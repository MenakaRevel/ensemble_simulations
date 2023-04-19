 #!/bin/sh

## CaMa-Flood: simulation map directory & simulation output dorectory
## Below is the example to prepare graphs for the validation purposes 
## for the result of sample simulation "test2"
## src/discharge_validation.py : discharge validation

CAMADIR="/cluster/data6/menaka/CaMa-Flood_v4"

# ##### For test2 simulation #####
# MAPDIR="../../map/conus_06min"
# OUTDIR="../../out/test2-conus_06min"
# OBSDIR="./obs_sample/discharge"
# LIST="./obs_sample/discharge/discharge_list_conus_06min.txt"  ## Discharge location list [as a text file]
# ## validation project tag
# TAG="conus"
# #select the output file type [netcdf/bin]
# OUTPUT="netcdf"
#################################

##### For test1 simulation #####
# MAPDIR="/cluster/data6/menaka/CaMa-Flood_v4/map/glb_15min"
MAPDIR="/cluster/data6/menaka/CaMa-Flood_v4/map/amz_06min"
# MAPDIR="/cluster/data6/menaka/CaMa-Flood_v4/map/glb_06min"

# OUTDIR="../../out/test3-glb_15min"
# OUTDIR="../../out/test5-amz_06min"
# OUTDIR="../../out/test6-amz_06min_ERA5"
# OUTDIR="../../out/test7-amz_06min_ERA5nc"
# OUTDIR="../../out/test8-amz_06min_ERA5"
# OUTDIR="/cluster/data6/menaka/ensemble_org/CaMa_out/AMZERA5001"
# OUTDIR="/cluster/data6/menaka/ensemble_org/CaMa_out/AMZERA5001"
# OUTDIR="/cluster/data6/menaka/ensemble_simulations/CaMa_out/AMZE2O008"
# OUTDIR="/cluster/data6/menaka/ensemble_simulations/CaMa_out/AMZCAL049E2O022"
# OUTDIR="/cluster/data6/menaka/ensemble_simulations/CaMa_out/AMZ049E2O022"
# OUTDIR="/cluster/data6/menaka/ensemble_org/CaMa_out/GLBERA5001"
# OUTDIR="/cluster/data6/menaka/ensemble_org/CaMa_out/AMZE2O003"
# OUTDIR="/cluster/data6/menaka/ensemble_org/CaMa_out/AMZCALecmwf001"
# OUTDIR="/work/a04/julien/CaMa-Flood_v4/out/coupled-model2"
# OUTDIR="/cluster/data7/menaka/ensemble_simulations/CaMa_out/AMZ050049ECMWF001"
OUTDIR="/cluster/data7/menaka/ensemble_simulations/CaMa_out/AMZ000049ECMWF001"
# OUTDIR="../../out/test-dev_HanazakiDam"
# OBSDIR="../../obs_sample/discharge"
OBSDIR="/cluster/data6/menaka/GRDC_2019"
# LIST="../../obs_sample/discharge/discharge_list_glb_15min.txt"  ## Discharge location list [as a text file]
LIST=${MAPDIR}"/grdc_loc.txt"  ## Discharge location list [as a text file]
#### validation project tag
# TAG="glb_Dam"
# TAG="glb_E2O"
# TAG="amz_E2O"
# TAG="amz_ERA5_02"
# TAG="amz_ERA5_org"
# TAG="amz_E2O_cal"
# TAG="amz_49_E2O_cal"
# TAG="amz_49_E2O"
# TAG="glb_ERA5"
# TAG="amz_ecmwf_nocal"
# TAG="amz_ecmwf_cal"
# TAG="SWOTH08"
# TAG="ECMWF050049_001"
TAG="ECMWF000049_001"

###select the output file type [netcdf/bin]
OUTPUT="bin"
# OUTPUT="netcdf"
#################################

echo "MAPDIR, OUTDIR, OBSDIR= " $MAPDIR, $OUTDIR, $OBSDIR


## specify validation period
SYEAR=2000
SMON=1
SDAY=1
EYEAR=2014
EMON=12
EDAY=31

##########

rm -f map
rm -f out
rm -f obs
ln -sf $MAPDIR map
ln -sf $OUTDIR out
ln -sf $OBSDIR obs

rm -f  list.txt
ln -sf $LIST  list.txt

mkdir -p fig/discharge
mkdir -p txt/discharge
##########

# make validation figures for discharge
echo "### DISCHARGE VISUALIZATION"
python src/discharge.py $SYEAR $SMON $SDAY $EYEAR $EMON $EDAY $OUTPUT

##########

## figures
rm -rf   fig_${TAG}/discharge
mkdir -p fig_${TAG}
mv       fig/discharge    fig_${TAG}/discharge
rm -rf   fig
echo "\n### figures saved in directory: fig_${TAG}/discharge"

## validation data
rm -rf   txt_${TAG}/discharge
mkdir -p txt_${TAG}
mv       txt/discharge    txt_${TAG}/discharge
rm -rf   txt
echo "\n### validation data saved in directory: txt_${TAG}/discharge"