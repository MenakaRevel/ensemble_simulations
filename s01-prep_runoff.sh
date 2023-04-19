#!/bin/sh
#====================
# create statitcs (e.g., mean, std) from netCDF4 files
# dimesnion (nx, ny)
# Menaka@IIS
# 2020/05/29
#===========================
#*** PBS setting when needed
#PBS -q F20
#PBS -l select=1:ncpus=20:mem=10gb
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N prep_runoff
#===========================
# cd $PBS_O_WORKDIR
# cd "/cluster/data6/menaka/ensemble_simulations"
cd "/cluster/data7/menaka/ensemble_simulations"

#================================================
# import virtual environment
source ~/.bashrc
source ~/.bash_conda

source activate pydef

which python

#================================================
# OpenMP Thread number
NCPUS=20
export OMP_NUM_THREADS=$NCPUS

USER=`whoami`

# input settings
syear=`python -c "import params; print (params.starttime()[0])"`
smonth=`python -c "import params; print (params.starttime()[1])"`
sdate=`python -c "import params; print (params.starttime()[2])"`
eyear=`python -c "import params; print (params.endtime()[0])"`
emonth=`python -c "import params; print (params.endtime()[1])"`
edate=`python -c "import params; print (params.endtime()[2])"`
# for reading bin
ssyear=`python -c "import params; print (params.start_year())"`
eeyear=`python -c "import params; print (params.end_year())"`
echo $ssyear" to "$eeyear
# names
CAMADIR=`python -c "import params; print (params.CaMa_dir())"`
outdir="./" #`python -c "import params; print (params.out_dir())"`
cpunums=`python -c "import params; print (params.cpu_nums())"`
mapname=`python -c "import params; print (params.mapname())"`
expname=`python -c "import params; print (params.expname())"`
runname=`python -c "import params; print (params.runname())"`
ens_num=`python -c "import params; print (params.ens_mem())"`
#=================================================
# prepare runoff pertubation
rundir="/work/a02/menaka/ERA5/bin"
outdir="/work/a06/menaka/ensemble_simulations/CaMa_in"
python ./src/prep_runoff.py $syear $eyear $ens_num $runname $rundir $outdir $method $beta $E $alpha $distopen $diststd $NCPUS $CAMADIR

# link the folder to ./CaMa_in
# need only of runoff was saved in another folder
ln -sf $outdir/$runname ./CaMa_in/$runname 

wait

conda deactivate