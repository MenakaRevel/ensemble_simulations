#! /usr/bin/bash
########################
#
# this program run the whole program
#
########################
#*** PBS seeting when needed
#PBS -q F20
#PBS -l select=1:ncpus=20:mem=40gb
#PBS -j oe
#PBS -m bea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N Ensemble_Sim

#
NCPUS=10
export OMP_NUM_THREADS=$NCPUS

# mkdir -p "/cluster/data7/menaka/ensemble_simulations"
#-----------------
# working directory
# cd $PBS_O_WORKDIR
# cd "/cluster/data6/menaka/ensemble_simulations"
cd "/cluster/data7/menaka/ensemble_simulations"
# cd "/work/a06/menaka/ensemble_simulations"

orgdir="/cluster/data6/menaka/ensemble_simulations"
echo ${orgdir}

# # copy source files
# cp -r ${orgdir}/params.py        ./params.py 
# cp -r ${orgdir}/run.py           ./run.py
# cp -r ${orgdir}/main_code.py     ./main_code.py
# cp -r ${orgdir}/prep_runoff.py   ./prep_runoff.py

# import virtual environment
source ~/.bashrc
source ~/.bash_conda

source activate pydef

which python

# link source codes
rm -rf ./run.py
rm -rf ./main_code.py

ln -sf ${orgdir}/src/run.py          ./run.py
ln -sf ${orgdir}/src/main_code.py    ./main_code.py

# run the simulations
echo "running the simulations..."

python run.py 

wait

echo "simulations done."

conda deactivate