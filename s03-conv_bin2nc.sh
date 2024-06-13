#!/bin/sh
#====================
# convert binary files to netCDF4 
# all the time steps will be in one file
# dimesnion (time, nx, ny)
# Menaka@IIS
# 2020/05/29
#===========================
#*** PBS setting when needed
#PBS -q F40
#PBS -l select=1:ncpus=40:mem=10gb
#PBS -j oe
#PBS -m ea
#PBS -M menaka@rainbow.iis.u-tokyo.ac.jp
#PBS -V
#PBS -N bin2nc
#===========================
# cd $PBS_O_WORKDIR
# cd "/cluster/data6/menaka/ensemble_simulations"
cd "/cluster/data7/menaka/ensemble_simulations"
#================================================
# OpenMP Thread number
NCPUS=40
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
#--
N=`python src/calc_days.py $syear $smonth $sdate $eyear $emonth $edate`

# # ens=1
# # ens_char=`printf %02d$ens`
# # echo "ens char "${ens_char}

# water surface elevation
varname="sfcelv"
#=================================================
# # # python conv_bin2nc.py $N $varname 
# # inputname=${expname}${runname}${ens_char}
# # # inputname=${expname}${runname}"001"
# # echo ./src/bin2nc $N $ssyear $eeyear $varname $mapname $inputname $CAMADIR $outdir
# # time ./src/bin2nc $N $ssyear $eeyear $varname $mapname $inputname $CAMADIR $outdir 
ens=1
# while [ $ens -le $ens_num ];
for ens in $(seq -f "%03g" 1 $ens_num)
do
    # ens_char=`printf %02d$ens`
    ens_char=$ens
    inputname=${expname}${runname}${ens_char}
    echo ./src/bin2nc $N $ssyear $eeyear $varname $mapname $inputname $CAMADIR $outdir
    time ./src/bin2nc $N $ssyear $eeyear $varname $mapname $inputname $CAMADIR $outdir &
    ## for parallel computation using multiple CPUs 
    NUM=`ps aux -U $USER | grep /src/bin2nc | wc -l | awk '{print $1}'`
    # echo $USER $NUM
    while [ $NUM -gt $NCPUS ];
    do
        sleep 1
        NUM=`ps aux -U $USER | grep /src/bin2nc | wc -l | awk '{print $1}'`
        echo $USER $NUM
    done
    # ens=$(( $ens + 1 ))
done
# ###########
# # discharge
# varname="outflw"
# #=================================================
# #python conv_bin2nc.py $N $varname
# inputname=${expname}${runname}"001"
# echo ./src/bin2nc $N $ssyear $eeyear $varname $mapname $inputname $CAMADIR $outdir
# time ./src/bin2nc $N $ssyear $eeyear $varname $mapname $inputname $CAMADIR $outdir 
# #ens=1
# while [ $ens -le $ens_num ];
# do
#     ens_char=`printf %03d$ens`
#     inputname=$expname$runname$ens_char
#     echo ./src/bin2nc $N $ssyear $eeyear $varname $mapname $inputname $CAMADIR $outdir
#     time ./src/bin2nc $N $ssyear $eeyear $varname $mapname $inputname $CAMADIR $outdir &
#     ens=$(( $ens + 1 ))
#     wait
# done

wait