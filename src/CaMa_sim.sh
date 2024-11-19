
#!/bin/sh
#==========================================================
# CaMa-Flood sample go script (1) global 15min simulation
# -- Multi 1-year simulations (2000 spinup -> 2000 -> 2001)
# -- Daily runoff forcing (plain binary) at 1deg resolution
#
# (C) D.Yamazaki & E. Dutra  (U-Tokyo/FCUL)  Aug 2019
#
# Licensed under the Apache License, Version 2.0 (the "License");
#   You may not use this file except in compliance with the License.
#   You may obtain a copy of the License at: http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License is 
#  distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the License for the specific language governing permissions and limitations under the License.
#==========================================================

#*** PBS setting when needed
#PBS -q E20
#PBS -l select=1:ncpus=20:mem=10gb
#PBS -j oe
#PBS -m ea
#PBS -V
#================================================
#edited by Menaka@IIS for ensemble simulations 2020/07/16
#================================================
# input settings
orgDIR=`pwd`

syear=$1

eyear=$2

ens_num=$3

CAMADIR=$4

cpunums=$5

runname=$6

mapname=$7

expname=$8

spinup_flag=$9

CaMa_opt=${10}

#e2oname=${10}
#================================================
# (0) Basic Setting (for workstation)

#*** 0a. Set CaMa-Flood base directory
BASE=$CAMADIR
OUTBASE=`pwd`
INBASE=`pwd`

echo $BASE

#*** 0b. Set dynamic library if needed
export IFORTLIB="/opt/intel/lib:/opt/intel/mkl/lib"
export DYLD_LIBRARY_PATH="${IFORTLIB}:${DYLD_LIBRARY_PATH}"

#*** 0c. OpenMP thread number
export OMP_NUM_THREADS=$cpunums                    # OpenMP cpu num

#================================================
# (1) Experiment setting
# -- some non-default options can be modified in NAMELIST section 

#============================
#*** 1a. Experiment directory setting
EXP=$expname$runname$ens_num                # experiment name (output directory name)
RDIR=${OUTBASE}/CaMa_out/${EXP}             # directory to run CaMa-Flood
EXE="MAIN_cmf"                              # Execute file name
PROG=${BASE}/src/${EXE}                     # location of Fortran main program
NMLIST="./input_cmf.nam"                    # standard namelist
LOGOUT="./log_CaMa.txt"                     # standard log output

echo $EXP
#============================
#*** 1b. Model physics option
DT=86400                                    # base DT (modified in physics loop by LADPSTP)
LADPSTP=".TRUE."                            # .TRUE. for adaptive time step

LFPLAIN=".TRUE."                            # .TRUE. to activate floodplain storage
LKINE=".FALSE."                             # .TRUE. to use kinematic wave equation
LFLDOUT=".TRUE."                            # .TRUE. to activate floodplain discharge
LPTHOUT=".FALSE."                           # .TRUE. to activate bifurcation flow, mainly for delta simulation
LDAMOUT=".FALSE."                           # .TRUE. to activate reservoir operation (under development)
LDAMYBY=".FALSE."                           # .TRUE. to use Year-By-Year dam activation scheme. .False. for All-reservoirs-in scheme
LiVnorm=".FALSE."                           # .TRUE. to use Noemal Volume as initial reservoir storage. False for zero-additional storage.
if [ $CaMa_opt = "bif" ] || [ $CaMa_opt = "all" ];then
     LPTHOUT=".TRUE."                       # .TRUE. to activate bifurcation flow, mainly for delta simulation
fi
# LDAMOUT=".FALSE."                           # .TRUE. to activate reservoir operation (under development)
if [ $CaMa_opt = "dam" ] || [ $CaMa_opt = "all" ];then
     LDAMOUT=".TRUE."                       # .TRUE. to activate reservoir operation (under development)
fi

#============================
#*** 1c. simulation time
YSTA=$syear                                 # start year ( from YSTA / Jan  1st _ 00:00)
YEND=$eyear                                 # end   year (until YEND / Dec 31st _ 24:00)
SPINUP=$spinup_flag                         # [0]: zero-storage start, [1]: from restart file
NSP=1                                       # spinup repeat time


#============================
#*** 1d. spinup setting

#* input restart file
LRESTART="" # see (3) set each year         # TRUE. to use restart initial condition
CRESTSTO="" # see (3) set each year         # input restart FIle
LSTOONLY=".FALSE."                          # .TRUE. for storage only restart (for assimilation)

#* output restart file
CRESTDIR="./"                               # output restart file directory
CVNREST="restart"                           # output restart file prefix
LRESTCDF=".FALSE."                          # .TRUE. to use netCDF restart file
IFRQ_RST="0"                                # output restat frequency.
                                            # [0]: only at last time, [1,2,3,...,24] hourly restart, [30]: monthly restart


#============================
#*** 1e. forcing setting
IFRQ_INP="24"                               # input forcing frequency: [1,2,3,...,24] hour
DROFUNIT="86400000"   # [mm/day->m/s]       # runoff unit conversion
if [ $runname = "E2O" ];then
     DROFUNIT="86400000"   # [mm/day->m/s]  # runoff unit conversion
elif [ $runname = "ERA20CM" ];then
     DROFUNIT="1000"       # [mm/s->m/s]    # runoff unit conversion
elif [ $runname = "ELSE_KIM2009" ];then
     DROFUNIT="86400000"   # [mm/day->m/s]  # runoff unit conversion
elif [ $runname = "VIC_BC" ];then
     DROFUNIT="86400000"   # [mm/day->m/s]  # runoff unit conversion
elif [ $runname = "S14FD" ];then
     DROFUNIT="86400000"   # [mm/day->m/s]  # runoff unit conversion
elif [ $runname = "isimip3a" ];then
     DROFUNIT="1000"       # [mm/day->m/s]  # runoff unit conversion
elif [ $runname = "ERA5" ];then
     DROFUNIT="86400"      # [m/day->m/s]  # runoff unit conversion
fi

#----- for plain binary runoff forcing
LINPCDF=".FALSE."                            # true for netCDF runoff
LINTERP=".TRUE."                            # .TRUE. to interporlate with input matrix
LINTERPCDF=".FALSE."                        # .TRUE. to use netCDF input matrix
###CROFDIR="${INBASE}"                         # runoff directory
####CROFPRE="Roff____"                          # runoff prefix/suffix  
####CROFSUF=".one"                              #   $(CROFPRE)YYYYMMDD$(CROFSUF)
###CROFDIR="${INBASE}/$runname/Roff"    #   runoff director
####CROFDIR="${INBASE}/CaMa_in/$runname"    #   runoff director
###CROFPRE="Roff__"
###CROFSUF="${ens_num}.one"
echo "$runname"
#----- for plain binary runoff forcing
LINPCDF=".FALSE."                            # true for netCDF runoff
LINTERP=".TRUE."                            # .TRUE. to interporlate with input matrix
LINTERPCDF=".FALSE."                        # .TRUE. to use netCDF input matrix
CROFDIR="${INBASE}/CaMa_in/${runname}/Roff" # runoff directory
CROFPRE="Roff__"                          # runoff prefix/suffix  
CROFSUF="${ens_num}.one"                    # runoff prefix/suffix 
# if [ $runname = "VIC_BC" ];then
#   #----- for plain binary runoff forcing
#   LINPCDF=".FALSE."                            # true for netCDF runoff
#   LINTERP=".TRUE."                            # .TRUE. to interporlate with input matrix
#   LINTERPCDF=".FALSE."                        # .TRUE. to use netCDF input matrix
#   CROFDIR="${INBASE}/CaMa_in/${runname}/Roff" # runoff directory
#   CROFPRE="Roff__"                          # runoff prefix/suffix  
#   CROFSUF="${ens_num}.one"                    # runoff prefix/suffix 
# elif [ $runname = "E2O" ];then
#   #----- for plain binary runoff forcing
#   LINPCDF=".FALSE."                            # true for netCDF runoff
#   LINTERP=".TRUE."                            # .TRUE. to interporlate with input matrix
#   LINTERPCDF=".FALSE."                        # .TRUE. to use netCDF input matrix
#   CROFDIR="${INBASE}/CaMa_in/${runname}/Roff" # runoff directory
#   CROFPRE="Roff__"                          # runoff prefix/suffix  
#   CROFSUF="${ens_num}.one"                    # runoff prefix/suffix 
# else
#   #----- for plain binary runoff forcing
#   LINPCDF=".FALSE."                            # true for netCDF runoff
#   LINTERP=".TRUE."                            # .TRUE. to interporlate with input matrix
#   LINTERPCDF=".FALSE."                        # .TRUE. to use netCDF input matrix
#   CROFDIR="${INBASE}/CaMa_in/${runname}/Roff" # runoff directory
#   CROFPRE="Roff__"                          # runoff prefix/suffix  
#   CROFSUF="${ens_num}.one"                    # runoff prefix/suffix 
# fi
# # if [ $runname = "VIC_BC" ];then
# # 	#----- for netCDF runoff forcing ###
# # 	LINPCDF=".TRUE."                              # true for netCDF runoff
# # 	LINTERP=".TRUE."                              # .TRUE. to interporlate with input matrix
# # 	LINTERPCDF=".FALSE."                          # .TRUE. to use netCDF input matrix
# # 	CROFDIR="${INBASE}"                            # runoff directory
# # 	CROFPRE="${runname}_runoff_"                  # runoff prefix/suffix  
# # 	CROFCDF=""     # see (3) set each year        # netCDF runoff file
# # 	CVNROF="runoff"                               # netCDF runoff    variable name
# # 	#CVNSUB=""                                     # netCDF runoffsub variable name
# # 	#SYEARIN=""     # see (3) set each year        #   netCDF runoff file, start date
# # 	#SMONIN=""      # see (3) set each year
# # 	#SDAYIN=""      # see (3) set each year
# # 	#SHOURIN=""     # see (3) set each year
# # 	#----- for plain binary runoff forcing
# # 	# LINPCDF=".FALSE."                           # true for netCDF runoff
# # 	# LINTERP=".TRUE."                            # .TRUE. to interporlate with input matrix
# # 	# LINTERPCDF=".FALSE."                        # .TRUE. to use netCDF input matrix
# # 	# CROFDIR="${INBASE}"                         # runoff directory
# # 	# CROFPRE="Roff____"                          # runoff prefix/suffix  
# # 	# CROFSUF=".one"                              #   $(CROFPRE)YYYYMMDD$(CROFSUF)
# #   # CROFSUF=".qtr"                              #   $(CROFPRE)YYYYMMDD$(CROFSUF)
# # 	# ###CROFDIR="${INBASE}/$runname/Roff"    #   runoff director
# # 	# ####CROFDIR="${INBASE}/CaMa_in/$runname"    #   runoff director
# # 	# ###CROFPRE="Roff__"
# # 	# ###CROFSUF="${ens_num}.one"
# # fi 
# # # #----- for netCDF runoff forcing ###
# # # LINPCDF=".TRUE."                              # true for netCDF runoff
# # # LINTERP=".TRUE."                              # .TRUE. to interporlate with input matrix
# # # LINTERPCDF=".FALSE."                          # .TRUE. to use netCDF input matrix
# # # CROFDIR="${INBASE}"                            # runoff directory
# # # CROFPRE="e2o_${e2oname}_wrr2_glob15_day_Runoff_"   # runoff prefix/suffix  
# # # CROFCDF=""     # see (3) set each year        # netCDF runoff file
# # # CVNROF="Runoff"                               # netCDF runoff    variable name
# # # #CVNSUB=""                                     # netCDF runoffsub variable name
# # # #SYEARIN=""     # see (3) set each year        #   netCDF runoff file, start date
# # # #SMONIN=""      # see (3) set each year
# # # #SDAYIN=""      # see (3) set each year
# # # #SHOURIN=""     # see (3) set each year

# # if [ $runname = "E2O" ];then
# # 	#----- for netCDF runoff forcing ###
# # 	LINPCDF=".TRUE."                              # true for netCDF runoff
# # 	LINTERP=".TRUE."                              # .TRUE. to interporlate with input matrix
# # 	LINTERPCDF=".FALSE."                          # .TRUE. to use netCDF input matrix
# # 	CROFDIR="${INBASE}"                            # runoff directory
# # 	CROFPRE="e2o_${e2oname}_wrr2_glob15_day_Runoff_"   # runoff prefix/suffix  
# # 	CROFCDF=""     # see (3) set each year        # netCDF runoff file
# # 	CVNROF="Runoff"                               # netCDF runoff    variable name
# # 	#CVNSUB=""                                     # netCDF runoffsub variable name
# # 	#SYEARIN=""     # see (3) set each year        #   netCDF runoff file, start date
# # 	#SMONIN=""      # see (3) set each year
# # 	#SDAYIN=""      # see (3) set each year
# # 	#SHOURIN=""     # see (3) set each year
# # fi 

###** sub-surface runoff scheme (not available with plain binary runoff)
LROSPLIT=".FALSE."                          # .TRUE. for sub-surface runoff
###CSUBDIR="NONE"                              # sub-surface runoff directory
###CSUBPRE="NONE"                              # sub-surface runoff prefix/suffix  
###CSUBSUF="NONE"                              #   $(PREFIX)YYYYMMDD$(SUFFIX)

#============================
#*** 1f. river map & topography
FMAP="${BASE}/map/${mapname}"                 # map directory

# Dam Parameter File
CDAMFILE="${FMAP}/dam_param.csv"             # dam parameter list

CDIMINFO="${FMAP}/diminfo_test-1deg.txt"    # dimention information file
CINPMAT=${FMAP}/inpmat_test-1deg.bin        # runoff input matrix for interporlation
#CDIMINFO="${FMAP}/diminfo_test-15min_nc.txt" # dimention information file
#CINPMAT=${FMAP}/inpmat_test-15min_nc.bin     # runoff input matrix for interporlation
if [ ${runname} = "E2O" ] ; then
     CDIMINFO="${FMAP}/diminfo-15min.txt" # dimention information file
     CINPMAT="${FMAP}/inpmat-15min.bin"     # runoff input matrix for interporlation
elif [ ${runname} = "ERA20CM" ] ; then
     CDIMINFO="${FMAP}/diminfo-1deg.txt"  # dimention information file
     CINPMAT="${FMAP}/inpmat-1deg.bin"      # runoff input matrix for interporlation
elif [ ${runname} = "ELSE_KIM2009" ] ; then
     CDIMINFO="${FMAP}/diminfo-1deg.txt"  # dimention information file
     CINPMAT="${FMAP}/inpmat-1deg.bin"      # runoff input matrix for interporlation
elif [ ${runname} = "VIC_BC" ] ; then
     CDIMINFO="${FMAP}/diminfo-05min.txt"  # dimention information file
     CINPMAT="${FMAP}/inpmat-05min.bin"      # runoff input matrix for interporlation
elif [ ${runname} = "isimip3a" ] ; then
     CDIMINFO="${FMAP}/diminfo-30min.txt"  # dimention information file
     CINPMAT="${FMAP}/inpmat-30min.bin"      # runoff input matrix for interporlation
elif [ ${runname} = "ERA5" ] ; then
     CDIMINFO="${FMAP}/diminfo-06min.txt"  # dimention information file
     CINPMAT="${FMAP}/inpmat-06min.bin"      # runoff input matrix for interporlation
fi

#----- for plain binary map input
#** basic topography
LMAPCDF=".FALSE."                           # .TRUE. for netCDF map
CNEXTXY="${FMAP}/nextxy.bin"                # downstream xy (river network map)
CGRAREA="${FMAP}/ctmare.bin"                # unit-catchment area   [m2]
CELEVTN="${FMAP}/elevtn.bin"                # channel top elevation [m]
CNXTDST="${FMAP}/nxtdst.bin"                # downstream distance   [m]
CRIVLEN="${FMAP}/rivlen.bin"                # channel length        [m]
CFLDHGT="${FMAP}/fldhgt.bin"                # floodplain elevation profile (height above 'elevtn') [m]

#** channel parameter
###CRIVWTH=${FMAP}/rivwth.bin"              # channel width [m] (empirical power-low)
CRIVWTH="${FMAP}/rivwth_gwdlr.bin"          # channel width [m] (GWD-LR + filled with empirical)
CRIVHGT="${FMAP}/rivhgt.bin"                # channel depth [m] (empirical power-low)
# CRIVHGT="${FMAP}/rivhgt_Xudong.bin"         # channel depth [m] (Xudong et al 2021)
CRIVMAN="${FMAP}/rivman.bin"                # manning coefficient river (The one in flood plain is a global parameter; set $PMANFLD below.)

#** bifurcation channel info
CPTHOUT="${FMAP}/bifprm.txt"                #   bifurcation channel list

###** groundwater delay (not available in plain binary runoff/map)
LGDWDLY=".FALSE."                           # .TRUE. to actuvate groundwater delay
#CGDWDLY=""                                 # ground water delay map

###** mean sea level
LMEANSL=".FALSE."                           # .TRUE. to use mean sea level data
#CMEANSL=""                                 # mean sea level map

#----- for netCDF map input 
###LMAPCDF=".TRUE."                         # .TRUE. for netCDF map
###CRIVCLINC=""                             # netCDF topography map
###CRIVPARNC=""                             # netCDF river parameters
###CMEANSLNC=""                             # netCDF mean sea level


#============================
#*** 1g. Dynamic Boundary Sea Level (not default)
LSEALEV=".FALSE."                           # .TRUE. to activate dynamic sea level 
###LSEALEVCDF=".FALSE."
###CSEALEVDIR="NONE"                        # Sea level boundary DIRECTORY
###CSEALEVPRE="NONE"                        # Sea level boundary PREFIX
###CSEALEVSUF="NONE"                        # Sea level boundary SUFFIX
###CSEALEVCDF="NONE"                        # * Sea level netCDF file name
###CVNSEALEV="sealev"                       # * Sea Level netCDF variable name
###SYEARSL=1                                # * netCDF sea level start year
###SMONSL=1                                 # * netCDF sea level start year
###SDAYSL=1                                 # * netCDF sea level start year
###SHOURSL=0                                # * netCDF sea level start year
###NSTATIONS=1                              # sea level data points
###CSLMAP="NONE"                            # sea level sta->XY conversion table

#============================
#*** 1h. Output Settings 
LOUTPUT=".TRUE."                            # .TRUE. to use CaMa-Flood standard output
IFRQ_OUT=24                                 # output frequency: [1,2,3,...,24] hour

LOUTCDF=".FALSE."                           # .TRUE. netCDF output, .FALSE. plain binary output
COUTDIR="./"                                # output directory 
CVARSOUT="sfcelv,outflw" # list output variable (comma separated)
# CVARSOUT="rivout,rivsto,rivdph,rivvel,fldout,fldsto,flddph,fldfrc,fldare,sfcelv,outflw,storge,pthflw,pthout,maxsto,maxflw,maxdph" # list output variable (comma separated)
COUTTAG=""  # see (3) set each year         #   output tag $(COUTDIR)/$(VARNAME)$(OUTTAG).bin

##### Model Parameters ################
PMANRIV="0.03D0"                            # manning coefficient river
PMANFLD="0.10D0"                            # manning coefficient floodplain
PCADP="0.7"                                 # satety coefficient for CFL condition
PDSTMTH="10000.D0"                          # downstream distance at river mouth [m]



#================================================
# (2) Initial setting

#*** 2a. create running dir 
mkdir -p ${RDIR}
cd ${RDIR}

#*** 2b. for new simulation, remove old files in running directory

if [ ${SPINUP} -eq 0 ]; then
  rm -rf ${RDIR}/????-sp*
  rm -rf ${RDIR}/*.bin
  rm -rf ${RDIR}/*.pth
  rm -rf ${RDIR}/*.vec
  rm -rf ${RDIR}/*.nc
  rm -rf ${RDIR}/*.log
  rm -rf ${RDIR}/*.txt
  rm -rf ${RDIR}/restart*
else
  NSP=0  # restart, no spinup
fi



#================================================
# (3) For each simulation year, modify setting
#--  loop 1-year simulation from $YSTART to $YEND

ISP=1           ## spinup count
IYR=${YSTA}     ## curent year
while [ ${IYR} -le ${YEND} ];
do 
  CYR=`printf %04d ${IYR}`   ## update file name, bugfix in v3.96a

  #*** 3a. modify restart setting
  if [ ${SPINUP} -eq 0 ];then
    LRESTART=".FALSE."                  ## from zero storage
    CRESTSTO=""
  else
    LRESTART=".TRUE."
    CRESTSTO="${CVNREST}${CYR}010100.bin"    ## from restart file
#    CRESTSTO="${CVNREST}${CYR}010100.nc"    ## from restart file
#    echo $CRESTSTO
  fi

  #*** 3b. update start-end year
  SYEAR=$IYR
  SMON=1
  SDAY=1
  SHOUR=0

  EYEAR=`expr $SYEAR + 1`
  EMON=1
  EDAY=1
  EHOUR=0

  ln -sf $PROG $EXE

  #*** 3c. update input / output file data
  CSYEAR=`printf %04d ${SYEAR}`
  COUTTAG=${CSYEAR}                  # output file tag

  # CROFCDF="${CROFDIR}/${CROFPRE}${CSYEAR}.nc"  # input netCDF runoff file
  # SYEARIN=$IYR
  # SMONIN=1
  # SDAYIN=1
  # SHOURIN=0
  #echo $CROFCDF


#================================================
# (4) Create NAMELIST for simulation year
# it is OK to remove optional variables (set to default in CaMa-Flood)

rm -f ${NMLIST}

#*** 0. config
cat >> ${NMLIST} << EOF
&NRUNVER
&NRUNVER
LADPSTP  = ${LADPSTP}                  ! true: use adaptive time step
LPTHOUT  = ${LPTHOUT}                  ! true: activate bifurcation scheme
LDAMOUT  = ${LDAMOUT}                  ! true: activate dam operation (under development)
LRESTART = ${LRESTART}                 ! true: initial condition from restart file
LSTOONLY = ${LSTOONLY}                 ! true: storage only restart (mainly for data assimilation)
/
&NDIMTIME
CDIMINFO = "${CDIMINFO}"               ! text file for dimention information
DT       = ${DT}                       ! time step length (sec)
IFRQ_INP = ${IFRQ_INP}                 ! input forcing update frequency (hour)
/
&NPARAM
PMANRIV  = ${PMANRIV}                  ! manning coefficient river
PMANFLD  = ${PMANFLD}                  ! manning coefficient floodplain
PDSTMTH  = ${PDSTMTH}                  ! downstream distance at river mouth [m]
PCADP    = ${PCADP}                    ! CFL coefficient
/
EOF

#*** 1. time
cat >> ${NMLIST} << EOF
&NSIMTIME
SYEAR   = ${SYEAR}                     ! start year
SMON    = ${SMON}                      !  month 
SDAY    = ${SDAY}                      !  day 
SHOUR   = ${SHOUR}                     !  houe
EYEAR   = ${EYEAR}                     ! end year
EMON    = ${EMON}                      !  month 
EDAY    = ${EDAY}                      !  day 
EHOUR   = ${EHOUR}                     !  hour
/
EOF

#*** 2. map
cat >> ${NMLIST} << EOF
&NMAP
LMAPCDF    = ${LMAPCDF}                ! * true for netCDF map input
CNEXTXY    = "${CNEXTXY}"              ! river network nextxy
CGRAREA    = "${CGRAREA}"              ! catchment area
CELEVTN    = "${CELEVTN}"              ! bank top elevation
CNXTDST    = "${CNXTDST}"              ! distance to next outlet
CRIVLEN    = "${CRIVLEN}"              ! river channel length
CFLDHGT    = "${CFLDHGT}"              ! floodplain elevation profile
CRIVWTH    = "${CRIVWTH}"              ! channel width
CRIVHGT    = "${CRIVHGT}"              ! channel depth
CRIVMAN    = "${CRIVMAN}"              ! river manning coefficient
CPTHOUT    = "${CPTHOUT}"              ! bifurcation channel table
CGDWDLY    = "${CGDWDLY}"              ! Groundwater Delay Parameter
CMEANSL    = "${CMEANSL}"              ! mean sea level
CRIVCLINC  = "${CRIVCLINC}"            ! * river map netcdf
CRIVPARNC  = "${CRIVPARNC}"            ! * river parameter netcdf (width, height, manning, ground water delay)
CMEANSLNC  = "${CMEANSLNC}"            ! * mean sea level netCDF
/
EOF

#*** 3. restart
cat >> ${NMLIST} << EOF
&NRESTART
CRESTSTO = "${CRESTSTO}"               ! restart file
CRESTDIR = "${CRESTDIR}"               ! restart directory
CVNREST  = "${CVNREST}"                ! restart variable name
LRESTCDF = ${LRESTCDF}                 ! * true for netCDF restart file
IFRQ_RST = ${IFRQ_RST}                 ! restart write frequency (1-24: hour, 0:end of run)
/
EOF

#*** 4. forcing
if [ ${LINPCDF} = ".FALSE." ]; then
cat >> ${NMLIST} << EOF
&NFORCE
LINPCDF  = ${LINPCDF}                  ! true for netCDF runoff
LINTERP  = ${LINTERP}                  ! true for runoff interpolation using input matrix
CINPMAT  = "${CINPMAT}"                ! input matrix file name
DROFUNIT = ${DROFUNIT}                 ! runoff unit conversion
CROFDIR  = "${CROFDIR}"                ! runoff             input directory
CROFPRE  = "${CROFPRE}"                ! runoff             input prefix
CROFSUF  = "${CROFSUF}"                ! runoff             input suffix
/
EOF

elif [ ${LINPCDF} = ".TRUE." ]; then
cat >> ${NMLIST} << EOF
&NFORCE
LINPCDF  = ${LINPCDF}                  ! true for netCDF runoff
LINTERP  = ${LINTERP}                  ! true for runoff interpolation using input matrix
CINPMAT  = "${CINPMAT}"                ! input matrix file name
DROFUNIT = ${DROFUNIT}                 ! runoff unit conversion
CROFCDF  = "${CROFCDF}"                ! * netCDF input runoff file name
CVNROF   = "${CVNROF}"                 ! * netCDF input runoff variable name
SYEARIN  = ${SYEARIN}                  ! * netCDF input start year
SMONIN   = ${SMONIN}                   ! * netCDF input start year
SDAYIN   = ${SDAYIN}                   ! * netCDF input start year
SHOURIN  = ${SHOURIN}                  ! * netCDF input start year
/
EOF

fi # (if LINPCDF)

#*** 5. outputs
cat >> ${NMLIST} << EOF
&NOUTPUT
COUTDIR  = "${COUTDIR}"                ! OUTPUT DIRECTORY
CVARSOUT = "${CVARSOUT}"               ! Comma-separated list of output variables to save 
COUTTAG  = "${COUTTAG}"                ! Output Tag Name for each experiment
LOUTVEC  = .FALSE                      ! TRUE FOR VECTORIAL OUTPUT, FALSE FOR NX,NY OUTPUT
LOUTCDF  = ${LOUTCDF}                  ! * true for netcdf output false for binary
NDLEVEL  = 0                           ! * NETCDF DEFLATION LEVEL 
IFRQ_OUT = ${IFRQ_OUT}                 ! output data write frequency (hour)
/
EOF

#*** Opt. Reservoir Operation
cat >> ${NMLIST} << EOF
&NDAMOUT
CDAMFILE = "${CDAMFILE}"               ! Reservoir Parameter File
LDAMTXT  = .TRUE.                      ! True for text-based reservoir data output
LDAMH22  = .FALSE.                     ! True to use Hanazaki 2022 dam scheme. (False for Yamazaki&Funato scheme)
LDAMYBY  = ${LDAMYBY}                  ! .TRUE. to use Year-By-Year dam activation scheme. .False. for All-reservoirs-in scheme
LiVnorm  = ${LiVnorm}                  ! .TRUE. to use Normal Volume as initial reservoir storage. False for zero-additional storage.
/
EOF

#### 6. sea level (optional) 
#cat >> ${NMLIST} << EOF
#&NBOUND
#LSEALEVCDF =  ${LSEALEVCDF}            ! * true : netCDF sea level boundary
#CSEALEVDIR = "${CSEALEVDIR}"           ! Sea level boundary DIRECTORY
#CSEALEVPRE = "${CSEALEVPRE}"           ! Sea level boundary PREFIX
#CSEALEVSUF = "${CSEALEVSUF}"           ! Sea level boundary SUFFIX
#CSEALEVCDF = "${CSEALEVCDF}"           ! * Sea level netCDF file name
#CVNSEALEV  = "${CVNSEALEV}"            ! * Sea Level netCDF variable name
#SYEARSL    = ${SYEARSL}                ! * netCDF sea level start year
#SMONSL     = ${SMONSL}                 ! * netCDF sea level start year
#SDAYSL     = ${SDAYSL}                 ! * netCDF sea level start year
#SHOURSL    = ${SHOURSL}                ! * netCDF sea level start year
#NSTATIONS  = ${NSTATIONS}              ! sea level data points
#CSLMAP     = "${CSLMAP}                ! station to XY conversion table
#IFRQ_SL    = ${IFRQ_SL}                ! sea level boundary update frequency (min)
#/
#EOF

#================================================
# (5) Execute main program

echo "start: ${SYEAR}" `date`  >> log.txt
time ./${EXE}                  >> log.txt 
echo "end:   ${SYEAR}" `date`  >> log.txt

mv ${LOGOUT} log_CaMa-${CYR}.txt


#================================================
# (6) manage spin up

# if curent spinup time $ISP < required spinup time $NSP
#   copy the restart file restart$(IYR+1) to restart${IYR}
#   copy the outputs to directory "${IYR}-sp1"

SPINUP=1
if [ ${IYR} -eq ${YSTA} ];
then
  if [ ${ISP} -le ${NSP} ];
    then
    IYR1=`expr ${IYR} + 1`
    CYR1=`printf %04d ${IYR1}`
    cp -f ${CVNREST}${CYR1}010100.bin ${CVNREST}${CYR}010100.bin-sp${ISP}         2> /dev/null
    mv -f ${CVNREST}${CYR1}010100.bin ${CVNREST}${CYR}010100.bin                  2> /dev/null

    cp -f ${CVNREST}${CYR1}010100.bin.pth ${CVNREST}${CYR}010100.bin.pth-sp${ISP} 2> /dev/null
    mv -f ${CVNREST}${CYR1}010100.bin.pth ${CVNREST}${CYR}010100.bin.pth          2> /dev/null

    cp -f ${CVNREST}${CYR1}010100.nc ${CVNREST}${CYR}010100.nc-sp${ISP}           2> /dev/null
    mv -f ${CVNREST}${CYR1}010100.nc ${CVNREST}${CYR}010100.nc                    2> /dev/null

    mkdir -p ${CYR}-sp${ISP}
    mv -f ./${CVNREST}${CYR}010100.bin-sp${ISP}      ${CYR}-sp${ISP}  2> /dev/null
    mv -f ./${CVNREST}${CYR}010100.bin.pth-sp${ISP}  ${CYR}-sp${ISP}  2> /dev/null
    mv -f ./${CVNREST}${CYR}010100.nc-sp${ISP}       ${CYR}-sp${ISP}  2> /dev/null
    mv -f ./*${CYR}.bin                              ${CYR}-sp${ISP}  2> /dev/null
    mv -f ./*${CYR}.pth                              ${CYR}-sp${ISP}  2> /dev/null
    mv -f ./o_*${CYR}.nc                             ${CYR}-sp${ISP}  2> /dev/null
    mv -f ./*${CYR}.log                              ${CYR}-sp${ISP}  2> /dev/null
    mv -f ./log_CaMa-${CYR}.txt                      ${CYR}-sp${ISP}  2> /dev/null

    ISP=`expr ${ISP} + 1`
  else
    ISP=0
    IYR=`expr ${IYR} + 1`
  fi
else
  IYR=`expr ${IYR} + 1`
fi

#================================================
# (7) End of each year loop. Back to (3)

done # loop to next year simulation

exit 0
