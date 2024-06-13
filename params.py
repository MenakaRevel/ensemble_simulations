#!/opt/local/bin/python
# -*- coding: utf-8 -*-
# import numpy as np

########################
#
# parameters list
#
########################
# 1. experimental settings
def expname():
    # return "ECMWF049_050"
    # return "ECMWF049_150"
    # return "AMZ049"
    # return "AMZCAL049"
    # return "AMZ" #amazone
    # return "AMZCAL"
    return "GLB" #global
    # return "CONUS50" #conus

def mapname():
    # return "amz_06min"
    return "glb_15min"
    # return "conus_06min"

# 2. time settings
def timestep():
    return 86400 # outer timestep in seconds

def starttime():
    return (2000,1,1) # start date: [year,month,date]

def endtime():
    return (2002,1,1) # end date: [year,month,date]
                      # *note: this date is not included

def start_year():
    return 2000

def end_year():
    return 2001

# 3. input runoff forcing
def runname():
    # return "ERA5"
    return "isimip3a"
    # return "VIC_BC"

def rundir():
    return "/work/a04/julien/CaMa-Flood_v4/inp/isimip3a/runoff" #isimip3a
    # return "/work/a02/menaka/ERA5/bin" #ERA5
    # return "/work/a06/menaka/VIC_BC/bin/"

def inputdir():
    return"/work/a06/menaka/ensemble_simulations/CaMa_in"
    # return "/cluster/data6/menaka/HydroDA/inp"
    #return "/home/yamadai/data/Runoff/E2O/nc"
    #return "/cluster/data6/menaka/ensemble_org/CaMa_in/VIC_BC/Roff"
    # return "/cluster/data6/menaka/covariance/CaMa_in/VIC/Roff"
    # return "/cluster/data6/menaka/HydroDA/inp/"
    # return "/work/a02/menaka/"
    # return "./CaMa_in/"

def input():
    # return "ECMWF"
    # return "E2O"
    # return "ERA5"
    return "isimip3a"
    # return "VIC_BC"

# 4. parameters for perturbing runoff
def pertub_method():
    # return "simple"
    # return "normal"
    return "lognormal"

def outputdir():
    return "/work/a06/menaka/ensemble_simulations/CaMa_in" # folder to save perturbed runoff 

def val_alpha():
    return 1.0 - (1.0/150.0)

def val_beta():
    return 0.0

def val_E():
    return 0.30

def val_distopen():
    return 0.50 # not needed for ERA20CM
    # corrupted runoff's percentage
    # 0.75 for original Data Assimilation simulation (25% reduced)
    # 1.25 for 25% increased simulation
    # 1.00 for simulation using 1 year before runoff
    # *note: also editing and and re-compile of control_inp at CaMa-Flood is nessecesary

def val_diststd():
    return 0.25 # not needed for ERA20CM
    # noise to make runoff input to scatter ensembles

# 5. ensemble members
def ens_mem():
    return 50
    # number of ensemble members
     
def mode():
    return 4
    # parameter to change assimilation mode
    # 1: Earth2Obs, 2: ERA20CM, 3: -25% ELSE_KIM2009, 4: ERA5, 5: ECMWF

def run_flag():
    return 0
    # 0 run all simulations
    # 1 run only corrupted and assimilated simulations
    # 2 run only true and assimilated simulations
    # 3 run only assimilated simulation

# 6. CaMa-Flood related settings
def spinup_flag():
    return 0
    # 1: no spinup simulation simulation 
    # 0: do spinup simulation simulation

def CaMa_dir():
    # return "/cluster/data6/menaka/CaMa-Flood_v395b_20191030"
    # return "/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
    return "/cluster/data6/menaka/CaMa-Flood_v410"
    # directory of CaMa-Flood
    # indicate the directory of ./map or ./src and other folders

def org_dir():
    # return "/cluster/data6/menaka/ensemble_simulations"
    return "/cluster/data7/menaka/ensemble_simulations"

def para_nums():
    return 40
    # setting number of parallels to run CaMa-Flood Model
    # default is 6, but may change depending on your system

def cpu_nums():
    return 1
    # number of cpus used 

def version():
    return "v1.0.0 (updated 2023-04-26)\n Simulating Ensembles"
    # version for perturbation  
