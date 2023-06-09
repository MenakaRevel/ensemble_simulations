#!/opt/local/bin/python
# -*- coding: utf-8 -*-
# import numpy as np

########################
#
# parameters list
#
########################

def timestep():
    return 86400 # outer timestep in seconds

def starttime():
    return (2000,1,1) # start date: [year,month,date]

def endtime():
    return (2021,1,1) # end date: [year,month,date]
                      # *note: this date is not included

def start_year():
    return 2000

def end_year():
    return 2020

# Input runoff forcing
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
    return "ERA5"

def ens_mem():
    return 20
    # if num == 1:
    #     return 49

    # if num == 2:
    #     return 10

    # if num == 3:
    #     return 1

    # if num == 4:
    #     return 21
    
    # number of ensemble members

def expname():
    # return "ECMWF049_050"
    # return "ECMWF049_150"
    # return "AMZ049"
    # return "AMZCAL049"
    # return "AMZ" #amazone
    # return "AMZCAL"
    # return "GLB" #global
    return "CONUS" #conus

# def runoff_mem():
#     return 10 # number of runoff ensemble members

# def manning_mem():
#     return 2 # number of runoff ensemble members

# def max_lat():
#     return 80. # maximum latitude of assimilation
#                # *note: SWOT ovservation is not available beyond 80 degs. this should be less or equal to 80
#                ## modified 2018-06-05

def distopen():
    return 0.50 # not needed for ERA20CM
    # corrupted runoff's percentage
    # 0.75 for original Data Assimilation simulation (25% reduced)
    # 1.25 for 25% increased simulation
    # 1.00 for simulation using 1 year before runoff
    # *note: also editing and and re-compile of control_inp at CaMa-Flood is nessessary

def diststd():
    return 0.25 # not needed for ERA20CM
    # noise to make runoff input to scatter ensembles

# def assimS():
#     return -75
#     # data Assimilation's Region (South Edge at latitude)
#     # *note: should be larger or equal to -80

# def assimN():
#     return 75
#     # data Assimilation's Region (North Edge at latitude)
#     # *note: should be smaller or equal to 80

# def assimW():
#     return -170
#     #return -68.25 # use this for disabling west side of the Amazon basin's observation
#     # data Assimilation's Region (West Edge at latitude)
#     # *note: should be larger or equal to -170

# def assimE():
#     return 170
#     # data Assimilation's Region (East Edge at latitude)
#     # *note: should be smaller or equal to 170

# def patch_size():
#     return 100
#     # the size of the local patch of LETKF(Local ** EnKF)
#     # 0: only 1 pixel (the pixel itself) belongs to its local patch
#     # 100: empirical local patch

# def err_expansion():
#     return 1.0
#     # variance-covariance expansion
#     # works well with 0.04

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

def spinup_flag():
    return 1 # 0: do spinup simulation simulation

def runname():
    return "ERA5"
    # if num == 1:
    #     return "E2O"

    # if num == 2:
    #     return "ERA20CM"

    # if num == 3:
    #     return "ELSE_KIM2009"

    # if num == 4:
    #     return "ERA5"

    # if num == 5:
    #     return "ECMWF"

def CaMa_dir():
    # return "/cluster/data6/menaka/CaMa-Flood_v395b_20191030"
    # return "/cluster/data6/menaka/CaMa-Flood_v396a_20200514"
    return "/cluster/data6/menaka/CaMa-Flood_v4"
    # directory of CaMa-Flood
    # indicate the directory of ./map or ./src and other folders

def mapname():
    # return "amz_06min"
    # return "glb_15min"
    return "conus_06min"

def org_dir():
    # return "/cluster/data6/menaka/ensemble_simulations"
    return "/cluster/data7/menaka/ensemble_simulations"

# def spinup_mode():
#     return 0
#     # 0: do spinup simulation for both (corrupted and true) simulation
#     # 1: do spin up only at corrupted simulation
#     # 2: do spin up only at true simulation
#     # 3: no spinup simulation at all
#     ### if initial restart file is ready, spinup simulation is no need

# def spinup_end_year():
#     return 2003

# def spinup_end_month():
#     return 12

# def spinup_end_date():
#     return 31

# def ovs_err():
#     return 0.1
#     # size of SWOT observation error in meters
#     # should be at least 0.02
#     # hope to be below 0.10

# def thersold():
#     return 0.60 # not needed for ERA20CM
#     # thersold to define the local patch

# def initial_infl():
#     return 1.01
#     # initial inflation parameter

# def MKLdir():
#     return "/opt/intel/compilers_and_libraries_2016.3.170/mac/mkl"
#     # directory of Intel MKL files
#     # Intel MKL is needed for doing data assimilation
#     # Please Download and Instal it to your System before running
#     # for more information --> https://software.intel.com/en-us/qualify-for-free-software/academicresearcher

# def output_er():
#     return 0
#     # setting for saving or deleting intermediate files
#     # 0 for saving & 1 for deleting
#     # those files may be more than 400GB, so erasing is recommended if not necessary

# def make_log():
#     return 1
    # setting for making log files
    # 1 is for making and 0 is for not making

def para_nums():
    return 40
    # setting number of parallels to run CaMa-Flood Model
    # defualt is 6, but may change depending on your system

# def slack_notification():
#     return 0
#     # setting for validating slack notification
#     # 1 for valid and 0 for invalid
#     # 0 is a default if you are not familiar with slack
#     # if you turn it to 1, you need to edit sendslack.py
#     # for more information refer https://api.slack.com/incoming-webhooks

# def ens_at_non():
#     return 1
#     # * At Recent version, ensemble generating random number is constant for full simulation.
#     # (For example, when ensemble 001 is corrupted with -0.1 at day 1, ensemble 001 will be always corrupted with 0.1 for full simulation period.)
#     # Previously, ensemble mean was used as an assimilated value for non-observed location.
#     # In this version, this treatment has changed and non-observed location is given with an ensemble value.
#     # To enable this new feature, set the return of params.py method “ens_at_non()”, “1”(DEFAULT SETTING).
#     # If you don’t want to use this, set it to “0”.

# # functions for corrupting manning coeffcient ###################
# # this is for corrupting manning coefficient at Corrupted Simulation
# # manning coefficient will be corrupted with random numbers generated from following functions
# # the random number is generated for each ensemble member
# # random number is made by gaussian noise of average = corruptman_base(), stddev = corruptman_std()
# def corruptman_base():
#     return 0.03

# def corruptman_std():
#     return 0.015

# def corruptele_base():
#     return 0.5 # not needed for ERA20CM

# def corruptele_std():
#     return 0.25 # not needed for ERA20CM

# def non_hgt():
#     return 7.0 # not needed for ERA20CM 
#     # nominal water height 

def cpu_nums():
    return 1
    # number of cpus used 

def version():
    return "v1.0.0 (updated 2021-06-17)\n Simulating Ensembles"
    # version for WSE assimilation / observation localization  
