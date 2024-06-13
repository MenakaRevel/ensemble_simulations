#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import sys
import os
from multiprocessing import Pool

# import params as pm
# copy mean and std to HydroDA
#=============================
def copy(inputlist):
    input_name=inputlist[0]
    mapname   =inputlist[1]
    tagout    =inputlist[2]
    runname   =inputlist[3]
    prefix    =inputlist[4]
    odir      =inputlist[5]
    #========================
    ens_char=input_name[-3::]
    # runname=pm.runname(pm.mode()) #+"050"
    # prefix="sfcelv_49"
    # prefix="cal_sfcelv_49"
    # prefix="sfcelv_corrupt"
    # prefix="sfcelv_bias"
    # prefix="sfcelv_bias_corrupt"
    # mean
    oname=odir+"/mean_"+prefix+"_"+runname+"_"+mapname+"_"+tagout+"_"+ens_char+".bin"
    # oname="/cluster/data6/menaka/HydroDA/dat/mean_sfcelv_"+runname+"_"+mapname+"_"+tagout+"_"+ens_char+".bin"
    # oname="/cluster/data6/menaka/HydroDA/dat/mean_sfcelv_49_"+runname+"_"+mapname+"_"+tagout+"_"+ens_char+".bin"
    # oname="/cluster/data6/menaka/HydroDA/dat/mean_cal_sfcelv_49_"+runname+"_"+mapname+"_"+tagout+"_"+ens_char+".bin"
    # oname="/cluster/data6/menaka/HydroDA/dat/mean_sfcelv_cal_"+runname+"_"+mapname+"_"+tagout+"_"+ens_char+".bin"
    iname="./CaMa_out/"+input_name+"/sfcelv_mean"+tagout+".bin"

    # print ("cp "+iname+" "+oname)
    os.system("cp "+iname+" "+oname)

    # std
    oname=odir+"/std_"+prefix+"_"+runname+"_"+mapname+"_"+tagout+"_"+ens_char+".bin"
    # oname="/cluster/data6/menaka/HydroDA/dat/std_sfcelv_"+runname+"_"+mapname+"_"+tagout+"_"+ens_char+".bin"
    # oname="/cluster/data6/menaka/HydroDA/dat/std_sfcelv_49_"+runname+"_"+mapname+"_"+tagout+"_"+ens_char+".bin"
    # oname="/cluster/data6/menaka/HydroDA/dat/std_cal_sfcelv_49_"+runname+"_"+mapname+"_"+tagout+"_"+ens_char+".bin"
    # oname="/cluster/data6/menaka/HydroDA/dat/std_sfcelv_cal_"+runname+"_"+mapname+"_"+tagout+"_"+ens_char+".bin"
    iname="./CaMa_out/"+input_name+"/sfcelv_std"+tagout+".bin"

    # print ("cp "+iname+" "+oname)
    os.system("cp "+iname+" "+oname)

    print (ens_char, oname)
    return
#=============================
syear  =int(sys.argv[1]) #pm.start_year()
eyear  =int(sys.argv[2]) #pm.end_year()
expname=sys.argv[3] #pm.expname() #"AMZ050049B" #"AMZ000049O" #"AMZ050049" #"AMZCAL049" 
runname=sys.argv[4] #pm.runname()
mapname=sys.argv[5] #pm.mapname()
ens_mem=int(sys.argv[6]) #pm.ens_mem()
ncpus  =int(sys.argv[7]) #pm.para_nums()
odir   =sys.argv[8] #pm.outdir()
prefix =sys.argv[9] #"sfcelv"
#=============================
inputlist=[]
syyyy="%04d"%(syear) #2000)
eyyyy="%04d"%(eyear) #2014) #2010)
tag=syyyy+"-"+eyyyy
for ens in np.arange(1,ens_mem+1):
    inputname="%s%s%03d"%(expname,runname,ens)
    inputlist.append([inputname,mapname,tag,runname,prefix,odir])
#==============
para=20
p=Pool(para)
p.map(copy,inputlist)
p.terminate()
#map(stats,inputlist)