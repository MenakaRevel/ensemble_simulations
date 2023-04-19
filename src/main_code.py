#!/opt/local/bin/python
# -*- coding: utf-8 -*-

#libralies
import os
import itertools
import numpy as np
import sys
import errno
from multiprocessing import Pool
from multiprocessing import Process
import datetime
import functools
import numpy.random as rd
import os.path
import datetime as dt
import glob
import shutil
import scipy.linalg as spla

#external python codes
import params as pm
import prep_runoff as roff
#import src.letkf_lib as lb
#import src.calc_one_daybef as odb


########################
#
# main program
#
########################

#main code for LETKF
############
## main control function
############
def main_act():

    # make necessary directories 
    print ("initial")
    # creating output folders
    mkdir("CaMa_out")
    mkdir("CaMa_in")

    # copy settings to ./assim_out/
    #shutil.copy("params.py","./assim_out/")

    # prepare runoff ensemble
    # print ("prepare runoff pertubations")
    # roff.prepare_input()

    inputlist=[]
    syear="%04d"%(pm.start_year())
    eyear="%04d"%(pm.end_year())

    #--
    print ("run ensemble simulations")
    for ens_num in np.arange(1,pm.ens_mem()+1):
         inputlist.append([syear,eyear,"%03d"%(ens_num)])

    # Run CaMa-Flood Model (ensemble simulations)
    p=Pool(pm.para_nums())
    p.map(simulation,inputlist)
    p.terminate()

#    start_dt=datetime.date(start_year,start_month,start_date)
#    #start_dt=datetime.date(2008,8,1)#5,29)
#    end_dt=datetime.date(end_year,end_month,end_date)
#    days_count=(end_dt-start_dt).days
#    #days_count=(end_dt-datetime.date(1991,1,1)).days 
#
#    # run daily simulations
#    for day in np.arange(days_count):
#
#        running_dt=start_dt+datetime.timedelta(days=day)
#        yyyy='%04d' % (running_dt.year)
#        mm='%02d' % (running_dt.month)
#        dd='%02d' % (running_dt.day)
#
#        if dd=="01" and pm.slack_notification()==1:
#            os.system("source src/sendslack.sh python_notification DA in progress "+yyyy+" "+mm+" "+dd)
#
#        one_day_loop(yyyy,mm,dd,day)
#
#    # clean all intermediate files
#    if pm.output_er()==1:
#        os.system("rm -Rf ./CaMa_out/"+yyyy+"*")

################
## single loop program
################
def simulation(inputlist):
    syear=inputlist[0]
    eyear=inputlist[1]
    ens_num=inputlist[2]
    #--
    dir0    =pm.org_dir()
    mode    =pm.mode()
    run_name=pm.runname()
    cpunums =pm.cpu_nums()
    CaMa_dir=pm.CaMa_dir()
    mapname =pm.mapname()
    expname =pm.expname()
    spinup  =pm.spinup_flag()
    # inputname=pm.inputdir()
    #---
    print ("=============== start simulation: "+ens_num+" ===============")
    print (syear+" "+eyear+" "+ens_num+" "+CaMa_dir+" "+str(cpunums)+" "+str(run_name))
    os.system("source "+dir0+"/src/CaMa_sim.sh "+syear+" "+eyear+" "+ens_num+" "+CaMa_dir
    +" "+str(cpunums)+" "+str(run_name)+" "+mapname+" "+expname+" "+str(spinup))
################################
def one_day_loop(yyyy,mm,dd,day):
    print "================================ start loop of "+yyyy+" "+mm+" "+dd+" ========================================"
    #
#    # True Simulation ####################################################################
    if pm.run_flag() == 0 or pm.run_flag() == 2:
        one_day_sim([yyyy,mm,dd,"000","true"])

    # Corrupted Simulation (Open Loop) ###################################################
    if pm.run_flag() == 0 or pm.run_flag() == 1:
        ODM_inputlist=[]
        randlist=np.fromfile("./CaMa_out/randlist.bin",np.float32)

        # set for ensemble simulations
        ens_num=1
        for ens_num in np.arange(1,pm.ens_mem()+1):
            ODM_inputlist.append([yyyy,mm,dd,'%03d'%ens_num,"open"])

        # Run CaMa-Flood Model (ensemble simulations)
        p=Pool(pm.para_nums())
        p.map(one_day_sim,ODM_inputlist)
        p.terminate()

        # copy ensemble corrupted simulation forecasts from CaMa_out (xa)
        cprestart=[]
        for ens_num in np.arange(1,pm.ens_mem()+1):
            cprestart.append([yyyy,mm,dd,ens_num])
        p=Pool(pm.para_nums())
        p.map(copy_corrupted_sfcelv,cprestart)
        p.terminate()
        #copy_corrupted_sfcelv(yyyy,mm,dd,num) # changed to copy sfcelv

        # copy restart files/ no need for recalculation - MODIFIED @ Menaka 
        p=Pool(pm.para_nums())
        p.map(copy_corrupted_restart,cprestart)
        p.terminate()

    # Assimilated Simulation #############################################################
    ODM_inputlist=[]
    randlist=np.fromfile("./CaMa_out/randlist.bin",np.float32)

    # set for ensemble simulations
    for ens_num in np.arange(1,pm.ens_mem()+1):
        ODM_inputlist.append([yyyy,mm,dd,'%03d'%ens_num,"assim"])


    # Run CaMa-Flood Model (ensemble simulations)
    p=Pool(pm.para_nums())
    p.map(one_day_sim,ODM_inputlist)
    p.terminate()

    # Calculate Ensemble Mean
    #os.system("./src/make_nonassim "+yyyy+mm+dd+" "+str(pm.ens_mem())+" "+"A")

    # make forecasted value for assimilated simulation
    # do assimilation (LETKF)
    data_assim(yyyy,mm,dd,day)

    # make restart MODIFIED v.1.1.0
    mkrestart=[]
    for ens_num in np.arange(1,pm.ens_mem()+1):
        mkrestart.append([yyyy,mm,dd,"assim",'%03d'%ens_num])

    # Modify the restart file 
    p=Pool(pm.para_nums())
    p.map(make_restart,mkrestart)
    p.terminate()

    # store river variable files
    store_out(yyyy,mm,dd)

#    # make rivout @menaka
#    mkrivout=[]
#    for ens_num in np.arange(1,pm.ens_mem()+1):
#        mkrivout.append([yyyy,mm,dd,"assim",'%03d'%ens_num])
#    p=Pool(pm.para_nums())
#    p.map(make_rivout,mkrivout)
#    p.terminate()

    # clean files
    if pm.output_er()==1:
        bef_dt=datetime.date(int(yyyy),int(mm),int(dd))-datetime.timedelta(days=1)
        bef_yyyy='%04d' %bef_dt.year
        bef_mm='%02d' %bef_dt.month
        bef_dd='%02d' %bef_dt.day
        os.system("rm -Rf ./CaMa_out/"+bef_yyyy+bef_mm+bef_dd+"*")


#######################################################################################


############
## main program functions
############
def spin_up(): #used
    # run spin up simulation
    # 1 year spin up for calculating initial value
    # one simulation for true
    # ensmble simulation for open

#    if loop=="open":
#        if pm.spinup_mode()==2 or pm.spinup_mode()==3:
#            return 0
#
#    if loop=="true":
#        if pm.spinup_mode()==1 or pm.spinup_mode()==3:
#            return 0

    dir2=pm.CaMa_dir()
    cpunums = pm.cpu_nums()
    yyyy = "%04d"%(pm.spinup_end_year())
    print pm.spinup_mode()
    if pm.spinup_mode()==3:
      return 0

    inputlist=[] 

    if pm.spinup_mode()==0 or pm.spinup_mode()==2:
        inputlist.append([yyyy,"true",'000'])
        #spinup_loop(inputlist)

    if pm.spinup_mode()==0 or pm.spinup_mode()==1:
        for ens_num in np.arange(1,pm.ens_mem()+1):
             inputlist.append([yyyy,"open",'%03d'%ens_num])

    # Run spinup simulations
    p=Pool(pm.para_nums())
    p.map(spinup_loop,inputlist)
    p.terminate()

    print "======================= end spinup =========================="

    return 0
###########################
def spinup_loop(inputlist):
    # Run spinup simulation
    yyyy=inputlist[0]
    loop=inputlist[1]
    ens_num=inputlist[2]
    dir2=pm.CaMa_dir()
    cpunums=pm.cpu_nums()
    mode=pm.mode()
    run_name=pm.runname(mode)
    print  "%s for %03d"%(loop,int(ens_num))
    os.system("source src/spin_up.sh "+str(yyyy)+" "+str(loop)+" "+ens_num+" "+dir2+" "+str(cpunums)+" "+str(run_name))
    return 0
###########################
def one_day_sim(inputlist):
    yyyy=inputlist[0]
    mm=inputlist[1]
    dd=inputlist[2]
    ens_num=inputlist[3]
    looptype=inputlist[4]
    mode=pm.mode()
    run_name=pm.runname(mode)

    # program for running one day model

    bef_dt=datetime.date(int(yyyy),int(mm),int(dd))-datetime.timedelta(days=1)
    bef_yyyy='%04d' %bef_dt.year
    bef_mm='%02d' %bef_dt.month
    bef_dd='%02d' %bef_dt.day

    print "oneday loop for",yyyy,mm,dd,ens_num,looptype
    dir2=pm.CaMa_dir()
    if looptype=="true":
        distopen="1.0"
    else:
        distopen=str(pm.distopen())

    print yyyy+" "+mm+" "+dd+" "+ens_num+" "+dir2+" "+looptype
    cpunums = pm.cpu_nums()
    os.system("source src/oneday_sim.sh "+yyyy+" "+mm+" "+dd+" "+ens_num+" "+dir2+" "+looptype+" "+str(cpunums)+" "+str(run_name))

    if looptype=="true":
        # copying "restart file" to ./CaMa_in/
        thisday=datetime.date(int(yyyy),int(mm),int(dd))
        nxt_day=thisday+datetime.timedelta(days=1)
        orgrestf="CaMa_out/"+yyyy+mm+dd+"T"+ens_num+"/restart"+'%04d'%(nxt_day.year)+'%02d'%(nxt_day.month)+'%02d'%(nxt_day.day)+".bin"
        newrestf="CaMa_in/restart/"+looptype+"/restart"+'%04d'%(nxt_day.year)+'%02d'%(nxt_day.month)+'%02d'%(nxt_day.day)+"T000.bin"
        #os.system("cp "+orgrestf+" "+newrestf)
        copy_stoonly(orgrestf,newrestf)

        # copying "WSE" as "xa_m" in ./assim_out
        oldfname="CaMa_out/"+yyyy+mm+dd+"T"+ens_num+"/sfcelv"+yyyy+".bin"
        newfname="assim_out/xa_m/"+looptype+"/"+yyyy+mm+dd+"_xam.bin"
        os.system("cp "+oldfname+" "+newfname)

    return 0
########################### # modified to run paralle @Menaka 
def copy_corrupted_sfcelv(inputlist):
    yyyy = inputlist[0]
    mm   = inputlist[1] 
    dd   = inputlist[2]
    num  = inputlist[3]
    numch='%03d'%num
    fname="./CaMa_out/"+yyyy+mm+dd+"C"+numch+"/sfcelv"+yyyy+".bin"
    os.system("cp "+fname+" ./assim_out/ens_xa/open/"+yyyy+mm+dd+"_"+numch+"_xa.bin")
    return 0
########################### # modified not calculate restart again/ no chage in WSE in corrupted @Menaka
def copy_corrupted_restart(inputlist):
    yyyy = inputlist[0]
    mm   = inputlist[1]
    dd   = inputlist[2]
    num  = inputlist[3]
    nxt_day = datetime.date(int(yyyy),int(mm),int(dd)) + datetime.timedelta(days=1)
    n_yyyy='%04d' % (nxt_day.year)
    n_mm='%02d' % (nxt_day.month)
    n_dd='%02d' % (nxt_day.day)
    numch='%03d'%num
    fname="./CaMa_out/"+yyyy+mm+dd+"C"+numch+"/restart"+n_yyyy+n_mm+n_dd+".bin"
    #os.system("cp "+fname+" ./CaMa_in/restart/open/restart"+n_yyyy+n_mm+n_dd+"C"+numch+".bin")
    copy_stoonly(fname,"./CaMa_in/restart/open/restart"+n_yyyy+n_mm+n_dd+"C"+numch+".bin")
    print "copy restart",n_yyyy,n_mm,n_dd,"C"+numch
    return 0
###########################
def copy_stoonly(iname,oname): # for CaMa_Flood v395b
    org=np.fromfile(iname,np.float32).reshape(6,-1)
    org[0:2].tofile(oname)
    return 0
###########################
def assim_at_fort(yyyy,mm,dd,day): #previous --> used
    dir1=pm.CaMa_dir()+"/"
    thisday=datetime.date(int(yyyy),int(mm),int(dd))
    nxt_day=thisday+datetime.timedelta(days=1)
    os.system("src/data_assim "+str(pm.assimN())+" "+str(pm.assimS())+" "+str(pm.assimW())+" "+str(pm.assimE())+" "+yyyy+mm+dd+" "+str('%02d'%SWOT_day(yyyy,mm,dd))+" "+str(pm.patch_size())+" "+str(pm.ens_mem())+" "+str(day)+" "+str('%04d'%(nxt_day.year)+'%02d'%(nxt_day.month)+'%02d'%(nxt_day.day))+" "+str(pm.err_expansion())+" "+dir1)
    return 0
###########################
def data_assim(yyyy,mm,dd,day): # new data assimilation function (2017-06-30)
    errrand=np.random.normal(0,pm.ovs_err())
    dir1=pm.CaMa_dir()+"/"
    thisday=datetime.date(int(yyyy),int(mm),int(dd))
    nxt_day=thisday+datetime.timedelta(days=1)
    print '%02d'%(nxt_day.day)
    os.system("src/data_assim "+str(pm.assimN())+" "+str(pm.assimS())+" "+str(pm.assimW())+" "+str(pm.assimE())+" "+yyyy+mm+dd+" "+str('%02d'%SWOT_day(yyyy,mm,dd))+" "+str(pm.patch_size())+" "+str(pm.ens_mem())+" "+str(day)+" "+str('%04d'%(nxt_day.year)+'%02d'%(nxt_day.month)+'%02d'%(nxt_day.day))+" "+str(pm.err_expansion())+" "+dir1+" "+str(errrand)+" "+str(pm.ovs_err())+" "+str(pm.thersold()))
#    os.system("src/data_assim "+str(pm.assimN())+" "+str(pm.assimS())+" "+str(pm.assimW())+" "+str(pm.assimE())+" "+yyyy+mm+dd+" "+str('%02d'%SWOT_day(yyyy,mm,dd))+" "+str(pm.patch_size())+" "+str(pm.ens_mem())+" "+str(day)+" "+str('%04d'%(nxt_day.year)+'%02d'%(nxt_day.month)+'%02d'%(nxt_day.day))+" "+str(pm.err_expansion())+" "+dir1+" "+str(errrand)+" "+str(pm.ovs_err()))
#    os.system("src/data_assim_fld "+str(pm.assimN())+" "+str(pm.assimS())+" "+str(pm.assimW())+" "+str(pm.assimE())+" "+yyyy+mm+dd+" "+str('%02d'%SWOT_day(yyyy,mm,dd))+" "+str(pm.patch_size())+" "+str(pm.ens_mem())+" "+str(day)+" "+str('%04d'%(nxt_day.year)+'%02d'%(nxt_day.month)+'%02d'%(nxt_day.day))+" "+str(pm.err_expansion())+" "+dir1+" "+str(errrand)+" "+str(pm.ovs_err()))
    return 0
###########################
def make_init_storge():
    bef_yyyy='%04d' % (pm.spinup_end_year())
    bef_mm='%02d' % (pm.spinup_end_month())
    bef_dd='%02d' % (pm.spinup_end_date())
    os.system("./src/make_nonassim_init ./CaMa_out/spinup_open/storge"+str(pm.spinup_end_year())+".bin "+"./assim_out/nonassim/open/nonasmC"+bef_yyyy+bef_mm+bef_dd+".bin")
    os.system("./src/make_nonassim_init ./CaMa_out/spinup_open/storge"+str(pm.spinup_end_year())+".bin "+"./assim_out/nonassim/assim/nonasmA"+bef_yyyy+bef_mm+bef_dd+".bin")
    return 0
###########################
def make_initial_restart(): # updated the name
    start_year,start_month,start_date=pm.starttime()
    yyyy="%04d"%(start_year)
    mm="%02d"%(start_month)
    dd="%02d"%(start_date)
    spinup_true="%04d%2d%02dT000"%(pm.spinup_end_year(),pm.spinup_end_month(),pm.spinup_end_date()) 
    #os.system("cp ./CaMa_out/"+spinup_true+"/restart"+yyyy+mm+dd+".bin ./CaMa_in/restart/true/restart"+yyyy+mm+dd+"T000.bin")
    copy_stoonly("./CaMa_out/"+spinup_true+"/restart"+yyyy+mm+dd+".bin","./CaMa_in/restart/true/restart"+yyyy+mm+dd+"T000.bin")

    print "cp ./CaMa_out/"+spinup_true+"/restart"+yyyy+mm+dd+".bin ./CaMa_in/restart/true/restart"+yyyy+mm+dd+"T000.bin"
    for num in np.arange(1,pm.ens_mem()+1):
        numch='%03d'%num
        spinup_open="%04d%2d%02dC%03d"%(pm.spinup_end_year(),pm.spinup_end_month(),pm.spinup_end_date(),num) 
        #os.system("cp ./CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin ./CaMa_in/restart/open/restart"+yyyy+mm+dd+"C"+numch+".bin")
        #os.system("cp ./CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin ./CaMa_in/restart/assim/restart"+yyyy+mm+dd+"A"+numch+".bin")
        copy_stoonly("./CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin","./CaMa_in/restart/open/restart"+yyyy+mm+dd+"C"+numch+".bin")
        copy_stoonly("./CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin","./CaMa_in/restart/assim/restart"+yyyy+mm+dd+"A"+numch+".bin")
###########################
def make_initial_restart_one(): # updated the name
    # copy restartyyyymmddC001 as restart for all simulations 
    start_year,start_month,start_date=pm.starttime()
    yyyy="%04d"%(start_year)
    mm="%02d"%(start_month)
    dd="%02d"%(start_date)
    spinup_true="%04d%2d%02dT000"%(pm.spinup_end_year(),pm.spinup_end_month(),pm.spinup_end_date())
    #os.system("cp ./CaMa_out/"+spinup_true+"/restart"+yyyy+mm+dd+".bin ./CaMa_in/restart/true/restart"+yyyy+mm+dd+"T000.bin")
    copy_stoonly("./CaMa_out/"+spinup_true+"/restart"+yyyy+mm+dd+".bin","./CaMa_in/restart/true/restart"+yyyy+mm+dd+"T000.bin")
    print "cp ./CaMa_out/"+spinup_true+"/restart"+yyyy+mm+dd+".bin ./CaMa_in/restart/true/restart"+yyyy+mm+dd+"T000.bin"
    spinup_open=spinup_true
    #spinup_open="%04d%2d%02dC%03d"%(pm.spinup_end_year(),pm.spinup_end_month(),pm.spinup_end_date(),1)
    for num in np.arange(1,pm.ens_mem()+1):
        numch='%03d'%num
        #os.system("cp ./CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin ./CaMa_in/restart/open/restart"+yyyy+mm+dd+"C"+numch+".bin")
        #os.system("cp ./CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin ./CaMa_in/restart/assim/restart"+yyyy+mm+dd+"A"+numch+".bin")
        copy_stoonly("./CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin","./CaMa_in/restart/open/restart"+yyyy+mm+dd+"C"+numch+".bin")
        copy_stoonly("./CaMa_out/"+spinup_open+"/restart"+yyyy+mm+dd+".bin","./CaMa_in/restart/assim/restart"+yyyy+mm+dd+"A"+numch+".bin")
###########################
def copy_out(loop): #old
    os.system("./src/copy_out "+loop)
    return 0
###########################
def reset_loop(yyyy): #old 
    os.system("rm -f ./assim_out/nonassim/*")
    os.system("rm -fR ./CaMa_out/global_15min_"+yyyy+"*")
###########################
def make_rand(std): #used
    # prepare random numbers for ensemable simulation
    randlist=rd.normal(0,std,pm.ens_mem())
    randlist=randlist.astype(np.float32)
    randlist.tofile("./CaMa_out/randlist.bin")
###########################
def initial(): #used
    # program for initialization

    # creating output folders
    mkdir("CaMa_out")
    mkdir("CaMa_in")
    mkdir("CaMa_in/restart")
    mkdir("CaMa_in/restart/assim")
    mkdir("CaMa_in/restart/open")
    mkdir("CaMa_in/restart/true")
    mkdir("assim_out")
    mkdir("assim_out/xa_m")
    mkdir("assim_out/xa_m/assim")
    mkdir("assim_out/xa_m/open")
    mkdir("assim_out/xa_m/true")

    mkdir("assim_out/ens_xa")         # NEW v.1.1.0
    mkdir("assim_out/ens_xa/assim")   # NEW v.1.1.0
    mkdir("assim_out/ens_xa/open")    # NEW v.1.1.0

    mkdir("assim_out/nonassim")
    mkdir("assim_out/nonassim/open")
    mkdir("assim_out/nonassim/assim")

    mkdir("assim_out/rest_true")
    mkdir("assim_out/rivout")
    mkdir("assim_out/rivout/open")
    mkdir("assim_out/rivout/assim")
    mkdir("assim_out/rivout/true")
    mkdir("assim_out/fldout")
    mkdir("assim_out/fldout/open")
    mkdir("assim_out/fldout/assim")
    mkdir("assim_out/fldout/true")
    mkdir("assim_out/flddph")
    mkdir("assim_out/flddph/open")
    mkdir("assim_out/flddph/assim")
    mkdir("assim_out/flddph/true")
    mkdir("err_out")
    mkdir("assim_out/fldarea/")
    mkdir("assim_out/fldarea/open")
    mkdir("assim_out/fldarea/assim")
    mkdir("assim_out/fldarea/true")
    # inflation parameter
    mkdir("inflation")

    os.system("touch assim_out/__init__.py")
    mkdir("logout")

    #mkdir("assim_out/rivhgt")
    #mkdir("assim_out/rivdph")
    #mkdir("assim_out/rivhgt/assim")
    #mkdir("assim_out/rivdph/assim")
    #mkdir("assim_out/rivdph/open")
    #mkdir("assim_out/rivdph/true")

    # link input files
    #os.system("ln -s "+pm.CaMa_dir()+"/inp/ELSE_GPCC ./CaMa_in")
    #slink(pm.CaMa_dir()+"/inp/ELSE_GPCC","./CaMa_in/ELSE_GPCC")

    #mkdir("./CaMa_in/ELSE_GPCC/mean_month")
    #mkdir("./CaMa_in/ELSE_KIM2009/mean_month")

    # make ./CaMa_in/ELSe_GPCC/* runoff files
    #mkdir("./CaMa_in/ELSE_GPCC/Roff_TRUE")
    #mkdir("./CaMa_in/ELSE_GPCC/Roff_CORR")
    #mkdir("./CaMa_in/ELSE_KIM2009/Roff_TRUE")
    #mkdir("./CaMa_in/ELSE_KIM2009/Roff_CORR")

    
    return 0
###########################
def compile_func(): #used
    # program for compiling
    # activate ifort
    #os.system("source /opt/intel/parallel_studio_xe_2017/psxevars.sh intel64")
    os.system("ifort src/make_nonassim_init.f90 -o src/make_nonassim_init -O2 -assume byterecl")
    os.system("ifort src/make_nonassim.f90 -o src/make_nonassim -O2 -assume byterecl")
    os.system("ifort src/copy_out.f90 -o src/copy_out -O2 -assume byterecl")
    os.system("ifort src/make_restart.f90 -o src/make_restart -O2 -assume byterecl -heap-arrays -nogen-interfaces -free -g -traceback  -lpthread -parallel")
    #os.system("ifort src/calc_stoerr.f90 -o src/calc_stoerr -O2 -assume byterecl")
    os.system("ifort src/make_rivout.f90 -o src/make_rivout -O2 -assume byterecl") 
    print "compile data assimilation codes..."
#    os.system("source src/compileMKL.sh "+pm.MKLdir())
    os.system("ifort  src/make_corrupt_rivhgt.f90 -o src/make_corrupt_rivhgt -O3 -assume byterecl -heap-arrays -nogen-interfaces -free -mkl -mcmodel=large -shared-intel")
#    os.system("ifort  src/data_assim.f90 -o src/data_assim -O3 -assume byterecl -heap-arrays 10 -nogen-interfaces -free -mkl -check bounds -g -fp-stack-check -g -traceback -lpthread -openmp")
    os.system("ifort  src/data_assim.f90 -o src/data_assim -O3 -assume byterecl -heap-arrays -nogen-interfaces -free -mkl=parallel -g -traceback  -lpthread -parallel") #-openmp
#    os.system("ifort  src/data_assim_bathy_fld.f90 -o src/data_assim_fld -O3 -assume byterecl -heap-arrays -nogen-interfaces -free -mkl -traceback -qopenmp")
    os.system("ifort  src/make_covariance.f90 -o src/make_covariance -O3 -assume byterecl -heap-arrays -nogen-interfaces -free -g -traceback -check bounds")
    if pm.patch_size()==0:
         os.system("ifort  src/data_assim_0.f90 -o src/data_assim -O3 -assume byterecl -heap-arrays -nogen-interfaces -free -mkl=parallel -g -traceback  -lpthread -parallel")
    return 0
###########################
def store_out(yyyy,mm,dd):
    # program for storing data #
    
    looptype = "true"
    # storing rivout
    numch = "000" 
    shutil.copy("./CaMa_out/"+yyyy+mm+dd+"T"+numch+"/rivout"+yyyy+".bin","assim_out/rivout/"+looptype+"/rivout"+yyyy+mm+dd+".bin")

    # storing rivdph
    #shutil.copy("./CaMa_out/"+yyyy+mm+dd+"T"+numch+"/rivdph"+yyyy+".bin","assim_out/rivdph/"+looptype+"/rivdph"+yyyy+mm+dd+".bin")

    # storing fldout
    shutil.copy("./CaMa_out/"+yyyy+mm+dd+"T"+numch+"/fldout"+yyyy+".bin","assim_out/fldout/"+looptype+"/fldout"+yyyy+mm+dd+".bin")

    # storing flddph
    shutil.copy("./CaMa_out/"+yyyy+mm+dd+"T"+numch+"/flddph"+yyyy+".bin","assim_out/flddph/"+looptype+"/flddph"+yyyy+mm+dd+".bin")

    # storing fldarea
    shutil.copy("./CaMa_out/"+yyyy+mm+dd+"T"+numch+"/fldare"+yyyy+".bin","assim_out/fldarea/"+looptype+"/fldarea"+yyyy+mm+dd+".bin")


    for CA in ["C","A"]:
        if CA == "C":
            looptype = "open"
        if CA == "A":
            looptype = "assim"

#        if CA == "C": 
#        # storing rivout
#            for num in np.arange(1,pm.ens_mem()+1):
#                numch = '%03d' % num 
#                shutil.copy("./CaMa_out/"+yyyy+mm+dd+CA+numch+"/rivdph"+yyyy+".bin","assim_out/rivdph/"+looptype+"/rivdph"+yyyy+mm+dd+"_"+numch+".bin")

        # storing rivout
        for num in np.arange(1,pm.ens_mem()+1):
            numch = '%03d' % num 
            shutil.copy("./CaMa_out/"+yyyy+mm+dd+CA+numch+"/rivout"+yyyy+".bin","assim_out/rivout/"+looptype+"/rivout"+yyyy+mm+dd+"_"+numch+".bin")

        # storing fldout
        for num in np.arange(1,pm.ens_mem()+1):
            numch = '%03d' % num 
            shutil.copy("./CaMa_out/"+yyyy+mm+dd+CA+numch+"/fldout"+yyyy+".bin","assim_out/fldout/"+looptype+"/fldout"+yyyy+mm+dd+"_"+numch+".bin")

        # storing flddph
        for num in np.arange(1,pm.ens_mem()+1):
            numch = '%03d' % num 
            shutil.copy("./CaMa_out/"+yyyy+mm+dd+CA+numch+"/flddph"+yyyy+".bin","assim_out/flddph/"+looptype+"/flddph"+yyyy+mm+dd+"_"+numch+".bin")

        # storing fldarea
        for num in np.arange(1,pm.ens_mem()+1):
            numch = '%03d' % num 
            shutil.copy("./CaMa_out/"+yyyy+mm+dd+CA+numch+"/fldare"+yyyy+".bin","assim_out/fldarea/"+looptype+"/fldarea"+yyyy+mm+dd+"_"+numch+".bin")

    return 0
###########################    
def SWOT_day(yyyy,mm,dd):
    st_year,st_month,st_date=pm.starttime()
    start_time=datetime.date(st_year,st_month,st_date)
    this_time=datetime.date(int(yyyy),int(mm),int(dd))
    days=this_time-start_time
    days=days.days
    return days%21+1
###########################
def make_restart(inputlist):
    # new version
    # sfcelv >> restart
    yyyy=inputlist[0]
    mm=inputlist[1]
    dd=inputlist[2]
    loop=inputlist[3]
    numch=inputlist[4]

    # make restart file
    # all in 4 bytes float (small endian) 1440*720
    # items are placed in the following order:
    #       river storage, floodplain storage, river outflow
    #       floodplain outflow, river depth, flood plain storage

    # built in hold
    print "finish assimilating"
    print "built in hold"
    print "press enter"

    # get the date of one day before
    bef_y=odb.calc_odb(yyyy,mm,dd,"year")
    bef_m=odb.calc_odb(yyyy,mm,dd,"month")
    bef_d=odb.calc_odb(yyyy,mm,dd,"date")

    yyyy_b="%04d"%bef_y
    mm_b="%02d"%bef_m
    dd_b="%02d"%bef_d

    # get the date of one day after
    nowdate=dt.datetime(int(yyyy),int(mm),int(dd))
    nextdate=nowdate+dt.timedelta(days=1)
    yyyy_n="%04d"%nextdate.year
    mm_n="%02d"%nextdate.month
    dd_n="%02d"%nextdate.day
    
    # calculate other variables from water storage
    dir1=pm.CaMa_dir()+"/"
    print "dir1",dir1
    os.system("./src/make_restart "+yyyy+mm+dd+" "+yyyy_b+mm_b+dd_b+" "+yyyy_n+mm_n+dd_n+" "+loop+" "+dir1+" "+str(pm.ens_mem())+" "+numch)

    print "finish restarting",numch
###########################
def make_rivout(inputlist):
    # new version
    # sfcelv >> rivout
    yyyy=inputlist[0]
    mm=inputlist[1]
    dd=inputlist[2]
    loop=inputlist[3]
    numch=inputlist[4]

    #print "calculate insantaneous discharge using Manning's equation"
   
    # calculate other variables from water storage
    dir1=pm.CaMa_dir()+"/"
    print "dir1",dir1
    os.system("./src/make_rivout "+yyyy+mm+dd+" "+loop+" "+dir1+" "+str(pm.ens_mem())+" "+numch)

###########################
def prepare_input_old(start_year,start_month,start_date,end_year,end_month,end_date):
    #if(len(glob.glob("./CaMa_in/ELSE_GPCC/Roff_TRUE/Roff*"))==0):
    if(len(glob.glob("./CaMa_in/ELSE_KIM2009/Roff_TRUE/Roff*"))==0):
        # true input file is not ready
        # need preparation
        #for item in glob.glob("./CaMa_in/ELSE_GPCC/Roff/Roff*"):
        #    shutil.copy(item,"./CaMa_in/ELSE_GPCC/Roff_TRUE/")
        for item in glob.glob("./CaMa_in/ELSE_KIM2009/Roff/Roff____2008*"):
            shutil.copy(item,"./CaMa_in/ELSE_KIM2009/Roff_TRUE/")
        for item in glob.glob("./CaMa_in/ELSE_KIM2009/Roff/Roff____2007*"):
            shutil.copy(item,"./CaMa_in/ELSE_KIM2009/Roff_TRUE/")



    # corrupt input file need preparation
    #os.system("rm -Rf ./CaMa_in/ELSE_GPCC/Roff_CORR/Roff*")
    #for item in glob.glob("./CaMa_in/ELSE_GPCC/Roff/Roff*"):
    #    shutil.copy(item,"./CaMa_in/ELSE_GPCC/Roff_CORR/")

    if pm.mode()==1:
      os.system("rm -Rf ./CaMa_in/ELSE_KIM2009/Roff_CORR/Roff*")
      for item in glob.glob("./CaMa_in/ELSE_KIM2009/Roff/Roff____2008*"):
          shutil.copy(item,"./CaMa_in/ELSE_KIM2009/Roff_CORR/")
      for item in glob.glob("./CaMa_in/ELSE_KIM2009/Roff/Roff____2007*"):
          shutil.copy(item,"./CaMa_in/ELSE_KIM2009/Roff_CORR/")


    start_dt=datetime.date(start_year,start_month,start_date)
    end_dt=datetime.date(end_year,end_month,end_date)
     
    if pm.mode()==2:
        # one year before experiment
        for day in np.arange((end_dt-start_dt).days):
            running_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (running_dt.year)
            mm='%02d' % (running_dt.month)
            dd='%02d' % (running_dt.day)
            #os.system("cp -f ./CaMa_in/ELSE_GPCC/Roff_CORR/Roff____"+"1990"+str(mm)+str(dd)+".one "+"./CaMa_in/ELSE_GPCC/Roff_CORR/Roff____tmp"+str(mm)+str(dd)+".one")
            #os.system("cp -f ./CaMa_in/ELSE_GPCC/Roff_CORR/Roff____"+"1991"+str(mm)+str(dd)+".one "+"./CaMa_in/ELSE_GPCC/Roff_CORR/Roff____1990"+str(mm)+str(dd)+".one")
            #os.system("cp -f ./CaMa_in/ELSE_GPCC/Roff_CORR/Roff____"+"tmp"+str(mm)+str(dd)+".one "+"./CaMa_in/ELSE_GPCC/Roff_CORR/Roff____1991"+str(mm)+str(dd)+".one")

            #os.system("cp -f ./CaMa_in/ELSE_KIM2009/Roff_CORR/Roff____"+"2007"+str(mm)+str(dd)+".one "+"./CaMa_in/ELSE_KIM2009/Roff_CORR/Roff____tmp"+str(mm)+str(dd)+".one")
            #os.system("cp -f ./CaMa_in/ELSE_KIM2009/Roff_CORR/Roff____"+"2008"+str(mm)+str(dd)+".one "+"./CaMa_in/ELSE_KIM2009/Roff_CORR/Roff____2007"+str(mm)+str(dd)+".one")
            #os.system("cp -f ./CaMa_in/ELSE_KIM2009/Roff_CORR/Roff____"+"tmp"+str(mm)+str(dd)+".one "+"./CaMa_in/ELSE_KIM2009/Roff_CORR/Roff____2008"+str(mm)+str(dd)+".one")
            #if mm=="02" and dd=="29":
            #  os.system("cp -f ./CaMa_in/ELSE_KIM2009/Roff/Roff____"+"20070228"+".one "+"./CaMa_in/ELSE_KIM2009/Roff_CORR/Roff____2008"+str(mm)+str(dd)+".one")
            #else:
            os.system("cp -f ./CaMa_in/ELSE_KIM2009/Roff/Roff____"+"2004"+str(mm)+str(dd)+".one "+"./CaMa_in/ELSE_KIM2009/Roff_CORR/Roff____2008"+str(mm)+str(dd)+".one")
            if mm=="02" and dd=="29":
              continue 
            os.system("cp -f ./CaMa_in/ELSE_KIM2009/Roff/Roff____"+"2003"+str(mm)+str(dd)+".one "+"./CaMa_in/ELSE_KIM2009/Roff_CORR/Roff____2007"+str(mm)+str(dd)+".one")

    if pm.mode()==3:
        # H08 runoff  
        os.system("rm -Rf ./CaMa_in/ELSE_KIM2009/Roff_CORR/Roff*")
        for item in glob.glob("./CaMa_in/ELSE_H08/Roff/Roff____2008*"):
            shutil.copy(item,"./CaMa_in/ELSE_KIM2009/Roff_CORR/")
        for item in glob.glob("./CaMa_in/ELSE_H08/Roff/Roff____2007*"):
            shutil.copy(item,"./CaMa_in/ELSE_KIM2009/Roff_CORR/")
 
    # calc monthly mean value
    #os.system("rm -Rf ./CaMa_in/ELSE_GPCC/mean_month/*")
    os.system("rm -Rf ./CaMa_in/ELSE_KIM2009/mean_month/*")
    threshold=0.1

    start_dt=datetime.date(start_year,start_month,start_date)
    for month in np.arange(12):
        roff_mon=np.zeros(360*180).reshape([180,360])
        count=np.zeros(360*180).reshape([180,360])
        for day in np.arange(month*30,(month+1)*30):
            day_num=day-month*30
            running_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (running_dt.year)
            mm='%02d' % (running_dt.month)
            dd='%02d' % (running_dt.day)

            #roff=np.fromfile("./CaMa_in/ELSE_GPCC/Roff_CORR/Roff____"+str(yyyy)+str(mm)+str(dd)+".one",np.float32).reshape([180,360])
            roff=np.fromfile("./CaMa_in/ELSE_KIM2009/Roff_CORR/Roff____"+str(yyyy)+str(mm)+str(dd)+".one",np.float32).reshape([180,360]) 
            roff_mon=roff_mon+roff*(roff>threshold)
            count=count+(roff>threshold)
        roff_mean=roff_mon/(count+1e-20)
        roff_mean=roff_mean.astype(np.float32)
        roff_mean=roff_mean+threshold
        #roff_mean.tofile("./CaMa_in/ELSE_GPCC/mean_month/mean_"+"%02d"%(month+1)+".bin")
        roff_mean.tofile("./CaMa_in/ELSE_KIM2009/mean_month/mean_"+"%02d"%(month+1)+".bin")
###########################
def cal_monthly_mean(start_year,end_year,months=24):
    # calc monthly mean value for two years
    runname=pm.runname(pm.mode())
    if runname=="E2O":
        nx=1440
        ny=720
    else:
        nx=360
        ny=180
    #os.system("rm -Rf ./CaMa_in/ELSE_GPCC/mean_month/*")
    mkdir("./CaMa_in/"+runname+"/mean_month")
    #os.system("rm -Rf ./CaMa_in/ELSE_KIM2009/mean_month/*")
    threshold=0.1

    start_dt=datetime.date(start_year,1,1)
    end_dt=datetime.date(end_year,12,31)
    for month in np.arange(months):
        ynow=int(start_year+int(month/12))
        ychar="%04d"%(ynow)
        mchar="%02d"%((month%12)+1)
        #print ychar, mchar
        roff_mon=np.zeros([ny,nx],np.float32)#.reshape([180,360])
        count=np.zeros([ny,nx],np.float32)#.reshape([180,360])
        for day in np.arange(month*30,(month+1)*30):
            day_num=day-month*30
            running_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (running_dt.year)
            mm='%02d' % (running_dt.month)
            dd='%02d' % (running_dt.day)
            roff=np.fromfile("./CaMa_in/"+runname+"/Roff/Roff____"+str(yyyy)+str(mm)+str(dd)+".one",np.float32).reshape([ny,nx]) 
            roff_mon=roff_mon+roff*(roff>threshold)
            count=count+(roff>threshold)
        roff_mean=roff_mon/(count+1e-20)
        roff_mean=roff_mean.astype(np.float32)
        roff_mean=roff_mean+threshold
        roff_mean.tofile("./CaMa_in/"+runname+"/mean_month/mean_"+ychar+mchar+".bin")
###########################
def cal_monthly_mean_ens(ens_num): #start_year,end_year,ens_num,months=24):
    # calc monthly mean value for two years
    start_year=pm.start_year() #inputlist[0]
    end_year=pm.end_year() #inputlist[1]
    #ens_num=inputlist[2]
    months=24 #inputlist[3]
    #threshold=inputlist[4]
    runname=pm.runname(pm.mode())
    if runname=="E2O":
        nx=1440
        ny=720
        threshold=0.1
        ens_mem=7
    else:
        nx=360
        ny=180
        threshold=1.0e-8
        ens_mem=10
    #os.system("rm -Rf ./CaMa_in/ELSE_GPCC/mean_month/*")
    mkdir("./CaMa_in/"+runname+"/mean_month")
    #os.system("rm -Rf ./CaMa_in/ELSE_KIM2009/mean_month/*")
    #threshold=0.1

    start_dt=datetime.date(start_year,1,1)
    end_dt=datetime.date(end_year,12,31)
    ens_char="%03d"%(ens_num)
    for month in np.arange(months):
        ynow=int(start_year+int(month/12))
        ychar="%04d"%(ynow)
        mchar="%02d"%((month%12)+1)
        #print ychar, mchar
        roff_mon=np.zeros([ny,nx],np.float32)#.reshape([180,360])
        count=np.zeros([ny,nx],np.float32)#.reshape([180,360])
        for day in np.arange(month*30,(month+1)*30):
            day_num=day-month*30
            running_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (running_dt.year)
            mm='%02d' % (running_dt.month)
            dd='%02d' % (running_dt.day)
            roff=np.fromfile("./CaMa_in/"+runname+"/Roff/Roff__"+str(yyyy)+str(mm)+str(dd)+ens_char+".one",np.float32).reshape([ny,nx])
            roff_mon=roff_mon+roff*(roff>threshold)
            count=count+(roff>threshold)
        roff_mean=roff_mon/(count+1e-20)
        roff_mean=roff_mean.astype(np.float32)
        roff_mean=roff_mean+threshold
        roff_mean.tofile("./CaMa_in/"+runname+"/mean_month/mean_"+ychar+mchar+ens_char+".bin")
###########################
def cal_monthly_total(start_year,end_year,months=24,threshold=0.1):
    # calc monthly mean value for two years
    runname=pm.runname(pm.mode())
    mkdir("./CaMa_in/"+runname+"/total_month")
    if runname=="E2O":
        nx=1440
        ny=720
    else:
        nx=360
        ny=180
    #os.system("rm -Rf ./CaMa_in/ELSE_GPCC/mean_month/*")
    #os.system("rm -Rf ./CaMa_in/"+runname+"/total_month/*")
    #threshold=0.1
    start_dt=datetime.date(start_year,1,1)
    end_dt=datetime.date(end_year,12,31)
    for month in np.arange(months):
        ynow=int(start_year+int(month/12))
        ychar="%04d"%(ynow)
        mchar="%02d"%((month%12)+1)
        #print ychar, mchar
        roff_mon=np.zeros([ny,nx],np.float32)#.reshape([180,360])
        #count=np.zeros([ny,nx],np.float32)#.reshape([180,360])
        for day in np.arange(month*30,(month+1)*30):
            day_num=day-month*30
            running_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (running_dt.year)
            mm='%02d' % (running_dt.month)
            dd='%02d' % (running_dt.day)
            roff=np.fromfile("./CaMa_in/"+runname+"/Roff/Roff____"+str(yyyy)+str(mm)+str(dd)+".one",np.float32).reshape([ny,nx]) 
            roff_mon=roff_mon+roff*(roff>threshold)
            #count=count+(roff>threshold)
        roff_total=roff_mon #/(count+1e-20)
        roff_total=roff_total.astype(np.float32)
        roff_total=roff_total+threshold
        roff_total.tofile("./CaMa_in/"+runname+"/total_month/total_"+ychar+mchar+".bin")
###########################
def cal_ens_mean():
    # calc ensemble mean value for two years
    start_year=pm.start_year() #inputlist[0]
    end_year=pm.end_year() #inputlist[1]
    #ens_num=inputlist[2]
    months=24 #inputlist[3]
    #threshold=inputlist[4]
    runname=pm.runname(pm.mode())
    if runname=="E2O":
        nx=1440
        ny=720
        threshold=0.1
    else:
        nx=360
        ny=180
        threshold=1.0e-8
    #os.system("rm -Rf ./CaMa_in/ELSE_GPCC/mean_month/*")
    mkdir("./CaMa_in/"+runname+"/mean_month")
    #os.system("rm -Rf ./CaMa_in/ELSE_KIM2009/mean_month/*")
    #threshold=0.1

    start_dt=datetime.date(start_year,1,1)
    end_dt=datetime.date(end_year,12,31)
    ens_char="%03d"%(ens_num)
    for month in np.arange(months):
        ynow=int(start_year+int(month/12))
        ychar="%04d"%(ynow)
        mchar="%02d"%((month%12)+1)
        #print ychar, mchar
        roff_mon=np.zeros([ny,nx],np.float32)#.reshape([180,360])
        count=np.zeros([ny,nx],np.float32)#.reshape([180,360])

###########################
def prepare_input():
    # spinup start_year
    # simulation end_year
    start_year=pm.start_year()
    end_year=pm.end_year()
    start_dt=datetime.date(start_year,1,1)
    last_dt=datetime.date(end_year,12,31)
    start=0
    last=int((last_dt-start_dt).days)+1
    #--------------
    # E2O
    if pm.mode()==1: # Earth2Observe
        distopen=0.75 #pm.distopen()
        diststd=0.25  #pm.diststd()
        # copy for TRUE simulation
        if(len(glob.glob("./CaMa_in/E2O/Roff_TRUE/Roff*"))!=0):
            pass
        # true input file is not ready
        # need preparation
        # make directories
        mkdir("./CaMa_in/E2O/Roff_TRUE")

        inputlist=[]
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            ens_char="T000"
            iname="./CaMa_in/E2O/Roff/Roff__"+yyyy+mm+dd+"003.one" # ECMWF
            oname="./CaMa_in/E2O/Roff_TRUE/Roff__"+yyyy+mm+dd+ens_char+".one"
            inputlist.append([iname,oname,"1.00"])

        # do parallel
        p=Pool(pm.para_nums())
        p.map(copy_runoff,inputlist)
        p.terminate()

        # calculate monthly mean
        # do parallel
        p=Pool(pm.para_nums())
        p.map(cal_monthly_mean_ens,np.arange(1,7+1))
        p.terminate()

        # calculate ensemble mean
        # do parallel
        p=Pool(pm.para_nums())
        p.map(cal_ens_mean,np.arange(1,7+1))
        p.terminate()


        # calculate monthly total

        # make courrpted runoff
        #if(len(glob.glob("./CaMa_in/E2O/Roff_CORR/Roff*"))!=0):
        #    pass
        # corrupted input file is not ready
        # need preparation
        # make directories
        mkdir("./CaMa_in/E2O/Roff_CORR")

        inputlist=[]
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            ens_num=1
            for runens in np.arange(1,7+1):
                run_num="%03d"%(runens)
                iname="./CaMa_in/E2O/Roff/Roff__"+yyyy+mm+dd+run_num+".one"
                ifile=np.fromfile("./CaMa_in/E2O/Roff/Roff__"+yyyy+mm+dd+run_num+".one",np.float32).reshape(720,1440)
                roff_mean=np.fromfile("./CaMa_in/E2O/mean_month/mean_"+yyyy+mm+run_num+".bin",np.float32).reshape(180,360)
                #roff_total=np.fromfile("./CaMa_in/E2O/total_month/total_"+yyyy+mm+".bin",np.float32).reshape(180,360)
                distopen_range=rd.normal(0,diststd,3)
                if runens == 3:
                    distopen_range=rd.normal(0,diststd,2)
                distopen_range=distopen_range.astype(np.float32)
                #distopen_range=[0.75,1.00,1.25]
                for dist in distopen_range:
                    #if runens == 1 and ens_num==2: #distopen == 1.00:
                    #    pass
                    ens_char="C%03d"%(ens_num)
                    ofile=ifile + roff_mean*dist
                    ofile.tofile("./CaMa_in/E2O/Roff_CORR/Roff__"+yyyy+mm+dd+ens_char+".one")
                    #oname="./CaMa_in/E2O/Roff_CORR/Roff__"+yyyy+mm+dd+ens_char+".one"
                    #inputlist.append([iname,oname,str(1.0 - dist)])
                    ens_num=ens_num+1

        # do parallel
        #p=Pool(pm.para_nums())
        #p.map(copy_runoff,inputlist)
        #p.terminate()
    #---------
    # ERA20CM
    if pm.mode()==2: #ECMWF ERA20CM
        distopen=0.75 #pm.distopen()
        diststd=0.5 #0.25  #pm.diststd()
        # copy for TRUE simulation
        if(len(glob.glob("./CaMa_in/ERA20CM/Roff_TRUE/Roff*"))!=0):
            pass
        # true input file is not ready
        # need preparation
        # make directory
        mkdir("./CaMa_in/ERA20CM/Roff_TRUE")

        inputlist=[]
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            ens_char="T000"
            iname="./CaMa_in/ERA20CM/Roff/Roff__"+yyyy+mm+dd+"001.one"
            oname="./CaMa_in/ERA20CM/Roff_TRUE/Roff__"+yyyy+mm+dd+ens_char+".one"
            inputlist.append([iname,oname,"1.00"])

        # do parallel
        p=Pool(pm.para_nums())
        p.map(copy_runoff,inputlist)
        p.terminate()

        # calculate monthly mean
        # do parallel
        p=Pool(pm.para_nums())
        p.map(cal_monthly_mean_ens,np.arange(1,10+1))
        p.terminate()

        # calculate monthly total

        # make courrpted runoff
        #if(len(glob.glob("./CaMa_in/ERA20CM/Roff_CORR/Roff*"))!=0):
        #    pass
        # corrupted input file is not ready
        # need preparation
        # make directories
        mkdir("./CaMa_in/ERA20CM/Roff_CORR")

        #
        inputlist=[]
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            ens_num=1
            for runens in np.arange(1,10+1):
                run_num="%03d"%(runens)
                #iname="./CaMa_in/ERA20CM/Roff/Roff__"+yyyy+mm+dd+run_num+".one"
                ifile=np.fromfile("./CaMa_in/ERA20CM/Roff/Roff__"+yyyy+mm+dd+run_num+".one",np.float32).reshape(180,360)
                roff_mean=np.fromfile("./CaMa_in/ERA20CM/mean_month/mean_"+yyyy+mm+run_num+".bin",np.float32).reshape(180,360)
                #roff_total=np.fromfile("./CaMa_in/ERA20CM/total_month/total_"+yyyy+mm+".bin",np.float32).reshape(180,360)
                distopen_range=rd.normal(0,diststd,2)
                distopen_range=distopen_range.astype(np.float32)
                #distopen_range=[0.75,1.25]
                for dist in distopen_range:
                    ens_char="C%03d"%(ens_num)
                    ofile=ifile + roff_mean*dist
                    ofile.tofile("./CaMa_in/ERA20CM/Roff_CORR/Roff__"+yyyy+mm+dd+ens_char+".one")
                    #oname="./CaMa_in/ERA20CM/Roff_CORR/Roff__"+yyyy+mm+dd+ens_char+".one"
                    #inputlist.append([iname,oname,str(distopen)])
                    ens_num=ens_num+1

        # do parallel
        #p=Pool(pm.para_nums())
        #p.map(copy_runoff,inputlist)
        #p.terminate()
    #--------------
    # -25% biased runoff experiment
    if pm.mode()==3: # Earth2Observe
        distopen=0.75 #0.75 #pm.distopen()
        diststd=0.25  #pm.diststd()
        # copy for TRUE simulation
        if(len(glob.glob("./CaMa_in/ELSE_KIM2009/Roff_TRUE/Roff*"))!=0):
            pass
        # true input file is not ready
        # need preparation
        # make directories
        mkdir("./CaMa_in/ELSE_KIM2009/Roff_TRUE")

        inputlist=[]
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            ens_char="T000"
            iname="./CaMa_in/ELSE_KIM2009/Roff/Roff____"+yyyy+mm+dd+".one"
            oname="./CaMa_in/ELSE_KIM2009/Roff_TRUE/Roff__"+yyyy+mm+dd+ens_char+".one"
            inputlist.append([iname,oname,"1.00"])

        # do parallel
        p=Pool(pm.para_nums())
        p.map(copy_runoff,inputlist)
        p.terminate()

        # calculate mothly mean
        cal_monthly_mean(start_year,end_year,24)

        # calculate monthly total
        cal_monthly_total(start_year,end_year,24)

        # make courrpted runoff
        #if(len(glob.glob("./CaMa_in/ELSE_KIM2009/Roff_CORR/Roff*"))!=0):
        #    pass
        # corrupted input file is not ready
        # need preparation
        # make directories
        mkdir("./CaMa_in/ELSE_KIM2009/Roff_CORR")

        inputlist=[]
        #std=rd.normal(0,diststd,pm.ens_mem())
        #std=std.astype(np.float32)
        for day in np.arange(start,last):
            target_dt=start_dt+datetime.timedelta(days=day)
            yyyy='%04d' % (target_dt.year)
            mm='%02d' % (target_dt.month)
            dd='%02d' % (target_dt.day)
            # make random values
            #std=np.fromfile("./CaMa_out/randlist.bin",np.float32)
            std=rd.normal(0,diststd,pm.ens_mem())
            std=std.astype(np.float32)
            std=np.sort(std)
            #print std
            ifile=np.fromfile("./CaMa_in/ELSE_KIM2009/Roff/Roff____"+yyyy+mm+dd+".one",np.float32).reshape(180,360)
            roff_mean=np.fromfile("./CaMa_in/ELSE_KIM2009/mean_month/mean_"+yyyy+mm+".bin",np.float32).reshape(180,360)
            roff_total=np.fromfile("./CaMa_in/ELSE_KIM2009/total_month/total_"+yyyy+mm+".bin",np.float32).reshape(180,360)
            for ens in np.arange(1,pm.ens_mem()+1):
                ens_char="C%03d"%(ens)
                #print pm.distopen(),std[ens-1]
                ofile=ifile*distopen + roff_mean*std[ens-1]#*10.0
                #ofile=ifile*(distopen + std[ens-1])
                #ofile=ifile*distopen + roff_total*std[ens-1]
                #ofile=ifile*distopen + roff_total*std[ens-1]*0.75
                ofile.astype(np.float32)
                ofile.tofile("./CaMa_in/ELSE_KIM2009/Roff_CORR/Roff__"+yyyy+mm+dd+ens_char+".one")

        # do parallel
        #p=Pool(pm.para_nums())
        #p.map(copy_runoff,inputlist)
        #p.terminate()
    return 0
###########################
def copy_runoff(inputlist): #do it parallel
    iname=inputlist[0]
    oname=inputlist[1]
    distopen=float(inputlist[2])
    runoff=np.fromfile(iname,np.float32)*distopen
    runoff.tofile(oname)
    return 0
###########################
def make_corrupt_man_old():
    # prepare random numbers for ensemable simulation
    manrandlist=rd.normal(pm.corruptman_base(),pm.corruptman_std(),pm.ens_mem())
    manrandlist=manrandlist.astype(np.float32)

    f=open("./CaMa_out/manrandlist.txt","w")
    for i in np.arange(0,pm.ens_mem()):
        f.write(str(manrandlist[i]))
        f.write("\n")
    f.close()

    return 0
###########################
def make_corrpt_rivhgt():
    # prepare random numbers for ensemable simulation
    elerandlist=rd.normal(pm.corruptele_base(),pm.corruptele_std(),pm.ens_mem())
    elerandlist=elerandlist.astype(np.float32)

    f=open("./CaMa_out/elerandlist.txt","w")
    for i in np.arange(0,pm.ens_mem()):
        f.write(str(elerandlist[i]))
        f.write("\n")
    f.close()

    return 0
###########################
def intial_assim_rivhgt():
    print "initialize elevation"
    for ens in np.arange(1,pm.ens_mem()+1):
      #shutil.copy(pm.CaMa_dir()+"/map/global_15min/rivhgt_%03dC.bin"%(ens),pm.CaMa_dir()+"/map/global_15min/rivhgt19910101_%03dA.bin"%(ens))
      shutil.copy(pm.CaMa_dir()+"/map/global_15min/rivhgt_%03dC.bin"%(ens),pm.CaMa_dir()+"/map/global_15min/rivhgt_%03dA.bin"%(ens))

    return 0
###########################
def courrpt_rivhgt():
    print "make_corrupt_rivhgt.f90"
    os.system("./src/make_corrupt_rivhgt "+str('%4d'%(pm.spinup_end_year()))+" "+str(pm.non_hgt)+" "+str(pm.ens_mem())+" "+pm.CaMa_dir())
    return 0
###########################
def make_initial_infl():
    parm_infl=np.ones([720,1440],np.float32)*pm.initial_infl()
    start_year,start_month,start_date=pm.starttime() # Start year month date
    yyyy='%04d' % (start_year)
    mm='%02d' % (start_month)
    dd='%02d' % (start_date)
    parm_infl.tofile("./inflation/parm_infl"+yyyy+mm+dd+".bin")
###########################
def mkdir(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
###########################
def slink(src,dst):
    try:
        os.symlink(src,dst)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST:
            os.remove(dst)
            os.symlink(src,dst)
        else:
            raise
###########################
def multivariate_normal_sampler(mean,covariance,n_samples=1):
  # create multivarite samples 
  L=spla.cholesky(covariance)
  Z=np.random.normal(size=(n_samples,covariance.shape[0]))
  return Z.dot(L) + mean
###########################
def make_corrupt_man():
  # make multivariate normal distribution if corrupted manning
  rivnum = "data/rivnum.bin"
  rivnum = np.fromfile(rivnum,np.int32).reshape(720,1440)
#--
  corr_ens=np.ones([pm.manning_mem(),720,1440],np.float32)
  print np.max(rivnum),np.shape(corr_ens)
  for num in np.arange(1,np.max(rivnum)):
    print num
    index=np.where(rivnum==num)
    l=np.shape(index)[1]
    #--
    fname="temp.txt"
    f = open(fname,"w")
    #--
    for i in np.arange(0,l):
      #print index[0][i], index[1][i]
      line="%04d    %04d\n"%(index[1][i]+1, index[0][i]+1)
      f.write(line)
    f.close()
    # find covariances
    print "/src/make_covariance "+str(num)+" "+str(l)+" "+str(pm.corruptman_std())+" "+pm.CaMa_dir()+"/"
    os.system("./src/make_covariance "+str(num)+" "+str(l)+" "+str(pm.corruptman_std())+" "+pm.CaMa_dir()+"/")
    # multivariate covariance
    print "multivariate covariance"
    mean=np.ones([l],np.float32)*pm.corruptman_base()
    covariance=np.fromfile("cov.bin",np.float32).reshape(l,l)
    corr_ens[:,index[0],index[1]]=multivariate_normal_sampler(mean,covariance,pm.manning_mem())
    #corr_ens[:,index[0],index[1]]=np.random.multivariate_normal(mean,covariance,pm.manning_mem())
   #--
  print "make ensembles"
  for ens in np.arange(1,pm.manning_mem()+1):
    fname="CaMa_out/corruptman%03d.bin"%(ens)
    corr_ens[ens-1].tofile(fname) 
    fname=pm.CaMa_dir()+"/map/glb_15min/rivmanC%03d.bin"%(ens)
    corr_ens[ens-1].tofile(fname)
  return 0
###########################
def make_corrlated_man():
  # make spatially correalted manning
  rivnum = "data/rivnum.bin"
  rivnum = np.fromfile(rivnum,np.int32).reshape(720,1440)
#--
  corr_ens=np.ones([720,1440],np.float32)
  print np.max(rivnum),np.shape(corr_ens)
  for num in np.arange(1,np.max(rivnum)):
    print num
    index=np.where(rivnum==num)
    l=np.shape(index)[1]
    #--
    fname="temp.txt"
    f = open(fname,"w")
    #--
    for i in np.arange(0,l):
      #print index[0][i], index[1][i]
      line="%04d    %04d\n"%(index[1][i]+1, index[0][i]+1)
      f.write(line)
    f.close()
    # find covariances
    print "/src/make_covariance "+str(num)+" "+str(l)+" "+str(pm.corruptman_std())+" "+pm.CaMa_dir()+"/"
    os.system("./src/make_covariance "+str(num)+" "+str(l)+" "+str(pm.corruptman_std())+" "+pm.CaMa_dir()+"/")
    # calculate covariance
    print "calculate covariance"
    mean=np.ones([l],np.float32)*pm.corruptman_base()
    covariance=np.fromfile("cov.bin",np.float32).reshape(l,l)
    #--
    # Compute the Cholesky decomposition
    c=cholesky(covariance, lower=True)
    #--
    corr_ens[index[0],index[1]]=np.dot(c,mean)
    
  fname=pm.CaMa_dir()+"/map/glb_15min/rivmanTRUE.bin"
  corr_ens.tofile(fname)
  return 0
###########################
def make_corrupt_man_simple():
  # prepare random numbers for ensemable simulation
  manrandlist=rd.normal(pm.corruptman_base(),pm.corruptman_std(),pm.ens_mem())
  manrandlist=manrandlist.astype(np.float32)
  # ---
  for ens in np.arange(1,pm.manning_mem()+1):
     fname=pm.CaMa_dir()+"/map/glb_15min/rivmanC%03d.bin"%(ens)
     manrand=np.ones([720,1440],np.float32)*manrandlist[ens-1]
     manrand.tofile(fname)
  return 0
###########################





