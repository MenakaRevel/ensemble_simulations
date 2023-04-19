#!/opt/local/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as cm
import sys
import os
import math
import calendar 
from multiprocessing import Pool
from multiprocessing import Process
from multiprocessing import sharedctypes
from numpy import ma
import netCDF4 as nc
import re

# os.system("ln -sf ../params.py params.py")
sys.path.append('../img/')
import params as pm
import read_grdc as grdc
import cal_stat as stat
#from matplotlib.font_manager import FontProperties
#fp = FontProperties(fname="jap.ttc",size=15)

# argvs = sys.argv
#=================
syear=2015
eyear=2020
#=================
#ERA5
ne=20 #pm.ens_mem()
output="bin"
runname="ERA5"
expname="CONUSERA5"
mapname="conus_06min"
CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
ncpus=20
# for runoff file
# ERA5
nXX=3600
nYY=1800
dRUnit=1e3
indir=".."

#=================
#E2O
# ne=49 #pm.ens_mem()
# output="bin"
# runname="E2O"
# expname="AMZE2O"
# mapname="amz_06min"
# CaMa_dir="/cluster/data6/menaka/CaMa-Flood_v4"
# ncpus=40
# # E2O
# nXX=1440
# nYY=720
# dRUnit=1.0
# indir="/cluster/data6/menaka/ensemble_simulations"
#--
# noise=np.sort(np.abs(np.random.normal(0,1,[ne])))
#---
os.system("mkdir ../figures")
os.system("mkdir ../figures/"+expname)
os.system("mkdir ../figures/"+expname+"/runoff")
#----
def SWOT_day(yyyy,mm,dd):
  st_year,st_month,st_date=pm.starttime()
  start_time=datetime.date(st_year,st_month,st_date)
  this_time=datetime.date(int(yyyy),int(mm),int(dd))
  days=this_time-start_time
  days=days.days
  return days%21+1
#----
def mk_dir(sdir):
  try:
    os.makedirs(sdir)
  except:
    pass
#----

#year=int(argvs[1])
month=1
date=1
start_dt=datetime.date(syear,1,1)
end_dt=datetime.date(eyear,12,31)
size=60
last=(end_dt-start_dt).days + 1
start=0
N=int(last)

green2="greenyellow"
#----
fname=CaMa_dir+"/map/"+mapname+"/params.txt"
with open(fname,"r") as f:
    lines=f.readlines()
#-------
nx     = int(filter(None, re.split(" ",lines[0]))[0])
ny     = int(filter(None, re.split(" ",lines[1]))[0])
gsize  = float(filter(None, re.split(" ",lines[3]))[0])
#--------------
nextxy = CaMa_dir+"/map/"+mapname+"/nextxy.bin"
rivwth = CaMa_dir+"/map/"+mapname+"/rivwth.bin"
rivhgt = CaMa_dir+"/map/"+mapname+"/rivhgt.bin"
rivlen = CaMa_dir+"/map/"+mapname+"/rivlen.bin"
elevtn = CaMa_dir+"/map/"+mapname+"/elevtn.bin"
nextxy = np.fromfile(nextxy,np.int32).reshape(2,ny,nx)
rivwth = np.fromfile(rivwth,np.float32).reshape(ny,nx)
rivhgt = np.fromfile(rivhgt,np.float32).reshape(ny,nx)
rivlen = np.fromfile(rivlen,np.float32).reshape(ny,nx)
elevtn = np.fromfile(elevtn,np.float32).reshape(ny,nx)
#----
# runname="ELSE_KIM2009"
#runname="E2O"
#runname="ERA20CM"
#runname=pm.runname(pm.mode())

#--
pname=[]
xlist=[]
ylist=[]
river=[]
#--
# rivernames = grdc.grdc_river_name()
rivernames = ["MISSOURI","MISSISSIPPI","COLORADO"]#,"AMAZON"] #
# rivernames = ["AMAZON"]
for rivername in rivernames: 
    grdc_id,station_loc,x_list,y_list = grdc.get_grdc_loc_v396(rivername,mapname="glb_06min")
    print (rivername, station_loc, x_list,y_list)
    river.append([rivername]*len(station_loc))
    pname.append(station_loc)
    xlist.append(x_list)
    ylist.append(y_list)

river=([flatten for inner in river for flatten in inner])
pname=([flatten for inner in pname for flatten in inner])
xlist=([flatten for inner in xlist for flatten in inner])
ylist=([flatten for inner in ylist for flatten in inner])

print (len(pname), len(xlist))

pnum=len(pname)


# multiprocessing array
sim=np.ctypeslib.as_ctypes(np.zeros([N,pnum,ne],np.float32))
shared_array_sim  = sharedctypes.RawArray(sim._type_, sim)

# for parallel calcualtion
inputlist=[]
for day in np.arange(start,last):
    target_dt=start_dt+datetime.timedelta(days=day)
    yyyy='%04d' % (target_dt.year)
    mm='%02d' % (target_dt.month)
    dd='%02d' % (target_dt.day)
    for num in np.arange(1,pm.ens_mem()+1):
        numch='%03d'%num
        ifile=indir+"/CaMa_in/"+runname+"/Roff/Roff__"+yyyy+mm+dd+numch+".one"
        inputlist.append([str(day),ifile,numch])
        # print (yyyy,ifile,numch)

#==============================
#=== function for read data ===
#==============================
def read_data(inputlist):
    day   = int(inputlist[0])
    ifile = inputlist[1]
    num   = int(inputlist[2]) - 1
    #--
    # print ("read data -> ", ifile)
    #--
    tmp_sim  = np.ctypeslib.as_array(shared_array_sim)

    # # year, mon, day
    # year=int(yyyy)
    
    # if calendar.isleap(year):
    #     dt=366
    # else:
    #     dt=365

    # # timings
    # target_dt=datetime.date(year,1,1)
    # st=(target_dt-start_dt).days
    # et=st+dt
    # if et >= N:
    #     et=None
    
    # # print ("read discharge")
    # # simulated discharge
    # if output == "bin":
    #     # fname=indir+"/outflw"+yyyy+".bin"
    #     # fname=indir+"/rivout"+yyyy+".bin"
    #     simfile=np.fromfile(ifile,np.float32).reshape([dt,ny,nx])
    # else:
    #     # fname=indir+"/o_outflw"+yyyy+".nc"
    #     # print ( fname )
    #     with nc.Dataset(ifile,"r") as cdf:
    #         simfile=cdf.variables["outflw"][:]

    simfile=np.fromfile(ifile,np.float32).reshape([nYY,nXX])
    print ("-- reading simulation file:", ifile )
    #-------------
    for point in np.arange(pnum):
        # print ("-- reading point:", pname[point] )
        ix1,iy1,ix2,iy2=grdc.get_grdc_station_v396(pname[point])
        if ix2 == -9999 or iy2 == -9999:
            tmp_sim[day,point,num]=simfile[iy1,ix1]
        else:
            tmp_sim[day,point,num]=simfile[iy1,ix1]+simfile[iy2,ix2]

#======================
#--read data parallel--
p=Pool(ncpus)
res = list(p.map(read_data, inputlist))
sim = np.ctypeslib.as_array(shared_array_sim)
p.terminate()

# for inpi in np.arange(inpn):
# res = map(read_data,inputlist)
# sim = np.ctypeslib.as_array(shared_array_sim)
#

# multiprocessing array for orginal data
org=np.ctypeslib.as_ctypes(np.zeros([N,pnum],np.float32))
shared_array_org  = sharedctypes.RawArray(org._type_, org)
#
#  for parallel calcualtion
inputlist=[]
for day in np.arange(start,last):
    target_dt=start_dt+datetime.timedelta(days=day)
    yyyy='%04d' % (target_dt.year)
    mm='%02d' % (target_dt.month)
    dd='%02d' % (target_dt.day)
    ifile="/work/a02/menaka/ERA5/bin/Roff____"+yyyy+mm+dd+".sixmin"
    # ifile=indir+"/CaMa_in/"+runname+"/Roff/Roff__"+yyyy+mm+dd+"007.one"
    # ifile="../CaMa_in/"+runname+"/Roff/Roff__"+yyyy+mm+dd+numch+".one"
    inputlist.append([str(day),ifile])

#==============================
# === function for org data === 
#==============================
def read_org_data(inputlist):
    day   = int(inputlist[0])
    ifile = inputlist[1]

    # print ("read data -> ", ifile)
    #--
    tmp_org  = np.ctypeslib.as_array(shared_array_org)

    # # year, mon, day
    # year=int(yyyy)
    
    # if calendar.isleap(year):
    #     dt=366
    # else:
    #     dt=365

    # # timings
    # target_dt=datetime.date(year,1,1)
    # st=(target_dt-start_dt).days
    # et=st+dt
    # if et >= N:
    #     et=None
    
    # # print ("read discharge")
    # # simulated discharge
    # if output == "bin":
    #     # fname=indir+"/outflw"+yyyy+".bin"
    #     # fname=indir+"/rivout"+yyyy+".bin"
    #     simfile=np.fromfile(ifile,np.float32).reshape([dt,ny,nx])
    # else:
    #     # fname=indir+"/o_outflw"+yyyy+".nc"
    #     # print ( fname )
    #     with nc.Dataset(ifile,"r") as cdf:
    #         simfile=cdf.variables["outflw"][:]

    simfile=np.fromfile(ifile,np.float32).reshape([nYY,nXX])
    print ("-- reading simulation file:", ifile )
    #-------------
    for point in np.arange(pnum):
        ix1,iy1,ix2,iy2=grdc.get_grdc_station_v396(pname[point])
        if ix2 == -9999 or iy2 == -9999:
            tmp_org[day,point]=simfile[iy1,ix1]
        else:
            tmp_org[day,point]=simfile[iy1,ix1]+simfile[iy2,ix2]

#======================
#--read data parallel--
p=Pool(ncpus)
res = list(p.map(read_org_data, inputlist))
org = np.ctypeslib.as_array(shared_array_org)
p.terminate()


def make_fig(point):
    plt.close()
    labels=[runname,"perturbed","mean"]
    colors=["#34495e","grey","grey"]
    # read observation data
    org1=org[:,point]
    # print (org)
    fig, ax = plt.subplots()
    lines=[ax.plot(np.arange(start,last),ma.masked_less(org1,0.0)*dRUnit,label="ERA",color="#34495e",linewidth=2.0,zorder=105)[0]] #,marker = "o",markevery=swt[point])
    for num in np.arange(0,ne):
        # print (num)
        lines.append(ax.plot(np.arange(start,last),sim[:,point,num]*dRUnit,label="perturbed",color="grey",linewidth=0.1,alpha=0.8,zorder=102)[0])
    lines.append(ax.plot(np.arange(start,last),np.mean(sim[:,point,:],axis=1)*dRUnit,label="mean",color="grey",linewidth=0.5,alpha=1.0,zorder=106)[0])
    # print (np.mean(sim[:,point,:],axis=1))
    # Make the y-axis label, ticks and tick labels match the line color.
    ax.set_ylabel('Runoff (mm/day)', color='k')
    ax.set_xlim(xmin=0,xmax=last+1)
    ax.tick_params('y', colors='k')
    # scentific notaion
    ax.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
    ax.yaxis.major.formatter._useMathText=True 
    # xlabels
    if eyear-syear > 5:
        dtt=2
        dt=int(math.ceil(((eyear-syear)+1)/2.0))
    elif eyear-syear > 10:
        dtt=5
        dt=int(math.ceil(((eyear-syear)+1)/5.0))
    else:
        dtt=1
        dt=(eyear-syear)+1
    xxlist=np.linspace(0,N,dt,endpoint=True)
    #xxlab=[calendar.month_name[i][:3] for i in range(1,13)]
    xxlab=np.arange(syear,eyear+1,dtt)
    ax.set_xticks(xxlist)
    ax.set_xticklabels(xxlab,fontsize=10)
    plt.legend(lines,labels,ncol=3,loc='upper center', bbox_to_anchor=(0.5, -0.05)) #,transform=ax1.transAxes)
    station_loc_list0="".join(pname[point].split())
    station_loc_list=station_loc_list0.split("/")
    station_name="-".join(station_loc_list) 
    print ('-- save  -> ',river[point] , station_name)
    figname0="".join(station_name.split())
    plt.savefig("../figures/"+expname+"/runoff/"+river[point]+"-"+station_name+".jpg",dpi=500)
    return 0

p=Pool(ncpus)
p.map(make_fig,np.arange(pnum))
p.terminate()

# map(make_fig,np.arange(pnum))