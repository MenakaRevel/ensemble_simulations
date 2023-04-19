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
import calendar 
from multiprocessing import Pool
from multiprocessing import Process
from numpy import ma

os.system("ln -sf ../params.py params.py")
sys.path.append('../assim_out/')
import params as pm
import read_grdc as grdc
import cal_stat as stat
#from matplotlib.font_manager import FontProperties
#fp = FontProperties(fname="jap.ttc",size=15)

argvs = sys.argv

os.system("mkdir ../img")
os.system("mkdir ../img/runoff")
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
year=2004
month=1
date=1
start_dt=datetime.date(year,month,date)
size=60

start=0
last=366 #int(argvs[1])
if calendar.isleap(year):
    last=366
else:
    last=365
N=int(last)
green2="greenyellow"
#--------------
runname="ELSE_KIM2009"
#runname="E2O"
#runname="ERA20CM"
#DROFFUNIT=86400
DROFFUNIT=1.0
#--------------
nextxy = pm.CaMa_dir()+"/map/glb_15min/nextxy.bin"
rivwth = pm.CaMa_dir()+"/map/glb_15min/rivwth.bin"
rivhgt = pm.CaMa_dir()+"/map/glb_15min/rivhgt.bin"
rivlen = pm.CaMa_dir()+"/map/glb_15min/rivlen.bin"
elevtn = pm.CaMa_dir()+"/map/glb_15min/elevtn.bin"
nextxy = np.fromfile(nextxy,np.int32).reshape(2,720,1440)
rivwth = np.fromfile(rivwth,np.float32).reshape(720,1440)
rivhgt = np.fromfile(rivhgt,np.float32).reshape(720,1440)
rivlen = np.fromfile(rivlen,np.float32).reshape(720,1440)
elevtn = np.fromfile(elevtn,np.float32).reshape(720,1440)
#----

pname=[]
xlist=[]
ylist=[]
river=[]
#--
org=[]
opn=[]
for day in np.arange(start,last):
    target_dt=start_dt+datetime.timedelta(days=day)
    yyyy='%04d' % (target_dt.year)
    mm='%02d' % (target_dt.month)
    dd='%02d' % (target_dt.day)
    print yyyy,mm,dd

    # make org
    #fname="../assim_out/xa_m/true/"+yyyy+mm+dd+"_xam.bin"
    fname="../CaMa_in/"+runname+"/Roff_TRUE/Roff__"+yyyy+mm+dd+"T000.one"
    orgfile=np.fromfile(fname,np.float32)#.reshape([180,360])
    #print orgfile[0,0]
    #org.append(np.mean(ma.masked_where(orgfile==-9999.0,orgfile)))
    org.append(np.ma.mean(ma.masked_less_equal(orgfile,0.0)))

    opn_ens=[]
    for num in np.arange(1,int(pm.ens_mem())+1):
        numch='C%03d'%num

        #fname="../assim_out/ens_xa/open/"+yyyy+mm+dd+"_"+numch+"_xa.bin"
        #fname="../CaMa_out/"+yyyy+mm+dd+"C"+numch+"/sfcelv"+yyyy+".bin"
        fname="../CaMa_in/"+runname+"/Roff_CORR/Roff__"+yyyy+mm+dd+numch+".one"
        opnfile=np.fromfile(fname,np.float32)#.reshape([180,360])
        #opn_ens.append(np.mean(ma.masked_where(opnfile==-9999.0,opnfile)))
        opn_ens.append(np.ma.mean(ma.masked_less_equal(opnfile,0.0)))
    opn.append(opn_ens)
#
#print opn
org=np.array(org)#.filled(0.0))
opn=np.array(opn)#.filled(0.0))
# unit conversion
org=org*DROFFUNIT
opn=opn*DROFFUNIT
#--------
plt.close()
fig, ax1 = plt.subplots()
ax1.plot(np.arange(start,last),org,label="true",color="black",linewidth=0.7,zorder=103)
for num in np.arange(0,int(pm.ens_mem())):
    ax1.plot(np.arange(start,last),opn[:,num],label="corrupted",color="blue",linewidth=0.1,alpha=0.3,zorder=101)
ax1.plot(np.arange(start,last),np.mean(opn,axis=1),label="corrupted",color="blue",linewidth=0.7,alpha=0.5,zorder=102)
ens_spr=np.amax(opn,axis=1)-np.amin(opn,axis=1)
#print ens_spr
es_mean=np.nanmean(ens_spr)
mes="Ens Spr: %3.1f"%(es_mean)
ax1.text(0.02,0.9,mes,ha="left",va="center",transform=ax1.transAxes,fontsize=10)
ax2=ax1.twinx()
ax2.plot(np.arange(start,last),ens_spr,color="maroon",zorder=104,alpha=0.3)
ax2.set_ylabel('Ensemble Spread', color='maroon')
ax2.tick_params('y', colors='maroon')
ax1.set_ylabel("$Runoff (mm/day)$", color="k")
plt.savefig("../img/runoff/mean_runoff"+runname+".png",dpi=500)

