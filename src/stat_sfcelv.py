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
import errno
from numpy import ma 
import matplotlib.gridspec as gridspec
import string
from scipy.fftpack import fft, ifft, fftfreq
#from slacker import Slacker
from multiprocessing import Pool
from multiprocessing import Process
import xarray as xr

# sys.path.append("../")
# import params as pm
#--
#----
def mk_dir(sdir):
    try:
        os.makedirs(sdir)
    except:
        pass
#--
def stats(inputlist):
  input_name=inputlist[0]
  syear_in=int(inputlist[1])
  eyear_in=int(inputlist[2])
  syear_out=int(inputlist[3])
  eyear_out=int(inputlist[4])
  #mk_dir(pm.out_dir()+"/"+input_name)
  #--read netCDF4--
  # start_year=pm.start_year()
  # end_year=pm.end_year()
  tag="%04d-%04d"%(syear_in,eyear_in)
  tagout="%04d-%04d"%(syear_out,eyear_out)
  # sfcelv
  fname=outdir+"/CaMa_out/"+input_name+"/sfcelv"+tag+".nc"
  options = dict(scale = 1, xdim='lon', ydim='lat', tdim='time')
  tdim, xdim, ydim = options['tdim'], options['xdim'], options['ydim']
  with xr.open_dataset(fname) as nc:
    # nc = nc.rename({xdim: 'lon', ydim: 'lat', tdim:'time'})
    #
    start='%02d-%02d-%04d'%(1,1,syear_out)
    end='%02d-%02d-%04d'%(31,12,eyear_out)
    sfcelv=nc.sfcelv.sel(time = slice(start,end)).load()
    sfcelv_mean=sfcelv.mean(axis=0)
    sfcelv_std=sfcelv.std(axis=0)
    sfcelv_mean=sfcelv_mean.values
    sfcelv_std=sfcelv_std.values
    sfcelv_mean.tofile(outdir+"/CaMa_out/"+input_name+"/sfcelv_mean"+tagout+".bin")
    sfcelv_std.tofile(outdir+"/CaMa_out/"+input_name+"/sfcelv_std"+tagout+".bin")
    print (np.shape(sfcelv_mean))
  return 0
#-----
def spatial_slice(data,lllat,urlat,lllon,urlon,res=1.0):
  if res == 1.0:  
    x1 = int((lllon + 0.5 - (-179.5))/1.0)
    x2 = int((urlon - 0.5 - (-179.5))/1.0)
    y1 = int((lllat + 0.5 - (  89.5))/1.0)
    y2 = int((urlat - 0.5 - (  89.5))/1.0)
  elif res == 0.25:     
    x1 = int((lllon + 0.125 - (-179.875))/0.25)
    x2 = int((urlon - 0.125 - (-179.875))/0.25) 
    y1 = int(((89.875) - (urlat + 0.125))/0.25)
    y2 = int(((89.875) - (lllat - 0.125))/0.25)
  elif res == 0.1:     
    x1 = int((lllon + 0.05 - (-179.95))/0.1)
    x2 = int((urlon - 0.05 - (-179.95))/0.1) 
    y1 = int(((89.95) - (urlat - 0.05))/0.1)
    y2 = int(((89.95) - (lllat + 0.05))/0.1)
  #--   
  data_region = data[y1:y2 + 1, x1:x2 + 1]  
  return data_region
#================================================================
syear  =int(sys.argv[1]) #pm.start_year()
eyear  =int(sys.argv[2]) - 1#pm.end_year()
expname=sys.argv[3] #pm.expname() #"AMZ050049B" #"AMZ000049O" #"AMZ050049" #"AMZCAL049" 
runname=sys.argv[4] #pm.runname()
mapname=sys.argv[5] #pm.mapname()
ens_mem=int(sys.argv[6]) #pm.ens_mem()
ncpus  =int(sys.argv[7]) #pm.para_nums()
inyear =int(sys.argv[8]) #pm.start_year()
outyear=int(sys.argv[9]) #pm.end_year()
outdir =sys.argv[10]
#================================================================
syyyy="%04d"%(syear)
eyyyy="%04d"%(eyear)
syyyy_out="%04d"%(inyear)
eyyyy_out="%04d"%(outyear)
#================================================================
inputlist=[]
for ens in np.arange(1,ens_mem+1):
    inputname="%s%s%03d"%(expname,runname,ens)
    # print inputname
    inputlist.append([inputname,syyyy,eyyyy,syyyy_out,eyyyy_out])
#==============
p=Pool(ncpus)
p.map(stats,inputlist)
p.terminate()
# map(stats,inputlist)