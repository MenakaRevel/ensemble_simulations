import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
import re
import math
import datetime
import calendar
from numpy import ma
#==================
def mk_fig(year,mon,day,ens,expname,runname="ERA5",nx=700,ny=350,mapname="conus_06min",CaMadir="/cluster/data6/menaka/CaMa-Flood_v4"):
    yyyy="%04d"%(year)
    mm="%02d"%(mon)
    dd="%02d"%(day)
    ens_char="%03d"%(ens)
    if calendar.isleap(year):
        days=366
    else:
        days=365
    # fname=
    fname="../CaMa_out/"+expname+runname+ens_char+"/outflw"+yyyy+".bin"
    outflw=np.fromfile(fname,np.float32).reshape(days,ny,nx)
    print (np.max(ma.masked_greater_equal(outflw,1.0e19)),np.min(ma.masked_greater_equal(outflw,1.0e19)))
    nextxy = CaMadir+"/map/"+mapname+"/nextxy.bin"
    nextxy = np.fromfile(nextxy,np.int32).reshape(2,ny,nx)
    outflw=outflw[day]*(nextxy[0]>0.0)*(nextxy[0]!=1e20)*1.0
    plt.clf()
    plt.imshow(ma.masked_equal(ma.masked_greater_equal(outflw,1.0e19),0.0),interpolation="nearest", origin="upper", cmap="viridis_r", vmin=0.0,vmax=3000.0)
    plt.colorbar()
    plt.savefig("../CaMa_out/"+expname+runname+ens_char+"/outflw__"+yyyy+mm+dd+ens_char+".jpg")
    print("Saving -> outflw__"+yyyy+mm+dd+ens_char+".jpg")
    plt.close()
    return 0
#===========
if __name__=="__main__":
    # mk_fig(2000,1,1,20)
    # mk_fig(2000,5,25,5)
    # mk_fig(2000,12,12,17)
    mk_fig(2015,1,12,43,"CONUS50",runname="VIC_BC")