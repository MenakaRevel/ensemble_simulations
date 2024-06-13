import numpy as np
import scipy.signal
import matplotlib.pyplot as plt
import re
import math
import datetime
from numpy import ma
#==================
def mk_fig(year,mon,day,ens,runname="ERA5",nx=3600,ny=1800,CaMadir="/cluster/data6/menaka/CaMa-Flood_v410"):
    yyyy="%04d"%(year)
    mm="%02d"%(mon)
    dd="%02d"%(day)
    ens_char="%03d"%(ens)
    # fname="/cluster/data8/abdul.moiz/20230511_CaMa-DA/Ensemble_Simulations/CaMa_in/ERA5/Roff/Roff__"+yyyy+mm+dd+ens_char+".one"
    fname="../CaMa_in/"+runname+"/Roff/Roff__"+yyyy+mm+dd+ens_char+".one"
    print fname
    roff=np.fromfile(fname,np.float32).reshape(ny,nx)
    print (np.max(roff),np.min(roff))
    # nextxy = CaMadir+"/map/glb_06min/nextxy.bin"
    # nextxy = np.fromfile(nextxy,np.int32).reshape(2,ny,nx)
    # roff=roff #*(nextxy[0]>0.0)*(nextxy[0]!=1e20)*1.0
    plt.clf()
    plt.imshow(ma.masked_less_equal(roff,0.0)*1e3,interpolation="nearest", origin="upper",vmin=0.0,vmax=5.0)
    plt.colorbar()
    # plt.savefig("./CaMa_in/"+runname+"/Roff/Roff__"+yyyy+mm+dd+ens_char+"_Moiz.jpg")
    plt.savefig("../CaMa_in/Roff__"+yyyy+mm+dd+ens_char+".jpg")
    print("Saving -> Roff__"+yyyy+mm+dd+ens_char+".jpg")
    plt.close()
    return 0
#===========
if __name__=="__main__":
    mk_fig(2000,1,1,42,runname="VIC_BC",nx=7200,ny=3000)
    # mk_fig(2000,1,1,1)
    # mk_fig(2000,1,1,20)
    # mk_fig(2000,5,25,5)
    # mk_fig(2000,12,12,17)