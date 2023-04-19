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
from numpy import ma
import random
import re
import calendar
import math

distopen=1.0 #pm.distopen(1)
diststd=0.1 #pm.diststd(1)
# true_run=pm.true_run(1) # for true ensemble
runname="E2O"#pm.runname(1) # get runoff name
open_list=np.arange(1,7+1)
diststd_num=7
numname=49
# dist std for ensembles
distopen_ranges={}
fname="./ensemble_E2O_%02d.txt"%(numname)
with open(fname,"w") as f:
    ens_num=1
    for runens in open_list:
        distopen_range=rd.normal(1,diststd,diststd_num)
        distopen_range=np.sort(distopen_range)
        distopen_range=distopen_range.astype(np.float32)
        for distopen in distopen_range:
            ens_char="%03d"%(ens_num)
            # distopen_ranges[ens_char]=distopen_range#[0.25,1.00,1.25]
            line="%s  %3.2f\n"%(ens_char,distopen)
            f.write(line)
            ens_num=ens_num+1