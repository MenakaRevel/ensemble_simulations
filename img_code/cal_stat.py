import numpy as np
import params as pm
from numpy import ma
#--
def rmse(assimilations,truths,add=1.0e-20):
  rms = []
  assimilations = np.array(assimilations)
  truths     = np.array(truths)
  # assimilations is 2d arry
  for i in range(0,len(truths)): 
    rms.append(np.sqrt(((assimilations[:,i]-truths[i]*np.ones([pm.ens_mem()],np.float32))**2).mean()))
  return rms#(rms/abs(truths+add))*100.0 
#--
def AI(assimilations,corruptions,truths):
  AIval = []
  assimilations = np.array(assimilations)
  corruptions = np.array(corruptions)
  truths     = np.array(truths)
  #--
  error =abs(truths-np.mean(corruptions,axis=1))/(truths+1e-20) 
  AIval = 1. - abs((truths-np.mean(assimilations,axis=1))/(truths-np.mean(corruptions,axis=1))+1.0e-20)
  
  return ma.masked_where(error<0.1,AIval).filled(0.0)
#--
def AI_new(assimilations,corruptions,truths):
  AIval = []
  assimilations = np.array(assimilations)
  corruptions = np.array(corruptions)
  truths     = np.array(truths)
  #--
  error =abs(truths-np.mean(corruptions,axis=1))/(truths+1e-20)  
  AIval = 1. - abs((np.mean(assimilations,axis=1)-np.mean(corruptions,axis=1))/(truths-np.mean(corruptions,axis=1)+1.0e-20)-1.)
  
  #return ma.masked_where(error<0.1,AIval).filled(0.0) 
  return AIval,error
#--
def AI_new1(assimilations,corruptions,truths):
  AIval = []
  assimilations = np.array(assimilations)
  corruptions = np.array(corruptions)
  truths     = np.array(truths)
  #--
  error =abs(truths-np.mean(corruptions,axis=1))/(truths+1e-20)
#  AIval =  abs((truths-np.mean(assimilations,axis=1))/(truths-np.mean(corruptions,axis=1)+1.0e-20)-1.)
  AIval =abs(1. - (np.mean(assimilations,axis=1)-np.mean(corruptions,axis=1))/(truths-np.mean(corruptions,axis=1)+1.0e-20))
  
  return ma.masked_where(error<0.1,AIval).filled(0.0) 
#--
def pBias(assimilations,corruptions,truths):
  pB = []
  assimilations = np.array(assimilations)
  corruptions = np.array(corruptions)
  truths     = np.array(truths)
  
  pB   = 100.0 * np.nansum(np.nanmean(assimilations,axis=1)-truths)/np.nansum(truths)+1e-20
  pB_c = 100.0 * np.nansum(np.nanmean(corruptions,axis=1)-truths)/np.nansum(truths)+1e-20

  return pB,pB_c
