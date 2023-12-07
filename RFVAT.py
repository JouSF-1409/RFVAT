#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 16:37:40 2019

@author: jason
"""

import os
import numpy as np
from tqdm import tqdm
from obspy import read
from VSP import PlotVAT
import sys
import configparser
import getopt

def Usage():
   print('python RFAVT.py  -r paraRFAVT.cfg')
   
def getrfdata(nsample,path):
    os.chdir(path)
    sacfiles = os.listdir()
    rayp = np.zeros(len(sacfiles))
    rf = np.zeros([nsample,len(sacfiles)])
    for i in range(len(sacfiles)):
        rf_seis = read(sacfiles[i])[0]
        rayp[i] = rf_seis.stats.sac.user0
        rf[:,i] = rf_seis.data
    indexsort = rayp.argsort()
    rayp = rayp[indexsort]
    rf = rf[:,indexsort]
    return rf,rayp

def VAT(rayp,rf,dt,t_max,e_max):
    p = np.mat(rayp**2)
    ntr = len(rayp)
    vpvs = np.transpose(np.mat(np.arange(0,e_max,0.05)))
    times = np.arange(0,t_max,dt)  
    PsSpectrum = np.zeros((len(times),len(vpvs)))
    PpPsSpectrum = np.zeros((len(times),len(vpvs)))
    trace = np.arange(ntr)
    for i in tqdm(range(len(times)),desc='RFVAT ing'):
        t = i*dt
        if t < 100: 
            Tps = t + t * vpvs * p/2
            sample = (Tps + shift)/dt
            Ps_amp = rf[sample.astype(int),trace.astype(int)]
            summary = Ps_amp.sum(axis=1)
            summary[abs(summary)<0.2] = 0
            PsSpectrum[i,:] = summary    
        Tppps = t - t * vpvs * p/2
        sample = (Tppps + shift) /dt
        PpPs_amp = rf[sample.astype(int),trace.astype(int)]
        summary = PpPs_amp.sum(axis=1)
        summary[abs(summary)<0.2] = 0
        PpPsSpectrum[i,:] = summary  
    return PsSpectrum,PpPsSpectrum

def PsEtime(t1,t2,polarity,t_max,e_max):
    Ps_max = np.zeros([len(t1),4])
    vpvs = np.transpose(np.mat(np.arange(0,e_max,0.05)))
    times = np.arange(0,t_max,dt)  
    for i in range(len(t1)):
        idx1 = int(t1[i]/dt)
        idx2 = int(t2[i]/dt)
        if polarity[i]:   
            Ps_max[i,0] = times[np.where(PsSpectrum == np.max(PsSpectrum[idx1:idx2,:]))[0][0]] #Tps
            Ps_max[i,1] = vpvs[np.where(PsSpectrum == np.max(PsSpectrum[idx1:idx2,:]))[1][0]] #E
        else:
            Ps_max[i,0] = times[np.where(PsSpectrum == np.min(PsSpectrum[idx1:idx2,:]))[0][0]] #Tps
            Ps_max[i,1] = vpvs[np.where(PsSpectrum == np.min(PsSpectrum[idx1:idx2,:]))[1][0]] #E  
    return Ps_max

if __name__ == '__main__':
    
    try:
       opts, args = getopt.getopt(sys.argv[1:], "r")
    except:
        print('Arguments are not found!')
        Usage()
        sys.exit(1)
    if opts == []:
        Usage()
        sys.exit(1)
    config = configparser.ConfigParser()
    config.read(args[0])
    data_path = config.get('path', 'data_path')
    out_path = config.get('path', 'out_path')
    nsample = int(config.get('para', 'nsample'))
    shift = float(config.get('para', 'shift'))
    dt = float(config.get('para', 'dt'))
    t_max = float(config.get('para', 't_max'))
    e_max = float(config.get('para', 'e_max'))
    t1 = np.array((config.get('para', 't1')).split(','),dtype=float)
    t2 = np.array((config.get('para', 't2')).split(','),dtype=float)
    polarity = np.array((config.get('para', 'polarity')).split(','),dtype=float)
    k1 = float(config.get('para', 'k1'))
    k2 = float(config.get('para', 'k2'))
    c = float(config.get('para', 'c'))
    thre = float(config.get('para', 'thre'))
    path = os.getcwd()
    
    rf,rayp = getrfdata(nsample = nsample,path = data_path)
    ntr = len(rayp)
    
    
    PsSpectrum,PpPsSpectrum = VAT(rayp=rayp,rf = rf , dt = dt,t_max=t_max,e_max=e_max)
    Ps_max = PsEtime(t1,t2,polarity=polarity,t_max=t_max,e_max=e_max)
    os.chdir(path)
    os.chdir(out_path)
    
    
    
    
    PlotVAT(rf=rf,PsSpectrum=PsSpectrum,PpPsSpectrum=PpPsSpectrum,Ps_max=Ps_max,\
            polarity=polarity,t_max=t_max,e_max=e_max,k1=k1,k2=k2,c=c,thre=thre,rayp=rayp,dt=dt,nsample=nsample,shift=shift)
