#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 22 17:26:52 2020

@author: planet_science
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def PlotVAT(rf,PsSpectrum,PpPsSpectrum,Ps_max,polarity,t_max,e_max,k1,k2,c,thre,rayp,dt,nsample,shift):
    vpvs = np.arange(0,e_max,0.05)
    times = np.arange(0,t_max,dt)
    EstimatePpPs = np.zeros(len(polarity))
    ntr = len(rayp)
    for i in range(len(polarity)):
        if Ps_max[i,0] < c:
            EstimatePpPs[i] = Ps_max[i,0] * k1
        else:
            EstimatePpPs[i] = Ps_max[i,0] * k2   
            
    mpl.use("Agg")
    plt.figure(figsize=(15,10))
    rect1 = [0.05,0.1,0.25,0.8]
    rect2 = [0.35,0.1,0.25,0.8]
    rect3 = [0.65,0.3,0.3,0.4]
    ax1 = plt.axes(rect1)
    ax2 = plt.axes(rect2)
    ax3 = plt.axes(rect3)  
    
    plt.sca(ax1)
    plt.pcolormesh(vpvs,times,PsSpectrum,vmax=2.5,cmap='jet')
    plt.ylabel(r'$T_{Ps}$ $(p=0)$ (s)',fontsize=12)
    for i in range(len(polarity)):
        plt.plot(Ps_max[i,1],Ps_max[i,0],'k*')
        plt.text(Ps_max[i,1]+3,Ps_max[i,0]+0.5,'Ps({})'.format(i+1),fontsize=7)
    plt.text(-11,2.35,'Ps(20)',fontsize=7)
    plt.text(-11,3.92,'Ps(35)',fontsize=7)
    plt.text(-13,22,'Ps(210)',fontsize=7)
   
    plt.xlabel(r'$E_{Ps}$',fontsize=12)
    plt.ylim(0,30)
    ax1.invert_yaxis()
    plt.colorbar()

    plt.sca(ax2)
    plt.pcolormesh(vpvs,times,PpPsSpectrum,vmax=2.5,cmap='jet')
    plt.ylabel(r'$T_{PpPs}$ $(p=0)$ (s)',fontsize=12)
    for i in range(len(polarity)):
        if Ps_max[i,0] < c:
            k = 2
        else:
            k = 4
        if polarity[i]: 
            plt.axhline(y=EstimatePpPs[i]-k, c="orange", ls="-.", lw=1)
            plt.axhline(y=EstimatePpPs[i]+k, c="orange", ls="-.", lw=1)
        else:
            plt.axhline(y=EstimatePpPs[i]-k, c="deepskyblue", ls="-.", lw=1)
            plt.axhline(y=EstimatePpPs[i]+k, c="deepskyblue", ls="-.", lw=1)                
        plt.text(80,EstimatePpPs[i]+1,'Est.PpPs({})'.format(i+1),fontsize=7)
        index1 = int((EstimatePpPs[i] - k)/dt)
        index2 = int((EstimatePpPs[i] + k)/dt)
        if np.max(PpPsSpectrum[index1:index2,:]) > thre:
            Ps_max[i,2] = times[np.where(PpPsSpectrum == np.max(PpPsSpectrum[index1:index2,:]))[0][0]] #Tps
            Ps_max[i,3] = vpvs[np.where(PpPsSpectrum == np.max(PpPsSpectrum[index1:index2,:]))[1][0]] #Tps            
            plt.plot(Ps_max[i,3],Ps_max[i,2],'k*')
    plt.xlabel(r'$E_{PpPs}$',fontsize=12)
    plt.ylim(0,100)
    ax2.invert_yaxis()
    plt.colorbar()     
    
    plt.sca(ax3)
    x = np.arange(ntr)*2.5 + 40
    for i in range(len(polarity)):
        if Ps_max[i,0] < c:
            k = 2
        else:
            k = 4
        index1 = int((EstimatePpPs[i] - k)/dt) 
        index2 = int((EstimatePpPs[i] + k)/dt)
        if np.max(PpPsSpectrum[index1:index2,:]) > thre:
            Ps_max[i,2] = times[np.where(PpPsSpectrum == np.max(PpPsSpectrum[index1:index2,:]))[0][0]] #Tps 
            Ps_max[i,3] = vpvs[np.where(PpPsSpectrum == np.max(PpPsSpectrum[index1:index2,:]))[1][0]] #Tps            
            tps = Ps_max[i,0] + Ps_max[i,0]*(Ps_max[i,1]/2)*rayp**2
            tppps = Ps_max[i,2] - Ps_max[i,2]*(Ps_max[i,3]/2)*rayp**2
            plt.plot(tps,x,'r--')
            plt.plot(tppps,x,'m--')
    plt.plot(tps,x,'r--',label='Ps')
    plt.plot(tppps,x,'m--',label='PpPs')  
    plt.legend(loc='upper right')
    times = np.arange(nsample)*dt - shift
    for i in range(ntr):
        rr = rf[:,i] * 20 + rayp[i] * 1000   
        plt.plot(times,rr,'b',linewidth=1)      
    plt.xlim(-1,100)
    plt.xlabel('Time(s)')
    plt.ylabel('Ray parameter (s/km)')
    y = np.around(np.arange(ntr/2)*0.005 + 0.04,decimals=5)
    labels = list(y)
    y = np.arange(ntr/2)*5 + 40
    plt.yticks(y,labels) 
    plt.savefig('Example.png',bbox_inches='tight')
