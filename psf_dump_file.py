#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pyIrfLoader
import numpy as np
import yaml
import math

def main():
    #psfnumero, theta
    PSF=np.zeros([4,4,18]) #array 3D 
    # la prima componente indica il 'numero' della PSF, la seconda campiona gli angoli theta e la terza l'energia
    val_energy=[] #questi array verranno riempiti con i valori su cui
    val_theta=[] # verranno campionati theta e energia
    val_energy=np.geomspace(10,10**6, num=18).tolist()
#    for ind_en in range(18):
#        np.linspace()
#        val_energy.append()
        #campiono da 10 a 10^6 MeV con 10 campionature per decade
    
    for ind_ang in range(4):
        val_theta.append(float(10+20*ind_ang))
        #campiono da 5° a 75° ogni 10°
    
    for ind_psf in range(4):
        for ind_ang in range(4):        
            for ind_en in range(18):
                PSF[ind_psf,ind_ang,ind_en]=float(compute_psf(ind_psf, val_energy[ind_en], val_theta[ind_ang]))
    
    contenitore={} #per dump su file YAML, inserisco le variabili calcolate
    contenitore['PSF_campionata']=PSF.tolist()
    contenitore['camp_ENERGY']=val_energy
    contenitore['camp_THETA']=val_theta
    with open('dati_psf2.yaml', 'w') as file_yaml:
        yaml.dump(contenitore, file_yaml,default_flow_style=False)
    
    
    
def compute_psf(psf_ind,E,theta):       #dà l'angolo di errore associato ad un fotone
    irf_name='P8R3_SOURCE_V3::PSF'+str(psf_ind)
    pyIrfLoader.Loader_go() 
    irfFactory = pyIrfLoader.IrfsFactory_instance()
    irfs = irfFactory.create(irf_name)
    psf=irfs.psf()
    phi=0
    r = np.concatenate(([0],np.logspace(-3.0,np.log10(45.0),1000)))
    deltax = np.radians(r[1:] - r[:-1])
    xc = 0.5*(r[:-1]+r[1:])
    y = np.array(psf.value(xc,E,theta,0.0))
    cdf = 2*np.pi*np.sin(np.radians(xc))*y*deltax
    cdf = np.cumsum(cdf)
    cdf = np.concatenate(([0],cdf))
    cdf /= cdf[-1]
    frac = 0.68
    indx = np.searchsorted(cdf, frac) - 1
    cont_radius= float((frac - cdf[indx])/(cdf[indx+1] - cdf[indx])  *(r[indx+1] - r[indx]) + r[indx])
    return cont_radius


if __name__ == "__main__":
    main()