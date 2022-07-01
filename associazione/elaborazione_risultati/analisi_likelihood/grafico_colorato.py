#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from astropy.io import fits
import yaml
import numpy as np 
from astropy.coordinates import SkyCoord
import pyIrfLoader

def main():
    global psf_camp
    global EN_camp
    global THETA_camp
    with open('/home/tiziano/Documenti/Tirocinio/dati_psf2.yaml') as f: #file YAML i valori campionati della PSF
        dict_psf = yaml.safe_load(f)
    psf_camp=np.array(dict_psf['PSF_campionata']) # è un array 3D 
    psf_camp.shape=(4,4,18)  #indice PSF, indice THETA, indice ENERGY
    EN_camp=dict_psf['camp_ENERGY'] #i valori su cui è stata campionata l'energia
    THETA_camp=dict_psf['camp_THETA'] # stessa cosa con theta
    
    
    h2=fits.open('6fot/outsrc.fits') #in input file .fits dato in output dal grsrcprob
    ra=[]
    dec=[]
    energy=[]
    psflist=[]
    prob=[]
    isodiff=[]
    galdiff=[]
    source=[]
    frb=[29.501,65.714]
    #LETTURA DATI
    for i in range(len(h2[1].data['ra'])):
        if (h2[1].data['180916s'][i]>0.000001 and ang(h2[1].data['dec'][i],h2[1].data['ra'][i],frb[1],frb[0])<3*compute_psf(int(psf(i,h2[1].data)),float(h2[1].data['energy'][i]),float(h2[1].data['theta'][i])) and ang(h2[1].data['dec'][i],h2[1].data['ra'][i],frb[1],frb[0])<2.5): #si tengono solo fotoni notevoli, con una probabilità che non è troppo bassa           
            ra.append(h2[1].data['ra'][i])
            dec.append(h2[1].data['dec'][i])
            energy.append(h2[1].data['energy'][i])
            prob.append(h2[1].data['180916s'][i])
            isodiff.append(h2[1].data['isodiff'][i])
            galdiff.append(h2[1].data['galdiff'][i])
            source.append(h2[1].data["4FGL J0205.7+6449"][i])
            psflist.append(compute_psf(int(psf(i,h2[1].data)),float(h2[1].data['energy'][i]),float(h2[1].data['theta'][i])))
    #si stampano a schermo tutte le info dei fotoni
    print(psflist)
    print(energy)
    print(dec)
    print(ra)
    print(prob)
    print('isodiff',isodiff)
    print('galdiff',galdiff)
    print('source',source)
    
    #DA QUI IN POI SI FA IL GRAFICO
    #qua si regola la dimensione dei punti che appaiono nel grafico, è stato fatto un tuning ad occhio
    for j in range(len(psflist)):
        psflist[j]=(psflist[j]*28)**2
    colori=np.linspace(start=0, stop=len(ra),num=len(ra))
    plt.scatter(ra,dec,s=psflist, marker='o', alpha=0.3,c=colori,cmap=get_cmap(len(ra)))
    plt.scatter([29.501],[65.714], marker='+', s=400) #posizione del FRB
    for j in range(len(psflist)):
        st=str(energy[j])+ " MeV \n" + str(prob[j]) #scrittura di energia e probabilità nel grafico
    plt.xlabel("RA(°)")
    plt.ylabel("DEC (°)")
    plt.title('Fotoni relativi al FRB180916 (cluster di 6 fotoni)')
    plt.xlim(min(ra)-1.,max(ra)+1.)
    plt.ylim(min(dec)-1,max(dec)+1)
    #plt.show()
    
    plt.savefig(fname='6fot.png',dpi=300)
    

def psf(indice,dati): #restituisce un integer che indica la PSF
    infopsf=dati['EVENT_TYPE'][indice][-6:-2] 
    if(infopsf[0]==True):
        psf_out=3
    elif (infopsf[1]==True):
        psf_out=2
    elif (infopsf[2]==True):
        psf_out=1
    elif (infopsf[3]==True):
        psf_out=0
    return psf_out


def compute_psf(psf_ind,E,theta):
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


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
def get_cmap(n, name='hsv'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)

def ang(b1, l1, b2,l2):
    lat_alfa = np.deg2rad(b1)           #np.pi * b1 / 180
    lat_beta = np.deg2rad(b2)         #np.pi * b2 / 180
    lon_alfa = np.deg2rad(l1)         #np.pi * l1 / 180
    lon_beta = np.deg2rad(l2)         #np.pi * l2 / 180
    fi = abs(lon_alfa - lon_beta)
    c = np.sin(lat_beta) * np.sin(lat_alfa) + np.cos(lat_beta) * np.cos(lat_alfa) * np.cos(fi)
    if c>1 or c<-1:
        c = int(c)
    p1 = float(np.arccos(c))
    p = np.rad2deg(p1)
    return p


if __name__ == "__main__":
    main()
