#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.io import fits
import numpy as np
import yaml
import gt_apps as my_apps
from tqdm import tqdm
import random

def main():
    
    global psf_camp
    global EN_camp
    global THETA_camp
    with open('/gpfs/glast/Users/pauletto/dativari/dati_psf2.yaml') as f: #file YAML i valori campionati della PSF 
        dict_psf = yaml.safe_load(f)


    psf_camp=np.array(dict_psf['PSF_campionata']) # è un array 3D 
    psf_camp.shape=(4,4,18)  #indice PSF, indice THETA, indice ENERGY
    EN_camp=dict_psf['camp_ENERGY'] #i valori su cui è stata campionata l'energia
    THETA_camp=dict_psf['camp_THETA'] # stessa cosa con theta
    niterazioni=int(input("Inserire numero di posizioni random su cui effettuare il clustering: "))
    for l in range(niterazioni):
        
        ra=random.uniform(0.,360.)
        dec=random.uniform(-90.,90.)
        
        
        print("Clustering sulla posizione:")
        print("RA" , ra)
        print("DEC", dec)
        
        #gtselect sulla posizione random su tutti i dati disponibili e con un raggio di 3°
        my_apps.filter['evclass'] = 128
        my_apps.filter['evtype'] = 3
        my_apps.filter['ra'] = ra
        my_apps.filter['dec'] = dec
        my_apps.filter['rad'] = 3
        my_apps.filter['emin'] = 100
        my_apps.filter['emax'] = 1000000.0
        my_apps.filter['zmax'] = 105.
        my_apps.filter['tmin'] = 239557417 
        my_apps.filter['tmax'] = 649610000
        my_apps.filter['infile'] = '/gpfs/glast/Users/gprincip/Fermi_13years/ft1.lst'
        my_apps.filter['outfile'] = 'filefiltrato_cluster.fits'
        my_apps.filter.run()
        
        hdul = fits.open('filefiltrato_cluster.fits') #apertura del file generato dal gtselect
            
    
        
        dati={#per comodità variabile con il dizionario che contiene solo le informazioni che ci servono per l'analisi
              'TIME' : [], 
              'RA' : [],
              'DEC' : [],
              'THETA' : [],
              'ENERGY' : [],
              'EVENT_TYPE' : []
              }
        
        for key in dati: #riempimento dizionario dati
            dati[key]=hdul[1].data[key]
        
        #opereremo il clustering con diversi valori di deltaT (in secondi), i valori sono in questo array:
        lista_deltaT=[10,100,500] 
        
        for tempo in lista_deltaT:
            
            risultato_cluster=ricerca_cluster(float(tempo),ra,dec,dati) 
            res_dict={} #contenitore che conterrà l'output
            
            #queste due righe sono un modo per evitare un errore nel salvataggio di questo array in formato yaml
            tempi_inizio=np.array(risultato_cluster['tempo_inizio_MET']) 
            res_dict['tempo_inizio_MET']=tempi_inizio.tolist()
            #riempimento dizionario di output
            res_dict['n_fotoni']=risultato_cluster['n_fotoni']
            res_dict['durata']=risultato_cluster['durata']
            res_dict['lista_fotoni']=risultato_cluster['lista_fotoni']
            
            #salvataggio
            with open('cluster_RA_'+str(ra)+"_DEC_"+str(dec)+'_DT_'+str(tempo)+'.yaml', 'w') as file_yaml:
                yaml.dump(res_dict, file_yaml,default_flow_style=False)
            
#FINE MAIN

def ricerca_cluster(deltaT,RA,DEC,dati): 

    lista_fotoni={  #dizionario che conterrà tutti i fotoni di tutti i cluster trovati
            "MET": [], 
            "ra" : [],
            "dec" : [],
            "energia": [],
            "psf" : [],
            "psf_ind": []
            }
    eventi_compatibili={ #dizionario con gli eventi compatibili spazialmente
               "MET": [], 
               "ra" : [],
               "dec" : [],
               "energia": [],
               "psf" : [],
               "psf_ind": []
               }
    #ciclo su tutti gli eventi che verifica la compatibilità spaziale tra essi e la direzione scelta 
    for i in tqdm(range(len(dati['TIME'])),desc="Clustering"):
        psf_temp=compute_psf(int(psf(i,dati)),float(dati['ENERGY'][i]),float(dati['THETA'][i]))
        if(ang(DEC,RA,dati['DEC'][i],dati['RA'][i])<3*psf_temp):
            #tutti gli eventi compatiili vengono inseriti nel dizionario eventi_compatibili
            eventi_compatibili['MET'].append(dati['TIME'][i])
            eventi_compatibili['ra'].append(dati['RA'][i])
            eventi_compatibili['dec'].append(dati['DEC'][i])
            eventi_compatibili['psf_ind'].append(psf(i,dati))
            eventi_compatibili['energia'].append(dati['ENERGY'][i])
            eventi_compatibili['psf'].append(compute_psf(int(psf(i,dati)),float(dati['ENERGY'][i]),float(dati['THETA'][i])))
            
    i=0 #contatore su tutti i fotoni compatibili trovati
    
    dict_cluster={ #dizionario con tutti i cluster che vengono trovati
            'tempo_inizio_MET': [],
            'durata' : [],
            'n_fotoni' : [],
            'lista_fotoni' : []
            }
    
    while(i<len(eventi_compatibili['MET'])-1): #ciclo su tutti gli eventi con cui abbiamo riempito il dizionario precedente
        k=i
        n_fotoni1=1
        diff_t=eventi_compatibili['MET'][k+1]-eventi_compatibili['MET'][i]
        #controllo la differenza dei tempi tra il fotone i-esimo e il i+1 esimo, se è <deltaT si procede
        while(diff_t<deltaT):
            n_fotoni1+=1
            k+=1 
            #controllo la differenza di tempi tra il fotone i-esimo e il i+2, i+3, i+4,... eccetera
            #finchè la differenza è <deltaT incremento il contatore di fotoni
            #se è >deltaT fermo il while
            #questi try ed except sono necessari perchè altrimenti una volta giunti al termine dell'array si avrebbe un IndexError
            try:
                diff_t=eventi_compatibili['MET'][k+1]-eventi_compatibili['MET'][i]
            except IndexError: 
                diff_t=deltaT+1
                
        if(n_fotoni1>=3): #se abbiamo rivelato almeno 3 fotoni vicini temporalmente meno di deltaT allora abbiamo un cluster
            lista_fotoni_temp={ #è il dizionario di questo singolo cluster
                    "MET": [], 
                    "ra" : [],
                    "dec" : [],
                    "energia": [],
                    "psf" : [],
                    "psf_ind": []
                    }
            
            for p in range(i,k+1): #salvo le info di ogni fotone del cluster
                for key in lista_fotoni_temp:
                    lista_fotoni_temp[key].append(float(eventi_compatibili[key][p])) #riempiamo il dizionario dei fotoni di questo cluster
                
            for key in lista_fotoni:    
                lista_fotoni[key].append(lista_fotoni_temp[key]) # aggiungiamo questo cluster al dizionario con tutti i fotoni di tutti i cluster
                
            #salviamo su dict_cluster tutte le info pertinenti a questo cluster
            dict_cluster['tempo_inizio_MET'].append(eventi_compatibili['MET'][i])
            dict_cluster['n_fotoni'].append(n_fotoni1)
            dict_cluster['durata'].append(float(eventi_compatibili['MET'][k]-eventi_compatibili['MET'][i]))
        n_fotoni1=1
        i+=1
    dict_cluster['lista_fotoni']=lista_fotoni
        
    return dict_cluster

def psf(indice,dati): #restituisce un integer che indica l'indice di qualità della PSF
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


def compute_psf(psf_ind,E,theta):       #dà l'angolo di errore associato ad un fotone
    global psf_camp
    global EN_camp
    global THETA_camp
    ind_THETA=find_nearest(THETA_camp,theta)
    ind_EN=find_nearest(EN_camp,E)
    return psf_camp[psf_ind,ind_THETA,ind_EN]


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

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

#qui si fa girare il main()
if __name__ == "__main__":
    main()
