#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from astropy.io import fits
import numpy as np
import yaml
import gt_apps as my_apps
from tqdm import tqdm

def main():
    global psf_camp #conterrà i valori campionati della PSF, array 3D (indice_PSF; THETA; ENERGY)
    global EN_camp
    global THETA_camp
    global dictft1 #dizionario che contiene le info sui file ft1
    global dictft2 #dizionario che contiene le info sui file ft1
    with open('/gpfs/glast/Users/pauletto/dativari/listadati/fileft1.yaml') as f: #file YAML i valori campionati della PSF
        dictft1 = yaml.safe_load(f)
    with open('/gpfs/glast/Users/pauletto/dativari/listadati/fileft2.yaml') as f: #file YAML i valori campionati della PSF
        dictft2 = yaml.safe_load(f)
    with open('/gpfs/glast/Users/pauletto/dativari/dati_psf2.yaml') as f: #file YAML i valori campionati della PSF
        dict_psf = yaml.safe_load(f)
    psf_camp=np.array(dict_psf['PSF_campionata']) # è un array 3D 
    psf_camp.shape=(4,4,18)  #indice PSF, indice THETA, indice ENERGY
    EN_camp=dict_psf['camp_ENERGY'] #i valori su cui è stata campionata l'energia
    THETA_camp=dict_psf['camp_THETA'] #i valori su cui è stata campionato THETA
    
    
    with open('/gpfs/glast/Users/pauletto/dativari/list_frbcat_chime_frb121102_SGR1935_vDec6.yaml') as f: #file YAML con la lista di FRB
        #per il file random usare /gpfs/glast/Users/pauletto/dativari/list_random_sources.yaml 
        dict_FRBs = yaml.safe_load(f)

    
    #preparazione del dizionario che verrà dato in output nel file .yaml
    res_dict={'FRB': [],
              'ev_associati': []
            }
    
    
    for l in tqdm(range(len(dict_FRBs['name']))): #CICLO SU TUTTI I FRB
        
        FRB={ #dizionario che contiene le info sul FRB preso in considerazione
         "name" : dict_FRBs['name'][l], 
         "ra" : dict_FRBs['ra'][l],
         "dec" : dict_FRBs['dec'][l],
         "met": dict_FRBs['met'][l]
         }
        
        res_dict['FRB'].append(FRB)
        risultato_associazione=associazione(FRB) #qua richiamamo la funzione di associazione
        res_dict['ev_associati'].append(risultato_associazione)
        
    
    with open('associazione_DT_510.yaml', 'w') as file_yaml:
        yaml.dump(res_dict, file_yaml,default_flow_style=False)
            
#FINE MAIN

def associazione(FRB):
    global dictft1 #i file sono stati aperti nel mine, usiamo variabili globali per passare le informazioni sui file ft1, ft2
    global dictft2
    
    indiceft2=find_nearest(dictft2["tstart"],FRB["met"])  #ricerca del file spacecraft relativo al FRB preso in considerazione
    if (dictft2["tstart"][indiceft2]>FRB["met"]):
        indiceft2=indiceft2-1
    ft2_file=dictft2["name"][indiceft2] #nome del file ft2 trovato
    
    in_fov=inside_fov(FRB['ra'],FRB['dec'],ft2_file,FRB['met']) #verifichiamo se la posizione di arrivo del FRB era nel FOV di Fermi
    if (not in_fov ):
        return "NOT IN FOV"
    #se non era nel FOV si ferma la funzione, altrimenti andiamo avanti    
    
    #Ricerca del file ft1 che contiene i dati di Fermi raccolti al tempo di arrivo del FRB
    indiceft1=find_nearest(dictft1["tstart"],FRB['met'])
    if (dictft1["tstart"][indiceft1]>FRB['met']):
        indiceft1=indiceft1-1
    ft1_file=dictft1["name"][indiceft1] #nome del file ft1 trovato
    
    
    #setup del gtselect sulla zona circostante il FRB
    my_apps.filter['evclass'] = 128
    my_apps.filter['evtype'] = 3
    my_apps.filter['ra'] = FRB['ra']
    my_apps.filter['dec'] = FRB['dec']
    my_apps.filter['rad'] = 4
    my_apps.filter['emin'] = 100
    my_apps.filter['emax'] = 1000000.0
    my_apps.filter['zmax'] = 105.
    my_apps.filter['tmin'] = FRB['met']-10.
    my_apps.filter['tmax'] = FRB['met']+500.
    #my_apps.filter['infile'] = '/gpfs/glast/Users/gprincip/Fermi_13years/ft1.lst' volendo si possono andare a vedere tutti i file ft1, è un po' più lento ma funziona ugualmente
    #in quel caso evitiamo il rischio che un FRB sia accaduto proprio a cavallo di due file ft1 diversi
    my_apps.filter['infile'] = ft1_file
    my_apps.filter['outfile'] = 'filefiltrato.fits'
    my_apps.filter.run() #esecuzione del gtselect
    #una volta fatto il gtselect si apre il file prodotto e si leggono gli eventi
    hdul = fits.open('filefiltrato.fits') 
    
    dati={#per comodità dizionario che contiene tutte le info che ci interessano sui fotoni
          'TIME' : [], 
          'RA' : [],
          'DEC' : [],
          'THETA' : [],
          'ENERGY' : [],
          'EVENT_TYPE' : []
          }
    for key in dati: #riempimento del dizionario
       dati[key]=hdul[1].data[key]
 
    
    dict_ev_associati={  #preparazione del dizionario di output
               "tempo_MET": [], 
               "ra" : [],
               "dec" : [],
               "energia": [],
               "psf_ind" : [],
               "psf": [],
               "distanza": [],
               "theta": [],
               "tempo_relativo": []
               }
    
    
    for i in range(len(dati['TIME'])): #ciclo su tutti i fotoni
        
        candidato={ #info sul fotone i-esimo preso in considerazione
            "tempo_MET": float(dati['TIME'][i]),
            "ra": float(dati['RA'][i]),
            "dec": float(dati['DEC'][i]),
            "energia" : float(dati['ENERGY'][i]) ,
            "psf_ind" : int(psf(i,dati)),
            "psf" : 0., 
            "distanza" : 0., #calcolate qualche riga più in basso
            "theta" : float(dati['THETA'][i]),
            "tempo_relativo": float(dati['TIME'][i]-FRB['met'])
            }
        candidato["distanza"]=float(ang(FRB['dec'],FRB['ra'],candidato["dec"],candidato["ra"]))
        candidato['psf']=float(compute_psf(int(candidato['psf_ind']),float(candidato['energia']),float(dati['THETA'][i])))
        if (candidato["distanza"]<3*candidato["psf"] and candidato["distanza"]<=3.): #ora vedo compatibilità spaziale, se c'è allora metto l'evento dentro dict_ev_associati
            for y in dict_ev_associati:
                dict_ev_associati[y].append(candidato[y])
    return dict_ev_associati

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
        

def compute_psf(psf_ind,E,theta):       #dà l'angolo di errore associato ad un fotone
    global psf_camp
    global EN_camp
    global THETA_camp
    ind_THETA=find_nearest(THETA_camp,theta)
    ind_EN=find_nearest(EN_camp,E)
    return psf_camp[psf_ind,ind_THETA,ind_EN]

def ang(b1, l1, b2,l2): #dec1,ra1, dec2,ra2
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


def inside_fov(rat,dect,ft2_file, time):
    hdu = fits.open(ft2_file,memmap=True)
    t = hdu[1].data['START']
    insaa = hdu[1].data['IN_SAA']
    raz = hdu[1].data['RA_SCZ']
    decz = hdu[1].data['DEC_SCZ']
    i=find_nearest(t,time)
    if (ang( dect,rat, decz[i],raz[i])< 90 and not insaa[i]):
        return True
    else:
        return False


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

#qui si fa girare il main()
if __name__ == "__main__":
    main()
