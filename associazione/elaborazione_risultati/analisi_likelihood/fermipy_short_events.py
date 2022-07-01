# IMPORT LIBRARIES
from fermipy.gtanalysis import GTAnalysis
import numpy as np

# DEFINE THE NAME OF MAIN SOURCES AND OUTPUT DIRECTORY 
my_source = "180916"
galname = "galdiff"
isoname = "isodiff"
outdir =  "/gpfs/glast/Users/pauletto/likelihood/6fot/out"

config_file = "/gpfs/glast/Users/pauletto/likelihood/6fot/fermipy_config.yaml"

# INITIALISE YOUR ANALYSIS AND CREATES countsmap, expmap ...

gta = GTAnalysis(config_file,logging={'verbosity' : 3}, fileio={"outdir":str(outdir)}) 
gta.setup()


# OPTMIZE THE MODEL
gta.optimize()
gta.print_roi()


# FREE OBJECTS in the model 
gta.free_source(galname)
gta.free_source(isoname)
#gta.free_sources(minmax_ts=[10.0,None], free=True,pars='norm')
gta.free_sources(distance=3.0, free=True,pars='norm')
gta.free_sources(minmax_ts=[50,None],pars='norm')
gta.free_source(my_source)


# PRELIMINARY FIT
gta.fit(update=True) #, min_fit_quality = 3)


#VISUALIZE THE OUTPUT MODEL
gta.write_roi('fit_model_1fit')
gta.print_roi()

#from astropy.io import fits
#hdu=fits.open(outdir + 'fit_model_1fit.fits')

# CHECKS ON THE PRELIMINARY MODEL
# CREATE RESIDMAP and TSMAP for a model check)
gta.residmap('initial_fit_with',make_plots=True)
gta.residmap('initial_fit_exclude',exclude=my_source,make_plots=True)
gta.tsmap('initial_ts_map_exclude', exclude=my_source,make_plots=True)
gta.tsmap('initial_ts_map_with',make_plots=True)


#CUSTOMIZE THE MODEL
# ELIMINATE SOURCES WITH TS<1 (EXCLUDING YOUR SOURCE)
for src in gta.roi.get_sources():
	if (src['ts']<1 or np.isnan(src['ts'])) and src.name!=my_source and src.name!='galdiff' and src.name!='isodiff':
		gta.delete_source(src.name)
            
            
# SEARCH FOR NEW SOURCES
model = {'Index' : 2.0, 'SpatialModel' : 'PointSource'}
srcs = gta.find_sources(model=model, sqrt_ts_threshold=5.0,min_separation=0.5)
gta.write_roi('fit_new_sources')


# SECOND FIT
gta.fit(update=True)
gta.write_roi('fit_model_1fit_new_sources')


# LOCALIZATION FITS
src = gta.roi[my_source]
if src['ts']>10 and src.name==my_source:
	gta.localize(my_source, free_background=True, update=True, make_plots=True)
	gta.fit(update=True)


# SED FIT
gta.sed(my_source, prefix=my_source+"sed_cl95", ul_confidence=0.95, loge_bins=[2.0,2.5,3.0,3.5,4.0,4.5,5.0,6.0], free_background=True, free_radius=2, write_fits=True, outfile='sed.fits', make_plots=True)
gta.write_roi('fit_model_sed1_1fit')


# CHECKS ON THE FINAL MODEL
# CREATE RESIDMAP and TSMAP for a model check)
gta.residmap('fit_with',make_plots=True)
gta.residmap('fit_exclude',exclude=my_source,make_plots=True)
gta.tsmap('ts_map_exclude', exclude=my_source,make_plots=True)
gta.tsmap('ts_map_with',make_plots=True)


# LIGHTCURVE (bin 3 hours - 10800 s)
#lc = gta.lightcurve(my_source, binsz=10800, write_fits=True, make_plots=True)



