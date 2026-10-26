# set up data from modelling
require(bio.lobster)
require(devtools)
require(sf)
require(PBSmapping)
require(bio.utilities)
require(dplyr)
require(ggplot2)

la()

p = list()
p$yrs = 2026

lobster.db('logs.redo',pH=p)
lobster.db('process.logs.redo')

lobster.db('logs41.redo')
lobster.db('greyzone_logs.redo')
lobster.db('fsrs.redo')
lobster.db('atSea.redo')
lobster.db("survey.redo")
lobster.db('scallop.redo')
lobster.db('landings_per')
nefsc.db('odbc.dump.redo')
nefsc.db('clean.redo')
groundfish.db('odbc.redo',datayrs = 2023:2025)
groundfish.db('gs_trawl_conversions_redo')

#this builds the tow dist and net mensuartion stuff
ILTS_ITQ_All_Data(redo_set_data = T, redo_base_data = T, biomass = F)  
g = compileAbundPresAbs_vessel_corr(redo = F,size=F)

#glorys bias corrected model needs to be run.

lobster.db('temperature.data.redo')

#get the glorys data -> git/bio.lobster.glorys/inst/R_script_to_downloadData.r
#then run in the following order 
#git/bio.lobster.glorys/inst/Analysis/CreateBiasCorr/Optimized_code/ 
# MergeGLORYS_Data.r
# biasCorrModel.r
# predictions.r


