from CRABClient.UserUtilities import config#, getUsernameFromSiteDB                                                                                                                                      
config = config()

# To submit to crab:                                                         
# crab submit -c crabConfig_data.py                                          
# To check job status:                                                       
# crab status -d <config.General.workArea>/<config.General.requestName>      
# To resubmit jobs:                                                          
# crab resubmit -d <config.General.workArea>/<config.General.requestName>    

# Local job directory will be created in:
# <config.General.workArea>/<config.General.requestName>

date = 'Feb26_2023'
#config.General.requestName = 'test_'+date+'_AA_1GeV'
config.General.requestName = 'test_'+date+'_GJ_pT40toInf'
config.General.workArea = 'crab_run'
config.General.transferOutputs = True
config.General.transferLogs = True

# CMS cfg file goes here:         
config.JobType.psetName = 'PhotonClassifier/RecHitAnalyzer/python/RecHitAnalyzer_cfg.py'
config.JobType.pluginName = 'Analysis'
#config.JobType.maxMemoryMB = 1500
#config.JobType.maxMemoryMB = 3072
#config.JobType.maxJobRuntimeMin = 2750 # mins
config.JobType.maxMemoryMB = 2500

# Define units per job here:          

#config.Data.inputDataset = '/GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
#config.Data.inputDataset = '/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'

config.Data.inputDataset = '/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'


#config.Data.inputDataset = '/HAHMHToAA_AToGG_MA-1GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'
#config.Data.inputDataset ='/HAHMHToAA_AToGG_MA-1GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM'

#config.Data.inputDataset = '/HAHMHToAA_AToGG_MA-0p4GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
#config.Data.inputDataset = '/HAHMHToAA_AToGG_MA-0p1GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
#config.Data.inputDataset = '/HAHMHToAA_AToGG_MA-1GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'


#config.Data.inputDataset = '/DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
config.Data.splitting = 'Automatic'
#config.Data.splitting = 'EventAwareLumiBased'#'Automatic'#'EventBased'
#config.Data.unitsPerJob = 100 #500
config.Data.totalUnits = -1 #50000
#config.Data.totalUnits = 10 #50000

#config.Data.lumiMask = 'Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
#config.Data.unitsPerJob = 100000
#config.Data.totalUnits = -1 # -1 to process all units 
config.Data.publication = False

config.Data.ignoreLocality = True
config.JobType.allowUndistributedCMSSW = True


# Output files will be stored in config.Site.storageSite at directory:   
# <config.Data.outLFNDirBase>/<config.Data.outputPrimaryDataset>/<config.Data.outputDatasetTag>/

config.Data.outputDatasetTag = config.General.requestName
#config.Site.storageSite = 'T3_US_FNALLPC' 
config.Site.whitelist = [ 'T2_US_Nebraska' ]
#config.Site.whitelist = [ 'T2_US*' ]
#config.Site.ignoreGlobalBlacklist = True
config.Site.storageSite = 'T3_US_CMU'
config.Data.outLFNDirBase = '/store/user/kypark/Updated_ATagger_GJ_HLT_Iso_{}/'.format(date)+date

