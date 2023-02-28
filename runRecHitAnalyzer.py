import os

cfg='PhotonClassifier/RecHitAnalyzer/python/RecHitAnalyzer_cfg.py'
#cfg='RecHitAnalyzer/python/SCRegressor_cfg_data.py'
#inputFiles_='file:../step_full.root'

#inputFiles_='/store/mc/RunIIFall17MiniAODv2/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/288ED03F-6C42-E811-96EC-0CC47A4C8E98.root'
#inputFiles_='/store/mc/RunIIFall17MiniAODv2/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/00AE0E2A-6F42-E811-8EA2-0025905B85AA.root'
#inputFiles_='/store/mc/RunIISummer16MiniAODv3/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/90000/82035724-E425-E911-9E1D-1866DA879ED8.root'

#inputFiles_ = '/store/mc/RunIIFall17MiniAODv2/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/00AE0E2A-6F42-E811-8EA2-0025905B85AA.root'                                                                                                                                                                                          

#GJet dataset 1
#inputFiles_ = '/store/mc/RunIIFall17MiniAODv2/HAHMHToAA_AToGG_MA-1GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/'

#inputFiles_ = 'root://cmsdata.phys.cmu.edu//store/mc/RunIISummer16DR80Premix/HAHMHToAA_AToGG_MA-1GeV_TuneCUETP8M1_PSweights_13TeV-madgraph_pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/FAC1D5E0-45F1-EA11-A9D1-C4346BC8D390.root'
#inputFiles_ = 'root://cmsdata.phys.cmu.edu//store/mc/RunIISummer16DR80Premix/HAHMHToAA_AToGG_MA-1GeV_TuneCUETP8M1_PSweights_13TeV-madgraph_pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/110000/FAC1D5E0-45F1-EA11-A9D1-C4346BC8D390.root'

#inputFiles_='/store/mc/RunIIFall17MiniAODv2/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/00AE0E2A-6F42-E811-8EA2-0025905B85AA.root'
#inputFiles_1='/store/mc/RunIIFall17MiniAODv2/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/04339D42-6C42-E811-BAFF-0CC47A7C360E.root'
#inputFiles_2='/store/mc/RunIIFall17MiniAODv2/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/04699198-6942-E811-ABF5-0025905A6076.root'
#inputFiles_3='/store/mc/RunIIFall17MiniAODv2/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/0A62663F-6C42-E811-BF6C-0CC47A78A4BA.root'
#inputFiles_4='/store/mc/RunIIFall17MiniAODv2/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/0AA40649-6C42-E811-8A6E-0025905A610A.root'
#inputFiles_5='/store/mc/RunIIFall17MiniAODv2/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/0AE9C549-6C42-E811-87C0-0025905B858E.root'
#inputFiles_6='/store/mc/RunIIFall17MiniAODv2/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/0E3C4952-C141-E811-BC92-484D7E8DF107.root'
#inputFiles_7='/store/mc/RunIIFall17MiniAODv2/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/0EACD1AE-6642-E811-B72C-0025905AA9CC.root'


### for finding nan
inputFiles_='/store/mc/RunIIFall17MiniAODv2/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/40000/30C55624-FE79-E811-B942-FA163E0FC4E5.root'

#inputFiles_='/store/mc/RunIIFall17MiniAODv2/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
#inputFiles_='/store/mc/RunIIFall17MiniAODv2/DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'

#inputFiles_='/store/mc/RunIIFall17MiniAODv2/QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/90000/FEE7EA67-5942-E811-9C25-0CC47A4D75EE.root'
#inputFiles_='/store/mc/RunIIFall17MiniAODv2/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/04B84384-FF41-E811-ACA3-7845C4FC3B18.root'
#root://cmsxrootd-site.fnal.gov

#maxEvents_=-1
maxEvents_=1000
skipEvents_=0
#outputFile_='output.root'
#inputTag=inputFiles_.strip('file:').strip('_FEVTDEBUG.root')
#inputTag='TEST'

#cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
#for ievt in range(1):
#if not os.path.isdir(inputTag):
#    os.system('mkdir %s'%(inputTag))
#cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,ievt)
cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,skipEvents_)
#cmd="cmsRun %s inputFiles=%s inputFiles=%s inputFiles=%s inputFiles=%s inputFiles=%s inputFiles=%s inputFiles=%s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_0,inputFiles_1,inputFiles_2,inputFiles_3,inputFiles_4,inputFiles_5,inputFiles_6,inputFiles_7,maxEvents_,skipEvents_)
print ('%s'%cmd)
os.system(cmd)
#os.system('mv c*.eps %s/'%(inputTag))
#    os.system('mv cEB*_%d.eps %s/'%(ievt+1,inputTag))

#os.system('scram b -j8')
