import ROOT
import glob
import os

site = 'root://cmsdata.phys.cmu.edu/'

fileDir = 'ATagger_AA_1GeV_pThist/Dec4_2022/HAHMHToAA_AToGG_MA-1GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8/test_Dec4_2022_AA_1GeV/221204_075117/0000'
#fileDir = 'ATagger_GJ_40toInf_pThist/Dec4_2022/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/test_Dec4_2022_GJ_pT40toInf/221204_074523/0000'
#fileDir = 'ATagger_GJ_40toInf/Nov28_2022/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/test_Nov28_2022_GJ_pT40toInf/221129_050015/0000'
fileDir = '/store/user/kypark/' + fileDir

ls_cmd = 'xrdfs {} ls {} > fileList.txt'.format(site, fileDir)

#inDir = 'root://cmsdata.phys.cmu.edu//store/user/kypark/ATagger_GJ_40toInf/Nov28_2022/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/test_Nov28_2022_GJ_pT40toInf/221129_050015/0000/'

os.system(ls_cmd)
fileList = open('fileList.txt', 'r').read().split('\n')
fileList.remove('')
fileList.remove(fileDir + '/log')
print(fileList)



#fileList = glob.glob('{}*.root'.format(inDir))
#print(fileList)

#fileList = ['root://cmsdata.phys.cmu.edu//store/user/kypark/ATagger_GJ_40toInf/Nov28_2022/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/test_Nov28_2022_GJ_pT40toInf/221129_050015/0000/output_1-100.root']

nPho_Trigger = 0
nPho_Presel = 0
nPho_Ancestor = 0

nAnalyzed = 0
nPassTrigger = 0
nPassPresel = 0

idx = 0

for inFile in fileList:
    inFile = site + inFile 
    f = ROOT.TFile.Open(inFile, 'READ')

    fevt = f.Get('fevt')

    # Get Histograms
    hNpassed_hlt = fevt.Get('hNpassed_hlt')
    hNpassed_img = fevt.Get('hNpassed_img')
    hNRecoPho_noCut = fevt.Get('hNRecoPho_noCut')
    hNRecoPho_HLT = fevt.Get('hNRecoPho_HLT')
    hNRecoPho_hits = fevt.Get('hNRecoPho_hits')
    hNRecoPho_genMatching = fevt.Get('hNRecoPho_genMatching')

    hLeadingPhoPt = fevt.Get('hLeadingPhoPt') 
    hSubleadingPhoPt = fevt.Get('hSubleadingPhoPt') 

    # Get nEvents for cutflow table
    nAnalyzed += int(hNpassed_hlt.GetBinContent(1))
    nPassTrigger += int(hNpassed_img.GetBinContent(1))
    nPassPresel += int(hNpassed_img.GetBinContent(2))

    for i in range(1, hNRecoPho_HLT.GetNbinsX()+1):
    #    print('{} photons after trigger: {}'.format(i-1, hNRecoPho_HLT.GetBinContent(i)))
        nPho_Trigger += (i-1)*hNRecoPho_HLT.GetBinContent(i)

    for i in range(1, hNRecoPho_hits.GetNbinsX()+1):
    #    print('{} photons after pre-selection: {}'.format(i-1, hNRecoPho_hits.GetBinContent(i)))
        nPho_Presel += (i-1)*hNRecoPho_hits.GetBinContent(i)

    for i in range(1, hNRecoPho_genMatching.GetNbinsX()+1):
    #    print('{} photons after ancestor matching: {}'.format(i-1, hNRecoPho_genMatching.GetBinContent(i)))
        nPho_Ancestor += (i-1)*hNRecoPho_genMatching.GetBinContent(i)

    # output
    outFileName = 'outFile{}.root'.format(idx)
    outFile = ROOT.TFile.Open(outFileName, 'RECREATE')

    outFile.cd()
    hNRecoPho_noCut.Write()
    hNRecoPho_HLT.Write()
    hLeadingPhoPt.Write()
    hSubleadingPhoPt.Write()
    outFile.Close()

    idx += 1

print('Events analyzed: {}'.format(nAnalyzed))
print('After trigger: {}'.format(nPassTrigger))
print('After photon pre-selection: {}'.format(nPassPresel))

print('Photons after Trigger: {}'.format(nPho_Trigger))
print('Photons after Presel: {}'.format(nPho_Presel))
print('Photons after Ancestor Matching: {}'.format(nPho_Ancestor))


'''
c1 = ROOT.TCanvas()
hNRecoPho_HLT.Draw()
c1.SaveAs('GJ_hNRecoPho_HLT.pdf', 'pdf')

hNRecoPho_hits.Draw()
c1.SaveAs('GJ_hNRecoPho_hits.pdf', 'pdf')

hNRecoPho_genMatching.Draw();
c1.SaveAs("GJ_hNRecoPho_genMatching.pdf", "pdf");
'''
