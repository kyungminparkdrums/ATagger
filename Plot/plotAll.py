import ROOT
import glob

ROOT.gROOT.SetBatch(1)

#inDir = '/eos/uscms/store/user/kyungmip/ATagger_Nov27_AA/Nov27_2022/HAHMHToAA_AToGG_MA-1GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8/test_Nov27_2022_AA_1GeV_Full/221127_091231/0000/'

inDir = '/eos/uscms/store/user/kyungmip/TMP/'

#inDir = '/eos/uscms/store/user/kyungmip/ATagger_Nov26_GJ_AncestorMatching/Nov26_2022/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/test_Nov26_2022_GJ_pT40toInf_All/221126_233246/0000/'

fileList = glob.glob('{}*.root'.format(inDir))[0:1]

nPho_Trigger = 0
nPho_Presel = 0
nPho_Ancestor = 0

nAnalyzed = 0
nPassTrigger = 0
nPassPresel = 0

hNRecoPho_HLT = []
hNRecoPho_hits = []
hNRecoPho_genMatching = []

for idx, inFile in enumerate(fileList):
    f = ROOT.TFile.Open(inFile, 'READ')

    fevt = f.Get('fevt')
    #fevt.ls()

    ttree = fevt.Get('RHTree')
    ttree.Print()
    # Get Histograms
    #hPhoE = fevt.Get('RHTree').Get('hPhoE')
    hPhoE = ttree.GetBranch('pho_E')

    hNpassed_hlt = fevt.Get('hNpassed_hlt')
    hNpassed_img = fevt.Get('hNpassed_img')
    hNRecoPho_HLT.append(fevt.Get('hNRecoPho_HLT'))
    hNRecoPho_hits.append(fevt.Get('hNRecoPho_hits'))
    hNRecoPho_genMatching.append(fevt.Get('hNRecoPho_genMatching'))

    # Get nEvents for cutflow table
    nAnalyzed += int(hNpassed_hlt.GetBinContent(1))
    nPassTrigger += int(hNpassed_img.GetBinContent(1))
    nPassPresel += int(hNpassed_img.GetBinContent(2))

    for i in range(1, hNRecoPho_HLT[idx].GetNbinsX()+1):
    #    print('{} photons after trigger: {}'.format(i-1, hNRecoPho_HLT.GetBinContent(i)))
        nPho_Trigger += (i-1)*hNRecoPho_HLT[idx].GetBinContent(i)

    for i in range(1, hNRecoPho_hits[idx].GetNbinsX()+1):
    #    print('{} photons after pre-selection: {}'.format(i-1, hNRecoPho_hits.GetBinContent(i)))
        nPho_Presel += (i-1)*hNRecoPho_hits[idx].GetBinContent(i)

    for i in range(1, hNRecoPho_genMatching[idx].GetNbinsX()+1):
    #    print('{} photons after ancestor matching: {}'.format(i-1, hNRecoPho_genMatching.GetBinContent(i)))
        nPho_Ancestor += (i-1)*hNRecoPho_genMatching[idx].GetBinContent(i)


    c1 = ROOT.TCanvas()
    hNRecoPho_HLT[idx].Draw()
    c1.SaveAs('GJ_hNRecoPho_HLT_{}.pdf'.format(idx), 'pdf')

    hNRecoPho_hits[idx].Draw()
    c1.SaveAs('GJ_hNRecoPho_hits_{}.pdf'.format(idx), 'pdf')

    hNRecoPho_genMatching[idx].Draw()
    c1.SaveAs("GJ_hNRecoPho_genMatching_{}.pdf".format(idx), "pdf")

    # output
    outFileName = 'outFile{}.root'.format(idx)
    outFile = ROOT.TFile.Open(outFileName, 'RECREATE')

    outFile.cd()
    hPhoE.Write()
    hNRecoPho_HLT[idx].Write()
    hNRecoPho_hits[idx].Write()
    hNRecoPho_genMatching[idx].Write()
    outFile.Close()

print('Events analyzed: {}'.format(nAnalyzed))
print('After trigger: {}'.format(nPassTrigger))
print('After photon pre-selection: {}'.format(nPassPresel))

print('Photons after Trigger: {}'.format(nPho_Trigger))
print('Photons after Presel: {}'.format(nPho_Presel))
print('Photons after Ancestor Matching: {}'.format(nPho_Ancestor))
