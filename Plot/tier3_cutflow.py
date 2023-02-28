import ROOT
import glob
import os

# Read filelist from Tier3
mass = '1'
date = 'Feb26_2023'

site = 'root://cmsdata.phys.cmu.edu/'
user_area = '/store/user/kypark/'

output_name = 'Ng2_Test_ATagger_AA_HLT_Iso_Ng2_{}GeV_{}/{}/'.format(mass, date, date)
#output_name = 'ATagger_AA_HLT_Iso_Ng2_Feb25_{}GeV_{}/{}/'.format(mass, date,date)
#output_name = 'ATagger_GJ_HLT_Iso_Ng2_{}/'.format(date)
file_loc = user_area + output_name

dsetdct = {
          #'HToAA_MA-0p1GeV': 'HAHMHToAA_AToGG_MA-0p1GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8'
          #'HToAA_MA-0p2GeV': 'HAHMHToAA_AToGG_MA-0p2GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8'
          #'HToAA_MA-0p4GeV': 'HAHMHToAA_AToGG_MA-0p4GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8'
          #'HToAA_MA-0p6GeV': 'HAHMHToAA_AToGG_MA-0p6GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8'
          #'HToAA_MA-0p8GeV': 'HAHMHToAA_AToGG_MA-0p8GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8'
          'HToAA_MA-1GeV': 'HAHMHToAA_AToGG_MA-1GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8'
          #'GJ-pT40toInf': 'GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8'
          }


def complete_path(filename):

    os.system('xrdfs root://cmsdata.phys.cmu.edu/ ls '+filename+' > tmp.txt')
#    os.system('xrdfs root://cmseos.fnal.gov/ ls '+filename+' > tmp.txt')
    path = open('tmp.txt', 'r').read()

    if path == '/' or path == '':
        return '' # error

    is_output = False
    count_out = 0
    for i in range(len(path.split('\n')) - 1):
        if 'output_' in path.split('\n')[i]:
            is_output = True
            count_out += 1

    if is_output:
        return filename
    else:
        return complete_path(path.split('\n')[-2])



subDir = {  #'GJ_QCD': 'ATagger_GJ_40toInf/Nov28_2022/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/test_Nov28_2022_GJ_pT40toInf/221129_050015/0000',
            #'DiPhoton': 'ATagger_DiPhoton_Jan29_2023/Jan29_2023/DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_pythia8/test_Jan29_2023_DiPhoton/230129_174304/0000'
            'AA_0p1GeV': 'ATagger_AA_0p1GeV_Jan25_2023/Jan25_2023/HAHMHToAA_AToGG_MA-0p1GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8/test_Jan25_2023_AA_0p1GeV/230125_175547/0000',
            #'AA_0p2GeV': 'ATagger_AA_0p2GeV/Nov29_2022/HAHMHToAA_AToGG_MA-0p2GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8/test_Nov29_2022_AA_0p2GeV/221129_102245/0000',
            #'AA_0p4GeV': 'ATagger_AA_0p4GeV/Nov29_2022/HAHMHToAA_AToGG_MA-0p4GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8/test_Nov29_2022_AA_0p4GeV/221129_102119/0000',
            #'AA_0p6GeV': 'ATagger_AA_0p6GeV/Nov29_2022/HAHMHToAA_AToGG_MA-0p6GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8/test_Nov29_2022_AA_0p6GeV/221129_102041/0000',
            #'AA_0p8GeV': 'ATagger_AA_0p8GeV/Nov29_2022/HAHMHToAA_AToGG_MA-0p8GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8/test_Nov29_2022_AA_0p8GeV/221129_102005/0000',
            #'AA_1GeV': 'ATagger_AA_1GeV/Nov28_2022/HAHMHToAA_AToGG_MA-1GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8/test_Nov28_2022_AA_1GeV/221129_045858/0000'
          }

total_nPho_Trigger = 0
total_nPho_Presel = 0
total_nPho_Ancestor = 0

total_nAnalyzed = 0
total_nPassTrigger = 0
total_nPassPresel = 0

for process in dsetdct:
#for process in subDir:
#for process in ['AA_0p1GeV']:
    #fileDir = '/store/user/kypark/' + subDir[process]
    #fileDir = '/store/user/kypark/' + subDir['GJ_QCD']

    #ls_cmd = 'xrdfs {} ls {} > fileList.txt'.format(site, fileDir)


    #os.system(ls_cmd)
    #fileList = open('fileList.txt', 'r').read().split('\n')
    #fileList.remove('')
    #fileList.remove(fileDir + '/log')

    path = complete_path( file_loc + dsetdct[process] )
    #path, count = complete_path( file_loc + dsetdct[key] + '/pfml_'+date+'_ntup_'+key+'/' )
    print('PATH: ', path)

    #os.system('xrdfs root://cmseos.fnal.gov/ ls '+path+' > tmp.txt')
    os.system('xrdfs root://cmsdata.phys.cmu.edu/ ls '+path+' > tmp.txt')
    fileList = open('tmp.txt', 'r').read().split('\n')
    fileList.remove('')
    fileList.remove(path+'/log')

    #fileList.remove(path+'/output_1-103.root')
    #fileList.remove(path+'/output_1-22.root')
 
    # Get cutflow
    nPho_Trigger = 0
    nPho_Presel = 0
    nPho_Ancestor = 0

    nAnalyzed = 0
    nPassTrigger = 0
    nPassPresel = 0

    for inFile in fileList:
        inFile = site + inFile
        f = ROOT.TFile.Open(inFile, 'READ')

        fevt = f.Get('fevt')

        # Get Histograms
        hNpassed_hlt = fevt.Get('hNpassed_hlt')
        hNpassed_img = fevt.Get('hNpassed_img')
        hNRecoPho_HLT = fevt.Get('hNRecoPho_HLT')
        hNRecoPho_presel = fevt.Get('hNRecoPho_presel')
        hNRecoPho_hits = fevt.Get('hNRecoPho_hits')
        hNRecoPho_genMatching = fevt.Get('hNRecoPho_genMatching')

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

    total_nAnalyzed += nAnalyzed
    total_nPassTrigger += nPassTrigger
    total_nPassPresel += nPassPresel

    total_nPho_Trigger += nPho_Trigger
    total_nPho_Presel += nPho_Presel
    total_nPho_Ancestor += nPho_Ancestor

    print('\n'+process)
    print('Events analyzed: {}'.format(nAnalyzed))
    print('After trigger: {}'.format(nPassTrigger))
    print('After photon pre-selection: {}'.format(nPassPresel))

    print('Photons after Trigger: {}'.format(int(nPho_Trigger)))
    print('Photons after Presel: {}'.format(int(nPho_Presel)))
    print('Photons after Ancestor Matching: {}'.format(int(nPho_Ancestor)))

print('\nTOTAL')
print('Events analyzed: {}'.format(total_nAnalyzed))
print('After trigger: {}'.format(total_nPassTrigger))
print('After photon pre-selection: {}'.format(total_nPassPresel))
print('Photons after Trigger: {}'.format(int(total_nPho_Trigger)))
print('Photons after Presel: {}'.format(int(total_nPho_Presel)))
print('Photons after Ancestor Matching: {}'.format(int(total_nPho_Ancestor)))




'''
c1 = ROOT.TCanvas()
hNRecoPho_HLT.Draw()
c1.SaveAs('GJ_hNRecoPho_HLT.pdf', 'pdf')

hNRecoPho_hits.Draw()
c1.SaveAs('GJ_hNRecoPho_hits.pdf', 'pdf')

hNRecoPho_genMatching.Draw();
c1.SaveAs("GJ_hNRecoPho_genMatching.pdf", "pdf");
'''
