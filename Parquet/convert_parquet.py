import ROOT
import pyarrow.parquet as pq
import pyarrow as pa # pip install pyarrow==0.7.1
import numpy as np
import glob
import os, sys
from utils import *
import yaml
import math
import time
import argparse

# parse arguments
parser = argparse.ArgumentParser(description='Arguments')

parser.add_argument('-i', '--inDir', default='ATagger_TwoPho_AA_0p1_NoIso_Mar24_2023', type=str, help='Input root file directory from crab jobs')
parser.add_argument('-p', '--decay', default='HToAA_MA-0p1GeV', type=str, help='Decay process.')
parser.add_argument('-o', '--outDir', default='/eos/uscms/store/user/kyungmip/2023_Mar/AA_NoIso/', type=str, help='Output pq file dir.')
parser.add_argument('-k', '--energy', default='0p1', type=str, help='Mass point or pT (Kinematics)')
parser.add_argument('-n', '--nTotalPho', default='100000', type=int, help='Number of photons')
parser.add_argument('-d', '--date', default='Mar24_2023', type=str, help='Date')

args = parser.parse_args()

# Process arguments

dsetdict = {
          'HToAA_MA-0p1GeV': 'HAHMHToAA_AToGG_MA-0p1GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8',
          'HToAA_MA-0p2GeV': 'HAHMHToAA_AToGG_MA-0p2GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8',
          'HToAA_MA-0p4GeV': 'HAHMHToAA_AToGG_MA-0p4GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8',
          'HToAA_MA-0p6GeV': 'HAHMHToAA_AToGG_MA-0p6GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8',
          'HToAA_MA-0p8GeV': 'HAHMHToAA_AToGG_MA-0p8GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8',
          'HToAA_MA-1GeV': 'HAHMHToAA_AToGG_MA-1GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8',
          'GJ-pT20to40': 'GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8',
          'GJ-pT40toInf': 'GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8'
          }

file_loc = '/store/user/kypark/{}/{}/'.format(args.inDir, args.date)

# Get the list of input files from cmu tier3
rhTreeStr = args.decay

def complete_path(filename):
    
    os.system('xrdfs root://cmsdata.phys.cmu.edu/ ls '+filename+' > tmp.txt')
    #os.system('xrdfs root://cmseos.fnal.gov/ ls '+filename+' > tmp.txt')
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
        return filename, count_out
    else:
        return complete_path(path.split('\n')[-2])

print(" >> Input file:", rhTreeStr)
rhTree = ROOT.TChain("fevt/RHTree")

#dbloc = 'root://cmseos.fnal.gov/'
dbloc = 'root://cmsdata.phys.cmu.edu/' 
path, count = complete_path( file_loc + dsetdict[args.decay] )
print('PATH: ', path)

#os.system('xrdfs root://cmseos.fnal.gov/ ls '+path+' > tmp.txt')
os.system('xrdfs {} ls '.format(dbloc) + path + ' > tmp.txt')
    
outlist = open('tmp.txt', 'r').read().split('\n')
outlist.remove('')
outlist.remove(path+'/log')

print('outlist: ', outlist)

for out in outlist:
    rhTree.Add(dbloc+out)

# Start processing
start = time.time()

nEvts = rhTree.GetEntries()
assert nEvts > 0
print(" >> nEvts:",nEvts)

iStart = 0
iEnd = iStart + nEvts

print(" >> Processing entries: [",iStart,"->",iEnd,")")

outStr = '{}/{}'.format(args.outDir, args.decay)
print(" >> Output file: {}.parquet".format(outStr))

##### EVENT SELECTION START #####

sw = ROOT.TStopwatch()
sw.Start()

nAs, nFakes, nPrompts = 0, 0, 0
isA, isFake, isPrompt = False, False, False

dataset = dsetdict[args.decay]

print('>> Dataset: ', dataset)

wgt = 1.
total_weights = 0.
nPhos, nSkipped = 0, 0

data = {} # Arrays to be written to parquet should be saved to data dict

if 'AA' in args.decay:
    isA = True
elif 'GJ' in args.decay:
    isFake = True
elif 'DiPhoton' in args.decay:
    isPrompt = True   


for iEvt in range(iStart,iEnd):
    
    # Initialize event
    rhTree.GetEntry(iEvt)
    
    if iEvt % 100 == 0:
        print(" .. Processing entry",iEvt,"nPhos:",nPhos,'time:',time.time() - start)
    

    '''
    if isA:
        nAEvt = len(rhTree.A_mass)
        for j in range(nAEvt):
            data['A'] = [
                    rhTree.A_mass[j],
                    rhTree.vA_DR[j],
                    rhTree.vA_E[j],
                    rhTree.vA_pT[j],
                    rhTree.vA_eta[j],
                    rhTree.vA_phi[j]
            ]
        #data['mHgen'] = rhTree.mHgen
    '''

    nPhoEvt = len(rhTree.pho_pT)
    idx = [rhTree.runId, rhTree.lumiId, rhTree.eventId]

    for i in range(nPhoEvt):
        y = 1 if rhTree.PhoObj_genMatchedPdgId[i]==35 else 0    # Truth value: A is 1, prompt or fake is 0
       
        # count the number of As, fakes, prompts 
        if y == 1:
            nAs += 1
        elif y == 0:
            if rhTree.PhoObj_AncestorPdgId[i]==22:
                nPrompts += 1
            else:
                nFakes += 1

        if isA:
            if rhTree.PhoObj_genMatchedPdgId[i] != 35:
                nSkipped += 1
                continue
        elif isFake:
            if rhTree.PhoObj_AncestorPdgId[i] == 22:
                nSkipped += 1
                continue
        elif isPrompt:
            if rhTree.PhoObj_AncestorPdgId[i] != 22:
                nSkipped += 1
                continue

        # save 
        data['dataset'] = dataset
        data['wgt'] = wgt
        data['idx'] = idx + [i]
        
        data['E'] = rhTree.pho_E[i]
        data['pt'] = rhTree.pho_pT[i]
        data['phi'] = rhTree.pho_phi[i]
        data['eta'] = rhTree.pho_eta[i]

        data['ancestry'] = [ rhTree.PhoObj_recoIdx[i], rhTree.PhoObj_genMatchedPdgId[i], rhTree.PhoObj_AncestorPdgId[i] ]
        data['genId'] = [ rhTree.PhoObj_AncestorPdgId[i] ]

        data['y'] = [ y ]
        data['bdt'] = [ rhTree.pho_bdt[i] ]
      
        total_weights += wgt

        data['A_genId'] =  rhTree.PhoObj_genMatchedPdgId[i]        
        #data['A_status'] = rhTree.PhoObj_Status[i]     

        '''
        data['pho_id'] = [
                rhTree.pho_r9[i]
                ,rhTree.pho_sieie[i]
                ,rhTree.pho_phoIso[i]
                ,rhTree.pho_chgIso[i]
                ,rhTree.pho_chgIsoWrongVtx[i]
                ,rhTree.pho_Eraw[i]
                ,rhTree.pho_phiWidth[i]
                ,rhTree.pho_etaWidth[i]
                ,rhTree.pho_scEta[i]
                ,rhTree.pho_sieip[i]
                ,rhTree.pho_s4[i]
            ]
        data['pho_vars'] = [
                rhTree.pho_r9[i]
                ,rhTree.pho_HoE[i]
                ,rhTree.pho_hasPxlSeed[i]
                ,rhTree.pho_sieie[i]
                ,rhTree.pho_phoIso[i]
                ,rhTree.pho_trkIso[i]
                ,rhTree.pho_chgIsoCorr[i]
                ,rhTree.pho_neuIsoCorr[i]
                ,rhTree.pho_phoIsoCorr[i]
                ,rhTree.pho_bdt[i]
            ]
        '''

        data['pu'] = np.array(rhTree.mPU)
        data['tpu'] = np.array(rhTree.mPUTrue)
        data['puIdx'] = np.array(rhTree.puBX)
        data['puReco'] = np.array(rhTree.nPU)
        data['puTrue'] = np.array(rhTree.puTrue)

        data['eventId'] = np.array(rhTree.eventId)
        
        # IMAGE LAYERS
        data['X_RH_energy']    = np.array(rhTree.RHimg_energy[i]).reshape(1,32,32)
        data['X_RH_energyT']   = np.array(rhTree.RHimg_energyT[i]).reshape(1,32,32) 
        data['X_RH_energyZ']   = np.array(rhTree.RHimg_energyZ[i]).reshape(1,32,32)
        data['X_PFT_pT']       = np.array(rhTree.PFTimg_pT[i]).reshape(1,32,32)
        data['X_PFT_d0']       = np.array(rhTree.PFTimg_d0[i]).reshape(1,32,32)
        data['X_PFT_z0']       = np.array(rhTree.PFTimg_z0[i]).reshape(1,32,32)
        
        data['X_PFT_QpT']      = np.array(rhTree.PFTimg_QpT[i]).reshape(1,32,32)
        data['X_PFT_pT_PV']    = np.array(rhTree.PFTimg_pT_PV[i]).reshape(1,32,32)
        data['X_PFT_pT_nPV']   = np.array(rhTree.PFTimg_pT_nPV[i]).reshape(1,32,32)
       
        ''' 
        # GRAPH NODES
        rh_nodes_img = np.array(rhTree.RHgraph_nodes_img[i])
        #print('rh nodes img', len(rh_nodes_img), len(rh_nodes_img)/3.)
        data['RHgraph_nodes_img'] = rh_nodes_img.reshape(1, int(len(rh_nodes_img)/3), 3)
        rh_nodes_r04 = np.array(rhTree.RHgraph_nodes_r04[i])
        #print('rh nodes r04', len(rh_nodes_r04), len(rh_nodes_r04)/3.)
        data['RHgraph_nodes_r04'] = rh_nodes_r04.reshape(1, int(len(rh_nodes_r04)/3), 3)

        pft_nodes_img = np.array(rhTree.PFTgraph_nodes_img[i])
        #print('pft nodes img', len(pft_nodes_img), len(pft_nodes_img)/6.)
        if len(pft_nodes_img) == 0:
            pft_nodes_img = np.array([-1., 0., 0., 0., 0., 0.])
        data['PFTgraph_nodes_img'] = pft_nodes_img.reshape(1, int(len(pft_nodes_img)/6), 6)
        pft_nodes_r04 = np.array(rhTree.PFTgraph_nodes_r04[i])
        #print('pft nodes r04', len(pft_nodes_r04), len(pft_nodes_r04)/6.)
        if len(pft_nodes_r04) == 0:
            pft_nodes_r04 = np.array([-1., 0., 0., 0., 0., 0.])
        data['PFTgraph_nodes_r04'] = pft_nodes_r04.reshape(1, int(len(pft_nodes_r04)/6), 6)

        skip = False
        for node in data['RHgraph_nodes_r04'][0]:
            if math.isnan(node[0]) or math.isinf(node[0]):
                print('NaN/Inf rec hit energy: {}'.format(node))
                skip = True
        for node in data['PFTgraph_nodes_r04'][0]:
            if node[0] < 0. and node[0] != -1. and node[-1] != 0.:
                print('Negative track pt: {}'.format(node))
            if math.isnan(node[0]) or math.isinf(node[0]):
                print('NaN/Inf track pt: {}'.format(node))
                skip = True
            if math.isnan(node[3]) or math.isinf(node[3]):
                print('NaN/Inf track d0: {}'.format(node))
                skip = True
            if math.isnan(node[4]) or math.isinf(node[4]):
                print('NaN/Inf track z0: {}'.format(node))
                skip = True
        if skip:
            print('skipping event '+str(i)+' due to NaN or inf')
            continue

        if skip:
            print('DID NOT SKIP')
            raise RuntimeError('did not skip event with NaN or inf')
        '''

        # create table for this event
        pqdata = [pa.array([d]) if np.isscalar(d) or type(d) == list else pa.array([d.tolist()]) for d in list(data.values())]
        table = pa.Table.from_arrays(pqdata, list(data.keys()))

        if nPhos == 0:
            writer = pq.ParquetWriter('{}.parquet'.format(outStr), table.schema, compression='snappy')

        # write table
        writer.write_table(table)

        nPhos += 1
     
    if nPhos > args.nTotalPho:
        nTotalEvt = iEvt
        break
   
writer.close()

print(" >> # of photons = {}, # of events = {}, # of events skipped = {}".format(nPhos, nTotalEvt, nSkipped))
print(" >> Real time:",sw.RealTime()/60.,"minutes")
print(" >> CPU time: ",sw.CpuTime() /60.,"minutes")
print(" >> ======================================")

sw.Stop()
print(" >> # of A = {}, # of fakes = {}, # of prompt photons = {}".format(nAs, nFakes, nPrompts))
print(" >> Real time:",sw.RealTime()/60.,"minutes")
print(" >> CPU time: ",sw.CpuTime() /60.,"minutes")
print(" >> ======================================")

pqIn = pq.ParquetFile('{}.parquet'.format(outStr))
print(pqIn.metadata)
print(pqIn.schema)

os.system('mv {}.parquet {}_{}ev_{}pho.parquet'.format(outStr, outStr, nTotalEvt, nPhos))
