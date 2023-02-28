import pyarrow.parquet as pq
import pyarrow as pa # pip install pyarrow==0.7.1
import ROOT
import numpy as np
import glob, os
import sys
from utils import *
import yaml
import math
import time
#np.set_printoptions(threshold=sys.maxsize)

import argparse

nPhoCut = 200000
#nPhoCut = 100000


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-i', '--infile', default=['test_Nov29_2022_AA_0p6GeV'], nargs='+', type=str, help='Input root file.')
#parser.add_argument('-i', '--infile', default=['test_Nov20_2022_GJ_pT40toInf_All'], nargs='+', type=str, help='Input root file.')
#parser.add_argument('-i', '--infile', default=['pfml_'+date+'_pq'], nargs='+', type=str, help='Input root file.')
parser.add_argument('-o', '--outdir', default='/eos/uscms/store/user/kyungmip/PARQUET_nPho_200k/', type=str, help='Output pq file dir.')
parser.add_argument('-d', '--decay', default='HToAA_MA-0p6GeV', type=str, help='Decay name.')
parser.add_argument('-n', '--idx', default=0, type=int, help='Output root file index.')
parser.add_argument('-s', '--split', default=1, type=int, help='How many files to create.')
parser.add_argument('-t', '--start', default=0, type=int, help='Event start index.')
parser.add_argument('-e', '--nevents', default=1000000, type=int, help='Number of events.')
args = parser.parse_args()

start = time.time()
rhTreeStr = args.infile

# root file locations from crab jobs
file_loc = '/store/user/kypark/ATagger_AA_0p6GeV/Nov29_2022/'

#file_loc = '/store/user/kyungmip/ATagger_Nov27_AA/Nov27_2022/'
#file_loc = '/eos/uscms/store/user/kyungmip/ATagger_Nov20_GJ/Nov20_2022/'
#file_loc = '/store/user/kyungmip/ATagger_Nov24_NoDeltaR_AA_Full/Nov24_2022_NoDeltaR/'

#dsetdct = {'gj1': 'GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8', 
#           'gj2': 'GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8',
#           'gj3': 'GJet_Pt-20toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8',
#           'dp': 'DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_pythia8' }

dsetdct = { 
          'HToAA_MA-0p6GeV': 'HAHMHToAA_AToGG_MA-0p6GeV_TuneCP5_PSweights_13TeV-madgraph_pythia8'
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
        return filename, count_out
    else:
        return complete_path(path.split('\n')[-2])

print(" >> Input file:",rhTreeStr)
rhTree = ROOT.TChain("fevt/RHTree")

#dbloc = 'root://cmseos.fnal.gov/'
dbloc = 'root://cmsdata.phys.cmu.edu/'

if 'test' not in rhTreeStr[0]:
#if 'pfml' not in rhTreeStr[0]:
    for f in rhTreeStr:
        rhTree.Add(f)
else: #  "pfml_July2021_pq" for reading from cmsdata

    key = args.decay
    #key = args.decay.split('-')[0]
    path, count = complete_path( file_loc + dsetdct[key] )
    #path, count = complete_path( file_loc + dsetdct[key] + '/ATagger_AA_0p6GeV/' )
    #path, count = complete_path( file_loc + dsetdct[key] + '/pfml_'+date+'_ntup_'+key+'/' )
    print('PATH: ', path)

    #os.system('xrdfs root://cmseos.fnal.gov/ ls '+path+' > tmp.txt')
    os.system('xrdfs root://cmsdata.phys.cmu.edu/ ls '+path+' > tmp.txt')
    outlist = open('tmp.txt', 'r').read().split('\n')
    outlist.remove('')
    outlist.remove(path+'/log')

    print('outlist: ', outlist)
    assert len(outlist) == count

    for out in outlist:
        rhTree.Add(dbloc+out)

nEvts = rhTree.GetEntries()
assert nEvts > 0
print(" >> nEvts:",nEvts)

iStart = args.start
iEnd = iStart + args.nevents
if iEnd > nEvts:
    iEnd = nEvts

print(" >> Processing entries: [",iStart,"->",iEnd,")")
outStr = '%s/%s.%s.parquet.%d'%(args.outdir, args.infile[0], args.decay, args.idx)
print(" >> Output file:",outStr)


##### EVENT SELECTION START #####

sw = ROOT.TStopwatch()
sw.Start()

nreals = 0
nfakes = 0

dataset = ''
wgt = 1.

#superweights_loc = '/afs/cern.ch/work/a/acrobert/gnn_classifier/superweights.yaml'
#print('Loading cross sections and selection efficiencies from '+superweights_loc)
#superweights_file = open(superweights_loc,'r')
#superweights = yaml.load(superweights_file, Loader=yaml.Loader)

datasetdict = dsetdct
#datasetdict = {'gj1': 'GJet_Pt-20to40', 'gj2': 'GJet_Pt-40toInf', 'gj3': 'GJet_Pt-20toInf', 'dp': 'DiPhotonJets'}
key = args.decay
#key = args.decay.split('-')[0]

#if key not in ['gj1','gj2','gj3','dp']:
#    wgt = 1.
    #dataset = 'GJet_arbitrary'
#elif 'gj' in key:
#    wgt = superweights[date][key]['xs'] * superweights[date][key]['nsel'] / superweights[date][key]['ntot']
#    dataset = datasetdict[key]
#else:
wgt = 1.
dataset = datasetdict[key]

print('>> Dataset: ', dataset)

total_weights = 0.
nPhos = 0
data = {} # Arrays to be written to parquet should be saved to data dict
nSkipped = 0

for iEvt in range(iStart,iEnd):
    
    # Initialize event
    rhTree.GetEntry(iEvt)
    
    #print(rhTree.nPU)
    #print(rhTree.pho_bdt)
    #print(rhTree.PFTimg_pT)
    #break    

    if iEvt % 100 == 0:
        print(" .. Processing entry",iEvt,"nPhos:",nPhos,'time:',time.time() - start)
    
    #if rhTree.m0 < 100. or rhTree.m0 > 110.:
    #  continue

    nAEvt = len(rhTree.A_mass)
    for j in range(nAEvt):
        data['A'] = [
                rhTree.A_mass[j],
                rhTree.A_DR[j],
                rhTree.A_E[j],
                rhTree.A_pT[j],
                rhTree.A_eta[j],
                rhTree.A_phi[j]
        ]

    data['mHgen'] = rhTree.mHgen

    nPhoEvt = len(rhTree.pho_pT)
    idx = [rhTree.runId, rhTree.lumiId, rhTree.eventId]
    #Xk_full = np.array(rhTree.TracksPt_EB).reshape(1,170,360)

    for i in range(nPhoEvt):

        y = 1 if rhTree.PhoObj_genMatchedPdgId[i]==35 else 0
        #y = 1 if rhTree.PhoObj_AncestorPdgId[i]==22 else 0
        if y == 1:
            nreals += 1
        elif y == 0:
            nfakes += 1

        if 'fakes' in args.decay:
            if rhTree.PhoObj_AncestorPdgId[i] == 22:
                nSkipped += 1
                continue
        elif 'photons' in args.decay:
            if rhTree.PhoObj_AncestorPdgId[i] != 22:
                nSkipped += 1
                continue

        data['dataset'] = dataset
        data['wgt'] = wgt
        data['idx'] = idx + [i]
        
        data['E'] = rhTree.pho_E[i]
        data['pt'] = rhTree.pho_pT[i]
        data['phi'] = rhTree.pho_phi[i]
        data['eta'] = rhTree.pho_eta[i]

        if data['eta'] >= 170-16:
            continue

        data['ancestry'] = [ rhTree.PhoObj_recoIdx[i], rhTree.PhoObj_genMatchedPdgId[i], rhTree.PhoObj_AncestorPdgId[i] ]
        data['genId'] = [ rhTree.PhoObj_AncestorPdgId[i] ]

        data['y'] = [ y ]
        data['bdt'] = [ rhTree.pho_bdt[i] ]
      
        total_weights += wgt

        data['A_genId'] =  rhTree.PhoObj_genMatchedPdgId[i]        
        data['A_status'] = rhTree.PhoObj_Status[i]     

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

        data['pu'] = np.array(rhTree.mPU)
        data['tpu'] = np.array(rhTree.mPUTrue)
        data['puIdx'] = np.array(rhTree.puBX)
        data['puReco'] = np.array(rhTree.nPU)
        data['puTrue'] = np.array(rhTree.puTrue)
        
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

        #pqkeys = []
        #pqdata = []
        #for key in data.keys():
        #    pqkeys.append(key)
        #    d = data[key]
        #    if np.isscalar(d):
        #        pqdata.append(pa.array([d]))
        #    elif type(d) == list:
        #        if len(d) != 1:
        #            pqdata.append(pa.array([d]))
        #        else:
        #            pqdata.append(pa.array(d))
        #    else:
        #        pqdata.append(pa.array([d.tolist()]))
        #table = pa.Table.from_arrays(pqdata, names=pqkeys)

        #create table for this event
        pqdata = [pa.array([d]) if np.isscalar(d) or type(d) == list else pa.array([d.tolist()]) for d in list(data.values())]
        table = pa.Table.from_arrays(pqdata, list(data.keys()))

        if nPhos == 0:
            writer = pq.ParquetWriter(outStr, table.schema, compression='snappy')

        #write table
        writer.write_table(table)

        nPhos += 1
     
    if nPhos > nPhoCut:
        break
   
writer.close()

print(" >> nPhos =",nPhos,"total wgts =",total_weights,"skipped:",nSkipped)
print(" >> Real time:",sw.RealTime()/60.,"minutes")
print(" >> CPU time: ",sw.CpuTime() /60.,"minutes")
print(" >> ======================================")

#os.system('globus-url-copy %s gsiftp://cmseos-gridftp.fnal.gov//eos/uscms/store/user/acrobert/pfdata/%s/%s'%(outStr,date,outStr)) 
#os.system('cp %s /afs/cern.ch/work/a/acrobert/gnn_classifier/%s'%(outStr,outStr))

sw.Stop()
print(" >> Photons:",nreals,"Fakes:",nfakes,"Sum of weights:",total_weights)
print(" >> Real time:",sw.RealTime()/60.,"minutes")
print(" >> CPU time: ",sw.CpuTime() /60.,"minutes")
print(" >> ======================================")

pqIn = pq.ParquetFile(outStr)
print(pqIn.metadata)
print(pqIn.schema)

'''
superweights_loc = '/afs/cern.ch/work/a/acrobert/gnn_classifier/superweights.yaml'
print('Loading cross sections and selection efficiencies from '+superweights_loc)
superweights_file = open(superweights_loc,'r')
superweights = yaml.load(superweights_file, Loader=yaml.Loader)

print(superweights)
if date not in superweights.keys():
    superweights[date] = {}
if key not in superweights[date].keys():
    superweights[date][key] = {'xs': superweights['June2022'][key]['xs']}
if 'nfk' not in superweights[date][key].keys():
    superweights[date][key]['nfk'] = 0
if 'npho' not in superweights[date][key].keys():
    superweights[date][key]['npho'] = 0
superweights[date][key]['nfk'] += nfakes
superweights[date][key]['npho'] += nreals
print(superweights)

new_superweights = yaml.dump(superweights, Dumper=yaml.Dumper)
with open(superweights_loc,'w') as loc:
    loc.write(new_superweights)

ns = open(superweights_loc,'r')
print('ns:')
for line in ns.readlines():
    print(line.strip('\n'))
'''
