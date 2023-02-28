import os
import pickle
import numpy as np
#import ROOT as r
import pyarrow as pa
import pyarrow.parquet as pq
import math
import time
import torch
#import torch_geometric as geo
import random
import yaml

names = ['idx','E','pt','eta','phi','pho_id','pi0_p4','m','Xtz','ieta','iphi','p4','pho_p4','pho_vars','X','dR','y','bdt','ancestry','genId','pteta','wgt','Xk','Xk_full','Xtzk','X6','X4','X5','pu']

date = 'July2022'

class ParquetDataset:
    def __init__(self, filename,cols=None):
        self.parquet = pq.ParquetFile(filename)
        self.cols = cols # if None, read all columns
        #self.cols = ['X_jets.list.item.list.item.list.item','y'] 
    def __getitem__(self, index):
        data = self.parquet.read_row_group(index, columns=self.cols).to_pydict()

        iters = self.cols if self.cols != None else names
        for name in iters:
            if name == 'idx':
                data['idx'] = np.int64(data['idx'][0])
            else:
                try:
                    data[name] = np.float64(data[name][0]) # everything else is doubles
                except KeyError:
                    pass
                
        return dict(data)
    def __len__(self):
        return self.parquet.num_row_groups

def round_to_pixel(pos):
    dst = 0.01726
    return dst * round(pos/dst)

class ParquetDatasetGNN:
    
    def __init__(self, filename, layer=None, knn=3):

        self.filename = filename
        self.parquet = pq.ParquetFile(filename)
        self.layer = layer
        self.cols = ['RHgraph_nodes_img','PFTgraph_nodes_img','pt','y','bdt','wgt','eta','pu']
        self.knn = knn
        
    def __getitem__(self, index):

        data = self.parquet.read_row_group(index, columns=self.cols).to_pydict()

        data['pt'] =  torch.reshape(torch.Tensor(data['pt']).float() , (1,))
        data['y'] =   torch.reshape(torch.Tensor(data['y']).float()  , (1,))
        data['bdt'] = torch.reshape(torch.Tensor(data['bdt']).float(), (1,))
        data['wgt'] = torch.reshape(torch.Tensor(data['wgt']).float(), (1,))
        data['eta'] = torch.reshape(torch.Tensor(data['eta']).float(), (1,))
        data['pu'] =  torch.reshape(torch.Tensor(data['pu']).float()   , (1,))

        # make graph
        node_list = []
        #node_pos = []
        for i in range(len(data['RHgraph_nodes_img'][0][0])):
            if data['RHgraph_nodes_img'][0][0][i][0] != 0 :
                xphi = round_to_pixel(data['RHgraph_nodes_img'][0][0][i][1])
                xeta = round_to_pixel(data['RHgraph_nodes_img'][0][0][i][2])
                node_list.append([1., data['RHgraph_nodes_img'][0][0][i][0], xphi, xeta, 0., 0., 0.])
                #node_pos.append( [xphi, xeta] ) 
        for i in range(len(data['PFTgraph_nodes_img'][0][0])):
            node_list.append([2.] + list(data['PFTgraph_nodes_img'][0][0][i][0:3]) + list(data['PFTgraph_nodes_img'][0][0][i][3:6]))
            #node_pos.append( list(data['PFTgraph_nodes_img'][0][0][i][1:3]) )

        edge_list = [[], []]
        edge_wgt = []
        if self.knn > 0:
            for i in range(len(node_list)):
                closest = []
                for j in range(len(node_list)):
                    if i != j:
                        #print(node_list[i], node_list[j])
                        closest.append( [ j, find_dist(node_list[i], node_list[j]) ] )

                #closest.sort(key = keyfun) 
                random.shuffle(closest)
                
                for j in range(min(self.knn, len(node_list) - 1)):
                    try:
                        if closest[j][1] > 0.:
                            edge_list[0].append(i)
                            edge_list[1].append(closest[j][0])
                            edge_wgt.append([ closest[j][1] ])
                    except:
                        print(j, len(closest), len(node_list))
                        raise

            data['G'] = geo.data.Data(x=torch.tensor(node_list), edge_index=torch.tensor(edge_list), edge_attr=torch.tensor(edge_wgt)) #pos=node_pos)
        else:
            data['G'] = geo.data.Data(x=torch.tensor(node_list))
        del data['RHgraph_nodes_img']
        del data['PFTgraph_nodes_img']
        
        return dict(data)

    def __len__(self):
        return self.parquet.num_row_groups
        
def keyfun(element):
    return element[1]

def find_dist(node_1, node_2):
    return np.sqrt( (node_1[3] - node_2[3])**2 + (node_1[2] - node_2[2])**2 ) 

class ParquetDatasetCNN:

    def __init__(self, filename,layer=None):

        self.parquet = pq.ParquetFile(filename)
        self.layer = layer
        self.cols = ['pt','y','bdt','wgt','eta','pu','X_RH_energyT', 'X_RH_energyZ', 'X_PFT_pT', 'X_PFT_d0', 'X_PFT_z0']
        self.norms = [40., 20., 15., 6., 20.]

    def __getitem__(self, index):
        data = self.parquet.read_row_group(index, columns=self.cols).to_pydict()

        data['pt'] = np.float32(data['pt'])        
        data['y'] = np.int64(data['y'][0])
        data['bdt'] = np.float32(data['bdt'])
        data['wgt'] = np.float64(data['wgt'])
        data['eta'] = np.float64(data['eta'])
        data['pu'] = np.float64(data['pu'])
        
        keys = ['X_RH_energyT', 'X_RH_energyZ', 'X_PFT_pT', 'X_PFT_d0', 'X_PFT_z0']
        X = [0,0,0,0,0]
        if self.layer == None:
            for i in range(len(X)):
                X[i] = np.float32(data[keys[i]][0][0])/self.norms[i] # to get pixel intensities roughly within [0, 1]
                del data[keys[i]]
            data['X'] = np.array(X)
        else:
            data['X'] = np.float32(data[keys[self.layer]][0])/self.norms[self.layer] # to get pixel intensities roughly within [0, 1]

        data['X'][data['X'] < 1.e-3] = 0.
        
        return dict(data)
    
    def __len__(self):
        return self.parquet.num_row_groups

class ParquetToNumpy:
    def __init__(self, type, nfiles, fname, cols=None, nevts = -1):
        
        start = time.time()

        colstr = ''
        if cols != None:
            for col in cols:
                colstr += col

        objloc = "obj/"+type+"_"+fname.replace('/','')+"_"+str(nevts)+"_"+colstr+"_dict.pkl"
        print(objloc)
        
        try:
            self.datadict = pickle.load( open( objloc, "rb" ) )
            if self.datadict != None:
                return
        except IOError:
            pass
        except EOFError:
            pass

        self.datadict = {}

        for i in range(nfiles):
            if nfiles > 1:
                parquet = ParquetDataset(fname+str(i))
            else:
                parquet = ParquetDataset(fname)
                
            if nevts == -1:
                nevts = len(parquet)
            print('>> Initiating file {} - location {}, length {}'.format(i, fname+str(i), len(parquet)))

            if cols == None:
                iters = names
            else:
                iters = cols
                
            for j in range(nevts):
                if j % 5000 == 0:
                    print(">> Running %d - nevts: %d time elapsed: %f s"%(i, j, time.time()-start))
                 
                row = parquet[j]
   
                for name in iters:
                    if name not in list(self.datadict.keys()):
                        self.datadict[name] = []
                                
                    self.datadict[name].append(row[name])
                    
                if i == 0 and j == 1:
                    print('>> Datadict at i == 0, j == 1:', self.datadict)

        for name in iters:
            try:
                self.datadict[name] = np.stack(self.datadict[name], axis=0)
            except KeyError:
                pass

        pickle.dump( self.datadict, open( objloc, "wb" ) )

    def __getitem__(self, name):
        return self.datadict[name]

    def getDict(self):
        return self.datadict

class timer:
    def __init__(self):
        self.start = time.time()

    def time(self):
        return time.time() - self.start

class obj_weights:

    def __init__(self):
        #with open('/uscms/home/acrobert/nobackup/gnn_classifier/pteta_ratios.yaml','r') as file:
        with open('./pteta_ratios.yaml','r') as file:
            self.wgt_dct = yaml.safe_load(file)
        self.ptbins = self.wgt_dct['ptbins']
        self.etabins = self.wgt_dct['etabins']
        self.nptbins = len(self.ptbins) - 1
        self.netabins = len(self.etabins) - 1

    def get_weights_torch(self,y,pt,eta,wgt):
        
        ly = len(y)
        wv = torch.zeros([ly,1],dtype=torch.float32).cuda()
        for i in range(ly):
            yi = int(y[i].item())
            if yi == 1:
                wv[i][0] = wgt[i][0]
            elif yi == 0:
                pti = pt[i].item()
                eti = eta[i].item()
                
                ipt = None
                for j in range(self.nptbins):
                    if pti >= self.ptbins[j] and pti < self.ptbins[j+1]:
                        ipt = j
                        break

                iet = None
                for j in range(self.netabins):
                    if eti >= self.etabins[j] and eti < self.etabins[j+1]:
                        iet = j
                        break
                    
                if ipt == None or iet == None:
                    print(pti, eti, ipt, iet)
                wv[i][0] = wgt[i][0] * self.wgt_dct[ipt][iet]

        return wv

    def get_weights_np(self,y,pt,eta,wgt):

        if y == 1:
            return wgt
        else:

            ipt = None
            for i in range(self.nptbins):
                if pt >= self.ptbins[i] and pt < self.ptbins[i+1]:
                    ipt = i
                    break
                    
            iet = None
            for i in range(self.netabins):
                if eta >= self.etabins[i] and eta < self.etabins[i+1]:
                    iet = i
                    break
                    
            wgt_adj = self.wgt_dct[ipt][iet]
            return wgt * wgt_adj
    
def imgplot(data,name,ind,indivs=False,ratio=False,ds2=None):
    ld = len(data)
    hist = r.TH2F('data '+name,name+' (n='+str(ld/1000)+'k); Relative i#phi; Relative i#eta; E [GeV]',32,0.,32.,32,0.,32.)

    for i in range(32):
        for j in range(32):
            if not ratio:
                bin_val = float(sum([ data[k][ind][i][j] for k in range(ld) ])) / ld
            else:
                bin_p = float(sum([ data[k][ind][i][j] for k in range(ld) ])) / ld
                bin_f = float(sum([ ds2[k][ind][i][j] for k in range(ld) ])) / ld
                if bin_f == 0.:
                    bin_val = 0.
                else:
                    bin_val = bin_p/bin_f

            hist.Fill(j,i,bin_val)

    hist.SetStats(0)
    r.gStyle.SetNumberContours(100)

    c = r.TCanvas('layer '+name)
    c.SetRightMargin(0.15)
    c.SetLogz()
    hist.Draw('colz')
    c.Write()

    del c

    if indivs:

        for i in range(5):
            indnm = name+' '+str(i)
            hist_indiv = r.TH2F('data_'+indnm,indnm+'; Relative i#phi; Relative i#eta; E [GeV]',32,0.,32.,32,0.,32.)
            for j in range(32):
                for k in range(32):
                    if not ratio:
                        hist_indiv.Fill(k,j,float( data[i][ind][j][k] ))
                    else:
                        if float( ds2[i][ind][j][k] ) == 0.:
                            hist_indiv.Fill(k,j,0.)
                        else:
                            hist_indiv.Fill(k,j,float( data[i][ind][j][k] )/float( ds2[i][ind][j][k] ))

            hist_indiv.SetStats(0)
            cx = r.TCanvas("indiv_layer_"+indnm)
            cx.SetRightMargin(0.15)
            cx.SetLogz()
            hist_indiv.Draw('colz')
            cx.Write()
            
            del cx
            del hist_indiv

    return hist

class imgplot_dyn:
    def __init__(self,ld,name,resc,zax):
        r.gStyle.SetNumberContours(100)
        self.n = str(int(round(ld/1000)))+'k' if ld >= 1000 else str(int(round(ld)))
        self.hist = r.TH2F('data '+name,name+' (n='+self.n+'); Relative i#phi; Relative i#eta; '+zax,32,0.,32.,32,0.,32.)
        self.hist.SetStats(0)
        self.resc = resc

        self.name = name
        #self.c = r.TCanvas('layer '+name)
        #self.c.SetRightMargin(0.15)
        #self.c.SetLogz()

    def fill(self,X):
        for i in range(32):
            for j in range(32):
                if X[i][j] != 0:
                    self.hist.Fill(j,i,X[i][j]/float(self.resc))

    def ret(self):
        #self.hist.Draw('colz')
        #self.c.Write()
        return self.hist, self.name

def drawlayerhist(hist,name):
    c = r.TCanvas('layer '+name)
    c.SetRightMargin(0.15)
    c.SetLogz()
    hist.Draw('colz')

    mini = 1000000.
    maxi = 0.
    for i in range(32):
        for j in range(32):
            val = hist.GetBinContent(i,j)
            if val != 0:
                if val > maxi:
                    maxi = val
                if val < mini:
                    mini = val

    hist.SetMinimum(mini)
    hist.SetMaximum(maxi)
    c.Write()
    
    del c

def combine_parquet_files(input_folder, target_path, npqs_per_file=-1):

    outlist = []

    if npqs_per_file == -1:
        npqs_per_file = 100000000

    made_writer = False

    ldir = os.listdir(input_folder)
    ld = len(ldir)

    noutfiles = int(math.ceil(float(ld) / float(npqs_per_file))) 
    print(noutfiles)

    for i in range(noutfiles):
        for j in range(npqs_per_file):
            idx = i*npqs_per_file + j
            if idx >= ld:
                break
            file_name = ldir[idx]
            print(file_name)
            
            pqds = ParquetDataset(os.path.join(input_folder, file_name))
            for k in range(len(pqds)): 
                dct = pqds[k]            
                pqdata = [pa.array([d]) if np.isscalar(d) or type(d) == list else pa.array([d.tolist()]) for d in dct.values()]
                table = pa.Table.from_arrays(pqdata, dct.keys())
    
                if not made_writer:
                    
                    writer = pq.ParquetWriter(target_path+'.'+str(i), table.schema, compression='snappy')
                    outlist.append(target_path+'.'+str(i))
                    made_writer = True
        
                writer.write_table(table)
                                
            del pqds
            
        writer.close()
        if idx >= ld:
            break
        made_writer = False
        del writer

    return outlist



def combine_parquet_files_evts(input_folder, target_path, nevts_per_file=-1, n_out_files=-1):

    outlist = []

    if nevts_per_file == -1:
        nevts_per_file = 10000000000

    made_writer = False

    ldir = os.listdir(input_folder)
    ld = len(ldir)

    for i in range(n_out_files):
        for j in range(nevts_per_file):

            pqds = ParquetDataset(os.path.join(input_folder, file_name))
            for k in range(len(pqds)):
                dct = pqds[k]

                pqdata = [pa.array([d]) if np.isscalar(d) or type(d) == list else pa.array([d.tolist()]) for d in dct.values()]
                table = pa.Table.from_arrays(pqdata, dct.keys())

                if not made_writer:

                    writer = pq.ParquetWriter(target_path+'.'+str(i), table.schema, compression='snappy')
                    outlist.append(target_path+'.'+str(i))
                    made_writer = True

                writer.write_table(table)

            del pqds

        writer.close()
        if idx >= ld:
            break
        made_writer = False

    return outlist

def check_nan(tensor, idxs):
    try:
        size = list(tensor.size())
    except:
        size = list(tensor.shape)

    lst = []

    if len(size) > 1:
        for i in range(size[0]):
            check = check_nan(tensor[i], idxs+[i])
            if check != []:
                lst += check
    else:
        for i in range(size[0]):
            if tensor[i] != tensor[i]:
                lst += [(idxs+[i],tensor[i].item())]

    return lst

def checknan(tensor):
    return check_nan(tensor, [])

def check_inf(tensor, idxs):
    try:
        size = list(tensor.size())
    except:
        size = list(tensor.shape)

    lst = []

    if len(size) > 1:
        for i in range(size[0]):
            check = check_inf(tensor[i], idxs+[i])
            if check != []:
                lst += check
    else:
        for i in range(size[0]):
            if math.isinf(tensor[i]):
                lst += [(idxs+[i],tensor[i].item())]

    return lst

def checkinf(tensor):
    return check_inf(tensor, [])

