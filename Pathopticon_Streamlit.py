import streamlit as st
import streamlit.components.v1 as components
import numpy as np
import pandas as pd
import scipy.stats as sc
from tqdm import tqdm
import pickle
import base64
import os
import json
import networkx as nx
from pyvis.network import Network
from bokeh.io import show
from bokeh.layouts import column
from bokeh.plotting import figure, ColumnDataSource
from bokeh.models import CustomJS, Select, HoverTool


@st.cache(show_spinner=False, allow_output_mutation=True)
def generate_allcells_list():
	allcells_list = ['A673', 'AGS', 'ASC', 'BT20', 'CD34', 'CL34', 'CORL23', 'COV644', 'DV90', 'EFO27', 'FIBRNPC', 'H1299', 
	'HCC15', 'HCT116', 'HEC108', 'HEK293T', 'HL60', 'HS27A', 'HS578T', 'HT115', 'HUH7', 'JHUEM2', 'JURKAT', 'LOVO', 'MCF10A', 
	'MDAMB231', 'MDST8', 'NCIH1694', 'NCIH1836', 'NCIH2073', 'NCIH508', 'NCIH596', 'NEU', 'NKDBA', 'NOMO1', 'OV7', 'PHH', 'PL21', 
	'RKO', 'RMGI', 'RMUGS', 'SKB', 'SKBR3', 'SKLU1', 'SKM1', 'SKMEL1', 'SKMEL28', 'SNGM', 'SNU1040', 'SNUC4', 'SNUC5', 'SW480', 
	'SW620', 'SW948', 'T3M10', 'THP1', 'TYKNU', 'U266', 'U937', 'WSUDLCL2']
	return allcells_list

@st.cache(show_spinner=False, allow_output_mutation=True)
def default_input_geneset():
	default_geneset_name = 'Jha_SARS-CoV-2_DEgenes'
	default_up = 'ACSL1\nADD3\nALDH3A1\nAQP3\nATG9B\nBLNK\nCCDC14\nCDH10\nCENPF\nCLCA2\nCPA4\nCREG2\nCTDSPL\nCXCL14\nCYP2B7P1\nCYP4F3\n\
DEK\nDSG3\nDST\nDUSP10\nENDOD1\nEPHA4\nERF\nF3\nFAM178A\nFAT2\nFAT4\nFLRT3\nFOXO1\nFYB\nGBP6\nGCLC\nGPNMB\nGPSM2\n\
GULP1\nHSD17B3\nIFITM10\nIL6R\nKCNMA1\nLCA5\nLRMP\nLY6D\nMAF\nMAP7D2\nMCOLN3\nMETTL7A\nMN1\nMTUS1\nMUC20\nMXRA5\n\
MYLK\nNANOS1\nNID1\nNMRK1\nNPNT\nNR1D2\nNR2F2\nNRP1\nOLFML2A\nPADI3\nPALMD\nPCDH7\nPDCD4\nPPARGC1A\nPRKACB\nPROS1\n\
PRSS23\nPTPRZ1\nRAB30\nRBM20\nRPL9\nSCOC\nSCP2\nSEPP1\nSERPINB13\nSESN3\nSLAMF7\nSLC26A2\nSLITRK6\nSMEK2\nSOX6\nSPARC\n\
SPTLC3\nSTON1\nTCF4\nTHBD\nTIMP3\nTP53INP1\nTSC22D3\nVAV3\nVTCN1\nZNF488'
	default_dn = 'ACLY\nACTN1\nADRB2\nALDH1A3\nALOX15B\nANGPTL4\nANXA3\nASS1\nATP1B1\nB4GALT2\nBCL2A1\nBCL3\nBID\nBIRC3\nBMP2\nBST2\n\
C15orf48\nC15orf52\nC1QTNF1\nC1S\nC2orf54\nC8orf4\nCCDC24\nCCT5\nCD82\nCDC25B\nCDC42EP2\nCDCP1\nCES1P1\nCFB\nCOL12A1\n\
COL8A1\nCRCT1\nCSF2\nCSF3\nCSNK1E\nCSRP1\nCTPS1\nCXCL1\nCXCL16\nCXCL2\nCXCL3\nCXCL5\nCXCL6\nCYP27B1\nDDX58\nDRAM1\nDTX2\n\
DUSP1\nEDN1\nEFNA1\nELOVL1\nEPHA2\nERRFI1\nETS1\nFAM129A\nFAM167A\nFDPS\nFEZ1\nFLNB\nFOSL1\nFST\nG0S2\nGBP5\nGJB2\nGJB3\n\
HAS3\nHBEGF\nHCAR3\nHDGF\nHELZ2\nHEPHL1\nHERC6\nHK2\nHMOX2\nHSD11B1\nHYI\nICAM1\nIDI1\nIER3\nIFI44\nIFI44L\nIFI6\nIFIH1\n\
IFIT1\nIFIT3\nIFITM3\nIFNGR1\nIGFBP3\nIL1A\nIL1B\nIL32\nIL36G\nIL7R\nIRAK2\nIRAK3\nIRF7\nIRF9\nISG20L2\nITGA5\nITGB3\n\
IVL\nKCTD11\nKIAA0247\nKRT16\nKRT17\nKRT23\nKRT24\nKRT4\nKYNU\nLCN2\nLDLR\nLIF\nLOC100862671\nLTB\nLTBP2\nMAFF\nMAP3K8\n\
MFSD2A\nMICAL2\nMMP1\nMMP13\nMMP9\nMRGPRX3\nMTSS1\nMX1\nMX2\nMYC\nMYEOV\nNAMPT\nNCF2\nNEDD9\nNEURL1B\nNFKB2\nNFKBIZ\n\
NLRP1\nNRCAM\nOAS1\nOAS2\nOAS3\nOSBP2\nOTUB2\nP2RY6\nPARP12\nPARP9\nPDGFB\nPDZK1IP1\nPGLYRP4\nPI3\nPLAT\nPLAU\nPLAUR\n\
PLSCR1\nPOU2F2\nPPP1R15A\nPRDM1\nPRELID1\nPRSS3\nPSME2\nPTAFR\nPTGS2\nRAB11FIP1\nRELB\nRFTN1\nRHOF\nRHOV\nS100A12\n\
S100A7\nS100A8\nS100A9\nS100P\nSAA2\nSAMD4A\nSAMD9L\nSAMHD1\nSAT1\nSBNO2\nSDC4\nSEMA4B\nSEMA7A\nSERPINA3\nSERPINB1\n\
SERPINB2\nSERPINB8\nSGPP2\nSH3KBP1\nSH3PXD2B\nSLAMF9\nSLC11A2\nSLC25A28\nSLC25A37\nSLC43A3\nSLC6A14\nSMAD3\nSMTN\nSOCS3\n\
SOD2\nSPRR1A\nSPRR1B\nSPRR2A\nSPRR2D\nSPRR2E\nSQRDL\nSTAP2\nSTAT1\nSTK25\nTAP1\nTFPI2\nTGFA\nTGM2\nTLCD1\nTLR2\nTNC\nTNF\n\
TNFAIP2\nTNFAIP3\nTNFRSF10A\nTNFSF14\nTNIP1\nTRAF3\nTRIM16\nTRIM47\nTRIML2\nTUBA1B\nTUBA1C\nTUBA4A\nTUBB\nTUBB4B\nTUFM\n\
TYMP\nUBE2L6\nUFD1L\nUSP54\nVARS\nVEGFA\nVNN1\nVNN3\nWNT7A\nWWC1\nXAF1\nXDH\nZBED2\nZC3H12A\nZC3H12C\nZSWIM4'
	return default_geneset_name, default_up, default_dn
		
@st.cache(show_spinner=False, allow_output_mutation=True)
def pandas_read_csv_tocache(f):
	out = pd.read_csv(f)
	return out

@st.cache(show_spinner=False, allow_output_mutation=True)
def pickle_open_tocache(f):	
	with open(f, 'rb') as handle:
		out = pickle.load(handle)
	return out
	
def calculate_auroc(TPR, FPR):    
    dTPR = np.concatenate((np.ediff1d(TPR), [0]))
    dFPR = np.concatenate((np.ediff1d(FPR), [0]))   
    return sum(TPR * dFPR) + sum(dTPR * dFPR)/2

@st.cache(show_spinner=False, allow_output_mutation=True)
def input_paths():

    input_paths_dict = {'QRZC_cids_inchi_smiles_path': proj_path + 'QRZC_drug_CIDs_inchi_smiles.csv',
                        'tool_path': proj_path + 'QRZC_activityStats_nooutliers_df_besttool.csv',
                        'L1000_gene_info_path': proj_path + 'GSE92742_Broad_LINCS_gene_info.txt',
                        'L1000_cell_info_path': proj_path + 'GSE92742_Broad_LINCS_cell_info.txt',
                        'L1000_inst_info_path': proj_path + 'GSE92742_Broad_LINCS_inst_info.txt',
                        'L1000_pert_info_path': proj_path + 'GSE92742_Broad_LINCS_pert_info.txt',
                        'L1000_sig_info_path': proj_path + 'GSE92742_Broad_LINCS_sig_info.txt',

                        'pos_edges_dict_path': proj_path + 'cell_pos_edges_dict_75.pickle',
                        'neg_edges_dict_path': proj_path + 'cell_neg_edges_dict_75.pickle',
                        'drugs_dict_path': proj_path + 'cell_drugs_dict_75.pickle',

                        'cgp_dir': proj_path + 'c2.cgp.v7.1.symbols.gmt',

                        'Enrichr_GEO_up_path': proj_path + 'Disease_Perturbations_from_GEO_up.txt',
                        'Enrichr_GEO_dn_path': proj_path + 'Disease_Perturbations_from_GEO_down.txt',

                        'TTD_drugs_path': proj_path + 'P1-02-TTD_drug_download.txt',
                        'TTD_InChI2CID_path': proj_path + 'TTD_drugs_InChI2CID.txt',
                        'TTD_drug_target_path': proj_path + 'P1-07-Drug-TargetMapping.csv',
                        'TTD_target_path': proj_path + 'P1-01-TTD_target_download.txt',
                        'LCB_networks_path': proj_path + 'LCB_networks/',
                        'CD_networks_path': proj_path + 'CD_networks/'}
    
    return input_paths_dict

@st.cache(show_spinner=False, allow_output_mutation=True)    
def import_L1000_metadata(L1000_gene_info_path, L1000_cell_info_path, L1000_inst_info_path, L1000_pert_info_path, L1000_sig_info_path):

    L1000_gene_info = pd.read_csv(L1000_gene_info_path, sep='\t', low_memory=False)
    L1000_cell_info = pd.read_csv(L1000_cell_info_path, sep='\t', low_memory=False)
    L1000_inst_info = pd.read_csv(L1000_inst_info_path, sep='\t', low_memory=False)
    L1000_pert_info = pd.read_csv(L1000_pert_info_path, sep='\t', low_memory=False)
    L1000_sig_info = pd.read_csv(L1000_sig_info_path, sep='\t', low_memory=False)
    
    return L1000_gene_info, L1000_cell_info, L1000_inst_info, L1000_pert_info, L1000_sig_info  

@st.cache(show_spinner=False, allow_output_mutation=True)
def generate_diG_dict(nodelist_df_dict, edgelist_df_dict, allcells):

    diG_dict = {}
    for c in tqdm(allcells['Cell_type'], position=0, leave=True):
        diG_dict[c] = nx.DiGraph()
        for ix, node in enumerate(nodelist_df_dict[c]['Id'].values):
            if nodelist_df_dict[c].iloc[ix]['Type']=='Drug':
                diG_dict[c].add_node(node, node_type='Drug', node_color='lime', node_size=15)
            elif nodelist_df_dict[c].iloc[ix]['Type']=='Gene':
                diG_dict[c].add_node(node, node_type='Gene', node_color='yellow', node_size=5) 

        for edge in edgelist_df_dict[c].values:       
            if edge[2]=='Up':
                diG_dict[c].add_edge(edge[0], edge[1], edge_type='Up', edge_color='red')
            elif edge[2]=='Down':
                diG_dict[c].add_edge(edge[0], edge[1], edge_type='Down', edge_color='deepskyblue')
                
    return diG_dict
    
#@st.cache(show_spinner=False)
def generate_input_spoke(geneset_name, geneset_up, geneset_dn):

	input_spoke = nx.DiGraph()	
	input_spoke.add_node(geneset_name, node_type='input_geneset_name', node_color='orange', node_size=20)
	for node in geneset_up:
		input_spoke.add_node(node, node_type='Gene', node_color='yellow', node_size=5)
		input_spoke.add_edge(geneset_name, node, edge_type='Up', edge_color='red')
	for node in geneset_dn:
		input_spoke.add_node(node, node_type='Gene', node_color='yellow', node_size=5)
		input_spoke.add_edge(geneset_name, node, edge_type='Down', edge_color='deepskyblue')
		
	return input_spoke
	
#@st.cache(show_spinner=False)
def generate_disease_spoke(disease_name, disease_up, disease_dn):

	disease_spoke = nx.DiGraph()	
	disease_spoke.add_node(disease_name, node_type='disease_name', node_color='pink', node_size=20)
	for node in disease_up:
		disease_spoke.add_node(node, node_type='Gene', node_color='yellow', node_size=5)
		disease_spoke.add_edge(disease_name, node, edge_type='Up', edge_color='red')
	for node in disease_dn:
		disease_spoke.add_node(node, node_type='Gene', node_color='yellow', node_size=5)
		disease_spoke.add_edge(disease_name, node, edge_type='Down', edge_color='deepskyblue')
		
	return disease_spoke

@st.cache(show_spinner=False, allow_output_mutation=True)
def process_MSigDB_CGP(cgp_dir):
    
    cgp = {}
    f = open(cgp_dir)
    for line in f:
        line = line.rstrip()
        cgp[line.split('http://')[0].split('\t')[0]] = line.split('http://')[1].split('\t')[1:]
    f.close()
        
    # identify the indices of the duplicates (after removing _UP and DN) to determine the genesets that have both _UP and _DN
    cgp_ix = pd.DataFrame(index=cgp.keys())
    cgp_ix = cgp_ix[(cgp_ix.index.str.endswith('_UP')) | (cgp_ix.index.str.endswith('_DN'))]
    cgp_ix = cgp_ix[cgp_ix.index.str[:-3].duplicated(keep=False)]

    # CGP dictionary with only the genesets with UP and DN
    cgp_updn = {i: cgp[i] for i in cgp_ix.index.values}
        
    cgp_updn_allgenes = set([j for i in cgp_updn.values() for j in i])
    cgp_updn_labels = sorted(list(set(cgp_ix.index.str[:-3])))
   
    return cgp_updn, cgp_updn_labels, cgp_updn_allgenes
    
@st.cache(show_spinner=False, allow_output_mutation=True)
def import_Enrichr_GEO(Enrichr_GEO_up_path, Enrichr_GEO_dn_path):

    Enrichr_GEO_disease_human_up = {}
    f = open(Enrichr_GEO_up_path)
    for line in f:
        line = line.rstrip()
        d = line.split('\t')[0]
        g = line.split('\t')[2:]
        if 'human' in d:
            Enrichr_GEO_disease_human_up[d] = g       
    f.close()

    Enrichr_GEO_disease_human_dn = {}
    f = open(Enrichr_GEO_dn_path)
    for line in f:
        line = line.rstrip()
        d = line.split('\t')[0]
        g = line.split('\t')[2:]
        if 'human' in d:
            Enrichr_GEO_disease_human_dn[d] = g       
    f.close()
    
    return Enrichr_GEO_disease_human_up, Enrichr_GEO_disease_human_dn

    
@st.cache(show_spinner=False, allow_output_mutation=True)
def import_TDD(TTD_drugs_path, TTD_InChI2CID_path, TTD_drug_target_path, TTD_target_path):
    
    TTD_drugs = pd.read_csv(TTD_drugs_path, sep='\t', skiprows=29)
    TTD_InChI2CID = pd.read_csv(TTD_InChI2CID_path, sep='\t', header=None)
    TTD_drugs_CID = pd.merge(TTD_drugs[TTD_drugs['DRUG__ID']=='DRUGINKE'], TTD_InChI2CID, left_on='D00AAN.1', right_on=0)
    TTD_drug_target = pd.read_csv(TTD_drug_target_path)

    f = open(TTD_target_path, 'r')
    TTD_targID_dict = {}
    for line in f:
        line = line.split('\t')
        if len(line) == 3:
            if line[1] == 'GENENAME':
                TTD_targID_dict[line[2].strip('\n')] = line[0]
    f.close()
    
    return TTD_drugs, TTD_InChI2CID, TTD_drugs_CID, TTD_drug_target, TTD_targID_dict
    
@st.cache(show_spinner=False, allow_output_mutation=True)
def get_geneset_targets(geneset_up, geneset_dn, TTD_targID_dict, TTD_drug_target, TTD_drugs_CID, pert_iname2CID):

    # getting the direct targets of the input geneset from TTD
    geneset_targetIDs = [TTD_targID_dict[i] for i in list(set(TTD_targID_dict.keys()) & (geneset_dn | geneset_up))]
    geneset_drugID_targetID = TTD_drug_target[TTD_drug_target['TargetID'].isin(geneset_targetIDs)]
    geneset_drugCIDs = TTD_drugs_CID[TTD_drugs_CID['D00AAN'].isin(geneset_drugID_targetID['DrugID'])]
    geneset_drugID_targetCID = pd.merge(geneset_drugID_targetID, geneset_drugCIDs, left_on='DrugID', right_on='D00AAN', how='left')


    # Get the subset of drugs targeting the given pathway above that is within the QRZC drugs by
    # (1) first getting the lookup table between CIDs, InCHI Keys, SMILES and pert_inames of QRZC compounds;
    # (2) then subsetting the drugs targeting the given pathway by the ones that are in QRZC, in terms of their pert_inames.
    # The resulting number of compounds is naturally small since the drugs in QRZC networks are a small subset of all the drugs 
    # in the TTD database.

    geneset_targets_pert_iname = pd.merge(geneset_drugID_targetCID, pert_iname2CID, 
                                                    left_on=1, right_on='pubchem_cid_x', how='left')
    geneset_targets_pert_iname = geneset_targets_pert_iname[~pd.isnull(geneset_targets_pert_iname['Pert_iname'])]

    return set(geneset_targets_pert_iname['Pert_iname'])


@st.cache(show_spinner=False, allow_output_mutation=True)
def process_QRZC_output(pos_edges_dict_path, neg_edges_dict_path, drugs_dict_path, L1000_gene_info):
		
    print('Processing QUIZ-C networks...')
    with open(pos_edges_dict_path, 'rb') as handle:
        cell_pos_edges_dict = pickle.load(handle)
    with open(neg_edges_dict_path, 'rb') as handle:
        cell_neg_edges_dict = pickle.load(handle)
    with open(drugs_dict_path, 'rb') as handle:
        cell_drugs_dict = pickle.load(handle)     
    
    edgelist_df_dict = {}
    nodelist_df_dict = {}
    for c in tqdm(cell_drugs_dict.keys(), position=0, leave=True):

        pos_edges_df = pd.DataFrame(np.array([(i, j) for i, j in zip(cell_drugs_dict[c], cell_pos_edges_dict[c])], dtype='object')).rename(columns={0: 'Drug', 1: 'Target'}).explode('Target')
        neg_edges_df = pd.DataFrame(np.array([(i, j) for i, j in zip(cell_drugs_dict[c], cell_neg_edges_dict[c])], dtype='object')).rename(columns={0: 'Drug', 1: 'Target'}).explode('Target')

        pos_edges_df = pos_edges_df[~pd.isnull(pos_edges_df['Target'])]
        neg_edges_df = neg_edges_df[~pd.isnull(neg_edges_df['Target'])]

        pos_edges_df = pos_edges_df.reset_index(drop=True)
        neg_edges_df = neg_edges_df.reset_index(drop=True)

        pos_edges_df['Target'] = pos_edges_df['Target'].astype(int)
        neg_edges_df['Target'] = neg_edges_df['Target'].astype(int)

        pos_edges_df['Direction'] = 'Up'
        neg_edges_df['Direction'] = 'Down'

        pos_edges_df = pd.merge(pos_edges_df, L1000_gene_info, left_on='Target', 
                                right_on='pr_gene_id', how='left')[['Drug', 'pr_gene_symbol', 'Direction']]
        neg_edges_df = pd.merge(neg_edges_df, L1000_gene_info, left_on='Target', 
                                right_on='pr_gene_id', how='left')[['Drug', 'pr_gene_symbol', 'Direction']]

        all_edges_df = pd.concat([pos_edges_df, neg_edges_df]).rename(columns={'pr_gene_symbol': 'Target'})

        if len(all_edges_df) > 0:
            edgelist_df_dict[c] = all_edges_df 

            nodelist_df_dict_drug = pd.DataFrame(columns=['Id', 'Label', 'Type'])
            nodelist_df_dict_drug['Id'] = edgelist_df_dict[c]['Drug'].unique()
            nodelist_df_dict_drug['Label'] = edgelist_df_dict[c]['Drug'].unique()
            nodelist_df_dict_drug['Type'] = 'Drug'
            nodelist_df_dict_target = pd.DataFrame(columns=['Id', 'Label', 'Type'])
            nodelist_df_dict_target['Id'] = edgelist_df_dict[c]['Target'].unique()
            nodelist_df_dict_target['Label'] = edgelist_df_dict[c]['Target'].unique()
            nodelist_df_dict_target['Type'] = 'Gene'

            nodelist_df_dict[c] = pd.concat([nodelist_df_dict_drug, nodelist_df_dict_target])


    allcells = pd.DataFrame(sorted(list(nodelist_df_dict.keys())), columns=['Cell_type'])
    allgenes = pd.DataFrame(sorted(list(set(pd.concat(nodelist_df_dict)[pd.concat(nodelist_df_dict)['Type']=='Gene']['Id']))), columns=['Gene_symbol'])
    alldrugs = pd.DataFrame(sorted(list(set(pd.concat(nodelist_df_dict)[pd.concat(nodelist_df_dict)['Type']=='Drug']['Id']))), columns=['Pert_iname'])

    allnodes = pd.DataFrame(columns=['Node_name'])
    allnodes['Node_name'] = np.concatenate([allgenes['Gene_symbol'].values, alldrugs['Pert_iname'].values])
    allnodes['Node_ID'] = allnodes.index.values
    
    return edgelist_df_dict, nodelist_df_dict, allcells, allgenes, alldrugs, allnodes

@st.cache(show_spinner=False, allow_output_mutation=True)
def import_CD_MODZ_networks(path):
    nodelist_df_dict = {}
    edgelist_df_dict = {}
    for f in os.listdir(path):
        if f.endswith('_nodes.csv'):
            nodelist_df_dict[f.split('_')[1]] = pd.read_csv(path + f)
            nodelist_df_dict[f.split('_')[1]] = nodelist_df_dict[f.split('_')[1]].rename(columns={'ID': 'Id'})
        elif f.endswith('_edges_dir.csv'):
            edgelist_df_dict[f.split('_')[1]] = pd.read_csv(path + f)
            edgelist_df_dict[f.split('_')[1]] = edgelist_df_dict[f.split('_')[1]].rename(columns={'pert_iname': 'Drug', 'gene_name': 'Target'})[['Drug', 'Target', 'Direction']]


    allcells = pd.DataFrame(sorted(list(nodelist_df_dict.keys())), columns=['Cell_type'])
    allgenes = pd.DataFrame(sorted(list(set(pd.concat(nodelist_df_dict)[pd.concat(nodelist_df_dict)['Type']=='Gene']['Id']))), 
                                columns=['Gene_symbol'])
    allgenes = allgenes.drop(0) # remove the row with gene name "-666", this is NA in Broad Institute convention
    alldrugs = pd.DataFrame(sorted(list(set(pd.concat(nodelist_df_dict)[pd.concat(nodelist_df_dict)['Type']=='Drug']['Id']))), 
                                columns=['Pert_iname'])

    allnodes = pd.DataFrame(columns=['Node_name'])
    allnodes['Node_name'] = np.concatenate([allgenes['Gene_symbol'].values, alldrugs['Pert_iname'].values])
    allnodes['Node_ID'] = allnodes.index.values
    
    return edgelist_df_dict, nodelist_df_dict, allcells, allgenes, alldrugs, allnodes

# Calculate the Signature Congruity Score (SCS) between two genesets, typically an input signature (or perturbation) and a disease signature.
def SCS(geneset_up, geneset_dn, disease_up, disease_dn):
        
        if len((geneset_up | geneset_dn) & (disease_up | disease_dn)) > 0:
            SCS_g_d = 1.0 * ((len(disease_up & geneset_up) + len(disease_dn & geneset_dn) 
                       - len(disease_dn & geneset_up) - len(disease_up & geneset_dn)) / len(disease_up | disease_dn))
            return SCS_g_d
        else:
            return np.nan        
            
@st.cache(show_spinner=False, allow_output_mutation=True)            
def PCP_geneset(geneset_up, geneset_dn, disease_sig_up_dict, disease_sig_dn_dict, geneset_name='Geneset'):

    print('Calculating PCP for input gene signature...')
    PCP_geneset_df = pd.DataFrame(index=[geneset_name], columns=sorted(disease_sig_up_dict.keys()), dtype='float')
    
    for d in tqdm(disease_sig_up_dict.keys(), position=0, leave=True):

        disease_up = set(disease_sig_up_dict[d])
        disease_dn = set(disease_sig_dn_dict[d])    
        
        PCP_geneset_df.at[geneset_name, d] = SCS(geneset_up, geneset_dn, disease_up, disease_dn)
        
    return PCP_geneset_df
    
@st.cache(show_spinner=False, allow_output_mutation=True)    
def PCP_perturbation(perturbation, edgelist_df_dict, disease_sig_up_dict, disease_sig_dn_dict):

    PCP_perturbation_df = pd.DataFrame(index=allcells['Cell_type'], columns=sorted(disease_sig_up_dict.keys()), dtype='float')
    
    for c in allcells['Cell_type']:

        cell_drugs = sorted(list(set(edgelist_df_dict[c]['Drug'])))

        if perturbation in cell_drugs:

            cell_sorted_df = edgelist_df_dict[c].sort_values(['Drug', 'Target']) 
            drug_up = set(cell_sorted_df[(cell_sorted_df['Drug']==perturbation) & (cell_sorted_df['Direction']=='Up')]['Target'].values)
            drug_dn = set(cell_sorted_df[(cell_sorted_df['Drug']==perturbation) & (cell_sorted_df['Direction']=='Down')]['Target'].values)

            for d in disease_sig_up_dict.keys():
                disease_up = set(disease_sig_up_dict[d])
                disease_dn = set(disease_sig_dn_dict[d])    

                PCP_perturbation_df.at[c, d] = SCS(drug_up, drug_dn, disease_up, disease_dn)

        else:
            PCP_perturbation_df = PCP_perturbation_df.drop(c)
            
    return PCP_perturbation_df
    
@st.cache(show_spinner=False, allow_output_mutation=True)    
def PCP_perturbation_alldrugs(alldrugs, edgelist_df_dict, disease_sig_up_dict, disease_sig_dn_dict, 
                              write_pickle=False, pickle_dir=None, pickle_name=None):
    
    PCP_perturbation_df_dict = {}

    for drug in tqdm(alldrugs['Pert_iname'], position=0, leave=True):
               
        PCP_perturbation_df_dict[drug] = PCP_perturbation(drug, edgelist_df_dict, disease_sig_up_dict, disease_sig_dn_dict)
        
    if write_pickle == True:
        with open(pickle_dir + pickle_name + '.pickle', 'wb') as handle:
            pickle.dump(PCP_perturbation_df_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    return PCP_perturbation_df_dict
    
@st.cache(show_spinner=False, allow_output_mutation=True)    
def PACOS(PCP_geneset_df, PCP_perturbation_df_dict, alldrugs, allcells, threshold=10, geneset_name='Geneset'):

    print('Running PACOS...')
    PACOS_spearman_rho_df = pd.DataFrame(index=alldrugs['Pert_iname'], columns=allcells['Cell_type'], dtype='float')
    PACOS_spearman_pval_df = pd.DataFrame(index=alldrugs['Pert_iname'], columns=allcells['Cell_type'], dtype='float')

    for drug in tqdm(alldrugs['Pert_iname'], position=0, leave=True):

        for c in PCP_perturbation_df_dict[drug].index:
            common_diseases = set(PCP_geneset_df.loc[geneset_name][~pd.isnull(PCP_geneset_df.loc[geneset_name])].index) &\
                                set(PCP_perturbation_df_dict[drug].loc[c][~pd.isnull(PCP_perturbation_df_dict[drug].loc[c])].index)
            if len(common_diseases) >= threshold:
                PACOS_spearman_rho_df.at[drug, c], PACOS_spearman_pval_df.at[drug, c] = sc.spearmanr(PCP_geneset_df.loc[geneset_name][common_diseases], 
                                                                                         PCP_perturbation_df_dict[drug].loc[c][common_diseases])
                
    return PACOS_spearman_rho_df, PACOS_spearman_pval_df
    
@st.cache(show_spinner=False, allow_output_mutation=True)    
def rankall_PACOS(PACOS_spearman_rho_df, PACOS_spearman_pval_df):
    
    PACOS_spearman_rho_df_copy = PACOS_spearman_rho_df.copy()
    PACOS_spearman_rho_df_copy[pd.isnull(PACOS_spearman_rho_df_copy)] = -666
    sorted_ix = np.argsort(PACOS_spearman_rho_df_copy.values.flatten())
    nonan_len = len(PACOS_spearman_rho_df_copy.values.flatten()[PACOS_spearman_rho_df_copy.values.flatten()!=-666])
    nonan_ranked = PACOS_spearman_rho_df_copy.values.flatten()[np.flip(sorted_ix)[0:nonan_len]]
    nonan_ranked_ix = np.flip(sorted_ix)[0:nonan_len]
    nonan_ranked_i, nonan_ranked_j = np.unravel_index(nonan_ranked_ix, PACOS_spearman_rho_df_copy.shape)

    PACOS_spearman_rho_all_ranked = pd.DataFrame()
    PACOS_spearman_rho_all_ranked['Pert_iname'] = PACOS_spearman_rho_df_copy.index[nonan_ranked_i]
    PACOS_spearman_rho_all_ranked['Cell_type'] = PACOS_spearman_rho_df_copy.columns[nonan_ranked_j]
    PACOS_spearman_rho_all_ranked['PACOS_Spearman_rho'] = np.array([PACOS_spearman_rho_df_copy.at[PACOS_spearman_rho_df_copy.index[i], 
                                                                                             PACOS_spearman_rho_df_copy.columns[j]] 
                                                                    for i, j in zip(nonan_ranked_i, nonan_ranked_j)])
    PACOS_spearman_rho_all_ranked['PACOS_Spearman_pval'] = np.array([PACOS_spearman_pval_df.at[PACOS_spearman_rho_df_copy.index[i], 
                                                                                               PACOS_spearman_rho_df_copy.columns[j]] 
                                                                    for i, j in zip(nonan_ranked_i, nonan_ranked_j)]) 
    
    return PACOS_spearman_rho_all_ranked 
    
@st.cache(show_spinner=False, allow_output_mutation=True)    
def combine_PACOS_tool(PACOS_spearman_rho_all_ranked, tool_scores, r, write_csv=False, csv_dir=None, csv_name=None, 
                       model_name='', geneset_name=''):

    PACOS_tool_merged_df = pd.merge(PACOS_spearman_rho_all_ranked, tool_scores, left_on='Pert_iname', right_on='Pert_iname', how='left')
    PACOS_tool_merged_df['PACOS_Spearman_rho_reverse'] = -1.0*PACOS_tool_merged_df['PACOS_Spearman_rho']
    PACOS_tool_merged_df['tool_score_imputed'] = PACOS_tool_merged_df['tool score scaled'].fillna(tool_scores['tool score scaled'].median())   
    PACOS_tool_merged_df['PACOS_tool_combined'] = (r*PACOS_tool_merged_df['PACOS_Spearman_rho'] + 
                                                   PACOS_tool_merged_df['tool_score_imputed']) / (r + 1.0)
    PACOS_tool_merged_df['PACOS_tool_combined_reverse'] = (-1.0*r*PACOS_tool_merged_df['PACOS_Spearman_rho'] + 
                                                   PACOS_tool_merged_df['tool_score_imputed']) / (r + 1.0)

    if write_csv == True:
        PACOS_tool_merged_df.to_csv(csv_dir + csv_name + '_%s_%s.csv' % (model_name, geneset_name), index=False)
        
    return PACOS_tool_merged_df
    
@st.cache(show_spinner=False, allow_output_mutation=True)    
def AUROC_AUPRC(ranked_array, positives, auc_threshold=5, binsize=1):

    ranked_array_positives = set(ranked_array) & positives                 
    
    if len(ranked_array_positives) >= auc_threshold:
        
        bins = np.arange(1, len(ranked_array), binsize)      
        TPR = np.zeros(len(bins))
        FPR = np.zeros(len(bins))
        precision = np.zeros(len(bins))
        recall = np.zeros(len(bins))    

        for i, n in enumerate(bins):

            topN = ranked_array[0:n]

            overlap = set(topN) & ranked_array_positives

            TP = 1.0 * len(overlap)
            FP = len(topN) - TP
            FN = len(ranked_array_positives) - TP
            TN = len(ranked_array) - (TP + FP + FN)
            TPR[i] = TP / (TP + FN)
            FPR[i] = FP / (FP + TN)
            precision[i] = TP / (TP + FP)
            recall[i] = TP / (TP + FN)

        auroc = calculate_auroc(TPR, FPR) 
        auprc = calculate_auroc(precision, recall)
    
        return auroc, auprc
    
    else:
        return np.nan, np.nan
        
@st.cache(show_spinner=False, allow_output_mutation=True)        
def PACOS_cell_AUC(PACOS_tool_merged_df, allcells, positives, auc_params={'auc_threshold':5, 'binsize':1}, rand_params={'randomize':False, 'Nrand':1},
                   models=['PACOS_Spearman_rho', 'PACOS_Spearman_rho_reverse', 'tool_score_imputed', 
                           'PACOS_tool_combined', 'PACOS_tool_combined_reverse'], write_pickle=False, pickle_dir=None, pickle_name=None):

    model_auroc_df_dict = {}
    model_auprc_df_dict = {}
    
    for m in models:
        
        cell_geneset_auroc_df = pd.DataFrame(index=allcells['Cell_type'], columns=np.arange(rand_params['Nrand']))
        cell_geneset_auprc_df = pd.DataFrame(index=allcells['Cell_type'], columns=np.arange(rand_params['Nrand']))        

        allcells_sorted = PACOS_tool_merged_df.sort_values(['Cell_type', m], ascending=[True, False])
        
        for nrand in (tqdm(np.arange(rand_params['Nrand']), position=0, leave=True) if rand_params['randomize'] == True else np.arange(rand_params['Nrand'])):
            for c in allcells['Cell_type']:  
                cell_sorted = allcells_sorted[allcells_sorted['Cell_type']==c]['Pert_iname'].values
                if rand_params['randomize'] == True: np.random.shuffle(cell_sorted)
                cell_geneset_auroc_df.at[c, nrand], cell_geneset_auprc_df.at[c, nrand] = AUROC_AUPRC(cell_sorted, positives, 
                                                                                                     auc_threshold=auc_params['auc_threshold'], 
                                                                                                     binsize=auc_params['binsize'])
            
        model_auroc_df_dict[m] = cell_geneset_auroc_df
        model_auprc_df_dict[m] = cell_geneset_auprc_df 
    
    if write_pickle == True:
        with open(pickle_dir + pickle_name + '.pickle', 'wb') as handle:
            pickle.dump(model_auroc_df_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    return model_auroc_df_dict, model_auprc_df_dict
    
@st.cache(show_spinner=False, allow_output_mutation=True)    
def PACOS_nested_prioritization(pacos_tool_merged_df, model_auroc_df_dict, rand_model_auroc_df_dict, 
                                model, write_csv=False, csv_dir=None, csv_name=None, 
                                model_name='', geneset_name=''):
    
    emp_pval_df = pd.DataFrame(index=model_auroc_df_dict[model].index, columns=['AUROC_emp_pval'])
    nonan_cells = model_auroc_df_dict[model][0][~pd.isnull(model_auroc_df_dict[model][0])].index
    for c in nonan_cells:
        emp_pval_df.at[c, 'AUROC_emp_pval'] = ((rand_model_auroc_df_dict[model].loc[c].values > 
                                                         model_auroc_df_dict[model].loc[c].values).sum() /
                                                        float(np.shape(rand_model_auroc_df_dict[model])[1]))
                    
    emp_pval_df = emp_pval_df.sort_values('AUROC_emp_pval')
    
    model_auroc_pval_df = pd.merge(model_auroc_df_dict[model], emp_pval_df, left_index=True, right_index=True)
    model_auroc_pval_df = model_auroc_pval_df.rename(columns={0: 'AUROC'})
    


    PACOS_nested_df = pd.merge(pacos_tool_merged_df[['Pert_iname', 'Cell_type', model]], 
                            model_auroc_pval_df, 
                            left_on='Cell_type', 
                            right_index=True).dropna(subset=['AUROC']).sort_values(['AUROC_emp_pval', 'AUROC', model], 
                                                                                   ascending=[True, False, False]).reset_index(drop=True)        

    if write_csv == True:
        PACOS_nested_df.to_csv(csv_dir + csv_name + '_%s_%s.csv' % (model_name, geneset_name), index=False)        
        
    return emp_pval_df, PACOS_nested_df


#########################################################################################################################
proj_path = '/Users/ardahalu/Research/CICS/L1000_project/Pathopticon_Streamlit_Docker/'


st.set_page_config(layout='wide')

st.sidebar.title(
	'''
	:eyeglasses: **Pathopticon ** 
	Linking gene and perturbation signatures through their **pathop**heno**t**yp**ic** c**on**gruity
	***
	
	'''
)
side_radio = st.sidebar.radio(
    'Please choose analysis type',
    ('Start here!', 'Inspect Networks', 'Run Pathopticon', 'Inspect Pathopticon Results')
)

with st.container():
	'''
	# :eyeglasses: **Pathopticon ** 
	Linking gene and perturbation signatures through their **pathop**heno**t**yp**ic** c**on**gruity
	***
	
	'''

input_paths_dict = input_paths()


if side_radio == 'Start here!':
	
	intro_cont = st.container()
	with intro_cont:	
		st.write(
		'''
		### Welcome! 
		Pathopticon is a computational tool to prioritize drugs that 
		'''
		)
		
		st.write('##')
		st.write(
		'''
		In a nutshell, Pathopticon takes as input (i) a gene-perturbation (i.e. drug) network based on Connectivity Map (CMap) data and 
		(ii) a set of upregulated and downregulated genes. 
		The gene-perturbation network used by default is one that we developed called **QU**antile-based **I**nstance **Z**-score **C**onsensus (QUIZ-C).
		
		Here is a primer on how to run Pathopticon and interpret its results.
		 
		'''	
		)
		
		st.write('##')
		intro_cols = st.columns(3)
		with intro_cols[0]:
			'''
			Inspecting QUIZ-C networks
			
			'''
			st.image(proj_path + 'QUIZ-C_image.png')
			
		with intro_cols[1]:
			'''
			Entering input genesets
			
			'''

		with intro_cols[2]:
			'''
			Setting Pathopticon parameters
			
			* Model: PACOS/PACOS-Tool to indicate perturbation rankings with respect to pathophenotypic congruity score (PACOS) only or PACOS and Tool score combined.
			'''

		'''
		Citation: Halu A et al. 
		 
		'''	


elif side_radio == 'Inspect Networks':

	expander_quizc = st.expander(label='Inspect QUIZ-C networks')
	with expander_quizc:
		allcells_list = generate_allcells_list()
		chosen_cell = st.selectbox('Choose cell type to be displayed', allcells_list)
	
		HtmlFile = open(proj_path + 'PACOS_cell_networks/%s_viz.html' % chosen_cell, 'r', encoding='utf-8')
		source_code = HtmlFile.read() 
		components.html(source_code, height = 1200, width=1000)
		
	expander_Enrichr = st.expander(label='Inspect the Enrichr Disease-Gene network')
	with expander_Enrichr:
	
		HtmlFile = open(proj_path + 'Enrichr_GEO_disease-gene_network_viz.html', 'r', encoding='utf-8')
		source_code = HtmlFile.read() 
		components.html(source_code, height = 1200, width=1000)

elif side_radio == 'Run Pathopticon':

# 	expander_input_desc = st.expander(label='Start here! How to use.')
# 	with expander_input_desc:	
# 		st.write('* Model: PACOS/PACOS-Tool to indicate perturbation rankings with respect to pathophenotypic congruity score (PACOS) only or PACOS and Tool score combined.')
# 		st.write('')	
					
	with st.form(key='input_panel'):
	
		form_cols1 = st.columns(3)
		mm1 = form_cols1[0].selectbox(
			'Model',
			 ['PACOS-Tool', 'PACOS'])
		mm2 = form_cols1[1].selectbox(
			'Effect',
			 ['Repress', 'Enhance'])
		nn1 = form_cols1[2].selectbox(
			'Gene-perturbation network',
			 ['QUIZ-C', 'MODZ', 'CD']) 
		gpnet_dict = {'QUIZ-C': 'QUIZ-C', 'MODZ': 'LCB', 'CD': 'CD'} # kluge to fix the old convention of referring to the MODZ networks as LCB.  
	
		form_cols2 = st.columns(3)
		rr = form_cols2[0].number_input('r value (default=2.0)', value=2.0, step=0.01)
		nn2 = form_cols2[1].number_input('Number of randomizations (default=400)', 400, step=100) 
		default_geneset_name, default_up, default_dn = default_input_geneset()
		geneset_name = form_cols2[2].text_input('Input gene signature name', default_geneset_name)
	
		form_cols2 = st.columns(2)
		input_up = form_cols2[0].text_area('Input gene signature: Up genes', default_up)
		input_dn = form_cols2[1].text_area('Input gene signature: Down genes', default_dn)
		submit_button = st.form_submit_button(label='Submit job to Pathopticon')
		#path = st.text_input('CSV file path')

		if (mm1=='PACOS') & (mm2=='Enhance'):
			PACOS_model = 'PACOS_Spearman_rho'
		elif (mm1=='PACOS') & (mm2=='Repress'):
			PACOS_model = 'PACOS_Spearman_rho_reverse'
		if (mm1=='PACOS-Tool') & (mm2=='Enhance'):
			PACOS_model = 'PACOS_tool_combined'
		elif (mm1=='PACOS-Tool') & (mm2=='Repress'):
			PACOS_model = 'PACOS_tool_combined_reverse'
	
	if submit_button:

		expander_status = st.expander(label='Pathopticon progress status')
		with expander_status:
			status_msg = st.container()

		geneset_up = set(input_up.split('\n')) - set([''])
		geneset_dn = set(input_dn.split('\n')) - set([''])

		if len(geneset_name)==0:
			form.error('Please provide an input signature name.')
		elif (len(geneset_up)==0) | (len(geneset_dn)==0):
			form.error('Please provide both up and down gene signatures.')	
		else:	
			data_load_state = status_msg.text('Loading input data...')
			QRZC_cids_inchi_smiles = pandas_read_csv_tocache(input_paths_dict['QRZC_cids_inchi_smiles_path'])
			TTD_drugs, TTD_InChI2CID, TTD_drugs_CID, TTD_drug_target, TTD_targID_dict = import_TDD(input_paths_dict['TTD_drugs_path'], 
																								   input_paths_dict['TTD_InChI2CID_path'],
																								   input_paths_dict['TTD_drug_target_path'], 
																								   input_paths_dict['TTD_target_path'])			
			geneset_targets = get_geneset_targets(geneset_up, geneset_dn, TTD_targID_dict, TTD_drug_target, 
														 TTD_drugs_CID, QRZC_cids_inchi_smiles)			

			if len(geneset_targets) >= 5:
				
				QRZC_activityStats_nooutliers_df_besttool = pandas_read_csv_tocache(input_paths_dict['tool_path'])			
				L1000_gene_info, L1000_cell_info, L1000_inst_info, L1000_pert_info, L1000_sig_info = import_L1000_metadata(input_paths_dict['L1000_gene_info_path'], 
																														   input_paths_dict['L1000_cell_info_path'],
																														   input_paths_dict['L1000_inst_info_path'], 
																														   input_paths_dict['L1000_pert_info_path'], 
																														   input_paths_dict['L1000_sig_info_path'])


				Enrichr_GEO_disease_human_up, Enrichr_GEO_disease_human_dn = import_Enrichr_GEO(input_paths_dict['Enrichr_GEO_up_path'], 
																								input_paths_dict['Enrichr_GEO_dn_path'])

				data_load_state.text('Loading input data...done.')



				data_load_state = status_msg.text('Processing %s networks...' % nn1)
				if nn1 == 'QUIZ-C':
					edgelist_df_dict, nodelist_df_dict, allcells, allgenes, alldrugs, allnodes = process_QRZC_output(input_paths_dict['pos_edges_dict_path'], 
																												 input_paths_dict['neg_edges_dict_path'],
																												 input_paths_dict['drugs_dict_path'],
																												 L1000_gene_info)
				else:
					edgelist_df_dict, nodelist_df_dict, allcells, allgenes, alldrugs, allnodes = import_CD_MODZ_networks(input_paths_dict['%s_networks_path' % gpnet_dict[nn1]])				
				diG_dict = generate_diG_dict(nodelist_df_dict, edgelist_df_dict, allcells)
				data_load_state.text('Processing %s networks...done.' %nn1)


				data_load_state = status_msg.text('Reading precomputed PCPs (on %s networks) for all perturbations...' % nn1)
				pcp_perturbation_df_dict = pickle_open_tocache(proj_path + '%s_pcp_perturbation_df_dict.pickle' % gpnet_dict[nn1])
				data_load_state.text('Reading precomputed PCPs (on %s networks) for all perturbations...done.' % nn1)

				data_load_state = status_msg.text('Calculating PCP for input gene signature...')
				pcp_geneset_df = PCP_geneset(geneset_up, geneset_dn, Enrichr_GEO_disease_human_up, Enrichr_GEO_disease_human_dn, geneset_name=geneset_name)
				data_load_state.text('Calculating PCP for input gene signature...done.')

				data_load_state = status_msg.text('Running Pathopticon...(this may take a few minutes)')
				pacos_spearman_rho_df, pacos_spearman_pval_df = PACOS(pcp_geneset_df, pcp_perturbation_df_dict, alldrugs, allcells, 
																	  threshold=10, geneset_name=geneset_name)
				pacos_spearman_rho_all_ranked = rankall_PACOS(pacos_spearman_rho_df, pacos_spearman_pval_df)
				pacos_tool_merged_df = combine_PACOS_tool(pacos_spearman_rho_all_ranked, QRZC_activityStats_nooutliers_df_besttool, rr)
				model_auroc_df_dict, model_auprc_df_dict = PACOS_cell_AUC(pacos_tool_merged_df, allcells, geneset_targets, models=[PACOS_model])
				data_load_state.text('Running Pathopticon...done.')

				data_load_state = status_msg.text('Performing cell line-specific randomization...(this may take a few minutes)')
				PACOS_rand_params={'randomize': True, 'Nrand': nn2}
				PACOS_output_args={'PACOS_output_write_csv': False, 'PACOS_output_csv_dir':None, 'PACOS_output_csv_name':None}
				PCP_perturbation_args={'run_PCP_perturbation':False, 'PCP_perturbation_readpath':None, 
													 'PCP_perturbation_pickle_writedir':None, 'PCP_perturbation_pickle_writename':None}
				
				rand_model_auroc_df_dict, rand_model_auprc_df_dict = PACOS_cell_AUC(pacos_tool_merged_df, allcells, 
																					geneset_targets, models=[PACOS_model],
																				   rand_params={'randomize': PACOS_rand_params['randomize'], 
																								'Nrand': PACOS_rand_params['Nrand']})

				auroc_emp_pval_df, pacos_nested_df = PACOS_nested_prioritization(pacos_tool_merged_df, model_auroc_df_dict, rand_model_auroc_df_dict, 
																				 PACOS_model, 
																				 write_csv=PACOS_output_args['PACOS_output_write_csv'], 
																				 csv_dir=PACOS_output_args['PACOS_output_csv_dir'], 
																				 csv_name=PACOS_output_args['PACOS_output_csv_name'])   
				data_load_state.text('Performing cell line-specific randomization...done.')	


				expander_results = st.expander(label='Pathopticon nested prioritization results')
				with expander_results:
					results_cont = st.container()

				results_cont.write(pacos_nested_df)   		
				filename = '%s_Pathopticon.csv' % geneset_name					
				csv_b64 = base64.b64encode(pacos_nested_df.to_csv(index=False).encode()).decode()
				results_cont.markdown(f'<a href="data:file/csv;base64,{csv_b64}" download={filename}>Download Pathopticon Results</a>', unsafe_allow_html=True)
			
			
				st.session_state['geneset_name'] = geneset_name		
				st.session_state['pacos_nested_df'] = pacos_nested_df
				st.session_state['pcp_perturbation_df_dict'] = pcp_perturbation_df_dict
				st.session_state['pcp_geneset_df'] = pcp_geneset_df
				st.session_state['diG_dict'] = diG_dict
				st.session_state['geneset_up'] = geneset_up
				st.session_state['geneset_dn'] = geneset_dn
				st.session_state['Enrichr_GEO_disease_human_up'] = Enrichr_GEO_disease_human_up
				st.session_state['Enrichr_GEO_disease_human_dn'] = Enrichr_GEO_disease_human_dn
	
			else:
				st.error('The number of true positives (i.e. drugs targeting the input genes) is too small to reliably calculate AUROCs. \
							Please consider expanding the input genesets.')			
	
			
elif side_radio == 'Inspect Pathopticon Results':

	if 'pacos_nested_df' in st.session_state:
		
		expander1 = st.expander(label='Choosing a cell line:')
		with expander1:
			'Hover your mouse over bars to see the AUROC value and AUROC empirical p-value \
			corresponding to that cell line. Red bars indicate cell lines with an AUROC empirical p-value < 0.05.'
	
		### plot cell type bar graph ### 	
		pacos_auroc_sorted = st.session_state['pacos_nested_df'][['Cell_type', 'AUROC', 'AUROC_emp_pval']].drop_duplicates().sort_values('AUROC', ascending=False)

		colors = ['red' if o<=0.05 else 'mistyrose' for o in pacos_auroc_sorted['AUROC_emp_pval']]
		source_data2 = {'Cell_type': pacos_auroc_sorted['Cell_type'], 'AUROC': pacos_auroc_sorted['AUROC'], 
						'AUROC_emp_pval': pacos_auroc_sorted['AUROC_emp_pval'], 'colors': colors} 
		source_s2 = ColumnDataSource(data=source_data2)
		hover2 = HoverTool(tooltips=[('Cell type', '@Cell_type'), ('AUROC', '@AUROC'), ('AUROC emp. p-value', '@AUROC_emp_pval')])
		s2 = figure(x_range=pacos_auroc_sorted['Cell_type'], plot_width=800, plot_height=300, tools='save, wheel_zoom, pan, reset')
		s2.vbar(x='Cell_type', top='AUROC', source=source_s2, color='colors', width=0.8, alpha=0.8)
		s2.xaxis.major_label_orientation = "vertical"
		s2.add_tools(hover2) 
		s2.background_fill_alpha = 0
		s2.border_fill_alpha = 0
		s2.yaxis.axis_label = 'AUROC'
		#s2.xaxis.axis_label = 'Cell type'
		#s2.xaxis.axis_label_text_font_style = "normal"
		# st.bokeh_chart(column(s1, s2))
		st.bokeh_chart(s2, use_container_width=True)

		expander2 = st.expander(label='Accessing perturbation rankings:')
		with expander2:
			'Choose a cell type in the dropdown menu below or start typing the name of the cell line.  \n \
			:bulb: **Tip:** You can sort by any column by clicking on the name of the column.'
			

		allcells = np.sort(st.session_state['pacos_nested_df']['Cell_type'].unique())
		chosen_cell = st.selectbox('Filter by cell type', np.concatenate([['All'], allcells]))
	
		if chosen_cell == 'All':
			#st.subheader('PACOS nested')
			st.write(st.session_state['pacos_nested_df'])   		
		else:
			#st.subheader('PACOS nested')
			pacos_nested_df_filt = st.session_state['pacos_nested_df'][st.session_state['pacos_nested_df']['Cell_type']==chosen_cell]
			st.write(pacos_nested_df_filt)   


			expander3 = st.expander(label='Choosing a perturbation for further exploration:')
			with expander3:
				'Choose a drug in the dropdown menu below. Then hover over the red dots to see each pathophenotype with its corresponding \
				SCS_gd (geneset-disease) and SCS_cpd(perturbation-disease given the cell line) value. Check the "show fit" box to show the linear fit.  \n \
				:bulb: **Tip:** Perturbations are ordered by their PACOS or PACOS-Tool value. If you need to quickly access a drug that is not top-ranked, simply start typing the name of the drug.'

	
			chosen_drug = st.selectbox('Choose a drug to inspect', pacos_nested_df_filt['Pert_iname'])
		
			temp_common_diseases = list(set(st.session_state['pcp_geneset_df'].loc[st.session_state['geneset_name']][~pd.isnull(st.session_state['pcp_geneset_df'].loc[st.session_state['geneset_name']])].index) &\
										set(st.session_state['pcp_perturbation_df_dict'][chosen_drug].loc[chosen_cell][~pd.isnull(st.session_state['pcp_perturbation_df_dict'][chosen_drug].loc[chosen_cell])].index))
			#temp_common_diseases_short = [s.split(' human')[0] for s in temp_common_diseases]
			
			source_data3 = {'SCS_gd': st.session_state['pcp_geneset_df'].loc[st.session_state['geneset_name']][temp_common_diseases].values,
							'SCS_cpd': st.session_state['pcp_perturbation_df_dict'][chosen_drug].loc[chosen_cell][temp_common_diseases].values, 
							'Pathophenotype': temp_common_diseases}  
			source_s3 = ColumnDataSource(data=source_data3)
					
			hover3 = HoverTool(tooltips=[('SCS_gd', '@SCS_gd'), ('SCS_cpd', '@SCS_cpd'), ('Pathophenotype', '@Pathophenotype')])
			s3 = figure(plot_width=650, plot_height=500, tools='save, wheel_zoom, pan, reset')
			s3.circle(x='SCS_gd', y='SCS_cpd', source=source_s3, size=20, color='red', alpha=0.25)
			s3.add_tools(hover3) 
			s3.background_fill_alpha = 0
			s3.border_fill_alpha = 0
			s3.yaxis.axis_label = 'SCS(perturbation-disease)'
			s3.xaxis.axis_label = 'SCS(input geneset-disease)'
			if st.checkbox('Show fit'):
				slope, intercept = np.polyfit(source_data3['SCS_gd'], source_data3['SCS_cpd'], 1)
				linx = np.linspace(min(source_data3['SCS_gd']), max(source_data3['SCS_gd']), 2)
				liny = intercept + slope*linx
				s3.line(x=linx, y=liny, color='red')

			st.bokeh_chart(s3,  use_container_width=True)							
			
			chosen_disease = st.selectbox('Choose a pathophenotype to inspect', temp_common_diseases)
			disease_common_up = set(st.session_state['Enrichr_GEO_disease_human_up'][chosen_disease]) & \
									(set(list(st.session_state['diG_dict'][chosen_cell].successors(chosen_drug))) | \
									st.session_state['geneset_up'] | st.session_state['geneset_dn'])
			disease_common_dn = set(st.session_state['Enrichr_GEO_disease_human_dn'][chosen_disease]) & \
									(set(list(st.session_state['diG_dict'][chosen_cell].successors(chosen_drug))) | \
									st.session_state['geneset_up'] | st.session_state['geneset_dn'])
			input_common_up = st.session_state['geneset_up'] & (set(st.session_state['Enrichr_GEO_disease_human_up'][chosen_disease]) | \
																set(st.session_state['Enrichr_GEO_disease_human_dn'][chosen_disease]) | \
																set(list(st.session_state['diG_dict'][chosen_cell].successors(chosen_drug))))
			input_common_dn = st.session_state['geneset_dn'] & (set(st.session_state['Enrichr_GEO_disease_human_up'][chosen_disease]) | \
																set(st.session_state['Enrichr_GEO_disease_human_dn'][chosen_disease]) | \
																set(list(st.session_state['diG_dict'][chosen_cell].successors(chosen_drug))))			
									
			disease_spoke = generate_disease_spoke(chosen_disease, disease_common_up, disease_common_dn)
			input_spoke = generate_input_spoke(st.session_state['geneset_name'], input_common_up, input_common_dn)
					
			neighbors1 = list(st.session_state['diG_dict'][chosen_cell].successors(chosen_drug))
			neighbors2 = []
			for i in neighbors1:
				neighbors2.extend(st.session_state['diG_dict'][chosen_cell].predecessors(i))
			neighbors2 = np.unique(np.array(neighbors2))
			if st.checkbox('Include 2nd neighbors of selected drug'):
				subG = st.session_state['diG_dict'][chosen_cell].subgraph(np.concatenate([[chosen_drug], neighbors1, neighbors2]))
			else:
				subG = st.session_state['diG_dict'][chosen_cell].subgraph(np.concatenate([[chosen_drug], neighbors1]))
			
			nt = Network('650px', '650px', directed=True, bgcolor='k', font_color='m')

			for node in list(input_spoke.nodes()):
				if input_spoke.nodes[node]['node_type'] == 'input_geneset_name':
					nt.add_node(node, color=input_spoke.nodes[node]['node_color'], size=input_spoke.nodes[node]['node_size'], physics=False, y=0.0)				
				else:
					nt.add_node(node, color=input_spoke.nodes[node]['node_color'], size=input_spoke.nodes[node]['node_size'])
			for e1, e2 in list(input_spoke.edges()):
				nt.add_edge(e1, e2, color=input_spoke.edges[e1, e2]['edge_color'])			
			for node in list(disease_spoke.nodes()):
				if disease_spoke.nodes[node]['node_type'] == 'disease_name':
					nt.add_node(node, color=disease_spoke.nodes[node]['node_color'], size=disease_spoke.nodes[node]['node_size'], physics=False, y=500.0)				
				else:
					nt.add_node(node, color=disease_spoke.nodes[node]['node_color'], size=disease_spoke.nodes[node]['node_size'])
			for e1, e2 in list(disease_spoke.edges()):
				nt.add_edge(e1, e2, color=disease_spoke.edges[e1, e2]['edge_color'])				
			for node in list(subG.nodes()):
				if node == chosen_drug:
					nt.add_node(node, color='greenyellow', size=20, physics=False, x=0.0)	
				else:
					nt.add_node(node, color=subG.nodes[node]['node_color'], size=subG.nodes[node]['node_size'])			
			for e1, e2 in list(subG.edges()):
				nt.add_edge(e1, e2, color=subG.edges[e1, e2]['edge_color'])

			nt.show('nx.html')
			f_html = open('nx.html', 'r', encoding='utf-8')
			source_html =  f_html.read()
			components.html(source_html, height=900, width=900)
					
	else:
		st.error('Please run Pathopticon first.')
