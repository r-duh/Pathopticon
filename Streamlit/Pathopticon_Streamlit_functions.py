import pandas as pd
import pickle
from tqdm import tqdm
import numpy as np
import os
import networkx as nx
import scipy.stats as st


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

@st.cache(show_spinner=False, allow_output_mutation=True)
def input_paths(proj_path):

    input_paths_dict = {'QUIZC_cids_inchi_smiles_path': proj_path + 'Pathopticon_intermediary_outputs/' + 'QUIZC_drug_CIDs_inchi_smiles.csv',
                        'tool_path': proj_path + 'Pathopticon_intermediary_outputs/' + 'QUIZC_activityStats_nooutliers_df_besttool.csv',
                        'benchmark_genesets_path': proj_path + 'Pathopticon_intermediary_outputs/' + 'benchmark_genesets.csv',
                        
                        'L1000_gene_info_path': proj_path + 'Pathopticon_external_data/L1000/' + 'GSE92742_Broad_LINCS_gene_info.txt',
                        'L1000_cell_info_path': proj_path + 'Pathopticon_external_data/L1000/' + 'GSE92742_Broad_LINCS_cell_info.txt',
                        'L1000_inst_info_path': proj_path + 'Pathopticon_external_data/L1000/' + 'GSE92742_Broad_LINCS_inst_info.txt',
                        'L1000_pert_info_path': proj_path + 'Pathopticon_external_data/L1000/' + 'GSE92742_Broad_LINCS_pert_info.txt',
                        'L1000_sig_info_path': proj_path + 'Pathopticon_external_data/L1000/' + 'GSE92742_Broad_LINCS_sig_info.txt',

                        'pos_edges_dict_path': proj_path + 'QUIZC75/' + 'cell_pos_edges_dict_75.pickle',
                        'neg_edges_dict_path': proj_path + 'QUIZC75/' + 'cell_neg_edges_dict_75.pickle',
                        'drugs_dict_path': proj_path + 'QUIZC75/' + 'cell_drugs_dict_75.pickle',

                        'cgp_dir': proj_path + 'Pathopticon_external_data/MSigDB/' + 'c2.cgp.v7.1.symbols.gmt',

                        'Enrichr_GEO_up_path': proj_path + 'Pathopticon_external_data/Enrichr/' + 'Disease_Perturbations_from_GEO_up.txt',
                        'Enrichr_GEO_dn_path': proj_path + 'Pathopticon_external_data/Enrichr/' + 'Disease_Perturbations_from_GEO_down.txt',

                        'TTD_drugs_path': proj_path + 'Pathopticon_external_data/TTD/' + 'P1-02-TTD_drug_download.txt',
                        'TTD_InChI2CID_path': proj_path + 'Pathopticon_external_data/TTD/' + 'TTD_drugs_InChI2CID.txt',
                        'TTD_drug_target_path': proj_path + 'Pathopticon_external_data/TTD/' + 'P1-07-Drug-TargetMapping.csv',
                        'TTD_target_path': proj_path + 'Pathopticon_external_data/TTD/' + 'P1-01-TTD_target_download.txt',
                        'MODZ_networks_path': proj_path + 'MODZ_networks/',
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
def process_QUIZC_output(pos_edges_dict_path, neg_edges_dict_path, drugs_dict_path, L1000_gene_info):

    with open(pos_edges_dict_path, 'rb') as handle:
        cell_pos_edges_dict = pickle.load(handle)
    with open(neg_edges_dict_path, 'rb') as handle:
        cell_neg_edges_dict = pickle.load(handle)
    with open(drugs_dict_path, 'rb') as handle:
        cell_drugs_dict = pickle.load(handle)     
    
    edgelist_df_dict = {}
    nodelist_df_dict = {}
    for c in tqdm(cell_drugs_dict.keys(), position=0, leave=True):

        pos_edges_df = pd.DataFrame(np.array([(i, j) for i, j in zip(cell_drugs_dict[c], cell_pos_edges_dict[c])])).explode(1)
        neg_edges_df = pd.DataFrame(np.array([(i, j) for i, j in zip(cell_drugs_dict[c], cell_neg_edges_dict[c])])).explode(1)

        pos_edges_df = pos_edges_df[~pd.isnull(pos_edges_df[1])]
        neg_edges_df = neg_edges_df[~pd.isnull(neg_edges_df[1])]

        pos_edges_df = pos_edges_df.reset_index(drop=True).rename(columns={0: 'Drug', 1: 'Target'})
        neg_edges_df = neg_edges_df.reset_index(drop=True).rename(columns={0: 'Drug', 1: 'Target'})

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

@st.cache(show_spinner=False, allow_output_mutation=True)
def generate_diG_dict(nodelist_df_dict, edgelist_df_dict, allcells):

    diG_dict = {}
    for c in tqdm(allcells['Cell_type'], position=0, leave=True):
        diG_dict[c] = nx.DiGraph()
        for ix, node in enumerate(nodelist_df_dict[c]['Id'].values):
            if nodelist_df_dict[c].iloc[ix]['Type']=='Drug':
                diG_dict[c].add_node(node, node_type='Drug', node_color='limegreen', node_size=150)
            elif nodelist_df_dict[c].iloc[ix]['Type']=='Gene':
                diG_dict[c].add_node(node, node_type='Gene', node_color='darkslateblue', node_size=50) 

        for edge in edgelist_df_dict[c].values:       
            if edge[2]=='Up':
                diG_dict[c].add_edge(edge[0], edge[1], edge_type='Up', edge_color='red')
            elif edge[2]=='Down':
                diG_dict[c].add_edge(edge[0], edge[1], edge_type='Down', edge_color='deepskyblue')
                
    return diG_dict

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

    PCP_geneset_df = pd.DataFrame(index=[geneset_name], columns=sorted(disease_sig_up_dict.keys()), dtype='float')
    
    for d in disease_sig_up_dict.keys():
        disease_up = set(disease_sig_up_dict[d])
        disease_dn = set(disease_sig_dn_dict[d])    
        
        PCP_geneset_df.at[geneset_name, d] = SCS(geneset_up, geneset_dn, disease_up, disease_dn)
        
    return PCP_geneset_df

@st.cache(show_spinner=False, allow_output_mutation=True)
def PCP_perturbation(perturbation, allcells, edgelist_df_dict, disease_sig_up_dict, disease_sig_dn_dict):

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
def PCP_perturbation_alldrugs(alldrugs, allcells, edgelist_df_dict, disease_sig_up_dict, disease_sig_dn_dict, 
                              proj_path, method_name, return_output=False):
    
    if not os.path.exists(proj_path + method_name + '_pcp_perturbation_df_dict.pickle'):
        
        print('Calculating PCPs for all perturbations...', flush=True)
        PCP_perturbation_df_dict = {}

        for drug in tqdm(alldrugs['Pert_iname'], position=0, leave=True):

            PCP_perturbation_df_dict[drug] = PCP_perturbation(drug, allcells, edgelist_df_dict, disease_sig_up_dict, disease_sig_dn_dict)

        with open(proj_path + method_name + '_pcp_perturbation_df_dict.pickle', 'wb') as fp:
            pickle.dump(PCP_perturbation_df_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
            
        if return_output==True:  
            return PCP_perturbation_df_dict
                
    else:
        
        print('Reading PCPs for all perturbations...', flush=True)
        with open(proj_path + method_name + '_pcp_perturbation_df_dict.pickle', 'rb') as fp:
            PCP_perturbation_df_dict = pickle.load(fp)      
    
        if return_output==True:  
            return PCP_perturbation_df_dict

def calculate_auroc(TPR, FPR):
    
    dTPR = np.concatenate((np.ediff1d(TPR), [0]))
    dFPR = np.concatenate((np.ediff1d(FPR), [0]))
    
    return sum(TPR * dFPR) + sum(dTPR * dFPR)/2

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
def PACOS(PCP_geneset_df, PCP_perturbation_df_dict, alldrugs, allcells, tool_scores, 
          proj_path, method_name, geneset_name, r=2.0, threshold=10, tqdm_off=False, messages=False):

    if not os.path.exists(proj_path +  '%s_%s_r%s.csv' % (method_name, geneset_name, int(r))):
        
        if messages: print('Running Pathopticon...', flush=True)
        PACOS_spearman_rho_df = pd.DataFrame(index=alldrugs['Pert_iname'], columns=allcells['Cell_type'], dtype='float')
        PACOS_spearman_pval_df = pd.DataFrame(index=alldrugs['Pert_iname'], columns=allcells['Cell_type'], dtype='float')
        
        for drug in tqdm(alldrugs['Pert_iname'], position=0, leave=True, disable=tqdm_off):

            for c in PCP_perturbation_df_dict[drug].index:
                common_diseases = set(PCP_geneset_df.loc[geneset_name][~pd.isnull(PCP_geneset_df.loc[geneset_name])].index) &\
                                    set(PCP_perturbation_df_dict[drug].loc[c][~pd.isnull(PCP_perturbation_df_dict[drug].loc[c])].index)
                if len(common_diseases) >= threshold:
                    PACOS_spearman_rho_df.at[drug, c], PACOS_spearman_pval_df.at[drug, c] = st.spearmanr(PCP_geneset_df.loc[geneset_name][common_diseases], 
                                                                      PCP_perturbation_df_dict[drug].loc[c][common_diseases])

        PACOS_spearman_rho_df[pd.isnull(PACOS_spearman_rho_df)] = -666
        sorted_ix = np.argsort(PACOS_spearman_rho_df.values.flatten())
        nonan_len = len(PACOS_spearman_rho_df.values.flatten()[PACOS_spearman_rho_df.values.flatten()!=-666])
        nonan_ranked = PACOS_spearman_rho_df.values.flatten()[np.flip(sorted_ix)[0:nonan_len]]
        nonan_ranked_ix = np.flip(sorted_ix)[0:nonan_len]
        nonan_ranked_i, nonan_ranked_j = np.unravel_index(nonan_ranked_ix, PACOS_spearman_rho_df.shape)

        PACOS_spearman_rho_all_ranked = pd.DataFrame()
        PACOS_spearman_rho_all_ranked['Pert_iname'] = PACOS_spearman_rho_df.index[nonan_ranked_i]
        PACOS_spearman_rho_all_ranked['Cell_type'] = PACOS_spearman_rho_df.columns[nonan_ranked_j]
        PACOS_spearman_rho_all_ranked['PACOS_Spearman_rho'] = np.array([PACOS_spearman_rho_df.at[PACOS_spearman_rho_df.index[i], 
                                                                                                 PACOS_spearman_rho_df.columns[j]] 
                                                                        for i, j in zip(nonan_ranked_i, nonan_ranked_j)])
        PACOS_spearman_rho_all_ranked['PACOS_Spearman_pval'] = np.array([PACOS_spearman_pval_df.at[PACOS_spearman_rho_df.index[i], 
                                                                                                   PACOS_spearman_rho_df.columns[j]] 
                                                                        for i, j in zip(nonan_ranked_i, nonan_ranked_j)]) 


        PACOS_tool_merged_df = pd.merge(PACOS_spearman_rho_all_ranked, tool_scores, left_on='Pert_iname', right_on='Pert_iname', how='left')
        PACOS_tool_merged_df['PACOS_Spearman_rho_reverse'] = -1.0*PACOS_tool_merged_df['PACOS_Spearman_rho']
        PACOS_tool_merged_df['tool_score_imputed'] = PACOS_tool_merged_df['tool score scaled'].fillna(tool_scores['tool score scaled'].median())   
        PACOS_tool_merged_df['PACOS_tool_combined'] = (r*PACOS_tool_merged_df['PACOS_Spearman_rho'] + 
                                                       PACOS_tool_merged_df['tool_score_imputed']) / (r + 1.0)
        PACOS_tool_merged_df['PACOS_tool_combined_reverse'] = (-1.0*r*PACOS_tool_merged_df['PACOS_Spearman_rho'] + 
                                                       PACOS_tool_merged_df['tool_score_imputed']) / (r + 1.0)


        PACOS_tool_merged_df.to_csv(proj_path +  '%s_%s_r%s.csv' % (method_name, geneset_name, int(r)), index=False)
        
        return PACOS_tool_merged_df
        
    else:
        
        PACOS_tool_merged_df = pd.read_csv(proj_path +  '%s_%s_r%s.csv' % (method_name, geneset_name, int(r)))
        
        return PACOS_tool_merged_df

@st.cache(show_spinner=False, allow_output_mutation=True)
def PACOS_cell_AUC(PACOS_tool_merged_df, allcells, positives, 
                   models=['PACOS_Spearman_rho', 'PACOS_Spearman_rho_reverse', 'PACOS_tool_combined', 'PACOS_tool_combined_reverse'],
                   auc_params={'auc_threshold':5, 'binsize':1}):
    
    model_auroc_df = pd.DataFrame(index=allcells['Cell_type'], columns=['%s_AUROC' % m for m in models])
    model_auprc_df = pd.DataFrame(index=allcells['Cell_type'], columns=['%s_AUPRC' % m for m in models])
    
    for m in models:
        
        allcells_sorted = PACOS_tool_merged_df.sort_values(['Cell_type', m], ascending=[True, False])
        
        for c in allcells['Cell_type']:  
            cell_sorted = allcells_sorted[allcells_sorted['Cell_type']==c]['Pert_iname'].values
            model_auroc_df.at[c, '%s_AUROC' % m], model_auprc_df.at[c, '%s_AUPRC' % m] = AUROC_AUPRC(cell_sorted, positives, 
                                                                                                                   auc_threshold=auc_params['auc_threshold'], 
                                                                                                                   binsize=auc_params['binsize'])
        
    return model_auroc_df, model_auprc_df

@st.cache(show_spinner=False, allow_output_mutation=True)
def PACOS_cell_AUC_randomize(PACOS_tool_merged_df, allcells, positives, 
                             models=['PACOS_Spearman_rho', 'PACOS_Spearman_rho_reverse', 'PACOS_tool_combined', 'PACOS_tool_combined_reverse'],
                             auc_params={'auc_threshold':5, 'binsize':1}, 
                             Nrand=200, tqdm_off=False):

    rand_model_auroc_df_dict = {}
    rand_model_auprc_df_dict = {}
    
    for m in models:
        
        cell_geneset_auroc_df = pd.DataFrame(index=allcells['Cell_type'], columns=np.arange(Nrand))
        cell_geneset_auprc_df = pd.DataFrame(index=allcells['Cell_type'], columns=np.arange(Nrand))        

        allcells_sorted = PACOS_tool_merged_df.sort_values(['Cell_type', m], ascending=[True, False])
        
        for nrand in tqdm(np.arange(Nrand), position=0, leave=True, disable=tqdm_off):
            for c in allcells['Cell_type']:  
                cell_sorted = allcells_sorted[allcells_sorted['Cell_type']==c]['Pert_iname'].values
                np.random.shuffle(cell_sorted)
                cell_geneset_auroc_df.at[c, nrand], cell_geneset_auprc_df.at[c, nrand] = AUROC_AUPRC(cell_sorted, positives, 
                                                                                                     auc_threshold=auc_params['auc_threshold'], 
                                                                                                     binsize=auc_params['binsize'])
            
        rand_model_auroc_df_dict[m] = cell_geneset_auroc_df
        rand_model_auprc_df_dict[m] = cell_geneset_auprc_df 
        
    return rand_model_auroc_df_dict, rand_model_auprc_df_dict

@st.cache(show_spinner=False, allow_output_mutation=True)
def PACOS_nested_prioritization(PACOS_tool_merged_df, model_auroc_df, rand_model_auroc_df_dict):
    
    emp_pval_df = pd.DataFrame(index=model_auroc_df.index, columns=['%s_AUROC_emp_pval' % m for m in list(rand_model_auroc_df_dict.keys())])
    
    for m in rand_model_auroc_df_dict.keys():
    
        nonan_cells = model_auroc_df['%s_AUROC' % m][~pd.isnull(model_auroc_df['%s_AUROC' % m])].index
        for c in nonan_cells:
            emp_pval_df.at[c, '%s_AUROC_emp_pval' % m] = ((rand_model_auroc_df_dict[m].loc[c].values > 
                                                             model_auroc_df['%s_AUROC' % m].loc[c]).sum() /
                                                            float(np.shape(rand_model_auroc_df_dict[m])[1]))
    
    model_auroc_pval_df = pd.merge(model_auroc_df, emp_pval_df, left_index=True, right_index=True)
    PACOS_nested_df = pd.merge(PACOS_tool_merged_df, model_auroc_pval_df, left_on='Cell_type', 
                               right_index=True).dropna(subset=['%s_AUROC' % m for m in list(rand_model_auroc_df_dict.keys())])        
        
    return emp_pval_df, PACOS_nested_df

