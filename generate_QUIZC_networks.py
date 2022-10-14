import sys 
import os
import numpy as np
import pandas as pd
import random
import time
import pickle
from multiprocessing import cpu_count, Pool


if __name__ == "__main__":

	parser=argparse.ArgumentParser()
	
	parser.add_argument("--L1000_path", help="Path to the folder where the L1000 metadata is located.")
	parser.add_argument("--L1000_feather_path", help="Path to the folder where the .feather files for Level 4 L1000 data are located.")
	parser.add_argument("--quantile", help="The quantile to be used as the consensus threshold. Default is 0.75.")
	parser.add_argument("--ncores", type=int, nargs='?', default=cpu_count(), const=1, help="The number of cores to be used in the computation.")
	parser.add_argument("--output_path", help="Path to the folder in which the network files will be written.")
	
	args=parser.parse_args()
	
	
	nCores = args.ncores
	print('Number of threads on this machine: %s; using %s ' % (cpu_count(), nCores))

	L1000_gene_info = pd.read_csv(args.L1000_path + 'GSE92742_Broad_LINCS_gene_info.txt', sep='\t', low_memory=False)
	L1000_cell_info = pd.read_csv(args.L1000_path + 'GSE92742_Broad_LINCS_cell_info.txt', sep='\t', low_memory=False)
	L1000_inst_info = pd.read_csv(args.L1000_path + 'GSE92742_Broad_LINCS_inst_info.txt', sep='\t', low_memory=False)
	L1000_pert_info = pd.read_csv(args.L1000_path + 'GSE92742_Broad_LINCS_pert_info.txt', sep='\t', low_memory=False)
	L1000_sig_info = pd.read_csv(args.L1000_path + 'GSE92742_Broad_LINCS_sig_info.txt', sep='\t', low_memory=False)

	allcells = sorted([x.split('.')[0].split('_')[-1] for x in os.listdir(args.L1000_feather_path)])

	cell_pos_edges_dict = {}
	cell_neg_edges_dict = {}
	cell_drugs_dict = {}

	for cell_type in allcells:
		level4 = pd.read_feather(args.L1000_feather_path + 'Level4_all_parsed_%s.feather' % cell_type)
		level4_inst_info = pd.merge( L1000_inst_info, level4.set_index('rid').T.sort_index(), right_index=True, left_on='inst_id')
		level4_pert_iname = level4_inst_info.set_index('pert_iname')[level4['rid']]
		level4_pert_iname_z = ((level4_pert_iname - level4_pert_iname.mean(axis=0))) / level4_pert_iname.std(axis=0)
		level4_pert_iname_z = level4_pert_iname_z.reset_index()
		level4_pert_iname_z_grouped = level4_pert_iname_z.groupby('pert_iname')
	
		def QUIZC_pos(d):
			#print(d)
			pert_df = level4_pert_iname_z_grouped.get_group(d)

			part_a_pos = (pert_df[pert_df.columns[1:]] >= 1.96).sum(axis=0)
			part_b_pos = (pert_df[pert_df.columns[1:]] >= np.quantile(pert_df[pert_df.columns[1:]], args.quantile, axis=0)).sum(axis=0)

			if len((part_a_pos > part_b_pos)[(part_a_pos > part_b_pos)]) > 0:
				pos_targets = (part_a_pos > part_b_pos)[(part_a_pos > part_b_pos)].index.values

				return pos_targets
		
		def QUIZC_neg(d):
			pert_df = level4_pert_iname_z_grouped.get_group(d)

			part_a_neg = (pert_df[pert_df.columns[1:]] <= -1.96).sum(axis=0)
			part_b_neg = (pert_df[pert_df.columns[1:]] <= np.quantile(pert_df[pert_df.columns[1:]], 1.0-args.quantile, axis=0)).sum(axis=0)

			if len((part_a_neg > part_b_neg)[(part_a_neg > part_b_neg)]) > 0:
				neg_targets = (part_a_neg > part_b_neg)[(part_a_neg > part_b_neg)].index.values        

				return neg_targets
		
		cell_drugs = sorted(list(set(level4_pert_iname_z['pert_iname'])))
		start = time.time()
		p = Pool(nCores)
		cell_pos_edges_dict[cell_type] = p.map(QUIZC_pos, cell_drugs)
		cell_neg_edges_dict[cell_type] = p.map(QUIZC_neg, cell_drugs)
		cell_drugs_dict[cell_type] = cell_drugs
	
		print(cell_type, time.time() - start)    
	
	with open(args.output_path + 'cell_pos_edges_dict_%s.pickle' % args.quantile, 'wb') as fp:
		pickle.dump(cell_pos_edges_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
	with open(args.output_path + 'cell_neg_edges_dict_%s.pickle' % args.quantile, 'wb') as fp:
		pickle.dump(cell_neg_edges_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)
	with open(args.output_path + 'cell_drugs_dict_%s.pickle' % args.quantile, 'wb') as fp:
		pickle.dump(cell_drugs_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)