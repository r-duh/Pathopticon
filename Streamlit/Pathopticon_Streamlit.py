from Pathopticon_Streamlit_functions import *
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

input_paths_dict = input_paths(proj_path)


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
			QUIZC_cids_inchi_smiles = pandas_read_csv_tocache(input_paths_dict['QUIZC_cids_inchi_smiles_path'])
			TTD_drugs, TTD_InChI2CID, TTD_drugs_CID, TTD_drug_target, TTD_targID_dict = import_TDD(input_paths_dict['TTD_drugs_path'], 
																								   input_paths_dict['TTD_InChI2CID_path'],
																								   input_paths_dict['TTD_drug_target_path'], 
																								   input_paths_dict['TTD_target_path'])
			geneset_targets = get_geneset_targets(geneset_up, geneset_dn, TTD_targID_dict, TTD_drug_target, TTD_drugs_CID, QUIZC_cids_inchi_smiles)			

			if len(geneset_targets) >= 5:
				
				L1000_gene_info, L1000_cell_info, L1000_inst_info, L1000_pert_info, L1000_sig_info = import_L1000_metadata(input_paths_dict['L1000_gene_info_path'], 
																														   input_paths_dict['L1000_cell_info_path'],
																														   input_paths_dict['L1000_inst_info_path'], 
																														   input_paths_dict['L1000_pert_info_path'], 
																														   input_paths_dict['L1000_sig_info_path'])   				
				
				QUIZC_activityStats_nooutliers_df_besttool = pandas_read_csv_tocache(input_paths_dict['tool_path'])			



				Enrichr_GEO_disease_human_up, Enrichr_GEO_disease_human_dn = import_Enrichr_GEO(input_paths_dict['Enrichr_GEO_up_path'], 
																								input_paths_dict['Enrichr_GEO_dn_path'])

				data_load_state.text('Loading input data...done.')



				data_load_state = status_msg.text('Processing %s networks...' % nn1)
				if nn1 == 'QUIZ-C':
					edgelist_df_dict, nodelist_df_dict, allcells, allgenes, alldrugs, allnodes = process_QUIZC_output(input_paths_dict['pos_edges_dict_path'], 
																												 input_paths_dict['neg_edges_dict_path'],
																												 input_paths_dict['drugs_dict_path'],
																												 L1000_gene_info)
				else:
					edgelist_df_dict, nodelist_df_dict, allcells, allgenes, alldrugs, allnodes = import_CD_MODZ_networks(input_paths_dict['%s_networks_path' % nn1])				
				diG_dict = generate_diG_dict(nodelist_df_dict, edgelist_df_dict, allcells)
				data_load_state.text('Processing %s networks...done.' %nn1)


				data_load_state = status_msg.text('Reading precomputed PCPs (on %s networks) for all perturbations...' % nn1)
				pcp_perturbation_df_dict = PCP_perturbation_alldrugs(alldrugs, allcells, edgelist_df_dict, 
                                                         Enrichr_GEO_disease_human_up, Enrichr_GEO_disease_human_dn,
                                                         proj_path, method_name=nn1, return_output=True) 
				data_load_state.text('Reading precomputed PCPs (on %s networks) for all perturbations...done.' % nn1)

				data_load_state = status_msg.text('Calculating PCP for input gene signature...')
				pcp_geneset_df = PCP_geneset(geneset_up, geneset_dn, Enrichr_GEO_disease_human_up, Enrichr_GEO_disease_human_dn, geneset_name=geneset_name)
				data_load_state.text('Calculating PCP for input gene signature...done.')

				data_load_state = status_msg.text('Running Pathopticon...(this may take a few minutes)')
				pacos_tool_merged_df = PACOS(pcp_geneset_df, pcp_perturbation_df_dict, alldrugs, allcells,
											 QUIZC_activityStats_nooutliers_df_besttool, proj_path, method_name=nn1, geneset_name=geneset_name, 
											 r=rr, threshold=PACOS_Spearman_threshold, tqdm_off=True)
				model_auroc_df_dict, model_auprc_df_dict = PACOS_cell_AUC(pacos_tool_merged_df, allcells, geneset_targets, models=[PACOS_model])
				data_load_state.text('Running Pathopticon...done.')

				data_load_state = status_msg.text('Performing cell line-specific randomization...(this may take a few minutes)')
				rand_model_auroc_df_dict, rand_model_auprc_df_dict = PACOS_cell_AUC_randomize(pacos_tool_merged_df, allcells, 
																							  geneset_targets, 
																							  models=[PACOS_model],
																							  auc_params={'auc_threshold':5, 'binsize':1}, 
																							  Nrand=nn2, tqdm_off=True) 

				emp_pval_df, pacos_nested_df = PACOS_nested_prioritization(pacos_tool_merged_df, model_auroc_df, rand_model_auroc_df_dict)    
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

			nt.show_buttons(filter_=["physics"])
			nt.show('nx.html')
			f_html = open('nx.html', 'r', encoding='utf-8')
			source_html =  f_html.read()
			components.html(source_html, height=900, width=900)
			
					
	else:
		st.error('Please run Pathopticon first.')
