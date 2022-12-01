# :eyeglasses: Pathopticon: Linking gene and perturbation signatures through their **pathop**heno**t**yp**ic** c**on**gruity

**Citation:** Arda Halu, Julius Decano, Joan Matamalas, Mary Whelan, Takaharu Asano, Namitra Kalicharran, 
Sasha Singh, Joseph Loscalzo, Masanori Aikawa, "Integrating pharmacogenomics and cheminformatics with diverse 
disease phenotypes for cell type-guided drug discovery"

***

### Welcome! 
**What is Pathopticon?** Pathopticon is a drug discovery platform that integrates pharmacogenomic and cheminformatic 
data with diverse disease phenotypes to predict cell type-specific candidate drugs for a given disease context.

**What does Pathopticon do?** Below, you can find a figure summarizing Pathopticon. Briefly, Pathopticon first creates 
cell type-specific gene-perturbation networks using a statistical approach that we call the QUantile-based Instance 
Z-score Consensus (QUIZ-C) method. Pathopticon then integrates these cell type-specific gene-perturbation networks 
with an "input" gene signature provided by the user (typically differentially up- and down-regulated genes from various 
omics experiments). Specifically, it measures the agreement between input and perturbation signatures within a global 
network of diverse disease phenotypes using what we call the PAthophenotypic COngruity Score (PACOS). After combining 
PACOS with pharmacological activity data, Pathopticon performs a nested prioritization that identifies the most suitable 
cell lines for an input gene signature, followed by the perturbations whose transcriptomic response best aligns with 
that of the input signature within each cell line. 

**Can I see these cell type-specific networks?** Absolutely, and interactively, too! Just click on the "Inspect Networks" 
option in the left panel and it will take you to a page in which you can choose any of the 60 available cell types. You 
can then use the dropdown menu to choose drugs or genes of interest, which will be highlighted in the network. You can 
pan and zoom in and out of the network, and choose and highlight _multiple_ genes/drugs and their connections using the 
Shift key to create your own subnetworks.

**How do I run Pathopticon?** It is very simple -- if you select the "Run Pathopticon" option in the left panel and 
enter a set of up- and down-regulated genes, Pathopticon will give you a ranked list of drugs and cell lines, which you 
can download. More details on the individual parameters can be found under "Run Pathopticon."

**How do I interpret the results?** The "Inspect Pathopticon Results" option in the left panel will guide you on which cell 
lines were significant, and let you access the ranked list of drugs within each cell line. You can then choose a drug and 
inspect the diseases with similar gene signatures, and finally focus on a disease to see a subnetwork consisting of the input 
signature, the chosen disease's gene signature and the chosen drug.		

