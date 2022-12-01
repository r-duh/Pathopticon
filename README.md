# :eyeglasses: Pathopticon: Linking gene and perturbation signatures through their **pathop**heno**t**yp**ic** c**on**gruity
***

'''
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
'''

