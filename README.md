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

![Overview of the Pathopticon framework](https://github.com/r-duh/Pathopticon/blob/main/Pathopticon_overview_fig.png?raw=true)

***

### Description of the Pathopticon Parameters

**Model:** Tells Pathopticon which scoring model to use. The two options are (i) PACOS (pharmacogenomic data) only and 
(ii) PACOS combined with Tool scores (pharmacogenomic data combined with cheminformatic data). The default option is the 
latter.

**Effect:** Tells Pathopticon which kind of perturbations to look for. The two options are (i) Repress, which ranks drugs 
according to how well they repress, or reverse, the input gene signature, and (ii) Enhance, which ranks drugs according to how well they enhance, or mimic, the 
input gene signature. The default option is Repress.

**Gene-perturbation network:** Tells Pathopticon which type of gene-perturbation network to use as the underlying network 
on which to make predictions. The three options are (i) QUIZ-C, (ii) MODZ and (iii) CD. The details on these networks 
can be found our manuscript. The default is QUIZ-C, which is the method we propose in this study.

**r-value:** The weight Pathopticon uses to balance the effect of pharmacogenomic data and cheminformatic data when combining
the two. This option is required for the combined PACOS-Tool model. r < 1 biases the combined score in favor of Tool scores 
and r > 1 biases the combined score in favor of PACOS scores. Based on our sensitivity analyses, the default value is 2.

**Number of randomizations:** The number of randomizations with which to calculate empirical p-values. The default value is 
400, which results in an empirical p-value resolution of 1/400 = 0.0025. Higher values will increase p-value resolution, but 
will also increase computation time.

**Input gene signature name:** The name defined by the user to describe the input signature. This will also be the name of the
output .csv file, so we recommend using underscores ("_") instead of spaces.

**Input gene signature - Up/Down genes:** The list of genes to be used as the input signature, which will be provided by the user. These are 
typically up- and down-regulated genes from omics experiments.

***

## Installing and running the Pathopticon Streamit app

We have two options: 

#### 1) Run Pathopticon Streamlit app as a Docker container (requires [Docker installation](https://docs.docker.com/get-docker/))
This option is somewhat slower to run but has the advantage of not depending on the specific package environment and operating system of the user.

- Once Docker is installed, the Pathopticon image "rduh/pathopticon-streamlit:slim" can either be pulled from Docker Hub
```
docker pull rduh/pathopticon-streamlit:slim
```
or loaded from the [.tar file]() using
```
docker load --input pathopticon_streamlit_slim.tar
```

- [Download]() the Pathopticon folder to be mounted as a volume to the Docker container. This local folder (i.e., located in the user's machine), named /Pathopticon_Streamlit_Docker_mount/, will act as the main folder in which Pathopticon Streamlit app's container will read and write files.

- To run the Pathopticon Streamlit app as a Docker container, type in the below command in the terminal. /path/to/Pathopticon_Streamlit_Docker_mount/ is where the folder you downloaded above is located in your computer. 
```
docker run -it -v /path/to/Pathopticon_Streamlit_Docker_mount/:/Pathopticon_Streamlit_app/ -p 8501:8501 rduh/pathopticon-streamlit:slim
```

- Finally, to see the Streamlit app, go to your browser and enter the address that appears in your terminal (It looks like this: "You can now view your Streamlit app in your browser. URL: `http://0.0.0.0:8501` "). So, typically `http://0.0.0.0:8501`. If you have more than one Streamlit instance running, this can be `http://0.0.0.0:8502`, `http://0.0.0.0:8503`, and so on.


#### 2) Run the Pathopticon Streamlit app script directly (requires conda to be installed). 
This is the faster option but requires familiarity with creating environments and running scripts.
- If not already done so, install [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) on your system. Next, type the following commands in your command prompt, or terminal, in the following order to set up and run Pathopticon.
- Using Conda, create a new environment named "pathopticon_streamlit" (or any name of your choosing) for the Pathopticon Streamlit app and install in it the dependencies needed by the app using the Pathopticon_Streamlit_requirements.yml file in /path/to/Pathopticon_Streamlit_Docker_mount/ (also available in the project GitHub [page](https://github.com/r-duh/Pathopticon/blob/main/Streamlit/Pathopticon_Streamlit_environment.yml)): 
```
conda env create -n pathopticon_streamlit -f /path/to/Pathopticon_Streamlit_Docker_mount/Pathopticon_Streamlit_environment.yml
```
- Activate the newly created Conda environment:
```
conda activate pathopticon_streamlit
```
- You have now created a new conda environment and installed in it all the packages the Pathopticon Streamlit app needs to run. The only remaining step is to run it. First, navigate to /path/to/Pathopticon_Streamlit_Docker_mount/
```
cd /path/to/Pathopticon_Streamlit_Docker_mount/
```
Then, run the streamlit run by typing the below command (note that it has the additional --proj_path flag, which needs to be set to the Pathopticon Streamlit app directory
```
streamlit run Pathopticon_Streamlit.py -- --proj_path=/path/to/Pathopticon_Streamlit_Docker_mount/
```

- Finally, to see the Streamlit app, go to your browser and enter the address that appears in your terminal (It looks like this: "You can now view your Streamlit app in your browser. URL: `http://0.0.0.0:8501` "). So, typically `http://0.0.0.0:8501`. If you have more than one Streamlit instance running, this can be `http://0.0.0.0:8502`, `http://0.0.0.0:8503`, and so on.

***

# Contact:
Created and maintained by Arda Halu. For requests for assistance, email arda.halu@channing.harvard.edu.
