# Microarray gene expression analysis workflow in R

Gene expression analysis of microarray dataset (raw probe intensities) has following steps: 

* Downloading and importing the dataset. 
* Filtering poor signals
* Normalization
* Identification of differentially expressed genes
* Gene set enrichment analysis

# Preparing the microarray dataset downloaded from NASA GeneLab

Download the microarray dataset (GLDS-3) from the provided url path (https://genelab-data.ndc.nasa.gov/genelab/static/media/dataset/GLDS-3_microarray_E-GEOD-23880.raw.1.zip?version=1) and uncompress it in a local repository. Before starting the analysis, remove samples with accession series (GSM588931-GSM588936 and GSM588940-GSM588945) because these are samples are from larval stages of D. melanogaster.  

# "Always remember that the step for analysis of microarray data from NASA GeneLab is similar to the NCBI-GEO."  

# This workflow

Here, I am demonstrating a step by step pipeline for identification of differentially expressed genes (DEGs) in space flight vs ground control studies of adult D. melanogaster. You'll be required to download all the libraries for running this workflow. Furthermore, provide the correct path to the repository containing downloaded microarray dataset and run the pipeline step by step.

Best wishes!

# Suggestions:

abdussalam@precisionmedicine.pk

salam.nurzai@gmail.com
