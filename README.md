# Pipeline Usage

## Setup 

#### 1. Clone this repository on the BYU HPC cluster and navigate to the directory.
```
git clone --recurse-submodules git@github.com:caysonjh/alphafold_kinase_screening.git
```

#### 2. Ensure that the required python libraries are installed
```
pip install -r requirements.txt
```


## Stage 1 -- AlphaFold3 Screens

#### 1. Customize script to bait protein
Edit `run_jobs.sh` to reflect the correct fasta file for the bait protein you are screening against. Ensure that the fasta file for said bait protein is in this directory. 

#### 2. Submit AlphaFold3 jobs
Ensure you are on a login node on the HPC cluster, and run `./run_jobs.sh`. This will retrieve all the kinase fasta files for kinases listed in `kinase_notkl.csv`, format them according to AlphaFold input, and submit the AlphaFold jobs. The output will be stored in output directories in an `output_dirs` folder that will be created. 

## Stage 2 -- Prepping for Analysis

Two scripts must be run to prepare the raw AlphaFold3 output for the later analysis scripts. 

#### 1. Run `generate_paeplots.sh`

This script will generate the PAE plot for each of the AlphaFold3 runs using the `generate_pae.py` script in this directory. 

#### 2. Run `prepare_download.sh` 

This script will create a new directory entitled `out_dirs` that will contain the necessary files for the later analysis scripts without the large data files, so that you could mass-download them to your machine with `sftp` if you would like. 

## Stage 3 -- Analysis 

#### 1. Run `run_full_pipeline.py`

Inside the `analysis/` directory is a python script `run_full_pipeline.py` that will perform the following steps  

- Create a final result directory within `final_dirs` for each run 
- Move the `.cif` model file into the directory 
- Collect the **ipTM** score from the AlphaFold3 output 
- Run **IPSAE** analysis using the `ipsae.py` from the included submodule 
- Run iLIS analysis using code modified from the original iLIS module 
- Concatenate the scores for each run into a single `all_scores.csv` file in the `final_dirs` directory 
- Generate standard and interactive violin plots for the score distributions
- Generate standard and interactive viewable screen rankings based on performance across the scoring metrics


Run the script with the `--project-root` parameter to specify the directory where your original `run_jobs.sh` was submitted. 
```
python run_full_pipeline.py --project-root /path/to/alphafold_project_dir
```


## Stage 4 -- Download and Viewing 

Any interesting screens can be downloaded to local machine for viewing in PyMOL or ChimeraX, or any further analysis. 
