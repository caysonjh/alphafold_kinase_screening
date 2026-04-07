# Kinase Screening

## Dependencies and Setup

### AlphaFold Permissions

In order to access AlphaFold3 on the BYU HPC cluster, permission must be obtained from Google to get access to AlphaFold3's parameters. For information on accessing the HPC, see [here](https://rc.byu.edu/wiki/?id=Logging+In).

Instructions on how to get permission can be found on BYU's [AlphaFold3 Page](https://rc.byu.edu/wiki/?page=AlphaFold+3) under the "How to Join" section. It often takes a day or two to get permissions from Google. 

Once permissions have been granted from Google, you will forward their email to rcsupport@byu.edu. They will then give you permission to run `module load alphafold3/3.0.1+`. 

### Environment Set-Up

#### 1. Clone this repository on the BYU HPC cluster and navigate to the directory (can be done locally if using sample/pre-run AlphaFold data).
```
git clone --recurse-submodules https://github.com/caysonjh/alphafold_kinase_screening.git
```

#### 2. Ensure that the required python libraries are installed
```
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

### Sample Data

To investigate the post-AlphaFold scripting and functionality, a sample dataset is located [here](https://byu.box.com/s/1llegwooxapbt2yr9b1mrary4phwfrxz). Download the entire linked folder into this cloned repository, it will be usable for the future analysis scripts. Ensure that the directory is named **OUTPUT_DIRS** as can be seen via the link. 


### Input File Formatting

There are two essential files needed to run the screens:  
- Bait Protein Fasta File -- fasta file for the protein that will be screened against  
- Test Protein CSV File -- csv file containing each protein to be screened against the bait protein  
    - **NOTE:** The CSV file should contain a column titled `UniprotID` that contains either the UniprotID or the protein name so that the fasta file can be retrieved from the Uniprot website via the [API](https://www.uniprot.org/help/api_retrieve_entries)
    <br><br>

---
---

## Screening Phase (Must Have AlphaFold Parameter Permissions)

Ensure that you are on a **LOGIN** node on the BYU supercomputer, the jobs will be submitted automatically via slurm.  

This is the general format for running the jobs:

```
./SUBMISSION_SCRIPTS/run_jobs.sh -i path/to/test_proteins.csv -p path/to/bait_protein.fasta
```
Sample bait proteins are included in the `bait_proteins` directory  

#### This script will: 
- Retrieve fasta files for all the test proteins
- Convert the fasta files into json input for AlphaFold
- Submit the AlphaFold slurm jobs
- Direct the AlphaFold output to `output_dirs` with a directory named after the protein's UniprotID or Name
<br><br>

## Output Formatting

Two scripts must be run to prepare the raw AlphaFold3 output for the later analysis scripts. 

#### 1. Run `generate_paeplots.sh`

This script will generate the PAE plot for each of the AlphaFold3 runs using the `generate_pae.py` script in this directory. 
```
./ANALYSIS_SCRIPTS/generate_paeplots.sh
```

#### 2. Run `prepare_download.sh` 

This script will create a new directory entitled `DOWNLOAD_DIRS` that will contain the necessary files for the later analysis scripts without the large data files, so that you could mass-download them to your machine with `sftp` if you would like. 
```
./ANALYSIS_SCRIPTS/prepare_download.sh
```

## Analysis and Figure Generation

#### 1. Run `run_full_pipeline.py`

Inside the `ANALYSIS_SCRIPTS/` directory is a python script `run_full_pipeline.py` that will perform the following steps  

- Create a final result directory within `FINAL_DIRS` for each run 
- Move the `.cif` model file into the directory 
- Collect the **ipTM** score from the AlphaFold3 output 
- Run **IPSAE** analysis using the `ipsae.py` from the included submodule 
- Run **iLIS** and **LIS** analysis using code modified from the original iLIS module 
- Concatenate the scores for each run into a single `all_scores.csv` file in the `FINAL_DIRS` directory 

Run the script with the `--project-root` parameter to specify the directory where your original `run_jobs.sh` was submitted. 
```
python ANALYSIS_SCRIPTS/run_full_pipeline.py --project-root /path/to/alphafold_project_dir
```

The output from the pipeline will include interactive html diagrams including: 
- `FINAL_DIRS/RANKINGS/ranking_report.html` -- This file will display the sorted Top 30 AlphaFold scans using a composite score based on ipTM, IPSAE, LIS, and iLIS.
- `FINAL_DIRS/all_scores_interactive.html` -- This file shows the violin plots and scatter matrix (comparing whether the high scores for each metric also correlate to high scores on the other ones) for each of ipTM, IPSAE, LIS, and iLIS.

#### These html files can be opened in any browser for visualization, and can be downloaded from the supercomputer using `scp` or `sftp`

You will also find static figures generated in the following files: 
- `FINAL_DIRS/RANKINGS/ranking_report.pdf`
- `FINAL_DIRS/all_scores_violin.png`

There will be csv files that can be used for further, more specific analysis: 
- `FINAL_DIRS/RANKINGS/all_scores_ranked.csv`
- `FINAL_DIRS/RANKINGS/top_iLIS.csv` -- Top 30 for iLIS 
- `FINAL_DIRS/RANKINGS/top_IPSAE.csv` -- Top 30 for IPSAE
- `FINAL_DIRS/RANKINGS/top_ipTM.csv` -- Top 30 for ipTM
- `FINAL_DIRS/RANKINGS/top_LIS.csv` -- Top 30 for LIS
- `FINAL_DIRS/RANKINGS/top_overall.csv` -- Top 30 overall
- `FINAL_DIRS/{your_protein}_all_scores.csv`


## Visualization

### Figure 2
Run the all_baits_violin_figure.py file to calculate top scores across bait proteins. The current files in the script are for SMO C-term, GLI1 N-term, GLI2, and SUFU.
```
python all_baits_violin_figure.py
```

#### If the steps up to this point were done on the supercomputer, you will need to download the files of interest using `scp` or `sftp`

#### If the steps were done on your local machine with a sample dataset, that is not neccessary

To view the `.cif` models that are found within the AlphaFold output directories, you will need a molecular visualization software such as [ChimeraX](https://www.rbvi.ucsf.edu/chimerax/download.html). From there you can identify the interacting residues, view the interaction structure, and more. 

For more information on using ChimeraX, see these [tutorials](https://www.rbvi.ucsf.edu/chimerax/download.html)

---
---

## All Script Descriptions

### SUBMISSION_SCRIPTS
- `run_jobs.sh` -- main AlphaFold jobs submission scripts, calls `get_fasta_files_from_uniprot.py` and submits AlphaFold3 slurm jobs. 
- `get_fasta_files_from_uniprot.py` -- uses the Uniprot API to get fasta files for each of the test proteins. 
- `prepare_download.sh` -- gets all the critical files from the raw AlphaFold3 output to create minimal folders for easier downloading to local machine

### ANAYSIS_SCRIPTS
- `extract_ilis_batch.py` -- modified from the [iLIS ipynb scripts](https://github.com/flyark/AFM-LIS) to calculate iLIS scores
- `generate_pae.py` -- uses the raw output from AlphaFold3 to generate a pae plot to see interacting residue areas
- `generate_paeplots.sh` -- navigates to each output directory and runs the `generate_pae.py` script there
- `plot_all_scores_interactive.py` -- creates the interactive html for violin plot scores across metrics
- `plot_all_scores_violin.py` -- creates a static .png with violin plots across metrics
- `rank_score_ids.py` -- creates ranking based on composite scoring of the test proteins
- `run_full_pipeline.py` -- runs each of the above scripts and formats the output
