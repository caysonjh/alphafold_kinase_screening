# Clear module environment
module purge

# Ensure that python libraries are installed on the local python environment being used 
if ! python -m requests --version >/dev/null 2>&1; then
  echo "Installing requests module..."
  pip install requests
else
  echo "requests module is already installed."
fi
if ! python -m pandas --version >/dev/null 2>&1; then
  echo "Installing pandas module..."
  pip install pandas
else
  echo "pandas module is already installed."
fi

# Parse command line arguments
input=""
protein=""
while getopts "i:p:" opt; do
  case "$opt" in
  i) input="$OPTARG" ;;
  p) protein="$OPTARG" ;;
  \?) echo "Invalid option -$OPTARG" >&2 ;;
  esac
done

# Make fasta file directory if doesn't exist
if [ ! -d "FASTA_FILES" ]; then
  mkdir FASTA_FILES
fi

# Get the fasta files for each of the testing proteins
python ./SUBMISSION_SCRIPTS/get_fasta_files_from_uniprot.py --kinase_file $input --init_protein_fasta $protein

# Load AlphaFold module
module load alphafold3/3.0.1+

# Make directories if needed
if [ ! -d "JSON_FILES" ]; then
  mkdir JSON_FILES
fi
if [ ! -d "OUTPUT_DIRS" ]; then
  mkdir OUTPUT_DIRS
fi

# For each of the fasta files that were downloaded, submit an AlphaFold job
for fasta_file in FASTA_FILES/*; do
  if [ -f "$fasta_file" ]; then
    echo $fasta_file
    base_name="$(basename "$fasta_file")"
    echo $base_name
    mkdir OUTPUT_DIRS/"${base_name%.*}"
    f2j.py "$fasta_file" JSON_FILES/"${base_name%.*}".json
    alphafold3_pipeline.sh JSON_FILES/"${base_name%.*}".json OUTPUT_DIRS/"${base_name%.*}"
  fi
done
