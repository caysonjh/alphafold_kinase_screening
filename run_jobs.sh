module purge 
pip install requests 
pip install pandas

mkdir fasta_files
python get_fasta_files.py --kinase_file kinases_notkl.csv --init_protein_fasta smo.fasta

module load alphafold3/3.0.1+
mkdir json_files
mkdir out_dirs
for fasta_file in fasta_files/*; do
	if [ -f "$fasta_file" ]; then
		echo $fasta_file
		base_name="$(basename "$fasta_file")"
		echo $base_name
		mkdir out_dirs/"${base_name%.*}"	
		f2j.py "$fasta_file" json_files/"${base_name%.*}".json
		alphafold3_pipeline.sh json_files/"${base_name%.*}".json out_dirs/"${base_name%.*}"
	fi

done
