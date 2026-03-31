cd ../OUTPUT_DIRS
for dir in *; do 
	cd $dir
	cd "$(dirname "$(find . -name "TERMS_OF_USE.md" -type f | head -n1)")"
	pwd
	if [ -f "TERMS_OF_USE.md" ]; then 
		python /home/caysonjh/AlphaFold_Stuff/ANALYSIS_SCRIPTS/generate_pae.py
	fi
	cd /home/caysonjh/AlphaFold_Stuff/OUTPUT_DIRS/
done
cd ../
