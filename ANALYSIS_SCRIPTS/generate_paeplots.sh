cd OUTPUT_DIRS
for dir in *; do 
	cd $dir
	cd "$(dirname "$(find . -name "TERMS_OF_USE.md" -type f | head -n1)")"
	pwd
	if [ -f "TERMS_OF_USE.md" ]; then 
		python ../../../ANALYSIS_SCRIPTS/generate_pae.py
	fi
	cd ../../../OUTPUT_DIRS/
done
cd ../
