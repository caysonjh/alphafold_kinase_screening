cd out_dirs
for dir in *; do 
	cd $dir
	cd "$(dirname "$(find . -name "TERMS_OF_USE.md" -type f | head -n1)")"
	pwd
	if [ -f "TERMS_OF_USE.md" ]; then 
		python /home/caysonjh/AlphaFold_Stuff/generate_pae.py
	fi
	cd /home/caysonjh/AlphaFold_Stuff/out_dirs/
done
cd ../
