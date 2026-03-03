cd out_dirs
for dir in *; do 
	mkdir ../download_dirs/$dir
	find $dir/ -type f \( -name "pae.png" -o -name "tr_model.cif" -o -name "sp_model.cif" -o -name "tr_summary_confidences.json" -o -name "sp_summary_confidences.json" \) -exec cp {} ../download_dirs/$dir/ \; 
done 
cd ../
