mkdir DOWNLOAD_DIRS
cd OUTPUT_DIRS
for dir in *; do 
	mkdir ../DOWNLOAD_DIRS/$dir
	find $dir/ -type f \( -name "sp_confidences.json" -o -name "pae.png" -o -name "tr_model.cif" -o -name "sp_model.cif" -o -name "tr_summary_confidences.json" -o -name "sp_summary_confidences.json" \) -exec cp {} ../DOWNLOAD_DIRS/$dir/ \; 
done 
cd ../
