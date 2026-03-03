mkdir out_dirs
cd output_dirs
for dir in *; do 
	mkdir ../out_dirs/$dir
	find $dir/ -type f \( -name "sp_confidences.json" -o -name "pae.png" -o -name "tr_model.cif" -o -name "sp_model.cif" -o -name "tr_summary_confidences.json" -o -name "sp_summary_confidences.json" \) -exec cp {} ../out_dirs/$dir/ \; 
done 
cd ../
