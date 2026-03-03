find out_dirs -type f -name 'pae.png' \
  -exec bash -c '
    for f do
      rel=${f#./}           
      path=${rel%/pae.png}  # strip filename
      second_dir=${path#*/}          # strip first component and slash
      second_dir=${second_dir%%/*}   # strip everything after next slash
      new="${second_dir}_pae.png"
      cp "$f" "./pae_plots/$new"
    done
  ' _ {} +
