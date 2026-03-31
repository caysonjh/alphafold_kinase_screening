import json 
import numpy as np 
import matplotlib.pyplot as plt 
import os

if "sp_confidences.json" in os.listdir(): 
    prefix = "sp"
else: 
    prefix = "tr"


with open(f'{prefix}_confidences.json', 'r') as f: 
    confidence_data = json.load(f) 
    
pae = np.array(confidence_data['pae'])

with open(f'{prefix}_data.json', 'r') as f: 
    data = json.load(f) 

chains = data['sequences'] 
lengths = [] 

for chain in chains:
    if 'protein' in chain: 
        lengths.append(len(chain['protein']['sequence']))


boundaries = np.cumsum(lengths)


plt.figure(figsize=(6, 5))
plt.imshow(pae, cmap="bwr", interpolation="nearest", vmin=0, vmax=30)
plt.colorbar(label="Predicted aligned error (Å)")
plt.xlabel("Residue index")
plt.ylabel("Residue index")
plt.title("Predicted aligned error (PAE)")

N = pae.shape[0] 
for b in boundaries[:-1]: 
    plt.axhline(b - 0.5, color='black', linestyle='--', linewidth=0.8)
    plt.axvline(b - 0.5, color='black', linestyle='--', linewidth=0.8)


plt.tight_layout()
plt.savefig("pae.png", dpi=300)
plt.close()
