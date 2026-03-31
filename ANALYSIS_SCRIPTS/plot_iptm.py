import plotnine as p9
import pandas as pd 
import json 

for file in os.listdir("summary_files/"): 
    if file.endswith(".json"): 
        with open(os.path.join("summary_files/", file), "r") as f: 
            data = json.load(f) 
            df = pd.DataFrame(data[['iptm', 'ranking_score', 'ptm']])
            df['kinase_id'] = file.split("/")[1]
            df['bait_protein'] = 'SMO C-term'


plot = (
        p9.ggplot(df, p9.aes(x='bait_protein', y='iptm')) + 
        p9.geom_voilin(draw_quantiles=[0.25, 0.5, 0.75]) + 
        p9.theme_bw() + 
        p9.xlab("IPTM") +
)

plot.save("iptm_voilin_plot.png", dpi=300)
