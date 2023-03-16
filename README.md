# Helper scripts to execute HDStIM

HDStIM package is published on CRAN at https://cran.r-project.org/package=HDStIM. The documentation website is https://niaid.github.io/HDStIM/. Preprint with the conceptual framework https://www.biorxiv.org/content/10.1101/2022.11.14.516371v1.

These are a set of helper scripts to set up and execute HDStIM on the cell population - stimulation combinations derived either from automated clustering for example Robinson's pipeline or through manual gating.

1. Concatenate per population FSC files into a data frame
*Note: this step is not required if the data is exported from the Robinson's pipeline.*

Besides concatenating FCS files this step also does arcsinh transformation. Conventionally we have used a co-factor of 5 for CyTOf and 150 for Flow data. 
```
scr01_00_fcs_to_tsv.R
```

2. Run HDStIM, generate diagnostic plots, and run marker ranking function
```
scr02_00_run_hdstim_cytof.R
```

3. Generate heatmap from marker ranking
```
scr03_00_marker_ranking_boruta.R
```

4. Export HDStIM statistics 
```
scr04_00_export_hdstim_stats.R
```
