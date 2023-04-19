# Run CaMa-Flood simultaneously

## Steps
1. Prepare runoff ensemble  `sh s01-prep_runoff.sh`
2. Run CaMa-Flood `sh s02-run_CaMa.sh`
3. Convert sfcelv to netCDF `sh s03-conv_bin2nc.sh`
4. Making statistics `sh s04-mk_stat.sh`
5. Copy statistics to HydroDA/dat folder `sh s05-copy2HydroDA.sh`

### Prepare runoff ensemble
Runoff can be perturbed using simple, normal, or lognormal methods

### Parallel simulation CaMa-Flood
Run CaMa-Flood simulation parallel

### Making statistics
Statistics such as mean and standard deviation can be made

### Copy to HydroDA folder
For anomaly or normalized DA method the statistics need to be saved at ./HydroDA/dat

## Visualization (./img_code)
1. runoff_ens.py - visualize perturbed runoff
2. digraph_ens.py - visualize simulated river discharge using perturbed runoff