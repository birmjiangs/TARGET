import os
import sys

# usage: python run_TARGET.py <tumor dataset name> <normal dataset name> <Hi-C data resolution> <mode>

T = "H209"
N = "MRC5"
resolution = 20000
mode = 1  #0: get DEG enrichment heatmap near altered TAD boundaries; 1: get TARGET candidate genes 

if len(sys.argv) == 5:
	T = sys.argv[1]
	N = sys.argv[2]
	resolution = sys.argv[3]

os.system("python ./1-get_insulation_altered_regions.py %s %s" % (T, N))
os.system("python ./2-get_altered_boundaries.py %s %s %s" % (T, N, resolution))
os.system("python ./3-enrichment_near_diffexp_genes.py %s %s" % (T, N))