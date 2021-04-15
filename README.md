## Input files:
1. Insulation score of tumor cells: a HOMER output .bedGraph file at ./data/<tumor dataset name>.Insulation.bedGraph
2. Insulation score of normal cells: a HOMER output .bedGraph file at ./data/<normal dataset name>.Insulation.bedGraph
3. Genome coordinates of TADs in tumor cells: A HOMER output .bed file at ./data/<tumor dataset name>.tad.2D.bed
4. Genome coordinates of TADs in normal cells: A HOMER output .bed file at ./data/<normal dataset name>.tad.2D.bed
5. A three-column tab-seperated(.tsv) file that provides information of differentially expressed genes at ./data/differentially_expressed_genes.tsv 
	column1: gene name. 
	column2: expression log2 fold change (tumor/normal). 
	column3: FDR. 
6. A file specifing the length of each chromosome at ./data/chr.sizes.genome
7. Genome coordinates of transcription start sites(TSSs) of all genes to be investigated at ./data/PCG_TSS


## Usage:
1. To perform TARGET analysis with the example dataset, use the following command without specifing any parameters:
```bash
python run_TARGET.py
```
2. To perform TARGET analysis using customized datasets, replace example datasets in ./data/.
The file name can be modified in <tumor dataset name> and <normal dataset name>, other parts of the file name should not be changed.
If file names are modified, perform TARGET analysis using the following command, where all parameters are required to be specified:
```bash
python run_TARGET.py <tumor dataset name> <normal dataset name> <Hi-C data resolution>
```
Make sure <tumor dataset name> and <normal dataset name> are identical to the corresponding part of the input file name.
<Hi-C data resolution> specifies Hi-C data resolution (The default resolution is set to 20000 for example data).
The output file that containing all TARGET candidate genes will be at ./TARGET_candidate_genes.txt
