import os
import sys
import numpy as np


def mkdir(path):
    isExists=os.path.exists(path)
    if not isExists:
        os.makedirs(path)
        return True
    else:
        return False

def savetxt(filename,x):
    np.savetxt(filename,x,delimiter = '\t',fmt='%s')

def sort_list(list_in): #use when using '\t' to seperate items
    list_out = sorted(list_in, key=lambda items: int(items.split('\t')[1]))
    list_out = sorted(list_out, key=lambda items: items.split('\t')[0])
    return list_out

#expand altered boundary by the half_expand_len to test a larger range of effect in boundary alteration to gene expression
T = sys.argv[1]
N = sys.argv[2]


half_expand_len = 40000
FC_thres = 2
count_altered_boundary = 0
out = []
f = open("2-altered_IS_boundary_%s_vs_%s.txt" % (T,N))   		
lines=f.readlines() 
nrow = len(lines)					
for i in range(len(lines)):
    L = lines[i].strip().split('\t')
    L[1] = str(int(float(L[1]) - half_expand_len))
    L[2] = str(int(float(L[2]) + half_expand_len))
    if float(L[3]) + float(L[4]) in [1,2]:  ##### float(L[3]) + float(L[4]) = 1: HOMER detect boundary alteration;  = 2: HOMER detect both celltype with boundaries at this bin
        # delete inconsistent result of IS change and HOMER TAD identification
        if (float(L[3]) - float(L[4])) * (float(L[5]) - float(L[6])) < 0:
            continue

        out.append(L[0]+'\t'+L[1]+'\t'+L[2])

f.close()
mkdir("./3-boundary_for_overlap")
savetxt("./3-boundary_for_overlap/altered", out)

gene_2_TSS = {}
f = open("./data/PCG_TSS")           
lines=f.readlines() 
nrow = len(lines)                   
for i in range(len(lines)):
    L = lines[i].strip().split('\t')
    chrN = L[0]
    start = str(int(float(L[1])-1000))
    end = str(int(float(L[1])+1000))
    if L[3] == "-":
        start = str(int(float(L[2])-1000))
        end = str(int(float(L[2])+1000))
    gene_2_TSS[L[4]] = chrN+'\t'+start+'\t'+end
f.close()

out = []
f = open("./data/differentially_expressed_genes.tsv")           
lines=f.readlines() 
nrow = len(lines)                   
for i in range(len(lines)):
    L = lines[i].strip().split('\t')
    gene_name = L[0]
    if gene_name in gene_2_TSS:
        gene_promoter = gene_2_TSS[gene_name]
        sig = "no"
        if float(L[2]) < 0.05 and abs(float(L[1])) > FC_thres:
            sig = "yes"

        sig_only_qvalue = "no"
        if float(L[2]) < 0.05:
            sig_only_qvalue = "yes"
        out.append(gene_promoter+'\t'+gene_name+'\t'+L[1]+'\t'+sig)

mkdir("./3-promoter_gene_sig_diff")
savetxt("./3-promoter_gene_sig_diff/diff_sig", out)


os.system("python ./overlap.py ./3-boundary_for_overlap/altered ./3-promoter_gene_sig_diff 1 0 ./3-overlapped")



#get gene list of genes near altered TAD boundary
gene_out = []
no_sig_gene = []
f = open("./3-overlapped/altered_diff_sig.txt")           
lines=f.readlines() 
nrow = len(lines)                   
for i in range(len(lines)):
    L = lines[i].strip().split('\t')
    for j in range(3,len(L),6):
        if L[j+5] == "yes":
            gene_out.append(L[j+3]+'\t'+L[j+4])
        if L[j+5] == "no":
            no_sig_gene.append(L[j+3]+'\t'+L[j+4])
f.close()
gene_out = list(set(gene_out))
no_sig_gene = list(set(no_sig_gene))

print(len(gene_out),len(no_sig_gene))
savetxt("./TARGET_candidate_genes.txt", gene_out)

