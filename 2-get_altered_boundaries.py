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
  
################
T = sys.argv[1]
N = sys.argv[2]
resolution = int(float(sys.argv[3]))
####### 0. index genome
pos_2_index = {}
index_2_pos = {}
index = 0
f = open("./data/chr.sizes.genome")           
lines=f.readlines() 
nrow = len(lines)                   
for i in range(len(lines)):
    L = lines[i].strip().split('\t')
    chrN = L[0]
    max_pos = int(L[1])
    for j in range(resolution,max_pos,resolution):
        start = str(j-resolution)
        end = j
        pos = chrN + '\t' + str(start) + '\t' + str(end)
        pos_2_index[pos] = index 
        index_2_pos[index] = pos
        index += 1
f.close()

##### 1. get all boundaries and ISs
boundaries = [[0,0] for i in range(index)]


f = open("./data/"+T+".tad.2D.bed")        # read HOMER output TAD file
lines=f.readlines() 
nrow = len(lines)                   
for i in range(len(lines)):
    L = lines[i].strip().split('\t')
    boundary1_pos = L[0] + '\t' + str(round(int(L[1])/20000)*20000) + '\t' + str(round(int(L[1])/20000+1)*20000)
    boundary2_pos = L[0] + '\t' + str(round(int(L[2])/20000-1)*20000) + '\t' + str(round(int(L[2])/20000)*20000)

    boundaries[pos_2_index[boundary1_pos]][0] = 1
    boundaries[pos_2_index[boundary2_pos]][0] = 1
f.close()

f = open("./data/"+N+".tad.2D.bed")          
lines=f.readlines() 
nrow = len(lines)                   
for i in range(len(lines)):
    L = lines[i].strip().split('\t')
    boundary1_pos = L[0] + '\t' + str(round(int(L[1])/20000)*20000) + '\t' + str(round(int(L[1])/20000+1)*20000)
    boundary2_pos = L[0] + '\t' + str(round(int(L[2])/20000-1)*20000) + '\t' + str(round(int(L[2])/20000)*20000)

    boundaries[pos_2_index[boundary1_pos]][1] = 1
    boundaries[pos_2_index[boundary2_pos]][1] = 1
f.close()


out = []
for i in range(index):
    pos = index_2_pos[i]
    out.append(pos+'\t'+"\t".join(list(map(str,boundaries[i]))))
savetxt("./2-boundaries_"+T+"_"+N+"", out)


##### get altered boundaries (tumor vs normal(normal))
pos_2_altered_tumor = {}
f = open("./1-diff_IS_positions/"+T+"_"+N+"")           
lines=f.readlines() 
nrow = len(lines)                   
for i in range(len(lines)):
    L = lines[i].strip().split('\t')
    pos = L[0] + '\t' + L[1] + '\t' + L[2]
    tumor_normal_IS = L[3] + '\t' + L[4]
    pos_2_altered_tumor[pos] = tumor_normal_IS
f.close()

out_altered_IS_boundary = []
for i in range(len(out)):
    L = out[i].split("\t")
    if L[3] == "0" and L[4] == "0":  #not boundary in either tumor or normal
        continue

    pos = L[0] + '\t' + L[1] + '\t' + L[2]
    if pos in pos_2_altered_tumor:
        out_altered_IS_boundary.append(out[i]+'\t'+pos_2_altered_tumor[pos])  #out[i]:pos, tumor_is_boundary, normal(normal)_is_boundary

savetxt("./2-altered_IS_boundary_"+T+"_vs_"+N+".txt", out_altered_IS_boundary)




####### part2, preparing data for 6.py: mark all bins with boundary(0/1) and altered IS(0/1)
## 0. get IS of all bins
ISs = [[float("nan") for n in range(3)] for i in range(index)]

f = open("./data/"+T+".Insulation.bedGraph")          
lines=f.readlines() 
nrow = len(lines)                   
for i in range(len(lines)):
    if i == 0:
        continue
    L = lines[i].strip().split('\t')
    pos = L[0] + '\t' + L[1] + '\t' + L[2]
    if pos not in pos_2_index:
        continue
    ISs[pos_2_index[pos]][0] = float(L[3])
f.close()

f = open("./data/"+N+".Insulation.bedGraph")          
lines=f.readlines() 
nrow = len(lines)                   
for i in range(len(lines)):
    if i == 0:
        continue
    L = lines[i].strip().split('\t')
    pos = L[0] + '\t' + L[1] + '\t' + L[2]
    if pos not in pos_2_index:
        continue
    ISs[pos_2_index[pos]][1] = float(L[3])
f.close()


## 1. do tumor vs normal
mkdir("./2-loci_proporties_for_6")
out = ["chr\tstart\tend\tis_boundary_"+T+"\tis_boundary_"+N+"\tis_altered_IS\tIS_"+T+"_if_altered\tIS_"+N+"_if_altered\tIS_"+T+"\tIS_"+N+""]
for i in range(index):
    pos = index_2_pos[i]
    boundary_of_3celltypes = list(map(str,boundaries[i]))
    if pos in pos_2_altered_tumor:  #this position corresponds to altered IS in tumor
        out.append(pos+'\t'+ boundary_of_3celltypes[0] + '\t' + boundary_of_3celltypes[1] + '\t1\t' + pos_2_altered_tumor[pos]+'\t'+str(ISs[pos_2_index[pos]][0])+'\t'+str(ISs[pos_2_index[pos]][1]))
        continue
    out.append(pos+'\t'+ boundary_of_3celltypes[0] + '\t' + boundary_of_3celltypes[1] + '\t0\t' + "nan\tnan"+'\t'+str(ISs[pos_2_index[pos]][0])+'\t'+str(ISs[pos_2_index[pos]][1]))
savetxt("./2-loci_proporties_for_6/"+T+"_vs_"+N+".tsv", out)

