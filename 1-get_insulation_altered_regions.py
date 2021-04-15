import os
import sys
import numpy as np
import statsmodels.api as sm
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.rcParams['pdf.fonttype'] = 42



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

def sort_data(list_in): #use when using '\t' to seperate items
    list_out = sorted(list_in, key=lambda items: items[0])
    return list_out
  
################
T = sys.argv[1]
N = sys.argv[2]

celltypes = [T, N]  # celltypes 
annotations = ["T","N"]

pos_2_IS_ct1 = {}
pos_2_IS_ct2 = {}

########### get insulation score of each celltype ############
f = open("./data/"+celltypes[0]+".Insulation.bedGraph")         # HOEMR input insulation score file   
lines=f.readlines() 
nrow = len(lines)                   
for i in range(len(lines)):
    if i == 0:
        continue
    L = lines[i].strip().split('\t')
    if L[0] == "chrY":
        continue

    pos = L[0] + '\t' + L[1] + '\t' + L[2]
    IS = float(L[3])
    pos_2_IS_ct1[pos] = IS
f.close()


f = open("./data/"+celltypes[1]+".Insulation.bedGraph")          
lines=f.readlines() 
nrow = len(lines)                   
for i in range(len(lines)):
    if i == 0:
        continue
    L = lines[i].strip().split('\t')
    pos = L[0] + '\t' + L[1] + '\t' + L[2]
    IS = float(L[3])
    pos_2_IS_ct2[pos] = IS
f.close()



##### get common positions
common_pos = []
for pos in pos_2_IS_ct1:
    if pos in pos_2_IS_ct2:
        if pos_2_IS_ct1[pos] * pos_2_IS_ct2[pos] == 0:  #skip zero IS(no HiC contact)
            continue
        common_pos.append(pos)

common_pos = sort_list(common_pos)

##### get diff IS between two celltypes
def compare_IS(pos_2_IS_ct1,pos_2_IS_ct2,ct1,ct2):

    diff_IS = []
    mean_IS = []
    
    std_used_points = 20
    std_diff_IS = []  #standard variation of adjacent 20 points

    for i in range(len(common_pos)):
        pos = common_pos[i]
        IS1 = pos_2_IS_ct1[pos]
        IS2 = pos_2_IS_ct2[pos]

        diff_IS.append(abs(IS1-IS2))
        mean_IS.append((IS1+IS2)/2)

    data_for_lowess = [[mean_IS[i],diff_IS[i]] for i in range(len(diff_IS))]
    data_for_lowess = sort_data(data_for_lowess)

    ## calculate std of diffIS for each mean IS:
    #@@for i in range(len(common_pos)):
#@@
    #@@	mean_IS_this = data_for_lowess[i][0]  # not useful for now, if need to set std calculation range based on, like current mean IS += 0.1, then this value is useful
#@@
    #@@	diff_ISs_for_std_calculation = []
    #@@	for j in range(max(0,int(i-std_used_points)), min(i+std_used_points,len(common_pos)-1)):
    #@@		if i == 0:
    #@@			print(j)
    #@@		diff_ISs_for_std_calculation.append(data_for_lowess[j][1])
    #@@	
    #@@	std = np.std(diff_ISs_for_std_calculation)
    #@@	std_diff_IS.append(std)
#@@
    #@@	if i % 1000 == 0:
    #@@		print(i, len(common_pos))
#@@
    #@@data_for_loess_std = [[data_for_lowess[i][0], std_diff_IS[i]] for i in range(len(diff_IS))]

    data_for_lowess = np.array(data_for_lowess)
    lowess_mean = sm.nonparametric.lowess(data_for_lowess[:,1], data_for_lowess[:,0], frac=0.01, it = 0)

    #@@data_for_loess_std = np.array(data_for_loess_std)
    #@@lowess_std = sm.nonparametric.lowess(data_for_loess_std[:,1], data_for_loess_std[:,0], frac=0.1, it = 0)
    #### (get 2 x fitted variation) <- do not use, using mean + 2 * std now!
    lowess_cutoff = np.copy(lowess_mean)
    lowess_cutoff[:,1] = lowess_mean[:,1] * 2

    mean_IS_2_diff_cutoffs = dict(zip(lowess_cutoff[:,0].tolist(), lowess_cutoff[:,1].tolist()))


    #### get bins that exceeds the lowess cutoff:
    out = []
    y_high_var = []
    x_high_var = []

    y_others = []  # non significant
    x_others = []
    for i in range(len(common_pos)):

        if i % 1000 == 0:
            print(i, len(common_pos))
        
        pos = common_pos[i]
        IS1 = pos_2_IS_ct1[pos]
        IS2 = pos_2_IS_ct2[pos]
        diff_IS = abs(IS1-IS2)
        mean_IS = (IS1+IS2)/2
        lowess_cutoff_fitted = mean_IS_2_diff_cutoffs[mean_IS]
        if lowess_cutoff_fitted < diff_IS:
            out.append(pos+'\t'+str(IS1)+'\t'+str(IS2))
            y_high_var.append(diff_IS)
            x_high_var.append(mean_IS)
        else:
        	y_others.append(diff_IS)
        	x_others.append(mean_IS)

    mkdir("./1-diff_IS_scatterplots")
    plt.figure(figsize=(3, 2), dpi=600)
    plt.rc('font',family='Arial')
    plt.rc('font',size = 9)
    bp = plt.scatter(x_others ,y_others, s = 0.5, color = [[0.6,0.6,0.6] for i in range(len(x_others))], linewidth = 0, alpha = 0.1)
    bp = plt.scatter(x_high_var ,y_high_var, s = 1, color = [[0.85,0.0,0.0] for i in range(len(x_high_var))], linewidth = 0, alpha = 0.5)
    plt.savefig("./1-diff_IS_scatterplots/"+ct1+"_"+ct2+".png")
    plt.savefig("./1-diff_IS_scatterplots/"+ct1+"_"+ct2+".pdf")
    plt.close()

    #@@plt.figure(figsize=(4, 2), dpi=600)
    #@@plt.rc('font',family='Arial')
    #@@plt.rc('font',size = 9)
    #@@bp = plt.scatter(data_for_loess_std[:,0] ,data_for_loess_std[:,1], s = 0.5, color = [[0.6,0.6,0.6] for i in range(len(data_for_loess_std))], linewidth = 0, alpha = 0.5)
    #@@for i in range(len(data_for_loess_std[:,0])):
    #@@    if i > 25 or i + 25 > len(data_for_loess_std[:,0]):
    #@@        if i % 1000 != 0:  # do not draw too many lines
    #@@            continue
    #@@    plt.plot([lowess_std[max(0,i-1),0], lowess_std[i,0]] ,[lowess_std[max(0,i-1),1], lowess_std[i,1]], color = [0.0,0.0,0.0], linewidth = 0.5)
    #@@plt.savefig("./1-diff_IS_scatterplots/"+ct1+"_"+ct2+"_std_of_diff_IS_at_each_mean_IS.png")
    #@@plt.close()

    mkdir("./1-diff_IS_positions")
    savetxt("./1-diff_IS_positions/"+ct1+"_"+ct2, out)
    #@@savetxt("./1-diff_IS_positions"+ct1+"_"+ct2+"_lowess_std.tsv", lowess_std)  # lowess regressed standard variation of diff IS against mean IS(2nd col) at each mean IS(1st col)
    #@@savetxt("./1-diff_IS_positions"+ct1+"_"+ct2+"_std.tsv", data_for_loess_std)  #standard variation of diff(2nd col) IS at each mean IS(1st col)
    return out

out1 = compare_IS(pos_2_IS_ct1,pos_2_IS_ct2,celltypes[0],celltypes[1])

