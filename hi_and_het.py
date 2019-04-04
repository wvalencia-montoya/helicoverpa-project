
####Bash script to run HI and interspecific heterozygosity


cd /home/wv221/rcs/rcs-cj107-heliconius/wv221/GVCFs_raw/raw_gvcfs/test_holi/anc_prop
python exe_py_p_anc_1.py


#####Python script

import pandas as pd
import numpy as np
import sys
import getopt
import os

os.chdir("/home/wv221/rcs/rcs-cj107-heliconius/wv221/GVCFs_raw/raw_gvcfs/test_holi/anc_prop")

#####this table is the same that I did previously
joint_set = pd.read_table("joint_set.tsv", sep='\t') ###Simon's geno


########Here just to pay the individuals you want to analyse

list_ind=["110_N704_S504",
"142_N703_S503",
"6_1_1",
"6_1_10",
"6_1_2",
"6_1_20",
"6_1_22",
"6_1_3",
"6_1_4",
"6_1_5",
"6_1_7"]


##########################


def tab_prop_anc(list_ind):
    def prop_ance_ind(individual):
        ind=individual
        pos_gen = ['zea_gen','arm_gen','het_a','het_b']
        #zea_gen_name = ind + '_' + pos_gen[0]
        zea_hom = []
        #
        for i in range(len(joint_set[ind])):
            if joint_set[ind][i] == joint_set[pos_gen[0]][i]:
             zea_hom.append(1.0)
            elif joint_set[ind][i] == joint_set[pos_gen[2]][i]:
             zea_hom.append(0.5)
            elif joint_set[ind][i] == joint_set[pos_gen[3]][i]:
              zea_hom.append(0.5)
            else:
              zea_hom.append(0.0)
              #
              #
        a = (sum(zea_hom))/len(zea_hom)
        return a
    #
    def het_ind(individual):
        ind=individual
        pos_gen = ['zea_gen','arm_gen','het_a','het_b']
        #zea_gen_name = ind + '_' + pos_gen[0]
        het = []
        #
        for i in range(len(joint_set[ind])):
            if joint_set[ind][i] == joint_set[pos_gen[0]][i]:
             het.append(0.0)
            elif joint_set[ind][i] == joint_set[pos_gen[2]][i]:
             het.append(1)
            elif joint_set[ind][i] == joint_set[pos_gen[3]][i]:
              het.append(1)
            else:
              het.append(0.0)
              #
              #
        h = (sum(het))/len(het)
        return h
        #
        #
    ind=[]
    anc_prop_ind=[]
    ho_ind = []
    for i in list_ind:
        ind.append(i)
        anc_prop_ind.append(prop_ance_ind(i))
        ho_ind.append(het_ind(i))
    #
    ance_tab_ind = pd.DataFrame({'ind': ind,'zea_anc': anc_prop_ind, 'het_prop': ho_ind})
    return(ance_tab_ind)
    
ance_tab_ho=tab_prop_anc(list_ind)

ance_tab_ho.to_csv('ance_tab_ho1.1.tsv', sep='\t')



