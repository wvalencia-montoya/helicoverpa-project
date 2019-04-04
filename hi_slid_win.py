cd /home/wv221/rcs/rcs-cj107-heliconius/wv221/GVCFs_raw/raw_gvcfs/test_holi/anc_prop
python exe_py_win_anc_1.py


####Python scripts

import pandas as pd
import numpy as np
import sys
import getopt
import os

os.chdir("/home/wv221/rcs/rcs-cj107-heliconius/wv221/GVCFs_raw/raw_gvcfs/test_holi/anc_prop")

joint_set = pd.read_table("joint_set.tsv", sep='\t') ###Simon's geno

#############
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
#############

def tab_prop_ance_win(list_ind, win_size):
    ###Function to calculate ancestry
    #
    def prop_ance(individual):
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
        print a
        #
        n_chr = list(joint_set['CHR_x'])
        pos = list(joint_set['POS_x'])
        chr_pos = list(joint_set['CHR-POS'])
        zea_anc = zea_hom
        anc_tab = pd.DataFrame(
        {'chr': n_chr,
         'pos': pos,
         'chr_pos': chr_pos,
         'zea_anc': zea_hom,
          })
        #
        return (anc_tab)
        #
    #
    ##Function to calculate ancestry in sliding windows
    def prop_ance_windows(anc_tab_ind, win_size):
        anc_tab = anc_tab_ind
        win_size = win_size
        #bed file of chromosomes lenghts
        chr_len = [15622547, 6414787, 10372743, 12978399, 12389108, 11903839, 9462763, 11503474, 12160917, 12835424, 4645711, 11616497,12328174, 9898702, 13242999, 10058692, 11387182,
        10601134, 10244286, 9437340, 10314289, 12896724, 12267648, 3773128, 9468847, 8159704, 7782505, 7377967, 11265821, 3551048, 6553157]
        #ance_sli_win = pd.DataFrame(columns=['chr','win_ini','win_mid','win_end','n_snps','prop_ance'])
        #
        ance_com = pd.DataFrame(columns=['chr', 'n_snps' , 'prop_ance' ,  'win_end' ,  'win_ini' ,  'win_mid'])
        #
        for i in range(1,32):
            temp_tab = anc_tab.loc[(anc_tab['chr']== i)]
            win_id = range(1,((chr_len[i-1] / win_size) + 1 ))
            win_pos = range(1,chr_len[i-1],win_size)
            chrm = [i] * (len(win_id))
            #
            win_ini=[]
            win_mid=[]
            win_end=[]
            n_snps = []
            prop_ance_win=[]
            #
            for k in win_id:
                win_ini.append(win_pos[k-1])
                win_mid.append((win_pos[k-1]+win_pos[k])/2)
                win_end.append(win_pos[k])
            #
            #
            for m in range(len(win_id)):
                prop_sum_w = []
                prop_ance_temp =[]
                for n in range(len(temp_tab['chr'])):
                    if win_ini[m] < temp_tab.iloc[n,2] < win_end[m]:
                       prop_sum_w.append(temp_tab.iloc[n,3])
                #
                prop_ance_temp = sum(prop_sum_w)/len(prop_sum_w)
                n_snps.append(len(prop_sum_w))
                prop_ance_win.append(prop_ance_temp)
            #
            ance_chr_win = pd.DataFrame(
            {'chr': chrm,
            'win_ini': win_ini,
            'win_mid': win_mid,
            'win_end': win_end,
            'n_snps': n_snps,
            'prop_ance' : prop_ance_win
             })
            ance_com = ance_com.append(ance_chr_win, ignore_index=True)
        #
        return(ance_com)
    #
    #
    for i in list_ind:
        print('for  ' + i + ' ....')
        print('Calculating ancestry table')
        anc_tab_ind=prop_ance(i)
        print('Calculating ancestry in windows table')
        temp_tab_ind=prop_ance_windows(anc_tab_ind, win_size)
        print('Ancestry in windows done')
        ance_ind_temp=pd.DataFrame(temp_tab_ind['prop_ance'])
        if list_ind.index(i) == 0:
            prop_anc_ind=pd.merge(temp_tab_ind,ance_ind_temp,right_index=True, left_index=True)
        else:
            prop_anc_ind=pd.merge(prop_anc_ind,ance_ind_temp,right_index=True, left_index=True)
            prop_anc_ind=prop_anc_ind.rename(columns={ prop_anc_ind.columns[-1]: i })
        print('Next individual_____________')
    #
    #
    return(prop_anc_ind)

ance_tab_win=tab_prop_ance_win(list_ind, 250000)

ance_tab_win.to_csv('ance_tab_win_1.tsv', sep='\t')


###In python
###then make any fake table of ancestry
sub_win_0=anc_win_9.iloc[:,1:7]
##and drop the last line
sub_win_0 = sub_win_0.drop(columns="prop_ance_x")

### and then paste the forst column 
anc_win_com = pd.merge(sub_win_0,sub_win_1, right_index=True, left_index=True)


### and then save
anc_win_com.to_csv('anc_win_com_named.tsv', sep='\t')


