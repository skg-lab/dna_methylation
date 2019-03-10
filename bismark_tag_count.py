import sys
import pandas as pd
import numpy as np

input_bismark_file = sys.argv[1]
input_bed_file = sys.argv[2]
out_put_file=sys.argv[3]

df_bismark = pd.read_csv(input_bismark_file, sep='\t', names=('chr', 'start', 'end', 'methyl_per','methyl','de_methyl'))

drop_row = ['end', 'methyl_per']
df_bismark_drop=df_bismark.drop(drop_row, axis=1)

df_bed = pd.read_csv(input_bed_file, sep='\t',names=('chr','start','end'))


list_df_bed=[]
df_bismark_list=[]

list_chr=['chr1','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19',
'chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chrM','chrX','chrY']

for i in range(len(list_chr)):
    df_bed_chr=df_bed[df_bed['chr'] == list_chr[i]]
    list_df_bed.append(df_bed_chr)

for i in range(len(list_chr)):
    df_bismark_chr=df_bismark_drop[df_bismark_drop['chr'] == list_chr[i]]
    df_bismark_list.append(df_bismark_chr)

df_tag_bed_list=[]

for i in range(len(list_chr)):
    bed=list_df_bed[i]
    bismark=df_bismark_list[i]
    bed_list=bed.values.tolist()

    for j in range(len(bed_list)):
        x=bed_list[j][1]
        y=bed_list[j][2]

        bismark_query=bismark.query("@x<= start<=@y")
        bed_list[j].append(bismark_query['methyl'].sum())
        bed_list[j].append(bismark_query['de_methyl'].sum())
    tag_bed = pd.DataFrame(bed_list,columns=('chr','start','end','methyl_sum','de_methyl_sum'))
    df_tag_bed_list.append(tag_bed)

all_tag_bed=pd.concat(df_tag_bed_list)
all_tag_bed.to_csv(out_put_file, mode='a',sep='\t',index=False)
