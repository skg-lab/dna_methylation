#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: yyasumizu
# @Date: 2019-03-15
# @Last Modified time: 2019-03-15

'''
usage : python hemi_methyl.py test/posi.txt test/nega.txt test/hemi.out.txt -t 8
'''

import subprocess
from io import StringIO
import os
import argparse    # 1. argparseをインポート

import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description='calculate p value for bins of positive and negative strand bismark.')    # 2. パーサを作る

# 3. parser.add_argumentで受け取る引数を追加していく
parser.add_argument('f_pos', help='positive strand')    # 必須の引数を追加
parser.add_argument('f_neg', help='negative strand')
parser.add_argument('f_out', help='output file')
parser.add_argument('-t', '--threads', type=int, default=4)
parser.add_argument('--ath', type=float, default=0.01, help='alpha for Bonferroni correction')

parser.add_argument('-p', '--path', default='methyl_diff/methyl_diff', help='path to methyl_diff')
parser.add_argument('--tmpdir', default='tmp', help='tmp dir')# よく使う引数なら省略形があると使う時に便利

args = parser.parse_args()    # 4. 引数を解析

f_pos = args.f_pos
f_neg = args.f_neg
f_out = args.f_out
threads = args.threads
a_th = args.ath
tmp_root = args.tmpdir
path_methyl_diff = args.path

print('input positive strand file : ', f_pos)
print('input negative strand file : ', f_neg)
print('output file : ', f_out)
print('threads : ', threads)

print('methyl_diff path : ', path_methyl_diff)
print('tmp dir : ', tmp_root)


df = pd.read_csv(f_pos, sep="\t", header=None)
df.columns=['chr', 'start', 'end', 'P_me', 'P_de', 'P_low', 'P_ratio', 'P_high']

df2= pd.read_csv(f_neg, sep="\t", header=None)
df2.columns=['chr', 'start', 'end', 'N_me', 'N_de', 'N_low', 'N_ratio', 'N_high']

df3= pd.merge(df, df2, on=["chr", "start", "end"])

print('loaded files')

index_split = np.array_split(df3.index, threads)

os.makedirs(tmp_root, exist_ok=True)

procs=[]
for i in range(threads):
    proc=subprocess.Popen(path_methyl_diff, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    procs.append(proc)


list_out = []
for i,proc in enumerate(procs):
    df3.loc[index_split[i], ['P_me', 'P_de', 'N_me', 'N_de']].to_csv(
        tmp_root+'/tmp_meth.{}.txt'.format(i), sep=' ', header=None, index=None)
    data = open(tmp_root+'/tmp_meth.{}.txt'.format(i), 'r').read()
    stdout_data, stderr_data = proc.communicate(data.encode('utf-8'))

    list_out.append(pd.read_csv(StringIO(stdout_data.decode('utf-8')), header=None))
    os.remove(tmp_root+'/tmp_meth.{}.txt'.format(i))

if os.listdir(tmp_root) == []:
    os.rmdir(tmp_root)

df3['Pvalue_PN'] = list(pd.concat(list_out).reset_index(drop=True)[0])
df3['Pvalue_NP'] = 1 - df3['Pvalue_PN']

flag =[]
thresh = a_th/df3.shape[0]

for pos,row in df3.iterrows():
    if row["Pvalue_PN"] < thresh:
        flag.append(1)
    elif row["Pvalue_NP"] < thresh:
        flag.append(2)
    else:
        flag.append(0)

df3["flag"] = flag

df3.to_csv(f_out, sep="\t", index=False)

print('completed')
