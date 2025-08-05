#!/usr/bin/env python3
# coding:utf-8

import pandas as pd
import time

df1=pd.read_excel('83.xlsx').astype(str)
print(df1)

df2=pd.read_excel('meta_sample_mb20240827.xlsx').astype(str)
print(df2)


df3 = pd.merge(df1,df2,left_on = 'sample',right_on = 'sample',how = 'left')
print(df3)
# df3=df3.dropna(axis = 0)
# print(df3)
time.sleep(300)
df3.to_csv('83mbR.xlsx',sep="\t",index=False)