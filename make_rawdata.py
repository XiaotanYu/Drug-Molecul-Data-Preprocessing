
# pands 数据框
import os 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import MACCSkeys
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

import deepchem as dc

try:
    os.mkdir("data")
except FileExistsError:
    pass

try:
    os.mkdir("figures")
except FileExistsError:
    pass

target = 'CDK'


# 读取excel
df = pd.read_csv("./excel/{0}.csv".format(target), sep=';')
columns = df.columns.tolist()
df2 = df.loc[:, ['Molecule ChEMBL ID',  
                 'Standard Type', 
                 'Standard Relation', 
                 'Standard Value',                  
                 'Standard Units', 
                 'Target Organism',
                 'Smiles',
                 ] 
             ]

# 排序
df2.sort_values(by=["Molecule ChEMBL ID", "Standard Value"], ascending=True)

# 去重
df2 = df2.drop_duplicates(subset=["Molecule ChEMBL ID", "Standard Value"], keep="first")
df2 = df2.drop_duplicates(subset=["Molecule ChEMBL ID"], keep="first")


# 定义计算函数，衍生新的变量 如: activity pIC50
def classifier(x):
    if x <= 1000:
        return 1
        
    elif x>1000:
        return 0

def calcpIC50(x):
    return -np.log10(x*10**(-9))


df2['activity'] = df2['Standard Value'].apply(lambda x: classifier(x))
df2['pIC50'] = df2['Standard Value'].apply(lambda x: calcpIC50(x))


df3 = df2.sort_values(by=["Standard Value"], ascending=True)
df3.index = range(len(df3))


# 输出做分类模型的数据
selections = ['Smiles', 'activity']
df_cls = df3[selections]
df_cls.to_excel('./excel/activity_smiles_{0}.xlsx'.format(target), index=False)


# 输出做回归模型的数据
selections = ['Smiles', 'pIC50']
df_reg = df3[selections]
df_reg.to_excel('./excel/pIC50_smiles_{0}.xlsx'.format(target), index=False)







