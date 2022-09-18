

# 设置工作目录
# pands 数据框
import os 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import MACCSkeys
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

# import deepchem as dc

try:
    os.mkdir("data")
except FileExistsError:
    pass

try:
    os.mkdir("figures")
except FileExistsError:
    pass

target = "CDK"

# 读取excel
df = pd.read_excel("./excel/activity_smiles_{0}.xlsx".format(target))

# 获取列名
columns = df.columns.tolist()

# 去重
df2 = df.drop_duplicates(subset=["Smiles"], keep="first")
print(df2.shape)

# 重置索引
df2.index = range(len(df2))
print(df2.shape)

# smiles = df2['Smiles'].tolist()
mol_list = [Chem.MolFromSmiles(i) for i in df2['Smiles'] ]


def calcMolDesAll(mol):
    # mol = Chem.MolFromSmiles('c1ccccc1C(=O)O')
    var_list = [x[0] for x in Descriptors._descList]
    calculator = MoleculeDescriptors.MolecularDescriptorCalculator(var_list)
    tuple1 = calculator.CalcDescriptors(mol)
    # dic = dict(zip(var_list,value_list))
    # 数据框
    df = pd.DataFrame(tuple1).T
    df.columns = var_list
    return df


# ECFP MACCS
def calcMolFP(mol_list):
    # mol = Chem.MolFromSmiles('c1ccccc1C(=O)O')
    # 计算分子指纹
    mol_list = [Chem.AddHs(i) for i in mol_list]
    ECFP = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for mol in mol_list]
    MACCS = [MACCSkeys.GenMACCSKeys(mol) for mol in mol_list]
    ECFP_array = np.asarray(ECFP, dtype=float)
    # ECFP_array.shape    
    MACCS_array = np.asarray(MACCS, dtype=float)
    # MACCS_array.shape    
    arr = np.concatenate([ECFP_array, MACCS_array],axis=1)    
    df = pd.DataFrame(arr)
    return df


var_list = [x[0] for x in Descriptors._descList]
DF = pd.DataFrame(columns=var_list)
for mol in mol_list:
    df = calcMolDesAll(mol)
    DF = pd.concat([DF, df], axis=0)


DF.index = range(len(DF))
print(DF.shape)


FP = calcMolFP(mol_list)
FP.index = range(len(FP))


var_names = ['ECFP' + str(i) for i in range(1,1025)] + [ 'MACCS' + str(i) for i in range(1,168)]
FP.columns = var_names
print(FP.shape)

final = pd.concat([df2, DF, FP], axis=1)
print(final.shape)


# 这里选择分子量小于1000的小分子，排除大分子量多肽, 
# 否则后续计算出来的Ipc参数过大导致模型计算出错
final = final.loc[final["MolWt"] <=1000, :]


final.to_excel("./data/final_analysis_{0}_0703.xlsx".format(target), index=False)


# 后期想构建自己的数据框，可以安装 MySQL
# from rdkit.ML.Descriptors import MoleculeDescriptors
# var_list = ["MaxEStateIndex"]
# calculator = MoleculeDescriptors.MolecularDescriptorCalculator(var_list)
# des_list = calculator.GetDescriptorSummaries()
# mol = mol_list[0]
# tuple1 = calculator.CalcDescriptors(mol)











