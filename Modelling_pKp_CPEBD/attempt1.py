#!/home/barbara/miniconda/envs/test/bin/python

"""
Predict properties

The used database is the "CPE DB"- the curated database of chemical penetration enhancers - http://intbio.org/cpedb/

We're going to predict the missing values of logKp
"""
import pandas as pd

#Import dataframe
raw=pd.read_csv('cpe_db.csv')

#Clean dataframe
tmp1=raw.drop('Unnamed: 0', axis=1)

#Rename tags

dic_tags=tmp1[['tag1', 'tag2', 'tag3', 'tag4', 'tag5', 'tag6', 'tag7']].iloc[1].to_dict()
tmp2=tmp1.rename(columns=dic_tags)

#Drop uninteresting columns
tmp3=tmp2.drop(['DRUGBANK','IUPAC NAME', 'CAS','CAS (cactus)', 'Drug status', 'DOI', 'SCAFFOLD name', 'Scaffold SMILES', 'Reference'], axis=1)

#Get dummies 
tmp4=pd.get_dummies(tmp3[['amide', 'amine','alkyl-amine','aromatic amine', 'aromatic nitrogen', 'basic nitrogen', 'acidic oxygen']])
tmp3=tmp3.drop(['amide', 'amine','alkyl-amine','aromatic amine', 'aromatic nitrogen', 'basic nitrogen', 'acidic oxygen'], axis=1)
tmp5=pd.concat([tmp3,tmp4],axis=1)
tmp6=pd.get_dummies(tmp5['CPE CLASS'])
tmp5=tmp5.drop('CPE CLASS', axis=1)
tmp7=pd.concat([tmp5,tmp6],axis=1)

#Clean smiles
tmp8=tmp7['SMILES'].apply(lambda x:x.strip('Copy'))
tmp7=tmp7.drop('SMILES',axis=1)
tmp9=pd.concat([tmp8, tmp7],axis=1)

c1.to_csv('training_set.csv')

#Select lines with logKp values and cleaning
tmp10=tmp9[tmp9['logKp'].notna()]
tmp11=tmp10['logKp'].apply(lambda x: float(x.replace(',','.')))
tmp10=tmp10.drop('logKp', axis=1)



clean1=pd.concat([tmp11,tmp10],axis=1)
clean1.to_csv('training_set.csv')

##########Manually clean smiles####falta ver tautomeros e aromaticos
#why the morganfingerprint
from rdkit.Chem  import  AllChem
from rdkit  import  Chem

new=pd.read_csv('training_set.csv')
new.rename(columns={'Unnamed: 0':'DBMolID'})
tmps=new['SMILES'].apply(lambda x: Chem.MolFromSmiles(x))
tmps2=tmps.apply(lambda x: AllChem.GetMorganFingerprintAsBitVect(x,2,nBits=1024))
tmps2.name='Vectorized_Smiles'
training=pd.concat([new, tmps2],axis=1)



clean_test1=tmp9[tmp9['logKp'].isna()].drop('logKp',axis=1)

##########Manually clean smiles


#Take care of smiles
