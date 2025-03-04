

import pandas as pd




# note
# that is the same with above, so no need to check 
#epic_df=pd.read_csv('resources/GPL21145_MethylationEPIC_15073387_v-1-0.csv',skiprows=7)


# step 1:

# df=pd.read_csv('resources/GPL13534_HumanMethylation450_15017482_v.1.1.csv',skiprows=7,dtype={'CHR':str})
# df=df.loc[:,['IlmnID','CHR','MAPINFO']]
# df=df.dropna()
# df['MAPINFO']=df['MAPINFO'].astype(int)
# df['CHR']="chr"+df['CHR'].astype(str)
# df['start']=df['MAPINFO']-1
# df['end']=df['MAPINFO']
# df['score']=0
# df=df.loc[:,['CHR','start','end','score','IlmnID']]
# df['start']=df['start'].astype(int)
# df['end']=df['end'].astype(int)
# df.to_csv('hg19.coordinates.bed',index=False,header=None,sep='\t')


# step 2:
# run in bash
#ml crossmap
#crossmap bed resources/hg19ToHg38.over.chain.gz hg19.coordinates.bed hg19tohg38.liftover.coordinates.bed 


# step 3:
# get the overlaps with ref_atlas CpGs
df=pd.read_csv('resources/reference_atlas.csv')
cpg_list=df['CpGs'].values

df=pd.read_csv('resources/hg19tohg38.liftover.coordinates.bed',header=None,sep='\t')
df=df.loc[:,[4,0,2]]
df=df.loc[ (df[4].isin(cpg_list)),]
df.columns=['cgid','chromosome','pos']
df.to_csv('resources/ref_atlas_CpGs_hg38.tsv',sep='\t',index=False)
