import pandas as pd
import os
import sys

df = pd.read_csv(sys.argv[1], header=0, sep='\t')

df2 = df.loc[df['comp'] == sys.argv[2]]
df2['path'] = df2['sample']
df2['path'] = sys.argv[3] + "/CpG_bwa/" + df2['path'] + ".bm_pe.deduplicated_CpG.bedGraph"
df2.to_csv(sys.argv[4], sep="\t",index=False)

df3 = df.loc[df['comp'] == sys.argv[2]]
df3['path'] = df3['sample']
df3['path'] = sys.argv[3] + "/CpG_bismark/" + df3['path'] + ".bismark_bt2_pe.deduplicated_CpG.bedGraph"
df3.to_csv(sys.argv[5], sep="\t",index=False)
