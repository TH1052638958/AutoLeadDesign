import os.path

import pandas as pd
from tools.tools import get_fingerpoint

fragment_df=pd.read_csv('/home/th2024/lmlf/LMLF-main/6IGX_two_stage_without/fragment_99.csv')
fragment=list(fragment_df.columns)
method=['MACCS','Topological','Morgan','Avalon']
for m in method:

    fp=(get_fingerpoint(fragment,m))
    fragment_dict={}
    for i in range(len(fp)):
        fp_tem=fp[i]
        fragment_tem=fragment[i]
        if fp_tem in fragment_dict:
            fragment_dict[fp_tem].append(fragment_tem)
        else:
            fragment_dict[fp_tem]=[fragment_tem]
    df_fragment = pd.DataFrame.from_dict(fragment_dict, orient='index').transpose()
    df_fragment.to_csv(os.path.join('fragment_ansys','fragment_ansys_{}.csv'.format(m)), index=False)



