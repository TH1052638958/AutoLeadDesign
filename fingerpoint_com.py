smiles_list1=[
'Cc1ccc(cc1C(=O)N[C@H](C)c2cccc3c2cccc3)N',
'C[C@H](c1cccc2c1cccc2)N(C)Cc3ccc(cc3)O',
'C[C@H](c1cccc2c1cccc2)N(C)Cc3c[nH]c4c3cccc4',

'Cc1ccc(cc1C(=O)N[C@H](C)c2cc(cc(c2)c3cnn(c3)C)c4cnn(c4)C)OCCN(C)C',
'Cc1ccc(cc1C(=O)N[C@H](C)c2cc(cc(c2)c3ccc(s3)CN)c4cccs4)OCCN(C)C',
'Cc1ccc(cc1C(=O)N[C@H](C)c2cc(cc(c2)c3cncs3)c4cnn(c4)C)OCCN(C)C',
'Cc1ccc(cc1C(=O)N[C@H](C)c2cc(cc(c2)c3cnn(c3)C(F)F)c4cnn(c4)C)OCCN(C)C',
'Cc1ccc(cc1C(=O)N[C@H](C)c2cc(cc(c2)c3cnn(c3)C4CC4)c5cnn(c5)C)OCCN(C)C',
'Cc1ccc(cc1C(=O)N[C@H](C)c2cc(cc(c2)c3cnn(c3)C(C)C)c4cnn(c4)C)OCCN(C)C',
'Cc1ccc(cc1C(=O)N[C@H](C)c2cc(cc(c2)c3cnn(c3)C(F)F)c4cnn(c4)COC)OCCN(C)C',
'CCn1ccc(n1)c2cc(cc(c2)[C@@H](C)NC(=O)c3cc(ccc3C)OCCN(C)C)c4cnn(c4)C'
]
smiles_list2=[
'CC1=C(C=C(C=C1)N)C(=O)N[C@H](C)C2=CC=CC3=CC=CC=C32',
'C[C@H](C1=CC=CC2=CC=CC=C21)N(C)CC3=CC=C(C=C3)O',
'C[C@H](C1=CC=CC2=CC=CC=C21)N(C)CC3=CNC4=CC=CC=C43',

'CC1=C(C=C(C=C1)OCCN(C)C)C(=O)N[C@H](C)C2=CC(=CC(=C2)C3=CN(N=C3)C)C4=CN(N=C4)C',
'CC1=C(C=C(C=C1)OCCN(C)C)C(=O)N[C@H](C)C2=CC(=CC(=C2)C3=CC=C(S3)CN)C4=CC=CS4',
'CC1=C(C=C(C=C1)OCCN(C)C)C(=O)N[C@H](C)C2=CC(=CC(=C2)C3=CN=CS3)C4=CN(N=C4)C',
'CC1=C(C=C(C=C1)OCCN(C)C)C(=O)N[C@H](C)C2=CC(=CC(=C2)C3=CN(N=C3)C(F)F)C4=CN(N=C4)C',
'CC1=C(C=C(C=C1)OCCN(C)C)C(=O)N[C@H](C)C2=CC(=CC(=C2)C3=CN(N=C3)C4CC4)C5=CN(N=C5)C',
'CC1=C(C=C(C=C1)OCCN(C)C)C(=O)N[C@H](C)C2=CC(=CC(=C2)C3=CN(N=C3)C(C)C)C4=CN(N=C4)C',
'CC1=C(C=C(C=C1)OCCN(C)C)C(=O)N[C@H](C)C2=CC(=CC(=C2)C3=CN(N=C3)C(F)F)C4=CN(N=C4)COC',
'CCN1C=CC(=N1)C2=CC(=CC(=C2)[C@@H](C)NC(=O)C3=C(C=CC(=C3)OCCN(C)C)C)C4=CN(N=C4)C'
]

from tools.tools import get_fingerpoint
method=['MACCS','Topological','Morgan','Avalon']
for m in method:
    fp1 = (get_fingerpoint(smiles_list1, m))
    fp2 = (get_fingerpoint(smiles_list2, m))
    for i in range(len(fp1)):
        if fp1[i] !=fp2[i]:
            print(f'Smiles {smiles_list1[i]} and {smiles_list2[i]} have dif fp at {m}')

