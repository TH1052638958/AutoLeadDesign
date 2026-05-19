from rdkit import DataStructs,Chem
from rdkit.Chem import AllChem
import numpy as np
from rdkit.Chem import MACCSkeys


from tools.tools import get_fingerpoint
# method=['MACCS','Topological','Morgan','Avalon']
method='Morgan'

def get_fp_list(smiles_list):
    fp_list = []
    for s in smiles_list:
        try:
            fp_list.append(get_fingerpoint(s))
        except Exception as e:
            print(f'Cont gen fp : {s}',e)
            continue
    return fp_list
def get_smiles_list(file_path):
    import pandas as pd



    a = pd.read_csv(file_path)
    smiles=list(a.iloc[:,0])
    return smiles

if __name__ =='__main__':
    file_list=['TSNE_PRMT5_random+LLY_LMLF/result.csv','TSNE_PRMT5_random+LLY_ours/result.csv','TSNE_PRMT5_random+LLY_ours1/init.csv']
    #good=['NC(C(C=CC=C1)=C1N=N2)=C2C3=CC=CC(C4=NC5=C(N=N4)C6=C5C=CC=C6)=C3','NC1=CC(C(C(C2=CC=C(CC3=C(O)C(NC4=C5C=CC(O)=C4)=C5C=C3)C=C2)N6)=O)=C6C=N1']
    good = ['O=C1N[C@@H](C[C@H]2CCc3cc(Cl)ccc23)C[C@@H](C(=O)[C@H](O)C[C@H](C)N)C1','NC(C1=CC=CC(O)=C1)C(O2)C(N)C(Br)C2N3CCC4=C3C=NN=C4N']
    #good=['C1=CC=C(C=C1)[C@H]([C@@H]2[C@H]([C@H]([C@@H](O2)N3C=NC4=C(N=CN=C43)N)O)O)O']
    #good=[]
    fp=[]
    label=[]
    num=[]
    for i in range(len(file_list)):
        file_path_init = file_list[i]
        smiles_init = get_smiles_list(file_path_init)
        #debug begin
        # if i==2:
        #     smiles_init=smiles_init[:100]
        #debug end
        fp_init = get_fingerpoint(smiles_init, method)
        label_tem=i
        label_init = [label_tem for j in fp_init]
        num.append(len(label_init))
        fp=fp+fp_init
        label=label+label_init
    #add good
    fp_good=get_fingerpoint(good,method)
    fp=fp+fp_good
    fp_victor=[]

    for fp_tem in fp :
        fp_filted=[]
        for i in fp_tem:
            fp_filted.append(int(i))
        fp_victor.append(fp_filted)
    fp=fp_victor

    begin1=0
    begin2=num[0]
    begin3=num[0]+num[1]
    begin4=num[0]+num[1]+num[2]-1
    begin5=num[0]+num[1]+num[2]



    #pca
    # from sklearn.decomposition import PCA
    #
    # import matplotlib.pyplot as plt
    #
    # pca = PCA()
    # X=fp
    #
    # Xt = pca.fit_transform(X)
    # x_tsne=Xt


    # X_ours=Xt[begin1:begin2-1,0]
    # Xt_ours=Xt[begin1:begin2-1,1]
    # label_ours=['indigo' for i in X_ours]
    # X_init = Xt[begin2:begin3 - 1,0]
    # Xt_init = Xt[begin2:begin3 - 1, 1]
    # label_init = ['yellow' for i in X_init]
    # X_LMLF = Xt[begin3:,0]
    # Xt_LMLF = Xt[begin3:, 1]
    # label_LMLF = ['lightseagreen' for i in X_LMLF]
    # X1=np.concatenate((X_ours,X_init))
    # Xt1=np.concatenate((Xt_ours,Xt_init))
    # label1=label_ours+label_init
    # X2 = np.concatenate((X_LMLF, X_init))
    # Xt2 = np.concatenate((Xt_LMLF, Xt_init))
    # label2 = label_LMLF + label_init
    #
    #
    # #0:black 1:yellow 2:green
    # plot = plt.scatter(X1, Xt1,c=label1)
    # plt.legend(handles=plot.legend_elements()[0],
    #            )

    # plt.show()
    # plot = plt.scatter(X2, Xt2, c=label2)
    # plt.legend(handles=plot.legend_elements()[0],
    #            )

    # plt.show()
    #
    #
    # X_all=np.concatenate((X_ours,X_LMLF,X_init))
    # Xt_all=np.concatenate((Xt_ours,Xt_LMLF,Xt_init))
    # label_all=label_ours+label_LMLF+label_init
    #
    # plot = plt.scatter(X_all, Xt_all,c=label_all)
    # plt.legend(handles=plot.legend_elements()[0],
    #            )

    # plt.show()
    #tsne

    import matplotlib.pyplot as plt

    from sklearn.manifold import TSNE


    tsne = TSNE(n_components=2)
    x =np.array( fp )

    x_tsne = tsne.fit_transform(x)
    # #第一部分 lmlf
    X_LMLF=x_tsne[begin1:begin2-1,0]
    Xt_LMLF=x_tsne[begin1:begin2-1,1]

    label_LMLF=['orange' for i in X_LMLF]
    #第二部分 ours
    X_ours = x_tsne[begin2:begin3 - 1,0]
    Xt_ours = x_tsne[begin2:begin3 - 1, 1]
    label_ours = ['steelblue' for i in X_ours]
    #第三部分 init
    X_init = x_tsne[begin3:begin4+1,0]
    Xt_init =x_tsne[begin3:begin4+1, 1]

    label_init = ['darkred' for i in X_init]
    #第四部分 reference
    X_reference = x_tsne[begin4:begin4+1,0]
    Xt_reference =x_tsne[begin4:begin4+1, 1]
    label_reference = ['green' for i in X_reference]
    #第五部分 good
    X_good = x_tsne[begin4+1:begin4+2,0]
    Xt_good =x_tsne[begin4+1:begin4+2, 1]
    label_good = ['black' for i in X_good]
    # 第五部分 good
    X_bed = x_tsne[begin4 + 2:, 0]
    Xt_bed = x_tsne[begin4 + 2:, 1]
    label_bed = ['black' for i in X_good]

    X1=np.concatenate((X_ours,X_init))
    Xt1=np.concatenate((Xt_ours,Xt_init))
    label1=label_ours+label_init
    X2 = np.concatenate((X_LMLF, X_init))
    Xt2 = np.concatenate((Xt_LMLF, Xt_init))
    label2 = label_LMLF + label_init

    #分别展示
    # #0:black 1:yellow 2:green
    size=10
    # plot = plt.scatter(X1, Xt1,s=size,c=label1)
    # plt.legend(handles=plot.legend_elements()[0],
    #            )
    # #plt.xlim(-80,80)
    # #plt.ylim(-70,65)
    # plt.show()
    # plot = plt.scatter(X2, Xt2,s=size, c=label2)
    # plt.legend(handles=plot.legend_elements()[0],
    #            )
    #plt.xlim(-85,150)
    # #plt.ylim(-70,65)
    # plt.show()


    X_all=np.concatenate((X_LMLF,X_ours,X_init,X_reference,X_good))
    Xt_all=np.concatenate((Xt_LMLF,Xt_ours,Xt_init,Xt_reference,Xt_good))
    label_all=label_LMLF+label_ours+label_init+label_reference+label_good

    #plot = plt.scatter(X_all, Xt_all,s=size,c=label_all)
    LMLF=plt.scatter(X_LMLF,Xt_LMLF,s=size,c=label_LMLF,label='Designed by LMLF')
    ALD=plt.scatter(X_ours, Xt_ours, s=size, c=label_ours, label='Designed by AutoLeadDesign')
    Init=plt.scatter(X_init, Xt_init, s=size, c=label_init, label='Init')
    #reference
    LLY=plot = plt.scatter(X_reference, Xt_reference, s=150, c='red',marker='*',label='LLY')
    good=plot = plt.scatter(X_good, Xt_good, s=150, c='purple', marker='*',label='PRA001')
    bed = plot = plt.scatter(X_bed, Xt_bed, s=150, c='yellow', marker='*', label='PRA001')
    plt.ylabel('', fontdict={'family': 'Times New Roman', 'size': 16})
    plt.xlabel('', fontdict={'family': 'Times New Roman', 'size': 16})
    plt.legend('', frameon=False)
    # plt.legend(handles=plot.legend_elements()[0],
    #            )
    #plt.legend(handles=[Init,LLY,LMLF,ALD,good],labels=['Init','LLY','Designed (LMLF)','Designed (AutoLeadDesign)','PRA001'],loc=1, borderaxespad=0)
    #plt.xlim(-80,80)
    #plt.ylim(-70,65)
    #save
    plt.savefig('Fig_PRMT5_reference.png',dpi=1000)
    plt.show()