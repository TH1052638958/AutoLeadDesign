#! /bin/bash
IDX=$1
E=$(($IDX+1))
echo "IDX:$IDX,E:$E"
BEGIN=$(($IDX*5))
END=$(($E*5))
echo "from $BEGIN to $END"
RPATH="crossdocked-100"
DPATH="datasets"
for ((i=$BEGIN;i<$END;i++))
do
	echo "processing ID:$i"
	OPATH="init_""$i"
	mkdir $OPATH
	
	LPATH="$RPATH/$i""-ligand.sdf"
	PPATH="$RPATH/$i""-protein.pdb"
	NLPATH="$DPATH/$i""_chainA_ligand.sdf"
	NPPATH="$DPATH/$i""_chainA_protein.pdb"
	cp $LPATH $NLPATH
	cp $PPATH $NPPATH
	cp datasets/init.csv $OPATH
	python init_reference.py --out-path $OPATH --protein-name $i
  python init.py --exp-path $OPATH --protein-name $i --cfg-smina config/config_smina.yaml
done
