#! /bin/bash
IDX=$1
PROTEIN=$2
REFERENCE=$3
i=0
API="AAA"
while read line
do
	if [ "$IDX" -eq "$i" ];then
		API=$line
	fi
	i=$((i+1))
	
	

done < /home/th2024/api.txt
EXPPATH="ALD_""$REFERENCE""_""$PROTEIN"
INITPATH="init_""$PROTEIN"
mkdir $EXPPATH
cp -rf "$INITPATH"/* "$EXPPATH"

echo "Running ALD in $REFERENCE Targeting $PROTEIN"
python doublefeedback_fragment_two_stage.py --exp-path $EXPPATH --protein-name $PROTEIN --num-molecules 100 --num-generations 20 --temperature 1.5 --model-engine deepseek-chat --api-key $API --plantform-url https://api.deepseek.com/v1 --cfg-smina config/config_smina.yaml
echo "Running LMLF in $REFERENCE Targeting $PROTEIN"
EXPPATH="LMLF_""$REFERENCE""_""$PROTEIN"
mkdir $EXPPATH
cp -rf "$INITPATH"/* "$EXPPATH"
Threshold=$(python prepare_random.py --init-path $INITPATH --exp-path $EXPPATH --save 0)
python doublefeedback.py --exp-path $EXPPATH --protein-name $PROTEIN --num-molecules 100 --num-generations 20 --model-engine deepseek-chat --cfg-smina config/config_smina.yaml --api-key $API --threshold $Threshold
echo "Done"

