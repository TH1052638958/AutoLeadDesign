#!/usr/bin/env bash

set -Eeuo pipefail
#检查报错
usage() {
    echo "Useage: $0 <Begin Index> <End Index><API_KEY>" >&2
    echo "Example: $0 1 10 your api key" >&2
}

if [[ $# -ne 3 ]]; then
    usage
    exit 1
fi

start="$1"
end="$2"
api_key="$3"

if ! [[ "$start" =~ ^[0-9]+$ && "$end" =~ ^[0-9]+$ ]]; then
    echo "Error" >&2
    usage
    exit 1
fi



if (( start > end )); then
    echo "Error: Bagin>End" >&2
    exit 1
fi


#main
work_dir="$PWD"
crossdocked_dir="$work_dir/crossdocked-100"
datasets_dir="$work_dir/datasets"




#检查报错
for required_path in \
    "$datasets_dir" \
    "$crossdocked_dir"\
    "$work_dir/init.csv" \
    "$work_dir/init.py" \
    "$work_dir/run.py"; do
    if [[ ! -e "$required_path" ]]; then
        echo "Error：cannt find $required_path" >&2
        exit 1
    fi
done


#开始迭代
for (( number = start; number <= end; number++ )); do
    source_pdb="$crossdocked_dir/${number}-protein.pdb"
    source_ligand="$crossdocked_dir/${number}-ligand.sdf"
    renamed_pdb="$datasets_dir/${number}_chainA_protein.pdb"
    renamed_ligand="$datasets_dir/${number}_chainA_ligand.pdb"
    number_dir="$work_dir/random_${number}"
    reference_dir="$work_dir/reference_${number}"

    echo "[$number] random begin"
    cp -- "$source_pdb" "$renamed_pdb"
    cp -- "$source_ligand" "$renamed_ligand"

    mkdir -p -- "$reference_dir"
    mkdir -p -- "$number_dir"
    #cp -- "$work_dir/datasets/init.csv" "$reference_dir/init.csv"
    #将reference SMILES 加入到init.csv中，保存到reference dir中
    python init_reference.py --out-path "$reference_dir/init.csv" --protein-name $number
    #将init+reference全对接
    python init.py --exp-path $reference_dir --protein-name $number
    #将除reference以外的结果保存到random路径下
    python prepare_random.py --exp-path $number_dir --init-path $reference_dir


    Threshold=$(python threshold.py  --exp-path $number_dir )
    python doublefeedback_fragment_two_stage.py --exp-path $number_dir --protein-name $number --threshold $Threshold  --num-molecules 100 --num-generations 20 --temperature 1.5 --model-engine deepseek-v4-flash --api-key $api_key --plantform-url https://api.deepseek.com/v1 --cfg-smina config/config_smina.yaml

    echo "[$number] random finished"

    echo "[$number] reference begin"
    python doublefeedback_fragment_two_stage.py --exp-path $reference_dir --protein-name $number --threshold $Threshold  --num-molecules 100 --num-generations 20 --temperature 1.5 --model-engine deepseek-v4-flash --api-key $api_key --plantform-url https://api.deepseek.com/v1 --cfg-smina config/config_smina.yaml



done

echo "Done：$start to $end"