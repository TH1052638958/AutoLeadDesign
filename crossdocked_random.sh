#!/usr/bin/env bash

set -Eeuo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
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







#开始迭代
for (( number = start; number <= end; number++ )); do
    source_pdb="$crossdocked_dir/${number}-protein.pdb"
    source_ligand="$crossdocked_dir/${number}-ligand.sdf"
    renamed_pdb="$datasets_dir/${number}_chainA_protein.pdb"
    renamed_ligand="$datasets_dir/${number}_chainA_ligand.sdf"
    renamed_pdb_com="/home/th2024/compare/mol_opt/datasets//${number}_chainA_protein.pdb"
    renamed_ligand_com="/home/th2024/compare/mol_opt/datasets//${number}_chainA_ligand.sdf"
    number_dir="$work_dir/random_${number}"
    reference_dir="$work_dir/reference_${number}"

    echo "[$number] random begin"
    cp -- "$source_pdb" "$renamed_pdb"
    cp -- "$source_ligand" "$renamed_ligand"

    mkdir -p -- "$reference_dir"
    mkdir -p -- "$number_dir"

    #将reference SMILES 加入到init.csv中，保存到reference dir中
    #切换环境
    conda activate lmlf
    #python init_reference.py --out-path "$reference_dir" --protein-name $number
    #将init+reference全对接
    #python init.py --exp-path $reference_dir --protein-name $number
    #将除reference以外的结果保存到random路径下
    #python prepare_random.py --exp-path $number_dir --init-path $reference_dir


    Threshold=$(python threshold.py  --exp-path $number_dir )
    #python doublefeedback_fragment_two_stage.py --exp-path $number_dir --protein-name $number --threshold $Threshold  --num-molecules 100 --num-generations 20 --temperature 1.5 --model-engine deepseek-v4-flash --api-key $api_key --plantform-url https://api.deepseek.com/v1 --cfg-smina config/config_smina.yaml

    echo "[$number] random finished"

    echo "[$number] reference begin"
    #python doublefeedback_fragment_two_stage.py --exp-path $reference_dir --protein-name $number --threshold $Threshold  --num-molecules 100 --num-generations 20 --temperature 1.5 --model-engine deepseek-v4-flash --api-key $api_key --plantform-url https://api.deepseek.com/v1 --cfg-smina config/config_smina.yaml
    #切换环境
    conda activate com
    cp -- "$source_pdb" "$renamed_pdb_com"
    cp -- "$source_ligand" "$renamed_ligand_com"
    #运行reinvent
    python /home/th2024/compare/mol_opt/run.py --protein $number --method reinvent --output_dir "reinvent_${number}" --max_oracle_calls 2000 --oracles QED
    #运行SMILES GA random
    python /home/th2024/compare/mol_opt/run.py --protein $number --method smiles_ga --smi_file "${number_dir}/init.csv" --output_dir "smilesga_${number}_random" --max_oracle_calls 2000 --oracles QED
    #运行SMILES GA reference
    python /home/th2024/compare/mol_opt/run.py --protein $number --method smiles_ga --smi_file "${reference_dir}/init.csv" --output_dir "smilesga_${number}_reference" --max_oracle_calls 2000 --oracles QED




done

echo "Done：$start to $end"