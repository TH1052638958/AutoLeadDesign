# A pipeline for AutoLeadDesign
## Get started
1. Clone github repository   
```bash
git clone "repo-url"
```
2. Build environment
```bash
cd AutoLeadDesign/
```
```bash
conda env create --name lmlf --file=environment.yaml
```
```bash
conda activate lmlf
```
3. Prepare the target proteins and ligand to identy pocket(For example: datasets/8UOB_chainA_protein.pdb and 8UOB_chainA_ligand.pdb) and init molecules (For example: datasets/init.csv).
 Notise: The target protein and ligand must be named as xxx_chainA_protein.pdb and xxx_chainA_ligand.pdb when you are running for your own targets.

4. Prepare the docking configration file(For example: config/config_smina.yaml) 
   
   
 5. Make experiment dir
```bash
mkdir run_test
```  
6. Update the init SMILES to exp path
For example:run_test/init.csv
7. Init Docking
```bash
python init.py --exp-path run_test --protein-name 14GS --init init.csv --cfg-smina config/config_smina.yaml
```






8. Prepare your LLM API keys to communicate with LLMs and get your API_KEYS and PLANTFORM_URL


9. Run
```bash
 python doublefeedback_fragment_two_stage.py --exp-path run_test --protein-name 14GS --num-molecules 100 --num-generations 20 --threshold 7 --temperature 0.9 --model-engine deepseek-v4-flash --api-key xxx(your API keys) --plantform-url https://api.deepseek.com --cfg-smina config/config_smina.yaml```
```
