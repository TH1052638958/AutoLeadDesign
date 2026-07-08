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
conda env create --name AutoLeadDesign --file=environments.yml
```
```bash
conda activate AutoLeadDesign
```
3. Prepare the target proteins and ligand to identy pocket(For example: datasets/8UOB_chainA_protein.pdb and 8UOB_chainA_ligand.pdb) and init molecules (For example: datasets/init.csv).
   
   Notise: The target protein and ligand must be named as xxx_chainA_protein.pdb and xxx_chainA_ligand.pdb when you are running for your own targets.
   
4.Prepare the docking configration file(For example: config/config_smina.yaml) 
5. Prepare your LLM API keys to communicate with LLMs and get your API_KEYS and PLANTFORM_URL 
6. Run
```bash
python doublefeedback_fragment_two_stage.py --exp-path 8UOB_test --protein-name 8UOB --num-molecules 100 --num-generations 20 --api-key $API_KEYS$ --plantform-url $PLANTFORM_URL$ --cfg-smina config/config_smina.yaml
```
