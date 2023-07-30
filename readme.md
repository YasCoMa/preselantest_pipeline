# PreSelAnTest - Epitope prediction, selection and multiple model based antigenicity test

Python3 pipeline to predict t-cell epitopes or parse b-cell bepipred results, filter epitopes and perform multiple-model antigenicity test.

## Summary

We have developed a comprehensive pipeline for prediction, simple curation and antigenicity prediction using multiple models. This pipeline contains the following functions: 
(1) Prediction of t-cell epitopes with netmhcpanII, with desired hla alleles, filtering the input fasta to select aa sequences according to this tool specifications; it also parse the results to organize the input for the other steps of the pipeline. Alternatively, it parses the fasta files aleady predicted by the b-cell bepipred3 predictor for linear epitopes. 
(2) Curation of these epitopes according to rank percentile of binding affinity (only t-cell), promiscuity of mhc alleles (only t-cell), overlapping with iedb epitopes (only t-cell), overlapping in protegen database bacterial protein sequences, and overlapping of epitopes with human proteins.
(3) Antigenicity prediction for these epitopes according to [paprec workflow](https://github.com/YasCoMa/papc_pipeline.git)
            
## Requirements:
* Python packages needed:
	- pip3 install numpy
	- pip3 install sklearn
	- pip3 install pandas
	- pip3 install matplotlib
	- pip3 install statistics
	- pip3 install python-Levenshtein
	- pip3 install boruta
	- pip3 install joblib
    
## Usage Instructions
### Preparation:
1. ````git clone https://github.com/YasCoMa/preselantest_pipeline.git````
2. ````cd preselantest_pipeline````
3. Umcompress trained_models.tar.xz
4. The main input file is the configuration in json, there is an example for b-cell (config_carb.json) and for t-cell (config.json). The second parameter of each function file is the step used in (1) and (2). Only in the prediction_s1.py, there is a third line command argument corresponding to the nemhcpanII tool path, only required if you are using for t-cell epitope prediction

### (1) Run Prediction of epitopes or parsing of prediction results:
- Configuration variables:  
    - folder_in: path to the working directory
    - cell_type: cell type for the epitopes (t-cell or b-cell), value is b or t
    - hlas: it is used only for t-cell, list of mhc alleles
    - bepipred3_output: it is used only for b-cell, bepipred3 fasta output (must be in the working directory)
- Run all:
    - ````python3 prediction_s1.py config_carb.json 0 ````
- Run only prediction:
    - ````python3 prediction_s1.py config.json 1 /path/to/nemhcpanII ````
- Run bepipred3 parsing:
    - ````python3 prediction_s1.py config_carb.json 1 ````
- Run t-cell results parsing:
    - ````python3 prediction_s1.py config.json 2 /path/to/nemhcpanII ````

### (2) Run epitope curation:
1. Configuration variables:  
    - folder_in: path to the working directory
    - cell_type: cell type for the epitopes (t-cell or b-cell), value is b or t
    - threshold_sim_iedb: it is used only for t-cell, threshold of similarity between iedb and prediction epitopes 
    - threshold_alleles: it is used only for t-cell, threshold of promiscuity of a epitope ith affinity predicted for multiple mhc alleles 
    - threshold_rank: it is used only for t-cell, threshold of percentile rank to filter results table from netmhcpanII
2. Run selection:
    - ````python3 curation_s2.py config_carb.json 1 ````
3. Run vaxijen2 summary output table parsing (only check Summary Mode in web server, copy the results in the webpage and save as vaxijen_{epitope|protein}_results.txt in {folder_in}/), this step must be made after the selection process:
    - ````python3 curation_s2.py config_carb.json 2 ````

### (3) Run epitope antigenicity prediction:
- Configuration variables:  
    - folder_in: path to the working directory
    - cell_type: cell type for the epitopes (t-cell or b-cell), value is b or t
    - control_epitope: it is not mandatory, only use if you want to compare the results with vaxijen (if you used step 2 of curation_s2.py, the value should be vaxijen_epitope_prediction.tsv)
    - control_protein: it is not mandatory, only use if you want to compare the results with vaxign-ml
    - For the control files, you may compare with any other tool, you just have to be sure to format according to the model, using the columns:
        - pepid: identifier for epitope/protein
        - sequence: aa sequence
        - value: original score of the tool
        - class: classification in binary number (1 or 0)
        - class_name: antigen or non-antigen
- Run selection:
    - ````python3 antigenicity_prediction_s3.py config_carb.json ````

## Reference

## Bug Report
Please, use the [Issues](https://github.com/YasCoMa/preselantest_pipeline/issues) tab to report any bug.
