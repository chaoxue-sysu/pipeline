## imputation
base=/home/xc/tmp/scz/SCZ_TaoLi/SCZ_Sichuen/GWAS_data/SCZ_CaseCtl
ref=/public1/data/resources/db_tools/phasing_imputation/1000GP_Phase3
progam_base=/home/xc/local/project/QTMediate/code
nohup python3 ${progam_base}/expredix/coopration/pipline/PhaseAndImpute.py ${base}/chr_haps ${ref} ${ref} ${base}/chr_impute /public1/data/resources/ref_genome/GRCh37/human_g1k_v37.chr.length 40 >>${progam_base}/expredix/coopration/pipline/PhaseAndImpute.py.log 2>&1 &

## convert impute2 (.gen) format to predixcan (dosage) format
nohup python3 ${progam_base}/expredix/coopration/predixcan/convert_GEN2VCF.py ${base}/chr_impute ${base}/chr_impute_dosage_predixcan 30 >>${progam_base}/expredix/coopration/predixcan/convert_GEN2VCF.py.log 2>&1 &

## predict expression by predixcan and combine the expression profile
nohup python3 ${progam_base}/expredix/coopration/predixcan/predict.py ${base}/chr_impute_dosage_predixcan ${base}/expression_predixcan ${base}/expression_predixcan_combine/combine.tsv.gz ${base}/scz_id_sex.txt 0 0 20 >> ${progam_base}/expredix/coopration/predixcan/predict.py.log 2>&1 &
