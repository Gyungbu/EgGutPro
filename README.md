# EgGutPro : Calculate the index and percentile ranks of human - gut microbiome.

![Python](https://img.shields.io/badge/Python-v3.9.0-blue.svg?style=flat&logo=python)&nbsp;
![Pandas](https://img.shields.io/badge/pandas-v2.0.1-blue.svg?style=flat&logo=pandas)&nbsp;
![Numpy](https://img.shields.io/badge/NumPy-v1.24.3-blue.svg?style=flat&logo=numpy)&nbsp;
![Scipy](https://img.shields.io/badge/SciPy-v1.10.1-blue.svg?style=flat&logo=scipy)&nbsp;
![scikit-bio](https://img.shields.io/badge/scikit_bio-grey.svg?style=flat&logo=scikit-bio)&nbsp;

## Installation

You can install the EgGutPro with following command.
	
	git clone https://github.com/Gyungbu/EgGutPro.git
 
The list of required packages for `script` is shown in the `requirements.txt` file. When you run the command below, the required packages will be downloaded. (version : `python 3.9.0`)
	
	conda create -n env_EgGutPro
	conda activate env_EgGutPro
	conda install pip  
	conda install python=3.9.0
	pip install -r ./EgGutPro/requirements.txt 

# EGgutPro_update_mrs : (Optional) Update Reference Data
## How to use

### 1. Prepare Merged Proportion File
Place the csv file of your Merged Proportion File in the `./EgGutPro/input/` folder.

Caveats: 

1. The first column of the Merged Proportion File is in the order of taxa, diversity, observed, and microbiome list.
2. From the second column of the Merged Proportion File, the sample name, the diversity value of the sample, the number of species in the sample, and the relative abundance value of each microbiome should be listed in the order.


### 2. Run EGgutPro_update_mrs
To run EGgutPro_update_mrs,
 
Run the command below:
  
    python ./EgGutPro/EGgutPro_update_mrs.py {path_exp}
    ### ex) python EGgutPro_update_mrs.py "/home/kbkim/EgGutPro/input/EGgutPro_mirror_output_3175.csv"
   
    
    
When EGgutPro_update_mrs is executed as above, the file `EGgutPro_db_abundance.xlsx`, `EGgutPro_mrs_db.xlsx`, `EGgutPro_percentile_rank_db.csv` will be created or modified in the `./EgGutPro/input/` folder.
And the file `EGgutPro_mrs_hist.png` will be created or modified in the `./EgGutPro/output/` folder.


# EGgutPro_percentile_rank : Calculate the index and percentile ranks of human - gut microbiome.
## How to use

### 1. Prepare Input data
Place the csv file of your proportion file in the `./EgGutPro/input/` folder.
Caveats: 

1. The first column of the Merged Proportion File is in the order of taxa, diversity, observed, and microbiome list.
2. From the second column of the Merged Proportion File, the sample name, the diversity value of the sample, the number of species in the sample, and the relative abundance value of each microbiome should be listed in the order.
3. In EGgutPro_percentile_rank.py, under if __name__ == '__main__': enter the path of the proportion file you want to analyze in the 'path_exp' value and save it.

### 2. Run EGgutPro_percentile_rank
To run EGgutPro_percentile_rank,
 
Run the command below:

    python ./EgGutPro/EGgutPro_percentile_rank.py
    

When EGgutPro_percentile_rank is executed as above, the file `EGgutPro_eval.csv`, `EGgutPro_percentile_rank.csv`, `EGgutPro_harmful.csv`, `EGgutPro_harmful_tot.csv`, `EGgutPro_beneficial.csv` and `EGgutPro_probio_tot.csv` will be created or modified in the `./EgGutPro/output/` folder.


