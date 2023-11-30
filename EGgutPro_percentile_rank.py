import os, datetime
import pandas as pd
from scipy.stats import percentileofscore, pearsonr
import sys
import math
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import scipy
from skbio.stats.composition import multiplicative_replacement, clr


#-------------------------------------------------------
# Common Function
#-------------------------------------------------------
def WriteLog(functionname, msg, type='INFO', fplog=None):
    #strmsg = "[%s][%s][%s] %s\n" % (datetime.datetime.now(), type, functionname, msg)
    #if( DEBUGMODE ): 
    #    print(strmsg)
    #else:
    #    if( fplog != None ):
    #        fplog.write(strmsg)
    
    head = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    writestr = f"[{head}][{functionname}] {msg}\n"
    #if( DEBUGMODE ):
    if( True ):
        #print(message)
        writestr = f"[{functionname}] {msg}\n"
        print(writestr)
        
    if( fplog != None ):
        fplog.write(writestr)
        fplog.flush()

def filter_species(taxon_str_for_domain):
    splited_by_taxonomy = taxon_str_for_domain.split("|")
    species_name = splited_by_taxonomy[len(splited_by_taxonomy) - 1]  # 계문강목과속종 중 종만 필터링!!!
    return species_name


def nogada_x_y_graph_v2(final_type, xpos, ypos):
    #print("GRAPH x: ", xpos, "y: ", ypos)

    if (final_type == 'D'):
        if (xpos < 8):
            xpos = 8
        if( ypos<8 ):
            ypos = 8
        if(xpos >30):
            xpos = 30
        if(ypos >30):
            ypos = 27

    if (final_type == 'B'):
        if (xpos > 45):
            xpos = 45
        if (ypos > 50):
            ypos = 50

    if (final_type == 'I'):
        if(xpos > 65):
            xpos=65
        if (ypos > 65):
            ypos = 65


    if (final_type == 'E'):
        if(xpos>85):
            xpos = 85
             
        if (ypos > 74):
            ypos = 68

    return (xpos, ypos)

def starts_with_s(value):
    return value.startswith('s__')

###################################
# MainClass
###################################
class EgGutProAnalysis:
    def __init__(self, path_exp, outdir=None, fplog=None):
        """
        Initializes a EgGutProAnalysis object.

        Parameters:
        path_exp (str): Path of Merged Proportion file to analyze.
        """
        self.path_exp = path_exp
        self.__fplog=fplog

        ## Path of Reference files
        curdir = os.path.dirname(os.path.abspath(__file__))
        self.path_ref = f"{curdir}/input/EGgutPro_mircobe_list.xlsx"       
        self.path_healthy = f"{curdir}/input/20231127_FMTDonor_OnlySp.xlsx"
        self.path_mrs_db = f"{curdir}/input/EGgutPro_mrs_db.xlsx"
        self.path_percentile_rank_db = f"{curdir}/input/EGgutPro_percentile_rank_db.csv"
        self.path_db = f"{curdir}/input/EGgutPro_db_abundance.xlsx"
        
        
        ###output
        if( outdir is not None ):
            self.outdir = outdir
        else:
            self.outdir = f"{curdir}/output"        
        
        
        ## Path of output files       
        self.path_percentile_rank_output = f"{self.outdir}/EGgutPro_percentile_rank.csv"
        self.path_eval_output = f"{self.outdir}/EGgutPro_eval.csv"
        self.path_scatterplot_output = f"{self.outdir}/EGgutPro_scatterplot.png"
        self.path_harmful = f"{self.outdir}/EGgutPro_harmful_10.csv"
        self.path_beneficial = f"{self.outdir}/EGgutPro_beneficial_10.csv"
        self.path_harmful_tot = f"{self.outdir}/EGgutPro_harmful_30.csv"
        self.path_probio_tot = f"{self.outdir}/EGgutPro_probio_19.csv"
        self.path_eval_output_proumed = f"{self.outdir}/proumed_output.csv"

        ## Dataframe of Reference files
        self.df_beta = None
        self.df_dysbiosis = None
        self.df_probio = None             
        self.df_healthy = None
        self.df_exp = None
        self.df_mrs_db = None
        self.df_percentile_rank_db = None
        self.df_db = None
        
        ## Dataframe of output files to calculate
        self.df_mrs = None
        self.df_percentile_rank = None
        self.df_eval = None
        self.df_harmful_10 = None
        self.df_beneficial_10 = None
        self.df_harmful_30 = None
        self.df_probio_19 = None
        self.df_proumed = None
        
        ## Lists used for calculation
        self.li_diversity = None
        self.li_observed = None
        self.li_new_sample_name = None
        self.li_phenotype = None
        self.li_microbiome = None
        self.li_ncbi_name = None
        self.li_probio_microbiome = None  
        
        self.observed_mean = None
        self.threshold = 0.00001

    # Load the DB file
    # df_beta : Data frame of of Phenotype-Microbiome information
    # df_exp : Data frame of Experimental result information - Abundance    
    def ReadDB(self):
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try:
            self.df_beta = pd.read_excel(self.path_ref, sheet_name = f"disease_taxa")
            self.df_dysbiosis = pd.read_excel(self.path_ref, sheet_name = f"harmful_beneficial_taxa")
            self.df_probio = pd.read_excel(self.path_ref, sheet_name = f"probio_taxa")
            
            self.df_healthy = pd.read_excel(self.path_healthy)
            self.df_mrs_db = pd.read_excel(self.path_mrs_db, index_col=0) 
            self.df_percentile_rank_db = pd.read_csv(self.path_percentile_rank_db)
            self.df_db = pd.read_excel(self.path_db)

            self.df_beta.rename(columns = {"Disease": "phenotype", "NCBI name": "ncbi_name", "MIrROR name": "microbiome", "Health sign": "beta", "subtract": "microbiome_subtract"}, inplace=True)
            self.df_beta = self.df_beta[["phenotype", "ncbi_name", "microbiome", "beta", "microbiome_subtract"]]
            self.df_beta['beta'] = self.df_beta['beta'].replace({'유해': 1, '유익': -1})

            self.df_dysbiosis.rename(columns = {"NCBI name": "ncbi_name", "MIrROR name": "microbiome", "Health sign": "beta", "subtract": "microbiome_subtract", "Exclude": "exclude"}, inplace=True)
            self.df_dysbiosis = self.df_dysbiosis[["ncbi_name", "microbiome", "beta", "microbiome_subtract", "exclude"]]
            self.df_dysbiosis['beta'] = self.df_dysbiosis['beta'].replace({'유해': 1, '유익': -1})

            self.df_probio.rename(columns = {"NCBI name": "ncbi_name", "MIrROR name": "microbiome", "subtract": "microbiome_subtract"}, inplace=True)
            self.df_probio = self.df_probio[["ncbi_name", "microbiome", "microbiome_subtract"]]
                        
            if self.path_exp.split(".")[-1] == 'txt':
                self.df_exp = pd.read_csv(self.path_exp, sep='\t', header=None)

                if len(self.df_exp.columns) == 2:     
                    sample_name = self.path_exp.split("_")[-1].split(".")[0]
                    self.df_exp.columns=["taxa", sample_name]
                
                elif len(self.df_exp.columns) == 4:
                    sample_name = self.df_exp.iloc[0,0]
                    observed = self.df_exp.iloc[1,2]
                    diversity = self.df_exp.iloc[0,3]

                    self.df_exp = self.df_exp.iloc[0:,[1,3]]              
                    self.df_exp.columns=["taxa", sample_name]
                    self.df_exp.loc[self.df_exp["taxa"] == 'observed', sample_name] = observed
                    self.df_exp.loc[self.df_exp["taxa"] == 'diversity', sample_name] = diversity
                                       
                else:
                    print("Check the proportion input file!")
                
            else:                    
                try:
                    self.df_exp = pd.read_csv(self.path_exp)
                
                except:
                    print("Check the proportion input file!")
            
            # Delete CST type row if it exists - df_exp
            if (self.df_exp['taxa']=='subCST').any():
                self.df_exp.drop(self.df_exp[(self.df_exp['taxa'] == 'subCST')].index, inplace=True)
                    
            print(self.df_exp)
            
            # Delete the diversity, observed rows
            if (list(self.df_exp['taxa'][0:2]) == ['diversity', 'observed']) & (list(self.df_db['taxa'][0:2]) == ['diversity', 'observed']):
                self.li_diversity = list(self.df_exp.iloc[0,1:]) # li_diversity : Alpha-Diversity list 
                self.li_observed = list(self.df_exp.iloc[1,1:]) # li_observed : Number of Microbiome list
                
                self.observed_mean = self.df_db[(self.df_db.taxa == "observed")].mean(axis=1, numeric_only=True).values[0]
                self.df_exp = self.df_exp.iloc[2:,:]
                self.df_db = self.df_db.iloc[2:,:]

            else: 
                print("Check the proportion input file!")
                
            # li_new_sample_name : Sample name list 
            # li_phenotype : Phenotype list 
            self.li_new_sample_name = list(self.df_exp.columns)[1:]  
            self.li_phenotype = list(dict.fromkeys(self.df_beta['phenotype']))
            
            print(self.df_beta)
                            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            
        return rv, rvmsg
    
    def TrimInputData(self): 
        """
        Trim the input data df_exp.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try:          
                   
            for idx in range(len(self.li_new_sample_name)):           
                min_abundance = self.df_exp[self.df_exp[self.li_new_sample_name[idx]] != 0].min()[self.li_new_sample_name[idx]]
                num_min_microbiome = len(self.df_exp.loc[self.df_exp[self.li_new_sample_name[idx]] == min_abundance])
                self.df_exp.loc[self.df_exp[self.li_new_sample_name[idx]] == min_abundance, self.li_new_sample_name[idx]] = 0
                
                self.df_exp[self.li_new_sample_name[idx]] /= (1-min_abundance*num_min_microbiome)
                
                print(self.df_exp)




        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
    
        return rv, rvmsg

    def CalculateMRS(self): 
        """
        Calculate the MRS (Microbiome Risk Score).

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try:                
            # df_mrs : Data frame of MRS corresponding to specific phenotype and sample
            self.df_mrs = pd.DataFrame(index = self.li_new_sample_name, columns = self.li_phenotype)
            self.df_mrs = self.df_mrs.fillna(0) 

            for i in range(len(self.li_new_sample_name)):
                for j in range(len(self.li_phenotype)):
                    condition_phen = (self.df_beta.phenotype == self.li_phenotype[j])   
                    mrs = 0

                    for idx_beta, row_beta in self.df_beta[condition_phen].iterrows():
                        condition_micro = (self.df_exp.taxa == row_beta['microbiome'])
                        abundance = 0

                        if (len(self.df_exp[condition_micro]) > 0):      
                            abundance += self.df_exp[condition_micro][self.li_new_sample_name[i]].values[0] 
                        
                        li_micro_sub = []

                        if pd.isna(row_beta['microbiome_subtract']) is False:
                            li_micro_sub = row_beta['microbiome_subtract'].split('\n')

                            for micro_sub in li_micro_sub:
                                condition_sub = (self.df_exp.taxa == micro_sub)

                                if len(self.df_exp[condition_sub]) > 0:
                                    abundance -= self.df_exp[condition_sub][self.li_new_sample_name[i]].values[0]     
                            
                        mrs += row_beta['beta'] * math.log10(100*abundance + 1) 

                    mrs /= len(self.df_beta[condition_phen])       
                    self.df_mrs.loc[self.li_new_sample_name[i], self.li_phenotype[j]] = mrs

        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
    
        return rv, rvmsg

    def CalculateDysbiosis(self): 
        """
        Calculate the Dysbiosis.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """         
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 
            self.li_microbiome = list(dict.fromkeys(self.df_dysbiosis['microbiome']))
            
            for i in range(len(self.li_new_sample_name)):
                dysbiosis_harmful = 0
                dysbiosis_beneficial = 0
                
                for j in range(len(self.li_microbiome)):
                    condition_harmful = (self.df_dysbiosis.microbiome == self.li_microbiome[j]) & (self.df_dysbiosis.beta == 1) & (self.df_dysbiosis.exclude != "Y") 
                    condition_beneficial = (self.df_dysbiosis.microbiome == self.li_microbiome[j]) & (self.df_dysbiosis.beta == -1) & (self.df_dysbiosis.exclude != "Y")
                    
                    if (len(self.df_dysbiosis[condition_harmful]) >= 1) & (len(self.df_dysbiosis[condition_beneficial]) == 0):
                        condition_micro = (self.df_exp.taxa == self.li_microbiome[j])
                        abundance = 0

                        if (len(self.df_exp[condition_micro]) > 0):      
                            abundance += self.df_exp[condition_micro][self.li_new_sample_name[i]].values[0]    
                            li_micro_sub = []
                            if pd.isna(self.df_dysbiosis[condition_harmful]['microbiome_subtract'].values[0]) is False:
                                li_micro_sub = self.df_dysbiosis[condition_harmful]['microbiome_subtract'].values[0].split('\n')

                                for micro_sub in li_micro_sub:
                                    condition_sub = (self.df_exp.taxa == micro_sub)

                                    if len(self.df_exp[condition_sub]) > 0:
                                        abundance -= self.df_exp[condition_sub][self.li_new_sample_name[i]].values[0]                                  
                                        
                            dysbiosis_harmful += math.log10(100*abundance + 1)            
                            
                    elif (len(self.df_dysbiosis[condition_harmful]) == 0) & (len(self.df_dysbiosis[condition_beneficial]) >= 1):
                        condition_micro = (self.df_exp.taxa == self.li_microbiome[j])
                        abundance = 0

                        if (len(self.df_exp[condition_micro]) > 0):      
                            abundance += self.df_exp[condition_micro][self.li_new_sample_name[i]].values[0]  
                            li_micro_sub = []
                            if pd.isna(self.df_dysbiosis[condition_beneficial]['microbiome_subtract'].values[0]) is False:
                                li_micro_sub = self.df_dysbiosis[condition_beneficial]['microbiome_subtract'].values[0].split('\n')                     
                                
                                for micro_sub in li_micro_sub:
                                    condition_sub = (self.df_exp.taxa == micro_sub)

                                    if len(self.df_exp[condition_sub]) > 0:
                                        abundance -= self.df_exp[condition_sub][self.li_new_sample_name[i]].values[0]       
                                        
                            dysbiosis_beneficial -= math.log10(100*abundance + 1)      
                            
                self.df_mrs.loc[self.li_new_sample_name[i], 'DysbiosisHarmful'] = dysbiosis_harmful
                self.df_mrs.loc[self.li_new_sample_name[i], 'DysbiosisBeneficial'] = -dysbiosis_beneficial
                self.df_mrs.loc[self.li_new_sample_name[i], 'Dysbiosis'] = -dysbiosis_harmful-dysbiosis_beneficial                         
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit() 
    
        return rv, rvmsg        
 

    def CalculateHealthyDistance(self): 
        """
        Calculate the Healthy Distance.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """           
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 
            self.df_mrs['HealthyDistance'] = 0     
            
            np_healthy_abundance = self.df_healthy['RA'].to_numpy()
            np_healthy_abundance = np.append(np_healthy_abundance, 1-np_healthy_abundance.sum())
            np_healthy_abundance = multiplicative_replacement(np_healthy_abundance)
            np_healthy_abundance = clr(np_healthy_abundance)
            
            # Subtract the abundance - df_exp_healthy
            for idx in range(len(self.li_new_sample_name)): 
                df_exp_one = self.df_exp[['taxa', self.li_new_sample_name[idx]]]
                df_exp_one = df_exp_one[df_exp_one[self.li_new_sample_name[idx]] != 0]
                df_exp_one = df_exp_one[df_exp_one[['taxa']].applymap(starts_with_s).any(axis=1)]
                
                np_abundance = np.array([], dtype=np.float64).reshape(0,1)
                np_abundance_others = np.ones((1,1), dtype=float)                
                
                for idx_healthy, row_healthy in self.df_healthy.iterrows(): 
                    micro = row_healthy['microbiome']
                    np_abundance_temp = np.zeros((1,1), dtype=float)

                    condition_append = (df_exp_one.taxa == micro)

                    if len(df_exp_one[condition_append]) > 0:
                        np_abundance_temp += df_exp_one[condition_append].to_numpy()[:,1:].astype(np.float64)
                        np_abundance_others -= df_exp_one[condition_append].to_numpy()[:,1:].astype(np.float64)


                    np_abundance = np.concatenate((np_abundance,np_abundance_temp),axis=0)

                np_abundance = np.concatenate((np_abundance,np_abundance_others),axis=0)
                np_abundance = np_abundance.transpose()

                # Apply multiplicative replacement and CLR transformations
                np_abundance = multiplicative_replacement(np_abundance)
                np_abundance = clr(np_abundance)   
            
                # Calculate healthy distance for each new sample
                healthy_dist = np.linalg.norm(np_abundance - np_healthy_abundance)  
                
                self.df_mrs.loc[self.li_new_sample_name[idx], 'HealthyDistance'] = healthy_dist
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
    
        return rv, rvmsg    
        
    def CalculatePercentileRank(self):
        """
        Calculate the Percentile Rank and Save the Percentile Rank data as an Csv file.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
         
        rv = True
        rvmsg = "Success"
        
        try:      
            self.df_mrs['Diversity'] = self.li_diversity
            
            # Append the Dysbiosis, HealthyDistance, Diversity, TotalRiskScore to phenotype list
            self.li_phenotype += ['DysbiosisHarmful', 'DysbiosisBeneficial', 'Diversity', 'HealthyDistance', 'Dysbiosis']

            # Create an empty data frame with the same index and columns as the df_mrs data frame
            self.df_percentile_rank = pd.DataFrame(index = self.li_new_sample_name, columns = self.li_phenotype)
            # Loop through all samples and phenotypes and calculate the percentile rank
            for i in range(len(self.li_new_sample_name)):
                for j in range(len(self.li_phenotype)):
                    if self.li_phenotype[j] != 'HealthyDistance':
                        self.df_percentile_rank.loc[self.li_new_sample_name[i], self.li_phenotype[j]] = (percentileofscore(list(self.df_mrs_db[self.li_phenotype[j]]), self.df_mrs.loc[self.li_new_sample_name[i], self.li_phenotype[j]], kind='mean')).round(1)
                    
                    else:                     
                        distance = self.df_mrs.loc[self.li_new_sample_name[i], self.li_phenotype[j]]                                         
                        max_distance = self.df_mrs_db[self.li_phenotype[j]].max() 
                                                      

                        self.df_percentile_rank.loc[self.li_new_sample_name[i], self.li_phenotype[j]] =  ((1-distance/max_distance) * 100).round(1) 
                        
            # Outliers
            # Replace percentile ranks that are less than or equal to 5 with 5, and those that are greater than or equal to 95 with 95
            for i in range(len(self.li_phenotype)):
                self.df_percentile_rank.loc[self.df_percentile_rank[self.li_phenotype[i]]<=5, self.li_phenotype[i]] = 5
                self.df_percentile_rank.loc[self.df_percentile_rank[self.li_phenotype[i]]>=95, self.li_phenotype[i]] = 95   

            # Define a dictionary to map species to their specific category and corresponding phenotypes
            self.species_specific_categories =  {
                    '뇌질환': ['치매', '우울증', '파킨슨병', '불면증', '불안&공황장애', 'ADHD', '자폐 스펙트럼'],
                    '심혈관질환': ['고혈압', '동맥경화', '심근경색'],
                    '장질환': ['변비', '설사', '과민성 장증후군', '염증성 장질환'],
                    '간질환': ['지방간', '간염', '간경변'],
                    '자가면역질환': ["아토피 피부염", "건선", "류마티스 관절염"],
                    '대사질환': ['비만', '제 2형 당뇨병']
                }
            
            # Main Category
            for category, phenotypes in self.species_specific_categories.items():
                self.df_percentile_rank[category] = self.df_percentile_rank[phenotypes].mean(axis=1)                
                
            self.df_percentile_rank['GMHS'] = ((self.df_percentile_rank['Diversity']*2) + self.df_percentile_rank['DysbiosisBeneficial'] + (1.5*(100-self.df_percentile_rank['DysbiosisHarmful'])) + self.df_percentile_rank['HealthyDistance']*0.5)/5

            for col in self.df_percentile_rank:
                self.df_percentile_rank[col] = self.df_percentile_rank[col].astype(float).round()
                     
            # Replace missing values with the string 'None'    
            self.df_percentile_rank = self.df_percentile_rank.fillna('None')
                        
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
    
        return rv, rvmsg
     
    def EvaluatePercentileRank(self):
        """
        Evaluate based on percentile rank value 

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
         
        rv = True
        rvmsg = "Success"
        
        try:                 
            self.df_eval = pd.DataFrame(index=self.df_percentile_rank.index)
            
            li_positive_var = ['GMHS', 'Diversity', 'DysbiosisBeneficial']
            
            for col in self.df_percentile_rank:
                if col in li_positive_var:
                    # Define the conditions and corresponding values
                    conditions = [
                        self.df_percentile_rank[col] >= 95,
                        (self.df_percentile_rank[col] > 80) & (self.df_percentile_rank[col] < 95),
                        (self.df_percentile_rank[col] > 50) & (self.df_percentile_rank[col] <= 80),
                        (self.df_percentile_rank[col] > 20) & (self.df_percentile_rank[col] <= 50),
                        self.df_percentile_rank[col] <= 20
                    ]
                    
                    values = ['좋음', '보통', '주의', '나쁨', '아주 나쁨']     
                    
                    self.df_eval[col] = np.select(conditions, values)  
   
                elif col == 'HealthyDistance':
                    # Define the conditions and corresponding values
                    conditions = [
                        self.df_percentile_rank[col] >= 50,
                        (self.df_percentile_rank[col] > 40) & (self.df_percentile_rank[col] < 50),
                        (self.df_percentile_rank[col] > 30) & (self.df_percentile_rank[col] <= 40),
                        (self.df_percentile_rank[col] > 20) & (self.df_percentile_rank[col] <= 30),
                        self.df_percentile_rank[col] <= 20
                    ]
                    
                    values = ['좋음', '보통', '주의', '나쁨', '아주 나쁨']     
                    
                    self.df_eval[col] = np.select(conditions, values)             
                                                      
                else: # MRS, DysbiosisHarmful
                    # Define the conditions and corresponding values
                    conditions = [
                        self.df_percentile_rank[col] >= 80,
                        (self.df_percentile_rank[col] >= 50) & (self.df_percentile_rank[col] < 80),
                        (self.df_percentile_rank[col] >= 20) & (self.df_percentile_rank[col] < 50),
                        (self.df_percentile_rank[col] > 5) & (self.df_percentile_rank[col] < 20),
                        self.df_percentile_rank[col] <= 5
                    ]     
                    
                    values = ['아주 나쁨', '나쁨', '주의', '보통', '좋음']   
                    
                    self.df_eval[col] = np.select(conditions, values)    
                    
            self.df_eval = self.df_eval.loc[:,'DysbiosisHarmful':]
            
            # Type E, B, I, D
            conditions = [
                (self.df_percentile_rank['GMHS'] >= 0) & (self.df_percentile_rank['GMHS'] < 45),

                (self.df_percentile_rank['GMHS'] >= 45) & (self.df_percentile_rank['GMHS'] < 60),

                (self.df_percentile_rank['GMHS'] >= 60) & (self.df_percentile_rank['GMHS'] < 75),

                (self.df_percentile_rank['GMHS'] >= 75)
            ]
            values = ['D', 'B', 'I', 'E']
            
            self.df_eval['Type'] = np.select(conditions, values)
            
            
            ##############
            '''
            conditions = [
                (self.df_percentile_rank['Diversity'] >= 60) & (self.df_percentile_rank['Dysbiosis'] >= 60),
                
                (self.df_percentile_rank['Diversity'] < 60) & (self.df_percentile_rank['Dysbiosis'] >= 60),
                
                (self.df_percentile_rank['Diversity'] >= 60) & (self.df_percentile_rank['Dysbiosis'] < 60),
                
                (self.df_percentile_rank['Diversity'] < 60) & (self.df_percentile_rank['Dysbiosis'] < 60)
            ]
            values = ['E', 'B', 'I', 'D']

            self.df_eval['Type'] = np.select(conditions, values)
            '''
            
            '''
            # Print the EBID percentages of the samples
            E_data = self.df_eval[(self.df_eval['Type'] == 'E')]
            B_data = self.df_eval[(self.df_eval['Type'] == 'B')]
            D_data = self.df_eval[(self.df_eval['Type'] == 'D')]
            I_data = self.df_eval[(self.df_eval['Type'] == 'I')]

            E_percent = len(E_data) / len(self.df_eval) * 100
            B_percent = len(B_data) / len(self.df_eval) * 100
            D_percent = len(D_data) / len(self.df_eval) * 100
            I_percent = len(I_data) / len(self.df_eval) * 100
            
            print("Percentage of samples in E: ", E_percent, '%')
            print("Percentage of samples in B: ", B_percent, '%') 
            print("Percentage of samples in D: ", D_percent, '%')
            print("Percentage of samples in I: ", I_percent, '%')           
            '''
            ##############
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
    
        return rv, rvmsg    

      
    def CalculateMicrobiomeRatio(self): 
        """
        Calculate the Beneficial Microbiome Ratio & Harmful Microbiome Ratio

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """         
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 
            self.li_microbiome = list(dict.fromkeys(self.df_dysbiosis['microbiome']))
            
            for i in range(len(self.li_new_sample_name)):
                harmful_abundance = 0
                beneficial_abundance = 0
                
                harmful_number = 0
                beneficial_number = 0
                
                for j in range(len(self.li_microbiome)):
                    condition_harmful = (self.df_dysbiosis.microbiome == self.li_microbiome[j]) & (self.df_dysbiosis.beta == 1) & (self.df_dysbiosis.exclude != "Y")
                    condition_beneficial = (self.df_dysbiosis.microbiome == self.li_microbiome[j]) & (self.df_dysbiosis.beta == -1) & (self.df_dysbiosis.exclude != "Y")
                    
                    if (len(self.df_dysbiosis[condition_harmful]) >= 1) & (len(self.df_dysbiosis[condition_beneficial]) == 0):
                        condition_micro = (self.df_exp.taxa == self.li_microbiome[j])
                        abundance = 0

                        if (len(self.df_exp[condition_micro]) > 0):      
                            abundance += self.df_exp[condition_micro][self.li_new_sample_name[i]].values[0]    
                            li_micro_sub = []
                            if pd.isna(self.df_dysbiosis[condition_harmful]['microbiome_subtract'].values[0]) is False:
                                li_micro_sub = self.df_dysbiosis[condition_harmful]['microbiome_subtract'].values[0].split('\n')

                                for micro_sub in li_micro_sub:
                                    condition_sub = (self.df_exp.taxa == micro_sub)

                                    if len(self.df_exp[condition_sub]) > 0:
                                        abundance -= self.df_exp[condition_sub][self.li_new_sample_name[i]].values[0]                                  
                                        
                        harmful_abundance += abundance  

                        if abundance > 0:
                            harmful_number += 1
                            
                    elif (len(self.df_dysbiosis[condition_harmful]) == 0) & (len(self.df_dysbiosis[condition_beneficial]) >= 1):
                        condition_micro = (self.df_exp.taxa == self.li_microbiome[j])
                        abundance = 0

                        if (len(self.df_exp[condition_micro]) > 0):      
                            abundance += self.df_exp[condition_micro][self.li_new_sample_name[i]].values[0]  
                            li_micro_sub = []
                            if pd.isna(self.df_dysbiosis[condition_beneficial]['microbiome_subtract'].values[0]) is False:
                                li_micro_sub = self.df_dysbiosis[condition_beneficial]['microbiome_subtract'].values[0].split('\n')                     
                                
                                for micro_sub in li_micro_sub:
                                    condition_sub = (self.df_exp.taxa == micro_sub)

                                    if len(self.df_exp[condition_sub]) > 0:
                                        abundance -= self.df_exp[condition_sub][self.li_new_sample_name[i]].values[0]       
                                        
                        beneficial_abundance += abundance   

                        if abundance > 0:
                            beneficial_number += 1
                            
                self.df_eval.loc[self.li_new_sample_name[i], 'harmful_abundance[%]'] = harmful_abundance * 100
                self.df_eval.loc[self.li_new_sample_name[i], 'beneficial_abundance[%]'] = beneficial_abundance * 100
                self.df_eval.loc[self.li_new_sample_name[i], 'num_total_species'] = int(self.li_observed[i])

                              
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
            
        return rv, rvmsg  

    def CalculateAverageMicrobiomeRatio(self): 
        """
        Calculate the average Beneficial Microbiome Ratio & Harmful Microbiome Ratio

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """         
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 

            harmful_mean_abundance = 0
            beneficial_mean_abundance = 0

            for j in range(len(self.li_microbiome)):
                condition_harmful = (self.df_dysbiosis.microbiome == self.li_microbiome[j]) & (self.df_dysbiosis.beta == 1) & (self.df_dysbiosis.exclude != "Y")
                condition_beneficial = (self.df_dysbiosis.microbiome == self.li_microbiome[j]) & (self.df_dysbiosis.beta == -1) & (self.df_dysbiosis.exclude != "Y")

                if (len(self.df_dysbiosis[condition_harmful]) >= 1) & (len(self.df_dysbiosis[condition_beneficial]) == 0):
                    condition_micro = (self.df_db.taxa == self.li_microbiome[j])
                    abundance_mean = 0

                    if (len(self.df_db[condition_micro]) > 0):      
                        abundance_mean += self.df_db[condition_micro].mean(axis=1, numeric_only=True).values[0]      
                        li_micro_sub = []
                        if pd.isna(self.df_dysbiosis[condition_harmful]['microbiome_subtract'].values[0]) is False:
                            li_micro_sub = self.df_dysbiosis[condition_harmful]['microbiome_subtract'].values[0].split('\n')

                            for micro_sub in li_micro_sub:
                                condition_sub = (self.df_db.taxa == micro_sub)

                                if len(self.df_db[condition_sub]) > 0:
                                    abundance_mean -= self.df_db[condition_sub].mean(axis=1, numeric_only=True).values[0]       

                        harmful_mean_abundance += abundance_mean       

                elif (len(self.df_dysbiosis[condition_harmful]) == 0) & (len(self.df_dysbiosis[condition_beneficial]) >= 1):
                    condition_micro = (self.df_db.taxa == self.li_microbiome[j])
                    abundance_mean = 0

                    if (len(self.df_db[condition_micro]) > 0):      
                        abundance_mean += self.df_db[condition_micro].mean(axis=1, numeric_only=True).values[0]     
                        li_micro_sub = []
                        if pd.isna(self.df_dysbiosis[condition_beneficial]['microbiome_subtract'].values[0]) is False:
                            li_micro_sub = self.df_dysbiosis[condition_beneficial]['microbiome_subtract'].values[0].split('\n')                     

                            for micro_sub in li_micro_sub:
                                condition_sub = (self.df_db.taxa == micro_sub)

                                if len(self.df_db[condition_sub]) > 0:
                                    abundance_mean -= self.df_db[condition_sub].mean(axis=1, numeric_only=True).values[0]     

                        beneficial_mean_abundance += abundance_mean   

            self.df_eval.loc[:,'harmful_mean_abundance[%]'] = harmful_mean_abundance * 100
            self.df_eval.loc[:,'beneficial_mean_abundance[%]'] = beneficial_mean_abundance * 100
                                
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
            
        return rv, rvmsg      
    
    def CalculateHarmfulMicrobiomeAbundance(self):
        """
        Calculate specific harmful microbiome abundance and average harmful microbiome abundance.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """  
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try:     
            json_abundance = []
            self.li_ncbi_name = list(dict.fromkeys(self.df_dysbiosis['ncbi_name']))
            
            for i in range(len(self.li_new_sample_name)):
                for j in range(len(self.li_ncbi_name)):

                    condition_ncbi = (self.df_dysbiosis.beta == 1) & (self.df_dysbiosis.ncbi_name == self.li_ncbi_name[j]) 

                    abundance = 0 
                    abundance_mean = 0
                    for idx_dysbiosis, row_dysbiosis in self.df_dysbiosis[condition_ncbi].iterrows(): 
                        condition_exp = (self.df_exp.taxa == row_dysbiosis['microbiome'])
                        condition_db = (self.df_db.taxa == row_dysbiosis['microbiome'])
                        
                        if len(self.df_exp[condition_exp]) > 0:
                            abundance += self.df_exp[condition_exp][self.li_new_sample_name[i]].values[0]
                            
                        li_micro_sub = []

                        if pd.isna(row_dysbiosis['microbiome_subtract']) is False:
                            li_micro_sub = row_dysbiosis['microbiome_subtract'].split('\n')

                            for micro_sub in li_micro_sub:
                                condition_sub = (self.df_exp.taxa == micro_sub)

                                if len(self.df_exp[condition_sub]) > 0:
                                    abundance -= self.df_exp[condition_sub][self.li_new_sample_name[i]].values[0]   
                            
                        
                        if len(self.df_db[condition_db]) > 0:
                            abundance_mean += self.df_db[condition_db].mean(axis=1, numeric_only=True).values[0]
                            
                        li_micro_sub = []

                        if pd.isna(row_dysbiosis['microbiome_subtract']) is False:
                            li_micro_sub = row_dysbiosis['microbiome_subtract'].split('\n')

                            for micro_sub in li_micro_sub:
                                condition_sub = (self.df_db.taxa == micro_sub)

                                if len(self.df_db[condition_sub]) > 0:
                                    abundance_mean -= self.df_db[condition_sub].mean(axis=1, numeric_only=True).values[0]                           
                        json_abundance.append({"sample_name" : self.li_new_sample_name[i], "ncbi_name" : self.li_ncbi_name[j], "abundance" : abundance, "abundance_mean" : abundance_mean})

            df_abundance = pd.DataFrame.from_dict(json_abundance)   

            df_abundance = df_abundance.drop_duplicates(['sample_name', 'ncbi_name'], keep='last')
               
            self.df_harmful_10 = pd.DataFrame(columns = ["sample_name", "ncbi_name", "abundance", "abundance_mean"])
            self.df_harmful_30 = pd.DataFrame(columns = ["sample_name", "ncbi_name", "abundance", "abundance_mean"])

            for i in range(len(self.li_new_sample_name)):
                condition = (df_abundance.sample_name == self.li_new_sample_name[i])
                df_new = df_abundance[condition].sort_values(by=['abundance_mean'], ascending=False).head(10)
                self.df_harmful_10 = pd.concat([self.df_harmful_10,df_new])
                
                df_tot = df_abundance[condition].sort_values(by=['abundance_mean'], ascending=False)
                self.df_harmful_30 = pd.concat([self.df_harmful_30,df_tot])
                
            self.df_harmful_10 = self.df_harmful_10.set_index(keys=['sample_name'], inplace=False, drop=True)           
            self.df_harmful_10.to_csv(self.path_harmful)   

            self.df_harmful_30 = self.df_harmful_30.set_index(keys=['sample_name'], inplace=False, drop=True)           
            self.df_harmful_30.to_csv(self.path_harmful_tot)   
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")
            sys.exit()
    
        return rv, rvmsg  
    
    
    def CalculateBeneficialMicrobiomeAbundance(self):
        """
        Calculate specific beneficial microbiome abundance and average beneficial microbiome abundance.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """  
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try:     
            json_abundance = []
            
            for i in range(len(self.li_new_sample_name)):
                for j in range(len(self.li_ncbi_name)):

                    condition_ncbi = (self.df_dysbiosis.beta == -1) & (self.df_dysbiosis.ncbi_name == self.li_ncbi_name[j]) 

                    abundance = 0 
                    abundance_mean = 0
                    for idx_dysbiosis, row_dysbiosis in self.df_dysbiosis[condition_ncbi].iterrows(): 
                        condition_exp = (self.df_exp.taxa == row_dysbiosis['microbiome'])
                        condition_db = (self.df_db.taxa == row_dysbiosis['microbiome'])
                        
                        if len(self.df_exp[condition_exp]) > 0:
                            abundance += self.df_exp[condition_exp][self.li_new_sample_name[i]].values[0]
                            
                        li_micro_sub = []

                        if pd.isna(row_dysbiosis['microbiome_subtract']) is False:
                            li_micro_sub = row_dysbiosis['microbiome_subtract'].split('\n')

                            for micro_sub in li_micro_sub:
                                condition_sub = (self.df_exp.taxa == micro_sub)

                                if len(self.df_exp[condition_sub]) > 0:
                                    abundance -= self.df_exp[condition_sub][self.li_new_sample_name[i]].values[0]   
                            
                        
                        if len(self.df_db[condition_db]) > 0:
                            abundance_mean += self.df_db[condition_db].mean(axis=1, numeric_only=True).values[0]
                            
                        li_micro_sub = []

                        if pd.isna(row_dysbiosis['microbiome_subtract']) is False:
                            li_micro_sub = row_dysbiosis['microbiome_subtract'].split('\n')

                            for micro_sub in li_micro_sub:
                                condition_sub = (self.df_db.taxa == micro_sub)

                                if len(self.df_db[condition_sub]) > 0:
                                    abundance_mean -= self.df_db[condition_sub].mean(axis=1, numeric_only=True).values[0]                           
                        json_abundance.append({"sample_name" : self.li_new_sample_name[i], "ncbi_name" : self.li_ncbi_name[j], "abundance" : abundance, "abundance_mean" : abundance_mean})

            df_abundance = pd.DataFrame.from_dict(json_abundance)   

            df_abundance = df_abundance.drop_duplicates(['sample_name', 'ncbi_name'], keep='last')

            self.df_beneficial_10 = pd.DataFrame(columns = ["sample_name", "ncbi_name", "abundance", "abundance_mean"])

            for i in range(len(self.li_new_sample_name)):
                condition = (df_abundance.sample_name == self.li_new_sample_name[i])
                df_new = df_abundance[condition].sort_values(by=['abundance_mean'], ascending=False).head(10)
                self.df_beneficial_10 = pd.concat([self.df_beneficial_10,df_new])

            self.df_beneficial_10 = self.df_beneficial_10.set_index(keys=['sample_name'], inplace=False, drop=True)           
            self.df_beneficial_10.to_csv(self.path_beneficial)    
    
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")
            sys.exit()
    
        return rv, rvmsg       
    
    def CalculateTotalProbioRatio(self): 
        """
        Calculate total probiotic abundance and average probiotic abundance.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """         
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try: 
            self.li_probio_microbiome = list(dict.fromkeys(self.df_probio['microbiome']))
            
            for i in range(len(self.li_new_sample_name)):
                
                probio_abundance = 0
                probio_abundance_mean = 0
                
                for j in range(len(self.li_probio_microbiome)):                                    
                    condition_micro = (self.df_exp.taxa == self.li_probio_microbiome[j])
                    condition_micro_db = (self.df_db.taxa == self.li_probio_microbiome[j])
                    
                    abundance = 0
                    abundance_mean = 0

                    if (len(self.df_exp[condition_micro]) > 0):      
                        abundance += self.df_exp[condition_micro][self.li_new_sample_name[i]].values[0]  
                        
                    if (len(self.df_db[condition_micro_db]) > 0):                        
                        abundance_mean += self.df_db[condition_micro_db].mean(axis=1, numeric_only=True).values[0]    


                    probio_abundance += abundance  
                    probio_abundance_mean += abundance_mean
                                                        
                self.df_eval.loc[self.li_new_sample_name[i], 'probio_abundance[%]'] = probio_abundance * 100
                
            self.df_eval.loc[:,'probio_abundance_mean[%]'] = probio_abundance_mean * 100
            self.df_eval.loc[:,'observed_mean'] = int(self.observed_mean)
                                   
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
            
        return rv, rvmsg  
    
    def CalculateSpecificProbioRatio(self):
        """
        Calculate specific probiotic abundance.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """  
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try:     
            json_probio_abundance = []
            self.li_ncbi_name = list(dict.fromkeys(self.df_probio['ncbi_name']))
            
            for i in range(len(self.li_new_sample_name)):
                for j in range(len(self.li_ncbi_name)):
                    
                    condition_ncbi = (self.df_probio.ncbi_name == self.li_ncbi_name[j]) 
                    
                    abundance = 0
                    for idx_probio, row_probio in self.df_probio[condition_ncbi].iterrows(): 
                    
                        condition_probio = (self.df_exp.taxa == row_probio["microbiome"])
                        if (len(self.df_exp[condition_probio]) > 0):                        
                            abundance += self.df_exp[condition_probio][self.li_new_sample_name[i]].values[0]
                    
                        json_probio_abundance.append({"sample_name" : self.li_new_sample_name[i], "ncbi_name" : row_probio["ncbi_name"], "abundance" : abundance})
            df_probio_abundance = pd.DataFrame.from_dict(json_probio_abundance)   
            df_probio_abundance = df_probio_abundance.drop_duplicates(['sample_name', 'ncbi_name'], keep='last')
                           
            self.df_probio_19 = pd.DataFrame(columns = ["sample_name", "ncbi_name", "abundance"])

            for i in range(len(self.li_new_sample_name)):
                condition = (df_probio_abundance.sample_name == self.li_new_sample_name[i])           
                df_tot = df_probio_abundance[condition].sort_values(by=['abundance'], ascending=False)
                self.df_probio_19 = pd.concat([self.df_probio_19,df_tot])
                
            self.df_probio_19 = self.df_probio_19.set_index(keys=['sample_name'], inplace=False, drop=True)           
            self.df_probio_19.to_csv(self.path_probio_tot)   

            for i in range(len(self.li_new_sample_name)):
                self.df_eval.loc[self.li_new_sample_name[i],'num_detected_harmful_microbiome'] = int(len(self.df_harmful_30.loc[(self.df_harmful_30['abundance'] > self.threshold) & (self.df_harmful_30.index == self.li_new_sample_name[i])])) 
                self.df_eval.loc[self.li_new_sample_name[i],'num_detected_probio'] = int(len(self.df_probio_19.loc[(self.df_probio_19['abundance'] > self.threshold) & (self.df_probio_19.index == self.li_new_sample_name[i])])) 
                
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")
            sys.exit()
    
        return rv, rvmsg     

    def CalculateCoordinateLocation(self):
        """
        Calculate the Coordinate Location.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """  
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try:    
            for i in range(len(self.li_new_sample_name)):
                GMHS_type = self.df_eval.at[self.li_new_sample_name[i], 'Type']
                xpos = self.df_percentile_rank.at[self.li_new_sample_name[i], 'GMHS']
                ypos = xpos

                xpos_cal, ypos_cal = nogada_x_y_graph_v2(GMHS_type, xpos, ypos)
                
                self.df_eval.loc[self.li_new_sample_name[i], 'Xpos'] = xpos_cal            
                self.df_eval.loc[self.li_new_sample_name[i], 'Ypos'] = ypos_cal                    
                       
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")
            sys.exit()
    
        return rv, rvmsg         

    def DetectProteusMirabilis(self):
        """
        Detect the Proteus Mirabilis

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """  
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try:    
            self.df_eval['Proteus_mirabilis'] = '불검출'
            if 's__Proteus_mirabilis' in self.df_exp.taxa.to_list():
                
                for i in range(len(self.li_new_sample_name)):                           
                    condition_proteus = (self.df_exp.taxa == 's__Proteus_mirabilis')               
                    proteus_abundance = self.df_exp[condition_proteus][self.li_new_sample_name[i]].values[0]  
                    if proteus_abundance > 0:
                        self.df_eval.loc[self.li_new_sample_name[i], 'Proteus_mirabilis'] = '검출'

            # Delete the Dysbiosis column
            self.df_percentile_rank = self.df_percentile_rank.drop('Dysbiosis', axis=1)
            # Save the output file - Percentile Rank of the samples
            self.df_percentile_rank.to_csv(self.path_percentile_rank_output, encoding="utf-8-sig", index_label='serial_number')
            
            # Delete the Dysbiosis column
            self.df_eval = self.df_eval.drop('Dysbiosis', axis=1)
            # Save the output file - df_eval
            self.df_eval.to_csv(self.path_eval_output, encoding="utf-8-sig", index_label='serial_number')   
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")
            sys.exit()
    
        return rv, rvmsg  

    def SaveProumedOutput(self):
        """
        Save the output file for PROUMED

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """  
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try:    
            self.df_proumed = self.df_percentile_rank[['GMHS','고혈압', '제 2형 당뇨병', '심근경색', '동맥경화', '간경변', '비만']].copy()
            
            self.df_proumed.insert(0,'Type', self.df_eval['Type'])
            
            
            self.df_proumed.rename(columns = {"serial_number": "Sample ID", 
                                              "Type": "Intestinal environment type", 
                                              "GMHS": "Intestinal environment score",  
                                              "고혈압": "High blood pressure score",
                                              "제 2형 당뇨병": "Diabetes score",
                                              "심근경색": "Myocardial infarction score",
                                              "동맥경화": "Arteriosclerosis score", 
                                              "간경변": "Cirrhosis score", 
                                              "비만": "Obesity score"}, inplace=True)
            
            
            # Save the output file - df_proumed
            self.df_proumed.to_csv(self.path_eval_output_proumed, encoding="utf-8-sig", index_label='Sample ID')   
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")
            sys.exit()
    
        return rv, rvmsg 
    
#####################################
# main
#####################################
if __name__ == '__main__':
    
    #path_exp = "input/EGgutPro_mirror_output_3175.csv"    
    #path_exp = "input/EGgutPro_one_sample.csv"
    #path_exp = "input/EGgutPro_sample_input.txt"
    path_exp = "input/BC72_EG23-HU08-NGS-S4H5.txt"

    
    eggutanalysis = EgGutProAnalysis(path_exp)
    eggutanalysis.ReadDB()
    #eggutanalysis.TrimInputData()      
    eggutanalysis.CalculateMRS()    
    eggutanalysis.CalculateDysbiosis()    
    eggutanalysis.CalculateHealthyDistance()
    eggutanalysis.CalculatePercentileRank()
    eggutanalysis.EvaluatePercentileRank()           
    eggutanalysis.CalculateMicrobiomeRatio()
    eggutanalysis.CalculateAverageMicrobiomeRatio()
    eggutanalysis.CalculateHarmfulMicrobiomeAbundance()
    eggutanalysis.CalculateBeneficialMicrobiomeAbundance()
    eggutanalysis.CalculateTotalProbioRatio()
    eggutanalysis.CalculateSpecificProbioRatio()   
    eggutanalysis.CalculateCoordinateLocation() 
    eggutanalysis.DetectProteusMirabilis() 
    #eggutanalysis.SaveProumedOutput()
        
    print('Analysis Complete')
    
    
    
    
    