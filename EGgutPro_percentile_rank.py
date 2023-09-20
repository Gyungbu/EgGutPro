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

def nogada_x_y_graph(type, xpos, ypos):
    #print("GRAPH x: ", xpos, "y: ", ypos)
    if (type == 'E'):
        if (xpos < 82):  # E타입 x축 임계점이 81.몇임. 이거 이하로 넘어가면 I타입 되버려서.... 예외처리
            xpos = 83

    if (type == 'I'):
        if (xpos < 54):
            xpos = 54  # I타입 X값 최소
        if (ypos < 54):
            ypos = 54

        if(xpos>75):
            xpos=75
        if(ypos>68):
            ypos=68

    if (type == 'B'):
        if (xpos < 23):
            xpos = 27
        if (ypos < 27):
            ypos = 33
        if (xpos>48):
            xpos=48
        if(ypos>50):
            ypos=50

    return (xpos, ypos)


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
        self.path_healthy = f"{curdir}/input/EGgutPro_healthy_person_profile_v2.xlsx"
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
        
        ## Lists used for calculation
        self.li_diversity = None
        self.li_observed = None
        self.li_new_sample_name = None
        self.li_phenotype = None
        self.li_microbiome = None
        self.li_ncbi_name = None
        self.li_probio_microbiome = None  
        
        self.observed_mean = None

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
            
            self.df_healthy = pd.read_excel(self.path_healthy, sheet_name="RA")
            self.df_exp = pd.read_csv(self.path_exp)
            self.df_mrs_db = pd.read_excel(self.path_mrs_db, index_col=0) 
            self.df_exp = pd.read_csv(self.path_exp)
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
            
            self.df_healthy = self.df_healthy[self.df_healthy['Taxonomy'].str.contains('s__')]
            self.df_healthy['Taxonomy'] = self.df_healthy['Taxonomy'].apply(filter_species)
            self.df_healthy['Taxonomy'] = self.df_healthy['Taxonomy'].str.replace(' ','_')
            self.df_healthy = self.df_healthy.rename(columns={'Taxonomy': 'taxa'})
            
            print(self.df_exp)
            # Delete the diversity, observed rows
            if (list(self.df_exp['taxa'][0:2]) == ['diversity', 'observed']) & (list(self.df_db['taxa'][0:2]) == ['diversity', 'observed']):
                self.li_diversity = list(self.df_exp.iloc[0,1:]) # li_diversity : Alpha-Diversity list 
                self.li_observed = list(self.df_exp.iloc[1,1:]) # li_observed : Number of Microbiome list
                
                self.observed_mean = self.df_db[(self.df_db.taxa == "observed")].mean(axis=1, numeric_only=True).values[0]
                self.df_exp = self.df_exp.iloc[2:,:]
                self.df_db = self.df_db.iloc[2:,:]
                            
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
                       
            # Subtract the abundance - df_exp_healthy
            for idx in range(len(self.li_new_sample_name)): 
                healthy_dist = 0
                
                df_exp_one = self.df_exp[['taxa', self.li_new_sample_name[idx]]]
                df_exp_one = df_exp_one[df_exp_one[self.li_new_sample_name[idx]] != 0]
                           
                for donor in self.df_healthy.iloc[:,1:].columns:
                    df_healthy_one = self.df_healthy[['taxa', donor]]
                    df_healthy_one = df_healthy_one[df_healthy_one[donor] != 0]
                                      
                    df_merged = pd.merge(df_exp_one, df_healthy_one, how='outer',on='taxa')
                    df_merged = df_merged.fillna(0)
                    
                    np_abundance = df_merged.iloc[:, 1:3]
                    np_abundance = np.transpose(np_abundance)
                    np_abundance = multiplicative_replacement(np_abundance)                     
                    np_abundance = clr(np_abundance)                  
                    np_abundance = np.transpose(np_abundance)             
                    
                    healthy_dist += np.linalg.norm(np_abundance[:, 0] - np_abundance[:, 1])  
                
                # Calculate healthy distance for each new sample
                self.df_mrs.loc[self.li_new_sample_name[idx], 'HealthyDistance'] = -healthy_dist / 8                  
                
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
            self.li_phenotype += ['DysbiosisHarmful', 'DysbiosisBeneficial', 'Diversity', 'HealthyDistance']

            # Create an empty data frame with the same index and columns as the df_mrs data frame
            self.df_percentile_rank = pd.DataFrame(index = self.li_new_sample_name, columns = self.li_phenotype)
            # Loop through all samples and phenotypes and calculate the percentile rank
            for i in range(len(self.li_new_sample_name)):
                for j in range(len(self.li_phenotype)):
                    self.df_percentile_rank.loc[self.li_new_sample_name[i], self.li_phenotype[j]] = (percentileofscore(list(self.df_mrs_db[self.li_phenotype[j]]), self.df_mrs.loc[self.li_new_sample_name[i], self.li_phenotype[j]], kind='mean')).round(1)
                 
            # Outliers
            # Replace percentile ranks that are less than or equal to 5 with 5, and those that are greater than or equal to 95 with 95
            for i in range(len(self.li_phenotype)):
                self.df_percentile_rank.loc[self.df_percentile_rank[self.li_phenotype[i]]<=5, self.li_phenotype[i]] = 5.0
                self.df_percentile_rank.loc[self.df_percentile_rank[self.li_phenotype[i]]>=95, self.li_phenotype[i]] = 95.0      

            # Define a dictionary to map species to their specific category and corresponding phenotypes
            species_specific_categories =  {
                    '뇌질환': ['알츠하이머', '이명', '불안장애', '불면증', '인지기능장애', '자폐증', '파킨슨병', '우울증'],
                    '심혈관질환': ['고혈압', '심근경색', '동맥경화'],
                    '간질환': ['지방간', '간경변', '간염'],
                    '소화기질환': ['변비', '설사', '염증성장', '과민성장'],
                    '대사질환': ["비만", "당뇨", "혈당조절"],
                    '자가면역질환': ['아토피', '건선', '류머티스']
                }
            
            # Main Category
            for category, phenotypes in species_specific_categories.items():
                self.df_percentile_rank[category] = self.df_percentile_rank[phenotypes].mean(axis=1)                
                
            self.df_percentile_rank['GMHS'] = ((self.df_percentile_rank['Diversity']*2) + self.df_percentile_rank['DysbiosisBeneficial'] + (1.5*(100-self.df_percentile_rank['DysbiosisHarmful'])) + self.df_percentile_rank['HealthyDistance'])/5.5
            
            for col in self.df_percentile_rank:
                self.df_percentile_rank[col] = self.df_percentile_rank[col].astype(float).round(1)
                     
            # Replace missing values with the string 'None'    
            self.df_percentile_rank = self.df_percentile_rank.fillna('None')

            # Save the output file - Percentile Rank of the samples
            self.df_percentile_rank.to_csv(self.path_percentile_rank_output, encoding="utf-8-sig", index_label='serial_number')
            
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
            
            li_positive_var = ['GMHS', 'Diversity', 'DysbiosisBeneficial', 'HealthyDistance']
            
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
                    
                    values = ['G', 'N', 'W', 'B', 'VB']     
                    
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
                    
                    values = ['VB', 'B', 'W', 'N', 'G']   
                    
                    self.df_eval[col] = np.select(conditions, values)    
                    
            self.df_eval = self.df_eval.loc[:,'DysbiosisHarmful':]
            
            # Type E, B, I, D
            conditions = [
                (self.df_percentile_rank['GMHS'] >= 0) & (self.df_percentile_rank['GMHS'] < 20),
                
                (self.df_percentile_rank['GMHS'] >= 20) & (self.df_percentile_rank['GMHS'] < 55),
                
                (self.df_percentile_rank['GMHS'] >= 55) & (self.df_percentile_rank['GMHS'] < 90),
                
                (self.df_percentile_rank['GMHS'] >= 90)
            ]
            values = ['D', 'B', 'I', 'E']

            self.df_eval['Type'] = np.select(conditions, values)
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
    
        return rv, rvmsg    

    def EvaluateBrainDiseaseException(self):
        """
        Evaluate Brain Disease for Exceptional Cases

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
         
        rv = True
        rvmsg = "Success"
        
        try:    
            li_phenotype_brain = ['알츠하이머', '이명', '불안장애', '불면증', '인지기능장애', '자폐증', '파킨슨병', '우울증']
            
            for i in range(len(self.li_new_sample_name)):
                for phenotype_brain in li_phenotype_brain:
                    brain_score = self.df_percentile_rank.at[self.li_new_sample_name[i], phenotype_brain]
                    if brain_score >= 80:
                        self.df_eval.loc[self.li_new_sample_name[i], '뇌질환'] = 'VB'            
            
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
                self.df_eval.loc[self.li_new_sample_name[i], 'num_total_species'] = self.li_observed[i]

                              
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
            self.df_eval.loc[:,'observed_mean'] = round(self.observed_mean)
                                   
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
                self.df_eval.loc[self.li_new_sample_name[i],'num_detected_beneficial_microbiome'] = len(self.df_probio_19.loc[(self.df_probio_19['abundance'] > 0) & (self.df_probio_19.index == self.li_new_sample_name[i])]) 
                self.df_eval.loc[self.li_new_sample_name[i],'num_detected_harmful_microbiome'] = len(self.df_harmful_30.loc[(self.df_harmful_30['abundance'] > 0) & (self.df_harmful_30.index == self.li_new_sample_name[i])]) 
                self.df_eval.loc[self.li_new_sample_name[i],'num_detected_probio'] = len(self.df_probio_19.loc[(self.df_probio_19['abundance'] > 0) & (self.df_probio_19.index == self.li_new_sample_name[i])]) 
                
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

                xpos_cal, ypos_cal = nogada_x_y_graph(GMHS_type, xpos, ypos)
                
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
               
            
            # Save the output file - df_eval
            self.df_eval.to_csv(self.path_eval_output, encoding="utf-8-sig", index_label='serial_number')   
            
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
    
    path_exp = "input/EGgutPro_one_sample.csv"
    
    eggutanalysis = EgGutProAnalysis(path_exp)
    eggutanalysis.ReadDB()
    eggutanalysis.CalculateMRS()    
    eggutanalysis.CalculateDysbiosis()    
    eggutanalysis.CalculateHealthyDistance()
    eggutanalysis.CalculatePercentileRank()
    eggutanalysis.EvaluatePercentileRank()    
    eggutanalysis.EvaluateBrainDiseaseException()        
    eggutanalysis.CalculateMicrobiomeRatio()
    eggutanalysis.CalculateAverageMicrobiomeRatio()
    eggutanalysis.CalculateHarmfulMicrobiomeAbundance()
    eggutanalysis.CalculateBeneficialMicrobiomeAbundance()
    eggutanalysis.CalculateTotalProbioRatio()
    eggutanalysis.CalculateSpecificProbioRatio()   
    eggutanalysis.CalculateCoordinateLocation() 
    eggutanalysis.DetectProteusMirabilis() 
        
    print('Analysis Complete')
    
    
    
    
    