##<Usage: python Script.py {path_exp}>
### ex) python EGgutPro_update_mrs.py "/home/kbkim/EgGutPro/input/EGgutPro_mirror_output_3175.csv"

import os, datetime
import pandas as pd
import sys
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import percentileofscore
from skbio.stats.composition import multiplicative_replacement, clr

# Check if the script is being called with the correct arguments
if len(sys.argv) != 2:
    print("Usage: python Script.py <path_exp>")
    print("Example: python EGgutPro_update_mrs.py \"/home/kbkim/EgGutPro/input/EGgutPro_mirror_output_3175.csv\"")
    sys.exit(1)
    
# path_exp : Path of Merged Proportion file to analyze
path_exp = sys.argv[1]

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
        
# Histogram Plot - mrs 
def save_histograms_to_file(df, filename):
    num_rows = df.shape[1]
    fig, axs = plt.subplots(num_rows, 1, figsize=(8, 6*num_rows))
    
    for i in range(num_rows):
        axs[i].hist(df.iloc[:,i], bins=10)
        axs[i].set_title(df.columns.to_list()[i])
    
    plt.tight_layout()
    plt.savefig(filename)

def filter_species(taxon_str_for_domain):
    splited_by_taxonomy = taxon_str_for_domain.split("|")
    species_name = splited_by_taxonomy[len(splited_by_taxonomy) - 1]  # 계문강목과속종 중 종만 필터링!!!
    return species_name

###################################
# MainClass
###################################
class EgGutProUpdateMRS:
    def __init__(self, path_exp, fplog=None):
        """
        Initializes a EgGutProUpdateMRS object.

        Parameters:
        path_exp (str): Path of Merged Proportion file to analyze.  
        """   
        
        self.path_exp = path_exp
        self.__fplog=fplog
        
        curdir = os.path.abspath('')
        self.path_ref = f"{curdir}/input/EGgutPro_mircobe_list.xlsx"
        self.path_healthy = f"{curdir}/input/EGgutPro_healthy_person_profile_v2.xlsx"
        self.path_db = f"{curdir}/input/EGgutPro_db_abundance.xlsx"
        self.path_mrs_db = f"{curdir}/input/EGgutPro_mrs_db.xlsx"
        self.path_hist = f"{curdir}/output/EGgutPro_mrs_hist.png"
        self.path_percentile_rank_db = f"{curdir}/input/EGgutPro_percentile_rank_db.csv"

        self.df_beta = None
        self.df_db = None
        self.df_exp = None
        self.df_mrs = None
        self.df_db_rev = None
        self.df_dysbiosis = None
        self.df_percentile_rank = None

        self.li_diversity = None
        self.li_new_sample_name = None
        self.li_phenotype = None
        self.li_microbiome = None            

            
    # Load the DB file
    # df_beta : Data frame of of Phenotype-Microbiome information
    # df_db : Data frame of accumulated Experimental result information - Abundance
    # df_exp : Data frame of Experimental result information - Abundance    
    def ReadDB(self):   
        """
        Read the data.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        
        rv = True
        rvmsg = "Success"
        
        try:       
            self.df_beta = pd.read_excel(self.path_ref, sheet_name = f"disease_taxa")
            self.df_dysbiosis = pd.read_excel(self.path_ref, sheet_name = f"harmful_beneficial_taxa")
            self.df_probio = pd.read_excel(self.path_ref, sheet_name = f"probio_taxa")

            self.df_healthy = pd.read_excel(self.path_healthy, sheet_name="RA")
            self.df_db = pd.read_excel(self.path_db)
            self.df_exp = pd.read_csv(self.path_exp)
    
            self.df_beta.rename(columns = {"Disease": "phenotype", "NCBI name": "ncbi_name", "MIrROR name": "microbiome", "Health sign": "beta", "subtract": "microbiome_subtract"}, inplace=True)
            self.df_beta = self.df_beta[["phenotype", "ncbi_name", "microbiome", "beta", "microbiome_subtract"]]
            self.df_beta['beta'] = self.df_beta['beta'].replace({'유해': 1, '유익': -1})    

            self.df_dysbiosis.rename(columns = {"NCBI name": "ncbi_name", "MIrROR name": "microbiome", "Health sign": "beta", "subtract": "microbiome_subtract", "Exclude": "exclude"}, inplace=True)
            self.df_dysbiosis = self.df_dysbiosis[["ncbi_name", "microbiome", "beta", "microbiome_subtract", "exclude"]]
            self.df_dysbiosis['beta'] = self.df_dysbiosis['beta'].replace({'유해': 1, '유익': -1})
            
            self.df_healthy = self.df_healthy[self.df_healthy['Taxonomy'].str.contains('s__')]
            self.df_healthy['Taxonomy'] = self.df_healthy['Taxonomy'].apply(filter_species)
            self.df_healthy['Taxonomy'] = self.df_healthy['Taxonomy'].str.replace(' ','_')
            self.df_healthy = self.df_healthy.rename(columns={'Taxonomy': 'taxa'})            
            
            print(self.df_healthy)
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            sys.exit()    
            
        return rv, rvmsg


    # Insert data into DB - Merge the data frame df_db & df_exp
    def InsertDataDB(self): 
        """
        Inserts data into the database by merging the data frames df_db and df_exp.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """   
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        rv = True
        rvmsg = "Success"
        
        try: 
            self.df_db = pd.merge(self.df_db, self.df_exp, how='outer',on='taxa', suffixes=['', '_right']) 
            self.df_db = self.df_db.fillna(0)
            self.df_db = self.df_db.filter(regex='^(?!.*_right).*') # Eliminate duplicate columns

            # Update the data - Convert df_exp to df_db
            self.df_exp = self.df_db       
                      
            self.df_db_rev = self.df_db.set_index(keys=['taxa'], inplace=False, drop=True)    
            self.df_db_rev.to_excel(self.path_db)

            # Delete the diversity, observed rows
            if (list(self.df_exp['taxa'][0:2]) == ['diversity', 'observed']) & (list(self.df_db['taxa'][0:2]) == ['diversity', 'observed']):
                self.li_diversity = list(self.df_exp.iloc[0,1:]) # li_diversity : Alpha-Diversity list 
                self.df_exp = self.df_exp.iloc[2:,:]
                self.df_db = self.df_db.iloc[2:,:]

            # li_new_sample_name : Sample name list 
            # li_phenotype : Phenotype list 
            self.li_new_sample_name = list(self.df_exp.columns)[1:]  
            self.li_phenotype = list(dict.fromkeys(self.df_beta['phenotype']))
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
            
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
            sys.exit() 
    
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
                    self.df_percentile_rank.loc[self.li_new_sample_name[i], self.li_phenotype[j]] = (percentileofscore(list(self.df_mrs[self.li_phenotype[j]]), self.df_mrs.loc[self.li_new_sample_name[i], self.li_phenotype[j]], kind='mean')).round(1)
                 
            # Outliers
            # Replace percentile ranks that are less than or equal to 5 with 5, and those that are greater than or equal to 95 with 95
            for i in range(len(self.li_phenotype)):
                self.df_percentile_rank.loc[self.df_percentile_rank[self.li_phenotype[i]]<=5, self.li_phenotype[i]] = 5.0
                self.df_percentile_rank.loc[self.df_percentile_rank[self.li_phenotype[i]]>=95, self.li_phenotype[i]] = 95.0      
                
            #self.df_percentile_rank['TotalScore'] = (self.df_percentile_rank['Dysbiosis']*1.1 + self.df_percentile_rank['HealthyDistance']*1.1 + self.df_percentile_rank['Diversity']*0.8)/3
            
            #self.df_percentile_rank['TotalScore'] = self.df_percentile_rank['TotalScore'].astype(float).round(1)
                       
            # Replace missing values with the string 'None'    
            self.df_percentile_rank = self.df_percentile_rank.fillna('None')
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
    
        return rv, rvmsg
    
    def UpdateMRS(self): 
        """
        Save the MRS data as an Excel file & Save the histogram as an png file.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """         
        myNAME = self.__class__.__name__+"::"+sys._getframe().f_code.co_name
        WriteLog(myNAME, "In", type='INFO', fplog=self.__fplog)
        rv = True
        rvmsg = "Success"
        
        try:             
            
            # Save the output file - Percentile Rank of the samples
            self.df_percentile_rank.to_csv(self.path_percentile_rank_db, encoding="utf-8-sig", index_label='serial_number')            
            
            # Save the df_mrs
            self.df_mrs.to_excel(self.path_mrs_db)
   
            # Histogram Plot - mrs 
            save_histograms_to_file(self.df_mrs, self.path_hist)    

        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print(f"Error has occurred in the {myNAME} process")    
            sys.exit()
    
        return rv, rvmsg       
    
####################################
# main
####################################
if __name__ == '__main__':
    
    eggutpro = EgGutProUpdateMRS(path_exp)
    eggutpro.ReadDB()
    eggutpro.InsertDataDB()
    eggutpro.CalculateMRS()    
    eggutpro.CalculateDysbiosis()  
    eggutpro.CalculateHealthyDistance()     
    eggutpro.CalculatePercentileRank() 
    eggutpro.UpdateMRS() 
    
    print('Update Complete') 