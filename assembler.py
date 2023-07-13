import pandas as pd
import os
import numpy as np
from GeneticCode import CODON_TUPLE, CODON_PAIRS_TUPLE

def assembler(directory: str, ORF: int):
    """Function that assembles individual dataframes with codon and codon pair conservation rates 
       (from a specific reading frame) into a single Pandas dataframe

    Args:
        directory (str): The directory where the gene-specific dataframes are located at
        ORF (int): Reading frame. It can be ORF+0, ORF+1 or ORF+2

    Raises:
        Exception: Invalid reading frame passed as parameter
        Exception: Error occurred during dataframe assembly 

    Returns:
        pandas.DataFrame: dataframe with reference and conservation counts for every codon across all species
        pandas.DataFrame: dataframe with reference and conservation counts for every codon pair across all species
    """
    if ORF not in {0,1,2}:
        raise Exception("Error. ORF must take values in {0,1,2}, but " + str(ORF) + " was given.")

    try:
        codon_list = list(CODON_TUPLE)
        codon_pairs_list = list(CODON_PAIRS_TUPLE)

        zeros1 = np.zeros(61)
        zeros2 = np.zeros(3721)

        df_codons = pd.DataFrame({'reference' : zeros1, 'conservation' : zeros1}, index = codon_list)
        df_bicodons = pd.DataFrame({'reference' : zeros2 , 'conservation' : zeros2}, index = codon_pairs_list)
        
        contents = os.listdir(directory)

        for file in contents:
            if os.path.isfile(os.path.join(directory, file)) and file.startswith("history_codons") and file.endswith(f'{ORF}.csv'):
                dataframe = pd.read_csv(file, header = 0, index_col = 0)
                df_codons = df_codons + dataframe
            elif os.path.isfile(os.path.join(directory, file)) and file.startswith("history_bicodons") and file.endswith(f'{ORF}.csv'):
                dataframe = pd.read_csv(file, header = 0, index_col = 0)
                df_bicodons = df_bicodons + dataframe

        codon_data = df_codons.to_csv(f'codon_data_ORF+{ORF}.csv')
        codon_pairs_data = df_bicodons.to_csv(f'bicodon_data_ORF+{ORF}.csv')

        return codon_data, codon_pairs_data
    except:
        print(f"An unknown error occurred while assembling the dataframes from ORF+{ORF}")
        return None, None