import pandas as pd
from scipy.stats import zscore
from fnmatch import fnmatch

class NotADataframe(Exception):
    """Custom exception for when a function receives a file that is not a CSV dataframe"""
    message = "Error: the function must receive a dataframe in CSV format."
    def __init__(self):
        super().__init__(self.message)


def normalize_conservation_rate(file_path: str) -> pd.DataFrame:
    """Normalizes conservation rates by calculating the z-score for each observed conservation rate.

    Args:
        file (pd.DataFrame): Pandas dataframe containing conservation rates

    Returns:
        pd.DataFrame: Pandas dataframe with normalized conservation rates
    """
    if fnmatch(file_path, '*.csv') == False:
        raise NotADataframe()

    dataframe = pd.read_csv(file_path, delimiter=',', header=0, index_col=0)
    
    # Calculate z-score
    dataframe[['NormalizedConservationRate']] = zscore(dataframe[['ConservationRate']], axis=None)
    
    # IMPORTANT: Modify this function (or create a new one)
    # to handle conservation rate values of zero

    return dataframe


# ----------------------------------------------------------------
# This might be deleted later on...
"""
def clean_data(file_path: str) -> pd.DataFrame:

    if fnmatch(file_path, '*.csv') == False:
        raise NotADataframe()
    
    dataframe = pd.read_csv(file_path, delimiter=',', header=0)

    # Replace column name
    if fnmatch(file_path, '*bicodon*'):
        dataframe.columns = dataframe.columns.str.replace('Unnamed: 0', 'CodonPair')
    else:
        dataframe.columns = dataframe.columns.str.replace('Unnamed: 0', 'Codon')

    # Remove rows with reference counts of zero
    corrected_dataframe = dataframe.loc[dataframe['ReferenceCount'] > 0]

    return dataframe
"""
# ----------------------------------------------------------------


def split_codon_pairs(file_path: str) -> pd.DataFrame:
    """Generate two columns with each constituent codon from a codon pair

    Args:
        file_path (str): File which contains relevant 
    Returns:
        pd.DataFrame: Pandas dataframe with two new columns
    """
    if fnmatch(file_path, '*.csv') == False:
        raise NotADataframe()
    
    dataframe = pd.read_csv(file_path, delimiter=',', header=0)

    # Get codon pairs
    codon_pairs = dataframe[dataframe.columns[0]]
    
    # Generate lists with each constituent codon
    first_codon, second_codon = [], []
    for pair in codon_pairs:
        first_codon.append(pair[:3])
        second_codon.append(pair[3:])
    
    # Join the lists to the dataframe
    df_codons = pd.DataFrame({
        'FirstCodon': pd.Series(first_codon, index=None),
        'SecondCodon': pd.Series(second_codon, index=None)
        })
    dataframe = dataframe.join(df_codons)

    return dataframe


def expected_codon_pair_conservation_rate(file_path: str, codons: str) -> pd.DataFrame:
    
    codon_information = pd.read_csv(codons, delimiter=',', header=0, index_col=0)
    dataframe = pd.read_csv(file_path, delimiter=',', header=0, index_col=0)
    
    # Auxiliary lists to store normalized conservation values
    first_codon_cr, second_codon_cr, expected_product = [], [], []
    
    for codon in dataframe['FirstCodon']:
        for i in codon_information['Codon']:
            if codon == i:
                first_codon_cr.append(codon_information.filter(index=i, axis=0)[['ConservationRate']]) 

    for codon in dataframe['SecondCodon']:
        for i in codon_information['Codon']:
            if codon == i:
                second_codon_cr.append(codon_information.filter(index=i, axis=0)[['ConservationRate']])
    
    for index, item in enumerate(first_codon_cr):
        expected_product.append(first_codon_cr[index] * second_codon_cr[index])

    ExpectedCodonPairConservationRate = pd.Series(expected_product, index=None)
    df_expectedcr = pd.DataFrame({
        'ExpectedCodonPairConservationRate': ExpectedCodonPairConservationRate
    })

    dataframe = dataframe.join(df_expectedcr)

    return dataframe