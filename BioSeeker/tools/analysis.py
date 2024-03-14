import pandas as pd
from scipy.stats import zscore
from fnmatch import fnmatch
from BioSeeker.core import BioSeekerExceptions


def normalize_conservation_rate(file_path: str) -> pd.DataFrame:
    """Normalizes conservation rates by calculating the z-score for each observed conservation rate.

    Args:
        file_path (pd.DataFrame): Pandas dataframe containing conservation rates

    Returns:
        pd.DataFrame: Pandas dataframe with normalized conservation rates
    """
    if not fnmatch(file_path, '*.csv'):
        raise BioSeekerExceptions.NotADataframe()

    dataframe = pd.read_csv(file_path, delimiter=',', header=0, index_col=0)
    
    # Calculate z-score
    dataframe[['NormalizedConservationRate']] = zscore(dataframe[['ConservationRate']], axis=None)
    
    # IMPORTANT: Modify this function (or create a new one)
    # to handle conservation rate values of zero

    return dataframe


def split_codon_pairs(file_path: str) -> pd.DataFrame:
    """Generate two columns with each constituent codon from a codon pair

    Args:
        file_path (str): File which contains relevant
    Raises:
         NotADataframe: file_path does not correspond to a dataframe file
    Returns:
        pd.DataFrame: Pandas dataframe with two new columns
    """
    if not fnmatch(file_path, '*.csv'):
        raise BioSeekerExceptions.NotADataframe()
    
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

    expectedCodonPairConservationRate = pd.Series(expected_product, index=None)
    df_expected_conservation_rate = pd.DataFrame({
        'ExpectedCodonPairConservationRate': expectedCodonPairConservationRate
    })

    dataframe = dataframe.join(df_expected_conservation_rate)

    return dataframe


df = normalize_conservation_rate('codon_data_ORF+0.csv')
