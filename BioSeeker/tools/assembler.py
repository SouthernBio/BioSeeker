import pandas as pd
import os
import numpy as np
from fnmatch import fnmatch
from BioSeeker.utils.GeneticCode import CODON_TUPLE, CODON_PAIRS_TUPLE
from BioSeeker.core.BioSeekerExceptions import InvalidReadingFrame


class ConservationRateCalculator:
    directory: str
    contents: list[str]
    reading_frame: int
    components: int
    dataframe_index: list
    assembling_codons: bool
    df: pd.DataFrame

    def __init__(self, _directory: str, _reading_frame: int, _components: int, is_codon: bool = True):
        if _reading_frame not in {0, 1, 2}:
            raise InvalidReadingFrame()
        self.directory = _directory
        self.contents = os.listdir(_directory)
        self.reading_frame = _reading_frame
        self.components = _components
        if is_codon:
            self.dataframe_index = list(CODON_TUPLE)
        else:
            self.dataframe_index = list(CODON_PAIRS_TUPLE)
            self.assembling_codons = False

    def __prepare(self):
        dataframe = pd.DataFrame({'ReferenceCount': np.zeros(self.components),
                                  'ConservationCount': np.zeros(self.components)
                                  }, index=self.dataframe_index)

        for file in self.contents:
            match self.assembling_codons:
                case True:
                    if os.path.isfile(os.path.join(self.directory, file) and fnmatch(file, f'history_codons*{self.reading_frame}')):
                        df = pd.read_csv(file, header=0, index_col=0)
                        dataframe = dataframe + df
                case False:
                    if os.path.isfile(os.path.join(self.directory, file) and fnmatch(file, f'history_codon_p*{self.reading_frame}')):
                        df = pd.read_csv(file, header=0, index_col=0)
                        dataframe = dataframe + df

        self.df = dataframe

        return None

    def calculator(self, reading_frame: int, is_codon: bool = True) -> None:
        """
        Assemble individual gene dataframes into one sign

        Args:
            reading_frame (int): the reading frame from which codons or codon pairs are established
            is_codon (bool): determines whether the program is working with individual codons or codon pairs

        Raises:
            InvalidReadingFrame: the provided reading frame was out of range.

        Returns:
            None
        """

        if reading_frame not in {0, 1, 2}:
            raise InvalidReadingFrame()

        # Calculate conservation rate
        self.df["ConservationRate"] = self.df["ConservationCount"] / self.df["ReferenceCount"]

        # Save the file
        os.makedirs('dataframes', exist_ok=True)
        match is_codon:
            case True:
                self.df.to_csv(f'dataframes/codon_data_ORF+{reading_frame}.csv')
            case False:
                self.df.to_csv(f'dataframes/codon_pairs_data_ORF+{reading_frame}.csv')

        return None
