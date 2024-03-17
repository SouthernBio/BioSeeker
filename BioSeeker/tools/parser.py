import numpy as np
import pandas as pd
from BioSeeker.utils.GeneticCode import CODON_TUPLE, CODON_PAIRS_TUPLE
from BioSeeker.core import BioSeekerExceptions


class FASTAParser:
    """
    Fill this up.
    """
    reading_frame: int                          # {0, 1, 2}
    codons: list                                # Codons present in each genetic sequence
    codon_pairs: list                           # Codon pairs present in each genetic sequence
    codon_size: int                             # This attribute is created to avoid using a hardcoded value of 3
    sequences_array: list
    number_of_reference_codons: int             # len(codon_reference)
    number_of_reference_codon_pairs: int        # len(codon_pairs_reference)
    codon_history: np.ndarray                   # np.zeros(number_of_reference_codons)
    conserved_codon_history: np.ndarray         # np.zeros(number_of_reference_codons)
    codon_pair_history: np.ndarray              # np.zeros(number_of_reference_codon_pairs)
    conserved_codon_pair_history: np.ndarray    # np.zeros(number_of_reference_codon_pairs)
    codon_reference: np.ndarray                 # Genetic sequence used as reference for codons
    codon_pair_reference: np.ndarray            # Genetic sequence used as reference for codon pairs
    reference_codon_history: np.ndarray         # How many times a codon is present in a reference sequence
    reference_codon_pair_history: np.ndarray    # How many times a codon pair is present in a reference sequence
    codon_conservation: np.ndarray              # How many times said codon has been preserved
    codon_pair_conservation: np.ndarray         # How many times said codon pair has been preserved

    def __init__(self):
        self.reading_frame = 0
        self.codon_size = 0
        self.codons = []
        self.codon_pairs = []
        self.sequences_array = []
        self.reference_codon_history = np.zeros(61, dtype="i")
        self.reference_codon_pair_history = np.zeros(3721, dtype="i")
        self.codon_conservation = np.zeros(61, dtype="i")
        self.codon_pair_conservation = np.zeros(3721, dtype="i")

    def __create_arrays(self, sequences_array: list, reading_frame: int) -> None:
        if reading_frame not in {0, 1, 2}:
            raise BioSeekerExceptions.InvalidReadingFrame()
        else:
            self.reading_frame = reading_frame

        for _, k in sequences_array:
            for j in k:
                codon = [j[i:i + self.codon_size] for i in range(reading_frame, len(j), self.codon_size)]
                codon_pair = [j[i:i + 2 * self.codon_size] for i in range(self.reading_frame, len(j), self.codon_size)]
                self.codons.append(codon)
                self.codon_pairs.append(codon_pair)

        return None

    def __establish_references(self) -> None:
        # Numpy: VisibleDeprecationWarning, data type = object.
        codons, codon_pairs = np.asarray(self.codons), np.asarray(self.codon_pairs)
        self.codon_reference, self.codon_pairs_reference = codons[0], codon_pairs[0]

        return None

    # noinspection PyMethodMayBeStatic
    def __compute_reference(self, reference: np.ndarray, codons_or_pairs) -> np.ndarray:
        """
        Method for computing reference sequence values. IMPORTANT: the output of this function must be added to
        self.codon_history or self.codon_pair_history depending on the context.

        Args:
            reference (list): genetic sequence of reference
            codons_or_pairs (list): an alias for self.codons or self.codon_pairs

        Returns:
            history (np.ndarray): an array with the presence history of each codon (or codon pair) from the reference
            sequence.
        """
        history = np.zeros(len(reference), dtype='i')
        for count, i in enumerate(reference, start=0):
            slicer = codons_or_pairs[:, count]
            for j in slicer:
                if j == i:
                    history[count] += 1

        return history

    # noinspection PyMethodMayBeStatic
    def __compute_conservation(self, codons_or_pairs: list, history: np.ndarray, conserved_history: np.ndarray) -> None:
        """
        Method for computing conservation history

        Args:
            codons_or_pairs (list): an alias for self.codons or self.codon_pairs
            history (list): an alias for self.codon_history or self.codon_pair_history
            conserved_history (list): an alias for self.conserved_codon_history or self.conserved_codon_pair_history

        Returns:
            None
        """
        aux = len(codons_or_pairs[:, 0])
        for count, x in enumerate(history, start=0):
            if x / aux > 0.9:
                conserved_history[count] += 1

        return None

    # noinspection PyMethodMayBeStatic
    def __fill_columns(self,
                       reference: np.ndarray,
                       history: np.ndarray,
                       m_tuple: tuple,
                       conservation_array: np.ndarray,
                       reference_history: np.ndarray) -> None:
        """
        Method for filling the conservation matrix reference columns for codons and codon pairs

        Args:
            reference: an alias for self.codon_reference or self.codon_pair_reference
            history: an alias for self.conserved_codon_history or self.conserved_codon_pair_history
            m_tuple (tuple): these tuples are CODON_TUPLE and CODON_PAIR_TUPLE from GeneticCode.py
            conservation_array: an alias for self.codon_conservation and self.codon_pair_conservation
            reference_history: an alias for self.reference_codon_history or self.reference_codon_pair_history

        Returns:
            None
        """
        conservation_matrix_reference = np.column_stack((np.asarray(reference), np.asarray(history)))

        for count, x in enumerate(m_tuple, start=0):
            for i in conservation_matrix_reference:
                if (i[1] != '0') and (i[0] == x):
                    conservation_array[count] += 1
            for j in reference:
                if x == j:
                    reference_history[count] += 1

        return None

    # noinspection PyMethodMayBeStatic
    def __create_dataframes(self,
                            m_tuple: tuple,
                            reference_array: np.ndarray,
                            conservation_array: np.ndarray,
                            is_codon: bool,
                            index: int,
                            reading_frame: int):
        """
        Method to generate the final dataframes to perform the statistical analysis

        Args:
            m_tuple (tuple): these tuples are CODON_TUPLE and CODON_PAIR_TUPLE from GeneticCode.py
            reference_array (np.ndarray): an alias for self.codon_reference or self.codon_pair_reference
            conservation_array (np.ndarray): an alias for self.codon_conservation or self.codon_pair_conservation
            is_codon (bool): specifies whether the generated dataframe contains information for single codons
            index (int): a numeric value to organize files by the order in which they are parsed
            reading_frame (int): the reading frame from which codons and codon pairs are determined
        """
        # Check input
        if reading_frame not in {0, 1, 2}:
            raise BioSeekerExceptions.InvalidReadingFrame()

        m_array = np.column_stack((list(m_tuple), reference_array, conservation_array))

        dataframe = pd.DataFrame(m_array)
        match is_codon:
            case True:
                dataframe = dataframe.rename(columns={0: 'codon', 1: 'ReferenceCount', 2: 'ConservationCount'})
                dataframe.to_csv(f"history_codons_{str(index)}_{str(reading_frame)}.csv", index=False)
            case False:
                dataframe = dataframe.rename(columns={0: 'codon_pair', 1: 'ReferenceCount', 2: 'ConservationCount'})
                dataframe.to_csv(f"history_codon_pairs_{str(index)}_{str(reading_frame)}.csv", index=False)

        return None

    def calculate_rate(self, sequences_array: list, indicated_index: int, m_reading_frame: int) -> None:
        """
        Public method to calculate codon and codon pair reference values and conservation count.

        Args:
            sequences_array (list)
            indicated_index (int): an integer value used to order files by number
            m_reading_frame (int): the reading frame from which codons and codon pairs are determined

        Raises:
            InvalidReadingFrame: reading frame can only take values in {0, 1, 2}, there are no other reading frames

        Returns:
            None
        """
        if m_reading_frame not in {0, 1, 2}:
            raise BioSeekerExceptions.InvalidReadingFrame()

        self.__create_arrays(sequences_array, m_reading_frame)
        self.__establish_references()

        self.codon_history = self.__compute_reference(self.codon_reference, self.codons)
        self.codon_pair_history = self.__compute_reference(self.codon_pair_reference, self.codon_pairs)

        self.__compute_conservation(self.codons, self.codon_history, self.conserved_codon_history)
        self.__compute_conservation(self.codon_pairs, self.codon_pair_history, self.conserved_codon_pair_history)

        # For codons
        self.__fill_columns(self.codon_reference,
                            self.conserved_codon_history,
                            CODON_TUPLE,
                            self.codon_conservation,
                            self.reference_codon_history)

        # For codon pairs
        self.__fill_columns(self.codon_pair_reference,
                            self.conserved_codon_pair_history,
                            CODON_PAIRS_TUPLE,
                            self.codon_pair_conservation,
                            self.reference_codon_pair_history)

        # Create dataframe for codons
        self.__create_dataframes(CODON_TUPLE,
                                 self.codon_reference,
                                 self.codon_conservation,
                                 is_codon=True,
                                 index=indicated_index,
                                 reading_frame=m_reading_frame)

        # Create dataframe for codon pairs
        self.__create_dataframes(CODON_PAIRS_TUPLE,
                                 self.codon_pair_reference,
                                 self.codon_pair_conservation,
                                 is_codon=False,
                                 index=indicated_index,
                                 reading_frame=m_reading_frame)

        return None

    def parse_info(self):
        pass
