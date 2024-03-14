import numpy as np
import pandas as pd
from BioSeeker.utils.GeneticCode import CODON_TUPLE, CODON_PAIRS_TUPLE
from BioSeeker.core import BioSeekerExceptions


class ConservationRateCalculator:
    """
    Fill this up.
    """
    reading_frame: int                          # {0, 1, 2}

    codons: list                                # Codons present in each genetic sequence
    codon_pairs: list                           # Codon pairs present in each genetic sequence
    codon_size: int                             # This attribute is created to avoid using a hardcoded value of 3
    sequences_array: list  #

    number_of_reference_codons: int             # len(codon_reference)
    number_of_reference_codon_pairs: int        # len(codon_pairs_reference)

    codon_history: np.ndarray                   # np.zeros(number_of_reference_codons)
    conserved_codon_history: np.ndarray         # np.zeros(number_of_reference_codons)

    codon_pair_history: np.ndarray              # np.zeros(number_of_reference_codon_pairs)
    conserved_codon_pair_history: np.ndarray    # np.zeros(number_of_reference_codon_pairs)

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

    def __create_arrays(self, sequences_array: list, reading_frame: int):
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

    def __establish_references(self):
        # Numpy: VisibleDeprecationWarning, data type = object.
        codons, codon_pairs = np.asarray(self.codons), np.asarray(self.codon_pairs)
        codon_reference, codon_pairs_reference = codons[0], codon_pairs[0]

        return codon_reference, codon_pairs_reference

    def __compute_reference(self, reference, is_codon: bool):
        for count, i in enumerate(reference, start=0):
            match is_codon:
                case True:
                    self.codon_history = np.zeros(len(reference), dtype='i')
                    slicer = self.codons[:, count]
                    for j in slicer:
                        if j == i:
                            self.codon_history[count] += 1
                case False:
                    self.codon_pair_history = np.zeros(len(reference), dtype='i')
                    slicer = self.codon_pairs[:, count]
                    for k in slicer:
                        if k == i:  # Changed from 'j' to 'k', see if it works
                            self.codon_pair_history[count] += 1

        return None

    def __compute_conservation(self):
        pass

    def __fill_columns(self):
        pass

    def __create_dataframes(self):
        pass

    def calculate_rate(self):
        pass

    def parse_info(self):
        pass


def calculations(sequences_array: list, indicated_index: int, reading_frame: int):
    """Function to calculate codon and codon pair conservation rates

    Args:
        sequences_array (list): An array of sequences
        indicated_index (int): auxiliary index to identify files
        reading_frame (int): Reading frame. Can only take values in {0, 1, 2}

    Raises:
        InvalidReadingFrame: invalid reading frame.

    Returns:
        pandas.DataFrame: codon dataframe
        pandas.DataFrame: codon pair dataframe
    """
    if reading_frame not in {0, 1, 2}:
        raise BioSeekerExceptions.InvalidReadingFrame()

    codons = []
    bicodons = []
    codon_size = 3

    #   Using the obtained array in the previous function, we divide the sequences in fragments of 3 and 6 bases
    #   Groups of codons and bicodons are assembled with a "codon_size" step of 3 nucleotides
    print(f"Creating codon and codon pair arrays at reading+{reading_frame}...")
    for _, k in sequences_array:
        for j in k:
            codon = [j[i:i+codon_size] for i in range(reading_frame, len(j), codon_size)]
            bicodon = [j[i:i+2*codon_size] for i in range(reading_frame, len(j), codon_size)]
            codons.append(codon)
            bicodons.append(bicodon)

    #   Converting to Numpy array
    codons = np.asarray(codons)
    bicodons = np.asarray(bicodons)  # Numpy: VisibleDeprecationWarning, dtype = object.

    #   Reference lists
    print(f"Establishing reference sequences at reading+{reading_frame}...")
    codon_reference = codons[0]
    bicodon_reference = bicodons[0]

    """
    1. Crear numpy array de codones y pares de codones
    2. Establecer cuál es la referencia de codones y pares de codones (como codons[0])
    3. Calcular el largo de la lista de codones de referencia (len(codon_reference))
    3.1. Guardarlo en un atributo de clase?
     
    """
    #   Creating accessory lists to store the calculations
    number_of_reference_codons = len(codon_reference)
    number_of_reference_bicodons = len(bicodon_reference)

    codon_history = np.zeros(number_of_reference_codons, dtype="i")
    bicodon_history = np.zeros(number_of_reference_bicodons, dtype="i")

    reference_codon_history = np.zeros(61, dtype="i")
    reference_bicodon_history = np.zeros(3721, dtype="i")

    conserved_codon_history = np.zeros(number_of_reference_codons, dtype="i")
    conserved_bicodon_history = np.zeros(number_of_reference_bicodons, dtype="i")

    codon_conservation = np.zeros(61, dtype="i")
    bicodon_conservation = np.zeros(3721, dtype="i")

    #   Calculations
    print(f"Computing conservation for codon_history and bicodon_history at reading+{reading_frame}...")
    for count, i in enumerate(codon_reference, start=0):
        slicer = codons[:, count]
        for j in slicer:
            if j == i:
                codon_history[count] += 1

    for count, i in enumerate(bicodon_reference, start=0):
        slicer = bicodons[:, count]
        for j in slicer:
            if j == i:
                bicodon_history[count] += 1

    print("Computation successful.")
    print(f"Computing conservation history at reading+{reading_frame}...")

    aux1 = len(codons[:, 0])
    aux2 = len(bicodons[:, 0])

    for count, x in enumerate(codon_history, start=0):
        if x/aux1 > 0.9:
            conserved_codon_history[count] += 1

    for count, x in enumerate(bicodon_history, start=0):
        if x/aux2 > 0.9:
            conserved_bicodon_history[count] += 1

    print("Computation successful.")
    # This part of the code is computationally expensive. We have to find a way to optimize it
    print(f"Assembling conservation matrix(cod_ref, hisc) for reading+{reading_frame}...")

    conservation_matrix_reference_codons = np.column_stack((np.asarray(codon_reference),
                                                            np.asarray(conserved_codon_history)))
    conservation_matrix_reference_bicodons = np.column_stack((np.asarray(bicodon_reference),
                                                              np.asarray(conserved_bicodon_history)))

    for count, codon in enumerate(CODON_TUPLE, start=0):
        for i in conservation_matrix_reference_codons:
            if (i[1] != '0') and (i[0] == codon):
                codon_conservation[count] += 1

    for count, bicodon in enumerate(CODON_PAIRS_TUPLE, start=0):
        for i in conservation_matrix_reference_bicodons:
            if (i[1] != '0') and (i[0] == bicodon):
                bicodon_conservation[count] += 1

    for count, cod in enumerate(CODON_TUPLE, start=0):
        for i in codon_reference:
            if cod == i:
                reference_codon_history[count] += 1

    for count, bicod in enumerate(CODON_PAIRS_TUPLE, start=0):
        for i in bicodon_reference:
            if bicod == i:
                reference_bicodon_history[count] += 1

    codon_array = np.column_stack((list(CODON_TUPLE), reference_codon_history, codon_conservation))
    bicodon_array = np.column_stack((list(CODON_PAIRS_TUPLE), reference_bicodon_history, bicodon_conservation))

    print(f"The assembly has been successful. Creating dataframes for reading+{reading_frame}...")

    #   Creating dataframes and changing column names
    codon_dataframe = pd.DataFrame(codon_array)
    bicodon_dataframe = pd.DataFrame(bicodon_array)

    codon_dataframe = codon_dataframe.rename(columns={0: 'codon', 1: 'ReferenceCount', 2: 'ConservationCount'})
    bicodon_dataframe = bicodon_dataframe.rename(columns={0: 'codon_pair', 1: 'ReferenceCount', 2: 'ConservationCount'})

    #   Exporting dataframes to CSV files
    dataframe_codons = codon_dataframe.to_csv("history_codons" + "_" + str(indicated_index)
                                              + "_" + str(reading_frame) + ".csv", index=False)
    dataframe_bicodons = bicodon_dataframe.to_csv("history_bicodons" + "_" + str(indicated_index)
                                                  + "_" + str(reading_frame) + ".csv", index=False)

    return dataframe_codons, dataframe_bicodons
