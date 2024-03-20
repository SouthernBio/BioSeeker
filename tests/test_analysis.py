import unittest
import pandas as pd
from scipy.stats import zscore
from BioSeeker.core.BioSeekerExceptions import NotADataframe
from BioSeeker.tools.analysis import normalize_conservation_rate, split_codon_pairs


class TestNormalizeConservationRate(unittest.TestCase):
    def test_file_path_not_csv(self):
        with self.assertRaises(NotADataframe):
            normalize_conservation_rate('not_a_csv.txt')

    def test_normalization(self):
        # Creamos un DataFrame de ejemplo y lo guardamos en un archivo CSV
        data = {'ConservationRate': [1, 2, 3, 4, 5]}
        df = pd.DataFrame(data)
        csv_file_path = 'example_data.csv'
        df.to_csv(csv_file_path)

        # Normalizamos usando la función con la ruta del archivo CSV
        normalized_df = normalize_conservation_rate(csv_file_path)

        # Comprobamos que la columna de Conservación Normalizada esté presente
        self.assertIn('NormalizedConservationRate', normalized_df.columns)

        # Comprobamos que la normalización se haya realizado correctamente
        expected_normalized_values = zscore(df['ConservationRate'])
        for i, row in normalized_df.iterrows():
            self.assertAlmostEqual(row['NormalizedConservationRate'], expected_normalized_values[i])

        # Eliminamos el archivo CSV de ejemplo después de las pruebas
        import os
        os.remove(csv_file_path)


class TestSplitCodonPairs(unittest.TestCase):
    def test_file_path_not_csv(self):
        with self.assertRaises(NotADataframe):
            split_codon_pairs('not_a_csv.txt')

    def test_split_codon_pairs(self):
        # Creamos un DataFrame de ejemplo
        data = {'CodonPairs': ['ATGCCG', 'CTAGGA', 'TTAGGC']}
        df = pd.DataFrame(data)
        csv_file_path = 'example_codon_pairs.csv'
        df.to_csv(csv_file_path)

        # Ejecutamos la función con el archivo CSV de ejemplo
        result_df = split_codon_pairs(csv_file_path)

        # Verificamos que los nombres de las columnas sean correctos
        self.assertEqual(result_df.columns.tolist(), ['CodonPairs', 'FirstCodon', 'SecondCodon'])

        # Verificamos que los valores de las columnas FirstCodon y SecondCodon sean correctos
        expected_first_codons = ['ATG', 'CTA', 'TTA']
        expected_second_codons = ['CCG', 'GGA', 'GGC']
        for i, row in result_df.iterrows():
            self.assertEqual(row['FirstCodon'], expected_first_codons[i])
            self.assertEqual(row['SecondCodon'], expected_second_codons[i])

        # Eliminamos el archivo CSV de ejemplo después de las pruebas
        import os
        os.remove(csv_file_path)


if __name__ == '__main__':
    unittest.main()
