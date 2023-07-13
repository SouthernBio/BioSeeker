# Work in progress...
import pandas as pd

codon_orf_1 = pd.read_csv("codon_data_ORF+1.csv")
codon_orf_1["ConservationRate"] = codon_orf_1["conservation"] / codon_orf_1["reference"]

print(codon_orf_1)