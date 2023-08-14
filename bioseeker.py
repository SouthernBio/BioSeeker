#!/usr/bin/env python3

import tools.listfiles as lf
import tools.seqxtract as sq
import tools.conservationrate as cr
from tools.assembler import assembler
from tools.delete import delete_csv_files
from utils.intro import intro_message
from os import getcwd

def main():
    intro_message()
    #   Generating list of files on current working directory
    BASE_PATH = getcwd()
    files = lf.list_of_files(BASE_PATH)

    # List of unreadable files that will be stored after the procedure is executed

    for indicated_index, file in enumerate(files):
        sequences_array = sq.extract_sequences(file)
        for ORF in {0, 1, 2}:
            cr.calculations(sequences_array, indicated_index, ORF)

    #   Dataframe assembly
    for ORF in {0, 1, 2}:
        assembler(BASE_PATH, ORF)

    #   Delete unnecessary CSV files
    delete_csv_files()

    return


if __name__ == "__main__":
    main()
