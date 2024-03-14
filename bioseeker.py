#!/usr/bin/env python3

from os import getcwd
import BioSeeker.tools.listfiles as lf
from BioSeeker.tools import seqxtract as sq, conservationrate as cr
from BioSeeker.tools.assembler import assembler
from BioSeeker.tools.delete import delete_csv_files
from BioSeeker.utils.intro import intro_message


def main():
    intro_message()
    #   Generating list of files on current working directory
    BASE_PATH = getcwd()
    files = lf.list_of_files(BASE_PATH)

    # List of unreadable files that will be stored after the procedure is executed

    for indicated_index, file in enumerate(files):
        sequences_array = sq.extract_sequences(file)
        for reading_frame in {0, 1, 2}:
            cr.calculations(sequences_array, indicated_index, reading_frame)

    #   Dataframe assembly
    for reading_frame in {0, 1, 2}:
        assembler(BASE_PATH, reading_frame)

    #   Delete unnecessary CSV files
    delete_csv_files()

    return


if __name__ == "__main__":
    main()
