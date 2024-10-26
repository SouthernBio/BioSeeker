#!/usr/bin/env python3

import tools.listfiles as lf
import tools.seqxtract as sq
import tools.conservationrate as cr
from tools.assembler import assembler
from tools.delete import delete_csv_files
from utils.intro import intro_message
from os import getcwd
from tools.upload_icps import upload_icps


def main():
    intro_message()
    
    # Ask the user to upload an ICP file
    file_path = input("Please enter the path of the ICP file to upload: ")
    upload_message = upload_icps(file_path)
    print(upload_message)

    
    #   Generating list of files on current working directory
    BASE_PATH = getcwd()
    files = lf.list_of_files(BASE_PATH)

    # List of unreadable files that will be stored after the procedure is executed
    unreadable_files = []

    for indicated_index, file in enumerate(files):
        try:
            sequences_array = sq.extract_sequences(file)
            for ORF in {0, 1, 2}:
                cr.calculations(sequences_array, indicated_index, ORF)
        except:
            print(f"An error was found during the analysis of {file}. It will be added to unreadable.txt")
            unreadable_files.append(file)

    print("The analysis is over. The following files have been ommited:\n", unreadable_files)

    #   Unreadable MSA files
    textfile = open("unreadable.txt", "w")
    for file in unreadable_files:
        textfile.write(file + "\n")
    textfile.close()

    #   Dataframe assembly
    for ORF in {0, 1, 2}:
        assembler(BASE_PATH, ORF)

    #   Delete unnecessary CSV files
    delete_csv_files()

    return


if __name__ == "__main__":
    main()
