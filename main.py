#!/usr/bin/env python3
import listfiles as lf
import seqxtract as sq
import conservationrate as cr
from assembler import assembler
import os
import subprocess

def main():
    #   Generating list of files on current working directory
    BASE_PATH = os.getcwd()
    files = lf.list_of_files(BASE_PATH)

    # List of unreadable files that will be stored after the procedure is executed
    unreadable_files = []

    indicated_index = 0
    for file in files:
        try:
            sequences_array = sq.extract_sequences(file)
            cr.calculations(sequences_array, indicated_index, 0)
            cr.calculations(sequences_array, indicated_index, 1)
            cr.calculations(sequences_array, indicated_index, 2)
            indicated_index += 1
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
    assembler(BASE_PATH, 0)   # ORF+0
    assembler(BASE_PATH, 1)    # ORF+1
    assembler(BASE_PATH, 2)    # ORF+2

    #   Delete unnecessary CSV files
    try:
        if os.name ==  "posix":
                subprocess.call('./delete_files.sh')
        else:
            subprocess.call([BASE_PATH + r'\\delete_files.bat'])
    except:
        print("An error was found during the deletion of unnecessary CSV files")

    return 0

if __name__ == '__main__':
    main()