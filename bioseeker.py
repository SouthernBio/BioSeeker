#!/usr/bin/env python3

import tools.listfiles as lf
import tools.seqxtract as sq
import tools.conservationrate as cr
from tools.assembler import assembler
from utils.intro import intro_message
import os
import sh
import sys

def main():
    #   Generating list of files on current working directory
    BASE_PATH = os.getcwd()
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
    assembler(BASE_PATH, 0)  # ORF+0
    assembler(BASE_PATH, 1)  # ORF+1
    assembler(BASE_PATH, 2)  # ORF+2

    #   Delete unnecessary CSV files
    try:
        # using os.system for file removal/terminal clearing as
        # sh library does not support wildcard completion/terminal clearing
        if os.name == "posix":
            os.system("rm history_*.csv")
            os.system("clear")
        else:
            os.system("del history_*.csv")
            os.system("cls")
        sh.echo('ALL FILES ASSEMBLED.', _out=sys.stdout, _err=sys.stderr)
    except Exception as error:
        print(f"An error was found during the deletion of unnecessary CSV files ({error})")

    return 0


def run():
    intro_message()
    main()
    return 0

if __name__ == "__main__":
    run()
