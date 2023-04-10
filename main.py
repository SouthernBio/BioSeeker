import listfiles as lf
import seqxtract as sq
import conservationrate as cr
from assembler import assembler
import os
import subprocess

#   Generating list of files
files = lf.list_of_files(r"C:\Users\YourUsername\WorkingDirectory")

#lista_de_dataframes_cod3 = []                                                   ???
#lista_de_dataframes_cod6 = []                                                   ???
unreadable_files = []

indicated_index = 0
for file in files:
    try:
        sequences_array = sq.extract_sequences(file)
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
directorio = "C:\\Users\\YourUsername\\WorkingDirectory\\"
assembler(directorio, 0)   # ORF+0
assembler(directorio, 1)    # ORF+1
assembler(directorio, 2)    # ORF+2

#   Delete unnecessary CSV files
if os.name ==  "posix":
    subprocess.call('./delete_files.sh')
else:
    subprocess.call([r'C:\Users\YourUSername\WorkingDirectory\delete_files.bat'])