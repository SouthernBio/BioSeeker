import os
from BioSeeker.core import BioSeekerExceptions


def extract_sequences(file: str):
    """Function to extract sequences from FASTA files
        This function takes as an argument a file from the list created with listfiles.py
        to open it and make a sorting of the lines present in it.
    Args:
        file (str): file that contains sequences to be sorted

    Raises:
        NotAFastaFile: Invalid file extension. The function only runs with FASTA files.
        Exception: Unexpected error encountered when reading the file

    Returns:
        Array: two-dimensional array of strings representing the species names and the aligned sequences.
    """
    # Checking file extension
    if not file.endswith('.afa'):
        raise BioSeekerExceptions.NotAFastaFile()

    try:
        # Creating auxiliary lists and variables
        text_list = []
        arr_x = []
        arr_y = []
        auxvar = True

        #   Storing every line of the file at "text_list"
        print(f"Extracting alignments from file {file}")
        open_file = open(os.path.expanduser(file), "r")
        for line in open_file:
            stripped_line = line.strip()
            line_list = stripped_line.split()
            text_list.append(line_list)
        print(f"Extraction from file {file} has been successful.")

        #   Sorting línes through the boolean "auxvar"
        #   This boolean acts as some sort of switch, and sends every line to one list or another
        #   The species will go to arr_x and the sequence goes to arr_y
        print("Creating alignment array...")
        for k in text_list:
            if auxvar:
                arr_x.append(k)
                auxvar = False
            else:
                arr_y.append(k)
                auxvar = True

        #   Creating array from the previous lists
        arr_xy = []
        for i in range(len(arr_x)):
            arr_xy.append([arr_x[i], arr_y[i]])
        print(f"The creation of the array from file {file} has been successful.")
        
        return arr_xy
    except Exception:
        print(f"The creation of the array from file {file} was not successful.")
        return None
