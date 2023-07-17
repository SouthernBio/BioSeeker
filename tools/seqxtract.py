import os

def extract_sequences(file: str):
    """Function to extract sequences from FASTA files
        This function takes as an argument a file from the list created with listfiles.py
        to open it and make a sorting of the lines present in it.
    Args:
        file (str): file that contains sequences to be sorted

    Raises:
        Exception: Invalid file extension. The function only runs with FASTA files.
        Exception: Unexpected error encountered when reading the file

    Returns:
        Array: two dimensional array of strings representing the species names and the aligned sequences.
    """
    # Checking file extension
    if file.endswith('.afa') == False:
        raise Exception('Error: invalid file extension, it must end with ".afa". Have you modified "listfiles.py"?')
    
    try:
        # Creating auxiliary lists and variables
        text_list = []
        arrX = []
        arrY = []
        auxvar = True

        #   Storing every line of the file at "text_list"
        print(f"Extracting alignments from file {file}")
        open_file = open(os.path.expanduser(file), "r")
        for line in open_file:
            stripped_line = line.strip()
            line_list = stripped_line.split()
            text_list.append(line_list)
        print(f"Extraction from file {file} has been successful.")

        #   Sorting líneas through the boolean "auxvar"
        #   This boolean acts as some sort of switch, and sends every line to one list or another
        #   The species will go to arrX and the sequence goes to arrY
        print("Creating alignment array...")
        for k in text_list:
            if auxvar == True:
                arrX.append(k)
                auxvar = False    
            else:
                arrY.append(k)
                auxvar = True

        #   Creating array from the previous lists
        arrXY = []
        for i in range(len(arrX)):
            arrXY.append([arrX[i], arrY[i]])
        print(f"The creation of the array from file {file} has been successful.")
        
        return arrXY
    except:
        print(f"The creation of the array from file {file} was not successful.")
        return None