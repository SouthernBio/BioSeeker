#   1   -   Function to find and store FASTA files' names
#           The function takes a directory as argument, which contains the alignment files
import os

def list_of_files(working_directory: str):
    """Function to find and store FASTA files' names
        The function takes a directory as argument, which contains the alignment files

    Args:
        working_directory (str): The working directory to look for FASTA files

    Returns:
        List: A list with file names
    """
    try:
        contents = os.listdir(working_directory)

        #   Creating a list of FASTA files, which will be processed by another function
        sequences = []
        print("Program initiated. Creating list of FASTA files in directory " + working_directory)
        for file in contents:
            if os.path.isfile(os.path.join(working_directory, file)) and file.endswith('.afa'):
                sequences.append(file)
        print("List created successfully.")

        return sequences
    except:
        print("An error occurred while trying to create the list of FASTA files.")
        return None


    