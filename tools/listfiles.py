import os

class AlignmentFilesMissing(Exception):
    """Custom exception for when BioSeeker cannot find any FASTA files in the current directory"""
    def __init__(self, message):
        super().__init__(message)


def list_of_files(working_directory: str):
    """Function to find and store FASTA files' names
        The function takes a directory as argument, which contains the alignment files

    Args:
        working_directory (str): The working directory to look for FASTA files
    
    Raises:
        AlignmentFilesMissing: Custom exception for when BioSeeker cannot find any FASTA files in the current directory

    Returns:
        List: A list with file names
    """

    contents = os.listdir(working_directory)

    #   Creating a list of FASTA files, which will be processed by another function
    sequences = []
    print("Program initiated. Creating list of FASTA files in directory " + working_directory)
    for file in contents:
        if os.path.isfile(os.path.join(working_directory, file)) and file.endswith('.afa'):
            sequences.append(file)

    if len(sequences) == 0:
        raise AlignmentFilesMissing('Error: BioSeeker could not find any alignment files on the working directory')
    else:
        print("List created successfully.")

    return sequences

