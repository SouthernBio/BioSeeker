#   1   -   Function to find and store FASTA files' names
#           The function takes a directory as argument, which contains the alignment files
def list_of_files(working_directory):
    import os
    directory = working_directory
    contents = os.listdir(directory)

    #   Creating a list of FASTA files, which will be processed by another function
    sequences = []
    print("Program initiated. Creating list of FASTA files in directory " + working_directory)
    for file in contents:
        if os.path.isfile(os.path.join(directory, file)) and file.endswith('.afa'):
            sequences.append(file)
    print("List created successfully.")

    return sequences