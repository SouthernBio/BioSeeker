class NotADataframe(Exception):
    """Custom exception for when a function receives a file that is not a CSV dataframe"""
    message = "Error: the function must receive a dataframe in CSV format."

    def __init__(self):
        super().__init__(self.message)


class InvalidReadingFrame(Exception):
    """Custom exception for when a function receives a reading frame out of range"""
    message = "Error: 'reading' can only take values in {0, 1, 2}"

    def __init__(self):
        super().__init__(self.message)


class AlignmentFilesMissing(Exception):
    """Custom exception for when BioSeeker cannot find any FASTA files in the current directory"""
    message = 'Error: BioSeeker could not find any alignment files on the working directory'

    def __init__(self):
        super().__init__(self.message)


class NotAFastaFile(Exception):
    """Custom exception for when extract_sequences() receives an invalid input"""
    message = 'Error: invalid file extension, it must end with ".afa". Have you modified "listfiles.py"?'

    def __init__(self):
        super().__init__(self.message)
