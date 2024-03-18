def intro_message() -> None:
    """A message that displays the current version of BioSeeker and the author."""

    message = """
    BioSeeker v0.0.1-beta
    Author: SouthernBio, Open Source Bioinformatics Collective from Argentina
    Official repository: https://github.com/SouthernBio/BioSeeker
    License: GNU General Public License v3.0, see LICENSE
    Wiki: https://github.com/SouthernBio/BioSeeker/wiki
    ---------------------------------------------------------------
    
    Press any key to start the analysis on the current directory...
    """
    print(message)

    return None
