"""BioSeeker: DNA Sequence Conservation Analyzer."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("bioseeker")
except PackageNotFoundError:
    # package is not installed
    __version__ = "unknown"
