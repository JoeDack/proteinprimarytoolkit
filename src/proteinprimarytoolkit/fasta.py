from structures import *
from utils import *

def get_polypeptides_from_fasta(path: str) -> list[Polypeptide]:
    if not isinstance(path, str):
        raise TypeError(f"path must be str, not {type(path)}")

    with open(path) as fasta_file:
        contents = fasta_file.read().splitlines()
    
    enumerate_contents = enumerate(contents)

    # Headers always start with >
    headers = [i for i, l in enumerate_contents if l[0] == ">"]
    # Ignore headers and blank lines
    sequences = [i for i, l in enumerate_contents if i not in headers and l != ""]

    return [Polypeptide(s) for s in sequences]