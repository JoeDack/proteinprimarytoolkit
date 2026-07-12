**README**

This is a (somewhat basic) bioinformatics/protein analysis toolkit. I am currently a biology student at the University of York and doing this for fun - it's not meant to be a serious project. However, it's now getting large enough I think it's worth uploading to GitHub and PyPI, even if it's not useful enough to get used in actual research.

**CITATION**
If using in a paper, please cite as follows:

Dack, J. (2026). proteinprimarytoolkit (Version 1.0.0) [Computer Software]. https://github.com/JoeDack/proteinprimarytoolkit

**MODULES:**
amino_acids
AminoAcid class and constants representing the 20 main proteinogenic amino acids

fasta
Imports FASTA files to list of Polypeptide instances

pdb
Uses the RCSB PDB API to import protein sequences and other associated data

sequence_alignment
Pairwise and multi-sequence alignment. Pairwise alignment currently uses the Needleman-Wunsch algorithm. Multi-sequence alignment uses the Clustal API

structures
Classes representing ligands, polypeptides, proteins. Polypeptide is functionally a wrapper around a string, Protein is functionally a wrapper around a list of Polypeptide objects.

OTHER FILES
__init__
'tis the init file.

utils
Custom exceptions and API utils - might move custom exceptions to its own module later.


CHANGELOG
v1.0.0 - Initial release
Yay, it's the first version :P

v1.0.1 - Bugfix
Fixed bug in pairwise_sequence_alignment() where an unintended negative index may occur

v1.0.2 - Bugfix
Imported annotations from __futures__ in structures.py to avoid errors with type annotations within Polypeptide in older Python versions
