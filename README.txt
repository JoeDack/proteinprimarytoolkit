**README**
This is a (somewhat basic) bioinformatics/protein analysis toolkit. 
I am currently a biology student at the University of York and doing this for fun - it's not meant to be a serious project. 
However, it's now getting large enough I think it's worth uploading to GitHub and PyPI, even if it's not useful enough to get used in actual research.

**CITATION**
If using in a paper , please cite as follows:
Dack, J. (2026). proteinprimarytoolkit (Version 1.0.8) [Computer Software]. https://github.com/JoeDack/proteinprimarytoolkit

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

**OTHER FILES**
__init__
'tis the init file.

utils
Custom exceptions and API utils - might move custom exceptions to its own module later.


**CHANGELOG**
v1.0.0 - Initial release
Yay, it's the first version :P

v1.0.1 - Bugfix
Fixed bug in pairwise_sequence_alignment() where an unintended negative index may occur

v1.0.2 - Bugfix 2, electric boogaloo
Imported annotations from __futures__ in structures.py to avoid errors with type annotations within Polypeptide in older Python versions

v1.0.3 - Yet again, thine bugs hath been fixed
Fixed logic bug in pairwise_sequence_alignment()

v1.0.4 - Error handling
Better (mostly API-related) error handling

v1.0.5 - Unreachable Code
Removed a bit of unreacheable code in pairwise_sequence_alignment(). The code would run if 2 variables were both less than 0, which was impossible.

v1.0.6 - A bugfix, a bugfix, my kingdom for a bugfix
Fixed bug in multi_sequence_alignment() where code would loop forever if a status code other than "FINISHED" was returned by the Clustal Omega API

v1.0.7 - Custom errors
Added a couple new custom errors - APIError and InvalidPDBIDError.
Also fixed the changelog - the previous update was listed as v1.0.4 instead of v1.0.6

v1.0.8 - API changes
More robust API helpers and better API error handling
Added default timeout for all functions which use APIs (currently 30s)
Removed explicit session type checks from API helpers