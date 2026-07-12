from __future__ import annotations
from re import finditer
from collections import Counter
from dataclasses import dataclass, field
from typing import Any, Iterator

from amino_acids import *
from utils import *


@dataclass(frozen=True, slots=True)
class Ligand:
    component_id: str
    name: str
    smiles: str | None

@dataclass(frozen=True, slots=True)
class Polypeptide:
    sequence: str = ""

    def __post_init__(self) -> None:
        if not self.sequence:
            raise ValueError("No sequence given")
        
        invalid_characters = set(self.sequence) - VALID_CODES_SET
        if invalid_characters:
            raise ValueError(f"Invalid characters present in sequence: {invalid_characters}")
        
        counts = Counter(self.sequence)

        #.__setattr__ used because Polypeptide is immutable
        object.__setattr__(self, "amino_acids", [AMINO_ACID_CODES[a] for a in self.sequence])
        object.__setattr__(self, "frequencies", {a: counts.get(a, 0) for a in VALID_CODES}) # Frequency of all amino acids
        object.__setattr__(self, "abundances", {a: counts.get(a, 0) / len(self.sequence) for a in VALID_CODES}) # Proportion of all amino acids
        object.__setattr__(self, "unique_acids", list(dict.fromkeys(self.sequence))) # Removes duplicates but retains order
        object.__setattr__(self, "acidic_acids", {a: counts.get(a.code, 0) for a in ACIDIC_ACIDS})
        object.__setattr__(self, "basic_acids", {a: counts.get(a.code, 0) for a in BASIC_ACIDS})

    @staticmethod
    def _validate_sequence(sequence: str | Polypeptide | AminoAcid) -> None:
        if isinstance(sequence, Polypeptide):
            sequence = sequence.sequence
        elif isinstance(sequence, AminoAcid):
            sequence = sequence.code

        invalid = [(i, a)for i, a in enumerate(sequence) if a not in VALID_CODES_SET]

        if invalid:
            if len(sequence) == 1:
                raise InvalidSequenceError(f"{sequence} is not a valid amino acid")
            raise InvalidSequenceError(f"Invalid amino acids at following positions: {invalid}")

    def __repr__(self) -> str:
        return self.sequence
    
    def __str__(self) -> str:
        return self.sequence
    
    def __bool__(self) -> bool:
        return bool(self.sequence)
    
    def __len__(self) -> int:
        return len(self.sequence)
    
    def __getitem__(self, key: int | slice) -> Polypeptide:
        return Polypeptide(self.sequence[key])
    
    def __contains__(self, item: str | AminoAcid) -> bool:
        if isinstance(item, AminoAcid):
            return item in self.amino_acids
        return item in self.sequence
    
    def __add__(self, other: str | Polypeptide):
        if isinstance(other, Polypeptide):
            other = other.sequence
        return Polypeptide(self.sequence + other)
    
    def __radd__(self, other: str | Polypeptide):
        if isinstance(other, Polypeptide):
            other = other.sequence
        return Polypeptide(other + self.sequence)
    
    def __iter__(self) -> Iterator[str]:
        return iter(self.sequence)
    
    def count(self, acid: str | AminoAcid) -> int:
        """
        Returns the number of a specific amino acid in the polypeptide sequence
        """
        if not (isinstance(acid, (str, AminoAcid))):
            raise TypeError(f"Acid must be type str or AminoAcid, not {type(acid)}")
        
        if isinstance(acid, AminoAcid):
            acid = acid.code
        
        self._validate_sequence(acid)
        
        return self.frequencies[acid]
    
    def abundance(self, acid: str | AminoAcid) -> float:
        """
        Returns the proportional abundance of a specific amino acid in the polypeptide sequence
        """
        if not (isinstance(acid, (str, AminoAcid))):
            raise TypeError(f"Acid must be str or AminoAcid, not {type(acid)}")
        
        if isinstance(acid, AminoAcid):
            acid = acid.code
        
        self._validate_sequence(acid)
        
        return self.abundances[acid]
    

    def substitution(self, position: int, new_acid: str | AminoAcid) -> Polypeptide:
        """
        Returns the amino acid sequence with an amino acid at a given position for another amino acid
        """
        if not (isinstance(new_acid, (str, AminoAcid))):
            raise TypeError(f"Acid must be str or AminoAcid, not {type(new_acid)}")
        
        if isinstance(new_acid, AminoAcid):
            new_acid = new_acid.code
        
        self._validate_sequence(new_acid)
        
        return Polypeptide(self.sequence[:position] + new_acid + self.sequence[position + 1:])
    
    def addition(self, position: int, new_acid: str | AminoAcid) -> Polypeptide:
        """
        Returns the amino acid sequence with an amino acid added at a given position
        """
        if not (isinstance(new_acid, (str, AminoAcid))):
            raise TypeError(f"Acid must be str or AminoAcid, not {type(new_acid)}")
        
        if isinstance(new_acid, AminoAcid):
            new_acid = new_acid.code

        self._validate_sequence(new_acid)
        
        return Polypeptide(self.sequence[:position] + new_acid + self.sequence[position:])
    
    def deletion(self, position: int) -> Polypeptide:
        """
        Returns the amino acid sequence with an amino acid at a given position removed
        """
        if not isinstance(position, int):
            raise TypeError(f"Position must be int, not {type(position)}")
        
        return Polypeptide(self.sequence[:position] + self.sequence[position + 1:])

    def molecular_weight(self) -> float:
        """
        Molecular weight of a polypeptide, adjusted for water loss during formation of peptide bonds.
        For approximations, consider 110.0 as an average amino acid weight.
        """
        return sum(a.weight for a in self.amino_acids) - 18.01528 * (len(self) - 1)
    
    def net_charge(self, ph: float | int) -> float:
        """
        Net charge of a protein at a given pH using the Henderson-Hasselback equation. 
        Takes the N-terminus and C-terminus into account.
        Assumes standard pKa values, doesn't take into account environmental effects.
        """
        if not isinstance(ph, (int, float)):
            raise TypeError(f"pH must be int or float, not {type(ph)}")
        acids = self.amino_acids

        charge = 0.0
        for acid, count in self.acidic_acids.items():
            charge -= count * ((1 + 10 ** (acid.side_chain_pka - ph)) ** -1)
        for acid, count in self.basic_acids.items():
            charge += count * ((1 + 10 ** (ph - acid.side_chain_pka)) ** -1)

        # Cysteine and tyrosine are still ionisable, despite being considered polar/uncharged
        # Cysteine
        charge -= self.count("C") * ((1 + 10 ** (CYSTEINE.side_chain_pka - ph)) ** -1)
        # Tyrosine
        charge -= self.count("Y") * ((1 + 10 ** (TYROSINE.side_chain_pka - ph)) ** -1)

        # N-terminus
        charge += (1 + 10 ** (ph - acids[0].amino_pka)) ** -1
        # C-terminus
        charge -= (1 + 10 ** (acids[-1].carboxyl_pka - ph)) ** -1

        return charge

    def isoelectric_point(self, lower: float | int = 0.0, upper: float | int = 14.0, tolerance: float = 1e-6) -> float:
        """
        Calculates a polypeptide's approximate isoelectric point using a binary search algorithm, assuming the polypeptide is monoisotopic.
        Default starting values are 0.0 and 14.0, default tolerance is 10^-6.
        """
        if not isinstance(lower, (int, float)):
            raise TypeError(f"Lower bound must be int or float, not {type(lower)}")
        if not isinstance(upper, (int, float)):
            raise TypeError(f"Upper bound must be int or float, not {type(upper)}")
        
        if not isinstance(tolerance, float):
            raise TypeError(f"Tolerance must be float, not {type(tolerance)}")
        elif tolerance <= 0:
            raise ValueError(f"Tolerance must be positive, not 0 or negative")
        
        midpoint = (upper + lower) / 2

        while upper - lower > tolerance:
            midpoint = (upper + lower) / 2

            if self.net_charge(midpoint) > 0:
                lower = midpoint
            else:
                upper = midpoint
        
        return midpoint
    
    def extinction_coefficient(self) -> int:
        """
        Extinction coefficient of a protein.
        Dependent on abundance of tryptophan, tyrosine, cysteine.
        """
        return 5500*self.count("W") + 1490*self.count("Y") + 125*self.count("C")
    
    def gravy(self) -> float:
        """
        GRAVY (Grand average of hydropathy_index) of a polypeptide
        """
        return sum(a.hydropathy_index for a in self.amino_acids) / len(self)
    
    def aliphatic_index(self) -> float:
        """
        Aliphatic index of a polypeptide.
        Dependent on the abundance of alanine, valine, leucine and isoleucine
        """
        return self.abundances["A"] + 2.9*self.abundances["V"] + 3.9*(self.abundances["I"] + self.abundances["L"])
    
    def aromatic_index(self) -> float:
        """
        Proportion of aromatic residues in a polypeptide
        """
        return (self.frequencies["F"] + self.frequencies["W"] + self.frequencies["Y"]) / len(self)
    
    def average_helix_propensity(self) -> float:
        """
        Average a-helix propensity of amino acids within a polypeptide
        """
        return sum(a.helix_propensity for a in self.amino_acids) / len(self)
    
    def average_sheet_propensity(self) -> float:
        """
        Average b-sheet propensity of amino acids within a polypeptide
        """
        return sum(a.sheet_propensity for a in self.amino_acids) / len(self)
    
    def average_turn_propensity(self) -> float:
        """
        Average turn propensity of amino acids within a polypeptide
        """
        return sum(a.turn_propensity for a in self.amino_acids) / len(self)
    
    def motif_search(self, motif: str) -> list[int]:
        """
        Given a specific motif, returns list of indices of instances of the motif
        """
        return [m.start() for m in finditer(f"(?={motif})", self.sequence)]
    
    def _possible_structures(self, attr_name: str, min_length: int = 6, threshold: int = 100, proportion: float = 0.667, break_on_proline_or_glycine: bool = False) -> list[tuple[int, int]]:
        """
        Given a specific secondary structure, searches for potential instances of the structure using amino acid propensities.
        Very simplified, based on the Chou-Fasman method
        """
        if not isinstance(min_length, int):
            raise TypeError(f"minLength must be int, not {type(min_length)}")
        elif min_length <= 0:
            raise ValueError(f"minLength must be positive, not 0 or negative")
        structures = []
        full_length = len(self)
        for start in range(full_length - min_length + 1):

            current = None
            window_length = min_length
            end = start + window_length

            region = self.amino_acids[start:end]

            propensity = sum(getattr(acid, attr_name) >= threshold for acid in region)
            while (propensity > proportion * window_length):

                current = (start, end)

                window_length += 1
                end += 1

                if end >= len(self):
                    break

                region = self.amino_acids[start:end]

                if break_on_proline_or_glycine:
                    last4 = region[-4:]

                    if any(a.code in {"P", "G"} for a in last4):
                        break
                
                propensity = sum(getattr(acid, attr_name) >= threshold for acid in region)

            if current:
                structures.append(current)

        return structures
    
    def possible_helices(self, min_length: int = 6) -> list[tuple[int, int]]:
        """
        Searches for potential helices, default minimum length is 6 residues. Returns list of start and stop positions
        """
        return self._possible_structures(attr_name="helix_propensity", min_length=min_length, break_on_proline_or_glycine=True)
    
    def possible_sheets(self, min_length: int = 6) -> list[tuple[int, int]]:
        """
        Searches for possible sheets, default minimum length is 6 residues. Returns list of start and stop positions
        """
        return self._possible_structures(attr_name="sheet_propensity", min_length=min_length, break_on_proline_or_glycine=True)
    
    def possible_turns(self, min_length: int = 6) -> list[tuple[int, int]]:
        """
        Searches for possible turns, default minimum length is 6 residues. Returns list of start and stop positions
        """
        return self._possible_structures(attr_name="turn_propensity", min_length=min_length)

@dataclass(frozen=True, slots=True)
class Protein:
    classification: dict[str, Any] | None
    sequences: list[Polypeptide] = field(default_factory=list)
    ligands: list[Ligand] = field(default_factory=list)
