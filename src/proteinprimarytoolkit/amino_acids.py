from dataclasses import dataclass
from enum import Enum, auto # For categories


@dataclass(frozen=True, slots=True)
class AminoAcid:
    code: str # 1 letter code
    weight: float # Monoisotopic weight
    carboxyl_pka: float # -COOH pKa
    amino_pka: float # -NH2 pKa
    side_chain_pka: float | None # -R pKa
    isoelectric_point: float # pH where charge = 0
    physiological_charge: float # Net charge at physiological pH (7.4)
    hydropathy_index: float # hydropathy_index index, used for GRAVY calculation
    helix_propensity: int # Chou-Fasman propensity for helix formation
    sheet_propensity: int # Chou-Fasman propensity for sheet formation
    turn_propensity: int # Chou-Fasman Propensity for turn formation
    category: str # hydrophobic, polar, acidic, basic, or aromatic. Glycine and proline are considered hydrophobic here

    def most_likely_secondary_structure(self) -> str:
        max_propensity = max(self.helix_propensity, self.sheet_propensity, self.turn_propensity)

        if max_propensity == self.helix_propensity:
            return "helix"
        elif max_propensity == self.sheet_propensity:
            return "sheet"
        return "turn"
    

class AminoAcidCategory(Enum):
    HYDROPHOBIC = auto()
    POLAR = auto()
    BASIC = auto()
    ACIDIC = auto()
    AROMATIC = auto()

GLYCINE = AminoAcid(code="G",
                    weight=57.02,
                    carboxyl_pka=2.32,
                    amino_pka=9.60,
                    side_chain_pka=None,
                    isoelectric_point=5.97,
                    physiological_charge=0.0,
                    hydropathy_index=-0.4,
                    helix_propensity=57,
                    sheet_propensity=75,
                    turn_propensity=156,
                    category=AminoAcidCategory.HYDROPHOBIC
                    )

ALANINE = AminoAcid(code="A",
                    weight=89.09,
                    carboxyl_pka=2.34,
                    amino_pka=9.69,
                    side_chain_pka=None,
                    isoelectric_point=6.0,
                    physiological_charge = 0.0,
                    hydropathy_index=1.8,
                    helix_propensity=142,
                    sheet_propensity=83,
                    turn_propensity=66,
                    category=AminoAcidCategory.HYDROPHOBIC
                    )

ARGININE = AminoAcid(code="R",
                     weight=156.10,
                     carboxyl_pka=2.17,
                     amino_pka=9.04,
                     side_chain_pka=12.48,
                     isoelectric_point=10.78,
                     physiological_charge=1.0,
                     hydropathy_index=-4.5,
                     helix_propensity=98,
                     sheet_propensity=83,
                     turn_propensity=66,
                     category=AminoAcidCategory.BASIC
                     )

ASPARAGINE = AminoAcid(code="N",
                       weight=114.04,
                       carboxyl_pka=2.02,
                       amino_pka=8.80,
                       side_chain_pka=None,
                       isoelectric_point=5.41,
                       hydropathy_index=-3.5,
                       physiological_charge=0.0,
                       helix_propensity=67,
                       sheet_propensity=89,
                       turn_propensity=156,
                       category=AminoAcidCategory.POLAR
                       )

ASPARTIC_ACID = AminoAcid(code="D",
                          weight=115.03,
                          carboxyl_pka=1.88,
                          amino_pka=9.60,
                          side_chain_pka=3.65,
                          isoelectric_point=2.77,
                          physiological_charge=-1.0,
                          hydropathy_index=-3.5,
                          helix_propensity=101,
                          sheet_propensity=54,
                          turn_propensity=146,
                          category=AminoAcidCategory.ACIDIC
                          )

CYSTEINE = AminoAcid(code="C",
                     weight=103.01,
                     carboxyl_pka=1.96,
                     amino_pka=8.18,
                     side_chain_pka=8.3,
                     isoelectric_point=5.07,
                     physiological_charge=0.0,
                     hydropathy_index=2.5,
                     helix_propensity=70,
                     sheet_propensity=119,
                     turn_propensity=119,
                     category=AminoAcidCategory.POLAR
                     )

GLUTAMIC_ACID = AminoAcid(code="E",
                          weight=129.04,
                          carboxyl_pka=2.19,
                          amino_pka=9.67,
                          side_chain_pka=4.25,
                          isoelectric_point=3.22,
                          physiological_charge=-1.0,
                          hydropathy_index=-3.5,
                          helix_propensity=151,
                          sheet_propensity=37,
                          turn_propensity=74,
                          category=AminoAcidCategory.ACIDIC
                          )

GLUTAMINE = AminoAcid(code="Q",
                      weight=128.06,
                      carboxyl_pka=2.17,
                      amino_pka=9.13,
                      side_chain_pka=None,
                      isoelectric_point=5.65,
                      physiological_charge=0.0,
                      hydropathy_index=-3.5,
                      helix_propensity=111,
                      sheet_propensity=110,
                      turn_propensity=98,
                      category=AminoAcidCategory.POLAR
                      )

HISTIDINE = AminoAcid(code="H",
                      weight=137.06,
                      carboxyl_pka=1.82,
                      amino_pka=9.17,
                      side_chain_pka=6.00,
                      isoelectric_point=7.59,
                      physiological_charge=0.1, # Histidine is only partially charged at physiological pH
                      hydropathy_index=-3.2,
                      helix_propensity=100,
                      sheet_propensity=87,
                      turn_propensity=95,
                      category=AminoAcidCategory.BASIC
                      )

ISOLEUCINE = AminoAcid(code="I",
                       weight=113.08,
                       carboxyl_pka=2.36,
                       amino_pka=9.60,
                       side_chain_pka=None,
                       isoelectric_point=6.02,
                       hydropathy_index=4.5,
                       physiological_charge=0.0,
                       helix_propensity=121,
                       sheet_propensity=130,
                       turn_propensity=59,
                       category=AminoAcidCategory.HYDROPHOBIC
                       )

LEUCINE = AminoAcid(code="L",
                    weight=113.08,
                    carboxyl_pka=2.36,
                    amino_pka=9.60,
                    side_chain_pka=None,
                    isoelectric_point=5.98,
                    hydropathy_index=3.8,
                    physiological_charge=0.0,
                    helix_propensity=121,
                    sheet_propensity=130,
                    turn_propensity=59,
                    category=AminoAcidCategory.HYDROPHOBIC
                    )

LYSINE = AminoAcid(code="K",
                   weight=129.09,
                   carboxyl_pka=2.18,
                   amino_pka=8.95,
                   side_chain_pka=10.53,
                   isoelectric_point=9.74,
                   physiological_charge=1.0,
                   hydropathy_index=-3.9,
                   helix_propensity=114,
                   sheet_propensity=74,
                   turn_propensity=60,
                   category=AminoAcidCategory.POLAR
                   )

METHIONINE = AminoAcid(code="M",
                       weight=131.04,
                       carboxyl_pka=2.28,
                       amino_pka=9.21,
                       side_chain_pka=None,
                       isoelectric_point=5.74,
                       physiological_charge=0.0,
                       hydropathy_index=1.9,
                       helix_propensity=145,
                       sheet_propensity=105,
                       turn_propensity=60,
                       category=AminoAcidCategory.HYDROPHOBIC
                       )

PHENYLALANINE = AminoAcid(code="F",
                          weight=147.07,
                          carboxyl_pka=1.83,
                          amino_pka=9.13,
                          side_chain_pka=None,
                          isoelectric_point=5.48,
                          physiological_charge=0.0,
                          hydropathy_index=2.8,
                          helix_propensity=113,
                          sheet_propensity=138,
                          turn_propensity=60,
                          category=AminoAcidCategory.AROMATIC)

PROLINE = AminoAcid(code="P",
                    weight=97.05,
                    carboxyl_pka=1.99,
                    amino_pka=10.60,
                    side_chain_pka=None,
                    isoelectric_point=6.30,
                    physiological_charge=0.0,
                    hydropathy_index=-1.60,
                    helix_propensity=57,
                    sheet_propensity=75,
                    turn_propensity=143,
                    category=AminoAcidCategory.HYDROPHOBIC)

SERINE = AminoAcid(code="S",
                   weight=87.03,
                   carboxyl_pka=2.21,
                   amino_pka=9.15,
                   side_chain_pka=None,
                   isoelectric_point=5.68,
                   physiological_charge=0.0,
                   hydropathy_index=-0.8,
                   helix_propensity=77,
                   sheet_propensity=75,
                   turn_propensity=143,
                   category=AminoAcidCategory.POLAR)

THREONINE = AminoAcid(code="T",
                      weight=101.05,
                      carboxyl_pka=2.09,
                      amino_pka=9.10,
                      side_chain_pka=None,
                      isoelectric_point=5.60,
                      physiological_charge=0.0,
                      hydropathy_index=-0.7,
                      helix_propensity=83,
                      sheet_propensity=119,
                      turn_propensity=96,
                      category=AminoAcidCategory.POLAR)

TRYPTOPHAN=AminoAcid(code="W",
                     weight=186.08,
                     carboxyl_pka=2.83,
                     amino_pka=9.39,
                     side_chain_pka=None,
                     isoelectric_point=5.89,
                     physiological_charge=0.0,
                     hydropathy_index=-0.9,
                     helix_propensity=108,
                     sheet_propensity=137,
                     turn_propensity=96,
                     category=AminoAcidCategory.AROMATIC)

TYROSINE = AminoAcid(code="Y",
                     weight=163.06,
                     carboxyl_pka=2.20,
                     amino_pka=9.11,
                     side_chain_pka=10.1,
                     isoelectric_point=5.66,
                     physiological_charge=0.00,
                     hydropathy_index=-1.3,
                     helix_propensity=69,
                     sheet_propensity=147,
                     turn_propensity=114,
                     category=AminoAcidCategory.AROMATIC)

VALINE = AminoAcid(code="V",
                   weight=99.07,
                   carboxyl_pka=2.32,
                   amino_pka=9.62,
                   side_chain_pka=None,
                   isoelectric_point=5.96,
                   physiological_charge=0.0,
                   hydropathy_index=4.2,
                   helix_propensity=106,
                   sheet_propensity=170,
                   turn_propensity=50,
                   category=AminoAcidCategory.HYDROPHOBIC)

AMINO_ACIDS = (ALANINE, ARGININE, ASPARAGINE, ASPARTIC_ACID, CYSTEINE, GLUTAMIC_ACID, GLUTAMINE, GLYCINE, HISTIDINE, LEUCINE, ISOLEUCINE, LYSINE, METHIONINE, PHENYLALANINE, PROLINE, SERINE, THREONINE, TRYPTOPHAN, TYROSINE, VALINE)
ACIDIC_ACIDS = tuple(a for a in AMINO_ACIDS if a.category == AminoAcidCategory.ACIDIC)
BASIC_ACIDS = tuple(a for a in AMINO_ACIDS if a.category == AminoAcidCategory.BASIC)
POLAR_ACIDS = tuple(a for a in AMINO_ACIDS if a.category == AminoAcidCategory.POLAR)
HYDROPHOBIC_ACIDS = tuple(a for a in AMINO_ACIDS if a.category == AminoAcidCategory.HYDROPHOBIC)
AROMATIC_ACIDS = tuple(a for a in AMINO_ACIDS if a.category == AminoAcidCategory.AROMATIC)

AMINO_ACID_CODES = {a.code: a for a in AMINO_ACIDS}
VALID_CODES = tuple(AMINO_ACID_CODES.keys())
VALID_CODES_SET = set(VALID_CODES)

AMBIGUOUS_CODES = {"X": AMINO_ACIDS,
                   "B": (ASPARTIC_ACID, ASPARAGINE),
                   "Z": (GLUTAMIC_ACID, GLUTAMINE),
                   "J": (LEUCINE, ISOLEUCINE),
                   }