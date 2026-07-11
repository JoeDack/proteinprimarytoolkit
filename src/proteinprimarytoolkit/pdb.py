from typing import Any # For type hints
from concurrent.futures import ThreadPoolExecutor # For concurrent API requests

from structures import *
from utils import *

BASE_PDB_API_URL = r"https://data.rcsb.org/rest/v1/core"


def get_polymer_sequence(pdb_id: str, polymer_id: str, timeout: int | float) -> str | None:
    """
    Retrieve a sequence for a polymer entity.
    Returns None if the entity is not a polypeptide.
    """
    session = get_session()

    entity_url = f"{BASE_PDB_API_URL}/polymer_entity/{pdb_id}/{polymer_id}"
    polypeptide = get_json_from_api(session, entity_url, timeout)

    if polypeptide["entity_poly"]["type"].startswith("polypeptide"):
        return polypeptide["entity_poly"]["pdbx_seq_one_letter_code_can"]

    return None


def get_ligand(pdb_id: str, ligand_id: str, timeout: int | float) -> Ligand:
    """
    Retrieve component ID, name and SMILES for a non-polymer entry. Returns Ligand object
    """
    session = get_session()

    ligand_url = f"{BASE_PDB_API_URL}/nonpolymer_entity/{pdb_id}/{ligand_id}"
    ligand = get_json_from_api(session, ligand_url, timeout)

    component_id = ligand["pdbx_entity_nonpoly"]["comp_id"]

    chemical_url = f"{BASE_PDB_API_URL}/chemcomp/{component_id}"
    chemical = get_json_from_api(session, chemical_url, timeout)

    name = chemical["chem_comp"]["name"]
    descriptor = chemical.get("rcsb_chem_comp_descriptor", {})
    smiles = descriptor.get("SMILES")

    return Ligand(component_id, name, smiles)


def get_entry(pdb_id: str, timeout: int) -> dict[str, Any]:
    session = get_session()
    entry_url = f"{BASE_PDB_API_URL}/entry/{pdb_id}"
    return get_json_from_api(session, entry_url, timeout)

def get_polypeptides_from_pdb(pdb_id: str, timeout: int = 10) -> list[Polypeptide]:
    """
    Retrieve sequences from RCSB PDB, returns list of Polypeptide instances
    """
    if not (len(pdb_id) == 4 and pdb_id.isalnum()):
        raise ValueError("PDB ID must be alphanumeric and 4 characters long")

    pdb_id = pdb_id.upper()

    entry = get_entry(pdb_id, timeout)
    polymer_ids = entry["rcsb_entry_container_identifiers"]["polymer_entity_ids"]

    # Dynamic allocation of max workers
    workers = min(16, len(polymer_ids)) or 1

    # Concurrent requests, much faster for proteins with many polypeptides
    with ThreadPoolExecutor(max_workers=workers) as executor:
        sequences = [s for s in executor.map(lambda p: get_polymer_sequence(pdb_id, p, timeout), polymer_ids) if s is not None]

    # Remove duplicates while preserving order
    polypeptides = [Polypeptide(s) for s in sequences]
    return polypeptides

def get_protein_from_pdb(pdb_id: str, timeout: int = 10) -> Protein:
    """
    Retrieve sequence and ligand information using the RCSB PDB Data API.
    """
    if not (len(pdb_id) == 4 and pdb_id.isalnum()):
        raise ValueError("PDB ID must be alphanumeric and 4 characters long")

    pdb_id = pdb_id.upper()

    entry = get_entry(pdb_id, timeout)

    classification = entry.get("rcsb_polymer_entity_annotation")

    polymer_ids = entry["rcsb_entry_container_identifiers"]["polymer_entity_ids"]
    ligand_ids = entry["rcsb_entry_container_identifiers"].get("non_polymer_entity_ids", [])

    # Dynamic allocation of max workers
    workers = min(16, len(polymer_ids)) or 1

    # Concurrent requests, much faster for proteins with many polypeptides
    with ThreadPoolExecutor(max_workers=workers) as executor:
        sequences = [s for s in executor.map(lambda p: get_polymer_sequence(pdb_id, p, timeout), polymer_ids) if s is not None]
        ligands = list(executor.map(lambda l: get_ligand(pdb_id, l, timeout), ligand_ids))

    # Remove duplicates while preserving order
    sequences = list(dict.fromkeys(sequences))
    polypeptides = [Polypeptide(s) for s in sequences]
    ligands = list(dict.fromkeys(ligands))

    return Protein(classification, polypeptides, ligands)