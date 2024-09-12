"""
This module defines the Protein class, which handles loading a PDB file,
extracting structural information, and computing various properties of a protein.

Classes:
- Protein: Handles the loading, parsing, and property computation of a protein from a PDB file.
"""

import numpy as np
from Bio import PDB
from Bio.PDB.DSSP import DSSP
from geometry_and_axis import Point, Vector

class Protein:
    """
    Handles the loading, parsing, and property computation of a protein from a PDB file.

    Attributes:
    - pdb_file (str): The path to the PDB file.
    - structure (Structure): The loaded protein structure.
    - chain (Chain): The chain of residues.
    - residues (list): List of residues in the chain.
    - C_alpha_coords (list): List of coordinates for C-alpha atoms.
    - solvent_accessibility (dict): Dictionary of solvent accessibility for each residue.
    - hydrophobicity (dict): Dictionary of hydrophobicity values for each residue.
    - center_of_mass (numpy.ndarray): The center of mass of the protein.
    - best_positions (list): List of best positions found (not used in current code).
    """

    def __init__(self, pdb_file):
        self.pdb_file = pdb_file
        self.structure = self._load_structure()
        self.chain = self.structure[0]['A']
        self.residues = list(self.chain)
        self.C_alpha_coords = self._extract_C_alpha_coords()
        self.solvent_accessibility = self._compute_solvent_accessibility()
        self.hydrophobicity = self._compute_hydrophobicity()
        self.center_of_mass = self._compute_center_of_mass()
        self.best_positions = []

    def _load_structure(self):
        """
        Loads the protein structure from the PDB file.

        Returns:
        - Structure: The loaded protein structure.
        """
        parser = PDB.PDBParser(QUIET=True)
        return parser.get_structure('protein', self.pdb_file)

    def _extract_C_alpha_coords(self):
        """
        Extracts the coordinates of C-alpha atoms from the residues.

        Returns:
        - list: List of (x, y, z) coordinates for C-alpha atoms.
        """
        return [res['CA'].coord for res in self.residues if 'CA' in res]

    def _compute_solvent_accessibility(self):
        """
        Computes the solvent accessibility of residues using DSSP.

        Returns:
        - dict: Dictionary with residue IDs as keys and relative ASA as values.
        """
        model = self.structure[0]
        dssp_obj = DSSP(model, self.pdb_file)
        
        solvent_accessibility = {}
        for item in dssp_obj:
            dssp_index, amino_acid, secondary_structure, relative_asa, *rest = item
            res_id = f"{amino_acid}{dssp_index}"
            solvent_accessibility[res_id] = relative_asa
        
        return solvent_accessibility

    def _compute_hydrophobicity(self):
        """
        Computes the hydrophobicity of residues.

        Returns:
        - dict: Dictionary with residue IDs as keys and hydrophobicity (1 or 0) as values.
        """
        hydrophobic_residues = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP'}
        return {f"{residue.get_resname()}{residue.get_id()}": 1 if residue.get_resname() in hydrophobic_residues else 0
                for residue in self.residues}

    def _compute_center_of_mass(self):
        """
        Computes the center of mass of the protein.

        Returns:
        - numpy.ndarray: The (x, y, z) coordinates of the center of mass.
        """
        coords = np.array(self.C_alpha_coords)
        return np.mean(coords, axis=0)
