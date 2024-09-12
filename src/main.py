"""
Main module for running the membrane finding analysis.

This script processes a PDB file to find the best membrane axis in a protein structure by analyzing hydrophobicity between planes.
It then visualizes the results in PyMOL.

Usage:
    python main.py <pdb_file>

Arguments:
    pdb_file (str): Path to the PDB file containing the protein structure.

Steps:
1. Load the protein structure from the given PDB file.
2. Find the best membrane axis by analyzing hydrophobicity between planes.
3. Print the best hydrophobicity and axis details.
4. Visualize the protein structure and planes in PyMOL.

Modules:
- Protein: Handles protein data and calculations.
- MembraneFinder: Finds the best membrane axis by analyzing hydrophobicity.
- show_in_pymol: Visualizes the results in PyMOL.
- Plane, Axis: Classes used for geometrical representation of planes and axes.
"""

import sys
from protein import Protein
from membrane_finder import MembraneFinder
from visualization import show_in_pymol
from geometry_and_axis import Plane, Axis

def main(pdb_file):
    """
    Main function to find the best membrane axis and visualize it.

    Args:
    - pdb_file (str): Path to the PDB file.

    Steps:
    1. Create a Protein object from the PDB file.
    2. Initialize a MembraneFinder with the Protein object.
    3. Find the best membrane axis and hydrophobicity.
    4. Print the best hydrophobicity and axis details.
    5. Visualize the protein structure and planes in PyMOL if a suitable axis is found.
    """
    protein = Protein(pdb_file)
    membrane_finder = MembraneFinder(protein)
    best_axis, best_hydrophobicity = membrane_finder.find_best_membrane()

    print(f"Best hydrophobicity: {best_hydrophobicity}")
    if best_axis:
        print("Best axis found:", best_axis)
        show_in_pymol(best_axis.plane1, best_axis.plane2, pdb_file)
    else:
        print("No suitable axis found.")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python main.py <pdb_file>")
    else:
        pdb_file = sys.argv[1]
        main(pdb_file)
