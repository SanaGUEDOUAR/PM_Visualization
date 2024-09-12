"""
This module provides functions to visualize the protein structure and planes in PyMOL.

Functions:
- create_plane_object: Creates a CGO object for a plane.
- show_in_pymol: Displays the protein structure and planes in PyMOL.
"""

import os
import pymol
from pymol import cmd, cgo
import numpy as np
from geometry_and_axis import Plane

def create_plane_object(plane, color, scale=30):
    """
    Creates a CGO object representing a plane for visualization in PyMOL.

    Args:
    - plane (Plane): The Plane object representing the plane to be visualized.
    - color (str): The color of the plane. Supported colors are 'blue', 'green', 'red', 'yellow'.
    - scale (float): The scale of the plane visualization. Default is 30.

    Returns:
    - list: A list of CGO commands for creating the plane in PyMOL.
    """
    color_dict = {
        'blue': [0.0, 0.0, 1.0],
        'green': [0.0, 1.0, 0.0],
        'red': [1.0, 0.0, 0.0],
        'yellow': [1.0, 1.0, 0.0]
    }

    if color not in color_dict:
        raise ValueError(f"Unsupported color: {color}")

    color_rgb = color_dict[color]
    num_points = 100
    x = np.linspace(-scale, scale, num_points)
    y = np.linspace(-scale, scale, num_points)
    X, Y = np.meshgrid(x, y)
    Z = (-plane.a * X - plane.b * Y - plane.d) / plane.c

    cgo_object = [cgo.BEGIN, cgo.TRIANGLES, cgo.COLOR, *color_rgb]

    for i in range(num_points - 1):
        for j in range(num_points - 1):
            v1 = [X[i, j], Y[i, j], Z[i, j]]
            v2 = [X[i + 1, j], Y[i + 1, j], Z[i + 1, j]]
            v3 = [X[i, j + 1], Y[i, j + 1], Z[i, j + 1]]
            v4 = [X[i + 1, j + 1], Y[i + 1, j + 1], Z[i + 1, j + 1]]

            cgo_object.extend([cgo.VERTEX, *v1, cgo.VERTEX, *v2, cgo.VERTEX, *v3])
            cgo_object.extend([cgo.VERTEX, *v2, cgo.VERTEX, *v3, cgo.VERTEX, *v4])

    cgo_object.append(cgo.END)
    return cgo_object

def show_in_pymol(plane1, plane2, pdb_file):
    """
    Visualizes the protein structure and planes in PyMOL and saves the session to a .pse file.

    Args:
    - plane1 (Plane): The first Plane object to be visualized.
    - plane2 (Plane): The second Plane object to be visualized.
    - pdb_file (str): The path to the PDB file of the protein structure.

    The function saves the PyMOL session as a .pse file in the '../results' directory.
    """
    pymol.finish_launching()
    cmd.load(pdb_file, 'protein_structure')

    plane1_cgo = create_plane_object(plane1, color='blue', scale=15)
    plane2_cgo = create_plane_object(plane2, color='green', scale=15)

    cmd.load_cgo(plane1_cgo, 'plane1')
    cmd.load_cgo(plane2_cgo, 'plane2')

    # Debugging information
    print(f"Plane1 - Normal: ({plane1.a}, {plane1.b}, {plane1.c}), d={plane1.d}")
    print(f"Plane2 - Normal: ({plane2.a}, {plane2.b}, {plane2.c}), d={plane2.d}")

    cmd.show('surface', 'plane1')
    cmd.show('surface', 'plane2')
    cmd.color('blue', 'plane1')
    cmd.color('green', 'plane2')
    cmd.zoom('protein_structure')
    
    # Ensure the 'results' directory exists
    os.makedirs('../results', exist_ok=True)
    
    # Save PyMOL session in the 'results' directory
    session_file = os.path.join('../results', pdb_file.split('/')[-1].replace('.pdb', '.pse'))
    cmd.save(session_file)

    print(f"Visualization in PyMOL completed. Session saved to {session_file}.")
