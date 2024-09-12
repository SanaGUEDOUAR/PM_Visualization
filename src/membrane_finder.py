"""
This module defines the MembraneFinder class, which is responsible for finding the best membrane
axis in a protein structure by analyzing hydrophobicity between planes.

Classes:
- MembraneFinder: Finds the best membrane axis by analyzing hydrophobicity between planes.
"""

import numpy as np
from geometry_and_axis import Point, Plane, Axis, Vector
from pymol import cmd, cgo

class MembraneFinder:
    """
    Finds the best membrane axis in a protein structure by analyzing hydrophobicity between planes.

    Attributes:
    - protein (Protein): The protein object containing structural information.
    """

    def __init__(self, protein):
        self.protein = protein

    def _count_residues_between_planes(self, plane1, plane2):
        """
        Counts hydrophobic and hydrophilic residues between two planes.

        Args:
        - plane1 (Plane): The first plane.
        - plane2 (Plane): The second plane.

        Returns:
        - float: The hydrophobicity score between the two planes.
        """
        in_between_planes = set()
        n_total_hydrophobic = 0
        n_total_hydrophile = 0
        nb_hydrophile_out_of_plan = 0
        n_hydrophobe_in_plan = 0

        # Identifying residues between planes
        for residue in self.protein.residues:
            if 'CA' in residue:
                ca_coord = residue['CA'].coord
                point = Point(ca_coord[0], ca_coord[1], ca_coord[2])
                
                if (plane1.is_below(point) and plane2.is_above(point)) or \
                   (plane2.is_below(point) and plane1.is_above(point)):
                    in_between_planes.add(point)
                
                is_hydrophobic = self.protein.hydrophobicity.get(f"{residue.get_resname()}{residue.get_id()}", 0) == 1
                if is_hydrophobic:
                    n_total_hydrophobic += 1
                else:
                    n_total_hydrophile += 1

        # Counting hydrophilic residues out of planes
        for residue in self.protein.residues:
            if 'CA' in residue:
                ca_coord = residue['CA'].coord
                point = Point(ca_coord[0], ca_coord[1], ca_coord[2])
                if point not in in_between_planes:
                    if not self.protein.hydrophobicity.get(f"{residue.get_resname()}{residue.get_id()}", 0):
                        nb_hydrophile_out_of_plan += 1

        # Counting hydrophobic residues in between planes
        for residue in self.protein.residues:
            if 'CA' in residue:
                ca_coord = residue['CA'].coord
                point = Point(ca_coord[0], ca_coord[1], ca_coord[2])
                if point in in_between_planes:
                    if self.protein.hydrophobicity.get(f"{residue.get_resname()}{residue.get_id()}", 0) == 1:
                        n_hydrophobe_in_plan += 1

        if len(in_between_planes) == 0 or len(in_between_planes) >= len(self.protein.residues):
            return 0

        # Calculating hydrophobicity
        return (nb_hydrophile_out_of_plan / (n_total_hydrophile + 1e-10)) + \
               (n_hydrophobe_in_plan / (n_total_hydrophobic + 1e-10))

    def find_best_membrane(self):
        """
        Finds the best membrane axis by exploring possible planes and calculating hydrophobicity.

        Returns:
        - Axis: The best axis found.
        - float: The hydrophobicity of the best axis.
        """
        best_hydrophobicity = -float('inf')
        best_axis = None
        num_points = 10

        center_point = Point(*self.protein.center_of_mass)
        for point in self._generate_points_on_hemisphere(num_points, center_point):
            axis_vector = Vector(*(point.get_coordinates() - center_point.get_coordinates()))
            for distance in np.arange(10, 20, 1):  # Adjust the step and range if necessary
                plane1 = Plane(center_point, axis_vector)
                plane2 = plane1.complementary(-distance)
                hydrophobicity = self._count_residues_between_planes(plane1, plane2)
                
                if hydrophobicity > best_hydrophobicity:
                    best_hydrophobicity = hydrophobicity
                    best_axis = Axis(plane1, plane2)

        return best_axis, float(best_hydrophobicity)

    def _generate_points_on_hemisphere(self, num_points, center_point):
        """
        Generates points on a hemisphere centered around a given point.

        Args:
        - num_points (int): Number of points to generate.
        - center_point (Point): The center point of the hemisphere.

        Returns:
        - list: List of generated points on the hemisphere.
        """
        phi = np.linspace(0, np.pi / 2, num_points)  # Phi ranges from 0 to pi/2 for a hemisphere
        theta = np.linspace(0, 2 * np.pi, num_points)  # Theta ranges from 0 to 2pi to cover the horizon
        points = []

        for phi_val in phi:
            for theta_val in theta:
                x = center_point.get_x() + np.sin(phi_val) * np.cos(theta_val)
                y = center_point.get_y() + np.sin(phi_val) * np.sin(theta_val)
                z = center_point.get_z() + np.cos(phi_val)
                if z > center_point.get_z():  # Ensure the point is above the plane z=center_z
                    points.append(Point(x, y, z))
                    
        return points
