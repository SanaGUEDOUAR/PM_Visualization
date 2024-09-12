"""
This module defines geometric entities used in the analysis of protein structures.
It includes classes for vectors, points, planes, and axes.

Classes:
- Vector: Represents a 3D vector.
- Point: Represents a 3D point with coordinates.
- Plane: Defines a plane using a point and a normal vector.
- Axis: Defines an axis as a pair of parallel planes.
"""

import numpy as np
import copy

class Vector:
    """
    Represents a 3D vector.

    Attributes:
    - coordinates (numpy.ndarray): The (x, y, z) coordinates of the vector.
    """

    def __init__(self, x, y, z):
        self.coordinates = np.array([x, y, z])

    def get_coordinates(self):
        """
        Returns the coordinates of the vector.

        Returns:
        - numpy.ndarray: The (x, y, z) coordinates of the vector.
        """
        return self.coordinates

class Point:
    """
    Represents a 3D point with coordinates.

    Attributes:
    - coordinates (numpy.ndarray): The (x, y, z) coordinates of the point.
    """

    def __init__(self, x, y, z):
        self.coordinates = np.array([x, y, z])

    def get_x(self):
        """
        Returns the x-coordinate of the point.

        Returns:
        - float: The x-coordinate.
        """
        return self.coordinates[0]

    def get_y(self):
        """
        Returns the y-coordinate of the point.

        Returns:
        - float: The y-coordinate.
        """
        return self.coordinates[1]

    def get_z(self):
        """
        Returns the z-coordinate of the point.

        Returns:
        - float: The z-coordinate.
        """
        return self.coordinates[2]

    def get_coordinates(self):
        """
        Returns the coordinates of the point.

        Returns:
        - numpy.ndarray: The (x, y, z) coordinates of the point.
        """
        return self.coordinates

    def __eq__(self, other):
        """
        Checks equality with another point.

        Args:
        - other (Point): The point to compare with.

        Returns:
        - bool: True if the points are equal, False otherwise.
        """
        if isinstance(other, Point):
            return np.array_equal(self.coordinates, other.get_coordinates())
        return False

    def __hash__(self):
        """
        Returns the hash value of the point.

        Returns:
        - int: The hash value of the point.
        """
        return hash(tuple(self.coordinates))

class Plane:
    """
    Defines a plane using a point and a normal vector.

    Attributes:
    - point (Point): A point on the plane.
    - normal (Vector): The normal vector of the plane.
    - a, b, c (float): Coefficients of the plane equation Ax + By + Cz + D = 0.
    - d (float): Constant term in the plane equation.
    """

    def __init__(self, point: Point, normal: Vector):
        self.point = point
        self.normal = normal
        self.a, self.b, self.c = normal.get_coordinates()
        self.d = -(self.a * self.point.get_x() +
                   self.b * self.point.get_y() +
                   self.c * self.point.get_z())

    def complementary(self, gap):
        """
        Creates a complementary plane with a specified gap.

        Args:
        - gap (float): The gap to add to the plane's constant term.

        Returns:
        - Plane: The complementary plane.
        """
        complementary_plane = copy.deepcopy(self)
        complementary_plane.d += gap
        return complementary_plane

    def is_above(self, point: Point):
        """
        Checks if a point is above the plane.

        Args:
        - point (Point): The point to check.

        Returns:
        - bool: True if the point is above the plane, False otherwise.
        """
        return (self.a * point.get_x() +
                self.b * point.get_y() +
                self.c * point.get_z() + self.d) > 0

    def is_below(self, point: Point):
        """
        Checks if a point is below the plane.

        Args:
        - point (Point): The point to check.

        Returns:
        - bool: True if the point is below the plane, False otherwise.
        """
        return (self.a * point.get_x() +
                self.b * point.get_y() +
                self.c * point.get_z() + self.d) < 0

    def __str__(self):
        """
        Returns a string representation of the plane.

        Returns:
        - str: The string representation of the plane.
        """
        return f"Plane with normal ({self.a:.3f}, {self.b:.3f}, {self.c:.3f}) and d={self.d:.3f}"

class Axis:
    """
    Defines an axis as a pair of parallel planes.

    Attributes:
    - plane1 (Plane): The first plane of the axis.
    - plane2 (Plane): The second plane of the axis.
    """

    def __init__(self, plane1, plane2):
        self.plane1 = plane1
        self.plane2 = plane2

    def __str__(self):
        """
        Returns a string representation of the axis.

        Returns:
        - str: The string representation of the axis.
        """
        return f"Axis with Plane1: {self.plane1} and Plane2: {self.plane2}"
