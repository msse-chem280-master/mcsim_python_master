"""
Functions for system generation.
"""

import math
import random


def read_xyz(filepath):
    """
    Reads coordinates from an xyz file.

    Parameters
    ----------
    filepath : str
       The path to the xyz file to be processed.

    Returns
    -------
    atomic_coordinates : list
        A two dimensional list containing atomic coordinates
    """

    with open(filepath) as f:
        box_length = float(f.readline().split()[0])
        num_atoms = float(f.readline())
        coordinates = f.readlines()

    atomic_coordinates = []

    for atom in coordinates:
        split_atoms = atom.split()

        float_coords = []

        # We split this way to get rid of the atom label.
        for coord in split_atoms[1:]:
            float_coords.append(float(coord))

        atomic_coordinates.append(float_coords)

    return atomic_coordinates, box_length


def generate_random_coordinates(num_atoms, density):
    """
    Generate random coordinates in a box.

    Parameters
    ----------
    num_atoms : int
        The number of atoms to place
    density : float
        The target system density

    Returns
    -------
    coordinates : list
        The generated coordinates
    box_length : float
        The box length for the number of atoms and desired density.
    """

    box_length = math.pow(num_atoms / density, (1 / 3))
    coordinates = []

    for i in range(num_atoms):
        x_val = random.uniform(0, box_length)
        y_val = random.uniform(0, box_length)
        z_val = random.uniform(0, box_length)
        coordinates.append([x_val, y_val, z_val])

    return coordinates, box_length
