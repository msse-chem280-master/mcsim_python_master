"""
Functions for running Monte Carlo Simulation.
"""

import os
import random
import math

import numpy as np

import matplotlib.pyplot as plt

def calculate_LJ(r_ij):
    """
    Calculate the Lennard Jones energy based on a separation distance.
    
    Parameters
    ----------
    r_ij : float
        The distance between the two particles.
        
    Returns
    -------
    pairwise_energy : float
        The Lennard Jones potential energy.
    """
    
    r6_term = np.power(1/r_ij, 6)
    r12_term = np.power(r6_term, 2)
    
    pairwise_energy = 4 * (r12_term - r6_term)
    
    return pairwise_energy

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

    box_length = math.pow(num_atoms/density, (1/3))
    coordinates = []

    for i in range(num_atoms):
        x_val = random.uniform(0, box_length)
        y_val = random.uniform(0, box_length)
        z_val = random.uniform(0, box_legnth)
        coordinates.append([x_val, y_val, z_val])
    
    return coordinates, box_length

def calculate_distance(coord1, coord2, box_length=None):
    """
    Calculate the distance between two 3D coordinates.
    Parameters
    ----------
    coord1, coord2: np.ndarray
        The atomic coordinates
    box_length : float
        The box length to be used for minimum image distance.

    Returns
    -------
    distance: float
        The distance between the two points.
    """
    coord_dist = coord1 - coord2

    if box_length:
        coord_dist = coord_dist - box_length * np.round(coord_dist / box_length)

    if coord_dist.ndim < 2:
        coord_dist = coord_dist.reshape(1, -1)

    coord_dist = coord_dist ** 2

    coord_dist_sum = coord_dist.sum(axis=1)

    distance = np.sqrt(coord_dist_sum)

    return distance

def calculate_tail_correction(num_particles, box_length, cutoff):
    """
    Calculate the long range tail correction
    """
    
    const1 = (8 * math.pi * num_particles ** 2) / (3 * box_length ** 3)
    const2 = (1/3) * (1 / cutoff)**9 - (1 / cutoff) **3
    
    return const1 * const2

def calculate_total_energy(coordinates, box_length, cutoff):
    """
    Numpy version - rewrite (calculate_total_energy)
    """
    
    total_energy = 0
    num_atoms = len(coordinates)

    for i in range(num_atoms):
        calc_coords = coordinates[i:]
        total_energy += calculate_pair_energy(calc_coords, 0, box_length, cutoff)
    
    return total_energy


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

def accept_or_reject(delta_e, beta):
    """
    Accept or reject based on change in energy and temperature.
    """
    
    if delta_e <= 0:
        accept = True
    else:
        random_number = random.random()
        p_acc = math.exp(-beta*delta_e)
        
        if random_number < p_acc:
            accept = True
        else:
            accept = False
    
    return accept

def calculate_pair_energy(coordinates, i_particle, box_length, cutoff):
    """
    Calculate the interaction energy of a particle with its environment (all other particles in the system) - rewrite
    
    Parameters
    ----------
    coordinates : list
        The coordinates for all particles in the system
        
    i_particle : int
        The particle number for which to calculate the energy
        
    cutoff : float
        The simulation cutoff. Beyond this distance, interactions are not calculated.
    
    Returns
    -------
    e_total : float 
        The pairwise interaction energy of he i_th particle with all other particles in the system.
    """
    
    e_total = 0.0
    i_position = coordinates[i_particle]
    
    distance_array = calculate_distance(coordinates, i_position, box_length)
    
    # Just so we don't use it for calculation
    distance_array[i_particle] = cutoff*2
    
    less_than_cutoff = distance_array[distance_array < cutoff]
    
    interaction_energies = calculate_LJ(less_than_cutoff)
 
    e_total = np.sum(interaction_energies)

    return e_total

def run_simulation(coordinates, box_length, cutoff, reduced_temperature, num_steps,max_displacement=0.1, freq=1000):
    """
    Run a Monte Carlo simulation with the specified parameters.
    """
    # Reporting information
    steps = []
    energies = []
    all_coordinates = []

    # Calculated quantities
    beta = 1/reduced_temperature
    num_particles = len(coordinates)

    # Calculated based on simulation inputs
    total_energy = calculate_total_energy(coordinates, box_length, cutoff)
    total_energy += calculate_tail_correction(num_particles, box_length, cutoff)

    for step in range(num_steps):
        
        # 1. Randomly pick one of num_particles particles
        random_particle = random.randrange(num_particles)
        
        # 2. Calculate the interaction energy of the selected particle with the system. Store this value.
        current_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff)
        
        # 3. Generate a random x, y, z displacement range (-max_displacement, max_displacement) - uniform distribution
        x_rand = random.uniform(-max_displacement, max_displacement)
        y_rand = random.uniform(-max_displacement, max_displacement)
        z_rand = random.uniform(-max_displacement, max_displacement)
        
        # 4. Modify the coordinate of selected particle by generated displacements.
        coordinates[random_particle][0] += x_rand
        coordinates[random_particle][1] += y_rand
        coordinates[random_particle][2] += z_rand
        
        # 5. Calculate the new interaction energy of moved particle, store this value.
        proposed_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff)
        
        # 6. Calculate energy change and decide if we accept the move.
        delta_energy = proposed_energy - current_energy
        
        accept = accept_or_reject(delta_energy, beta)
        
        # 7. If accept, keep movement. If not revert to old position.
        if accept:
            total_energy += delta_energy
        else:
            # Move is not accepted, roll back coordinates
            coordinates[random_particle][0] -= x_rand
            coordinates[random_particle][1] -= y_rand
            coordinates[random_particle][2] -= z_rand
        
        # 8. Print the energy and store the coordinates at certain intervals
        if step % freq == 0:
            print(step, total_energy/num_particles)
            steps.append(step)
            energies.append(total_energy/num_particles)
            all_coordinates.append(coordinates)
    
    return coordinates
    