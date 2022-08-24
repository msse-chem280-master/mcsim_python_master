"""
Simulation Loop
"""

import math
import random

from copy import deepcopy


def accept_or_reject(delta_e, beta):
    """
    Accept or reject based on change in energy and temperature.
    """

    if delta_e <= 0:
        accept = True
    else:
        random_number = random.random()
        p_acc = math.exp(-beta * delta_e)

        if random_number < p_acc:
            accept = True
        else:
            accept = False

    return accept


def run_simulation(
    coordinates,
    box_length,
    cutoff,
    reduced_temperature,
    num_steps,
    max_displacement=0.1,
    freq=1000,
    use_numpy=True,
):
    """
    Run a Monte Carlo simulation with the specified parameters.
    """

    if use_numpy:
        import numpy as np
        import mcsim.monte_carlo_numpy as mc

        coordinates = np.array(coordinates)
    else:
        import mcsim.monte_carlo as mc

    # Reporting information
    steps = []
    energies = []
    all_coordinates = []
    pressures = []

    # Calculated quantities
    beta = 1 / reduced_temperature
    num_particles = len(coordinates)

    # Calculated based on simulation inputs
    total_energy = mc.calculate_total_energy(
        coordinates=coordinates, box_length=box_length, cutoff=cutoff
    )
    total_energy += mc.calculate_tail_correction(
        num_particles=num_particles, box_length=box_length, cutoff=cutoff
    )

    for step in range(num_steps):

        # 1. Randomly pick one of num_particles particles
        random_particle = random.randrange(num_particles)

        # 2. Calculate the interaction energy of the selected particle with the system. Store this value.
        current_energy = mc.calculate_pair_energy(
            coordinates=coordinates,
            i_particle=random_particle,
            box_length=box_length,
            cutoff=cutoff,
        )

        # 3. Generate a random x, y, z displacement range (-max_displacement, max_displacement) - uniform distribution
        x_rand = random.uniform(-max_displacement, max_displacement)
        y_rand = random.uniform(-max_displacement, max_displacement)
        z_rand = random.uniform(-max_displacement, max_displacement)

        # 4. Modify the coordinate of selected particle by generated displacements.
        coordinates[random_particle][0] += x_rand
        coordinates[random_particle][1] += y_rand
        coordinates[random_particle][2] += z_rand

        # 5. Calculate the new interaction energy of moved particle, store this value.
        proposed_energy = mc.calculate_pair_energy(
            coordinates=coordinates,
            i_particle=random_particle,
            box_length=box_length,
            cutoff=cutoff,
        )

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
            pressure = mc.calculate_total_pressure(
                coordinates, box_length, cutoff, reduced_temperature
            )
            print(step, total_energy / num_particles, pressure)
            steps.append(step)
            energies.append(total_energy / num_particles)
            pressures.append(pressure)
            all_coordinates.append(deepcopy(coordinates))

    return all_coordinates
