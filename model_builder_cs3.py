import numpy as np
import copy
import math
from dimod import ConstrainedQuadraticModel, QuadraticModel, BinaryQuadraticModel
from itertools import combinations
from math import radians as DegToRad


def GaussianWake(frame_coords, turb_diam):
    """Return each turbine's total loss due to wake from upstream turbines"""
    # Equations and values explained in <iea37-wakemodel.pdf>
    num_turb = len(frame_coords)

    # Constant thrust coefficient
    CT = 4.0*1./3.*(1.0-1./3.)
    # Constant, relating to a turbulence intensity of 0.075
    k = 0.0324555
    # Array holding the wake deficit seen at each turbine
    loss = np.zeros(num_turb)

    for i in range(num_turb):            # Looking at each turb (Primary)
        loss_array = np.zeros(num_turb)  # Calculate the loss from all others
        for j in range(num_turb):        # Looking at all other turbs (Target)
            x = frame_coords.x[i] - frame_coords.x[j]   # Calculate the x-dist
            y = frame_coords.y[i] - frame_coords.y[j]   # And the y-offset
            if x > 0.:                   # If Primary is downwind of the Target
                sigma = k*x + turb_diam/np.sqrt(8.)  # Calculate the wake loss
                # Simplified Bastankhah Gaussian wake model
                exponent = -0.5 * (y/sigma)**2
                radical = 1. - CT/(8.*sigma**2 / turb_diam**2)
                loss_array[j] = (1.-np.sqrt(radical)) * np.exp(exponent)
            # Note that if the Target is upstream, loss is defaulted to zero
        # Total wake losses from all upstream turbs, using sqrt of sum of sqrs
        loss[i] = np.sqrt(np.sum(loss_array**2))

    return loss


def build_cqm(wind_tower_labels, coordinate, wind_direction_bins, wind_speed_bins,
              rated_wind_speed, diameter, cut_in_wind_speed, cut_out_wind_speed,
              wind_turbine_power, wind_direction_probability, wind_speed_probability,
              DirPower, n_turbines, potential_nodes, minimum_distance):
    """
    Build a Constrained Quadratic Model (CQM) with the specified wind farm optimization objective.

    Parameters:
    - wind_tower_labels (list): List of labels for wind tower locations in Cartesian coords.
    - coordinate (dtype): Data type for coordinate values.
    - wind_direction_bins (list): List of wind direction bins.
    - wind_speed_bins (list): List of wind speed bins.
    - rated_wind_speed (float): Rated wind speed for the turbines.
    - diameter (float): Diameter of the wind turbines.
    - cut_in_wind_speed (float): Cut-in wind speed for the turbines.
    - cut_out_wind_speed (float): Cut-out wind speed for the turbines.
    - wind_turbine_power (float): Power output of a wind turbine.
    - wind_direction_probability (list): List of probabilities corresponding to wind direction bins.
    - wind_speed_probability (list): List of probabilities corresponding to wind speed bins.
    - DirPower (function): Function to calculate the power produced by each turbine.
    - n_turbines (int): Number of turbines to be selected.
    - potential_nodes (list): List of potential nodes for wind turbine locations.
    - minimum_distance (float): Minimum distance between any pair of turbines.

    Returns:
    - ConstrainedQuadraticModel: The constructed CQM model.
    """

    interactions = {}
    
    for node1 in wind_tower_labels:
        for node2 in wind_tower_labels:
            # Skip if it's the same node
            if node1 == node2:
                continue

            # Extract coordinates from node1
            x1, y1 = map(float, node1.split('_'))

            # Extract coordinates from node2
            x2, y2 = map(float, node2.split('_'))

            turb_coords = np.asarray([[x1, y1], [x2, y2]])
            
            coordinate_dtype = np.dtype([('x', 'f8'), ('y', 'f8')])
            frame_coords = np.recarray(len(turb_coords), coordinate)
            #frame_coords.x, frame_coords.y = turb_xc, turb_yc


            power_values = []
            for wind_direction_bin in wind_direction_bins:
                for wind_speed_bin in wind_speed_bins:
                    # Assuming wind_speed_bin is the actual wind speed at this point
                    power_value = DirPower(
                        frame_coords,
                        GaussianWake(frame_coords, diameter),
                        wind_speed_bin,
                        cut_in_wind_speed,
                        cut_out_wind_speed,
                        wind_turbine_power,
                        rated_wind_speed  # Adjust this line based on the correct arguments
                    )
                    power_values.append(power_value)

            # Store the results in the interactions dictionary
            interactions[(node1, node2)] = copy.deepcopy(power_values)

    refactored_interactions = {}

    for dir_index, wind_direction_bin in enumerate(wind_direction_bins):
        for speed_index, wind_speed_bin in enumerate(wind_speed_bins):
            index = dir_index * len(wind_speed_bins) + speed_index
            bin_key = index  # Using the bin index as the new key
            bin_values = {}

            for pair, actions in interactions.items():
                pair_key = pair  # Keeping the pair as a key in the nested dictionary
                action_value = actions[index]  # Extracting the action for the current bin
                #probability = wind_direction_probability[dir_index] * wind_speed_probability[speed_index]
                probability = wind_direction_probability[dir_index] * wind_speed_probability[speed_index][0]

                # Multiply the action by the probability
                action_value *= probability * 8760

                bin_values[pair_key] = action_value

            refactored_interactions[bin_key] = bin_values

    # Create a ConstrainedQuadraticModel
    cqm = ConstrainedQuadraticModel()

    objective = BinaryQuadraticModel(vartype='BINARY')

    # Iterate through the refactored_interactions dictionary
    for bin_index, bin_values in refactored_interactions.items():
        for dv_pair, action_value in bin_values.items():
            dv1, dv2 = dv_pair  # Assuming dv_pair is a tuple with two DV labels

            # Skip if it's the same DV
            if dv1 == dv2:
                continue

            # Add the quadratic term to the objective
            objective.add_quadratic(dv1, dv2, -action_value)

    cqm.set_objective(objective)

    # Constraint: Choose exactly n_turbines
    constraint1 = QuadraticModel()
    constraint1.add_linear_from(
        ((wind_tower_labels[i], 1) for i in range(len(wind_tower_labels))),
        default_vartype='BINARY'
    )
    cqm.add_constraint(constraint1, sense="==", rhs=n_turbines, label='n_turbines', copy=False)

    # Create a dictionary to associate each tuple of Cartesian coordinates with a BinaryQuadraticModel
    wind_towers = {}
    for x, y in potential_nodes:
        label = f'{x:.1f}_{y:.0f}'  # Adjust the label format
        key = (round(x, 1), round(y, 0))  # Round the values when creating the key
        wind_towers[key] = BinaryQuadraticModel({label: 1.0}, {}, 0.0, 'BINARY')

    pair_list = list(combinations(wind_towers.keys(), 2))

    # Calculate the minimum distance constraint
    for i, (a, b) in enumerate(pair_list):
        x1, y1 = a
        x2, y2 = b

        # Calculate the distance between two points in Cartesian coordinates
        dist = math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

        # Enforce the minimum distance constraint 
        if dist <= minimum_distance:
            cqm.add_constraint(wind_towers[a] + wind_towers[b] <= 1, label=f"min_distance_{a}_{b}")

    return cqm

# Example usage:
# build_cqm(wind_tower_labels, coordinate, wind_direction_bins, wind_speed_bins,
#           rated_wind_speed, diameter, cut_in_wind_speed, cut_out_wind_speed,
#           wind_turbine_power, wind_direction_probability, wind_speed_probability,
#           DirPower, n_turbines, potential_nodes
