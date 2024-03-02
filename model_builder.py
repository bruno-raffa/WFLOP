
import numpy as np
import copy, math
from dimod import ConstrainedQuadraticModel, Integer, Binary, quicksum, QuadraticModel, BinaryQuadraticModel
from itertools import combinations
from math import radians as DegToRad, cos, radians, sin, sqrt

def build_cqm(wind_tower_labels, coordinate, wind_speed_bins,rated_wind_speed, diameter, 
              cut_in_wind_speed, cut_out_wind_speed,wind_turbine_power, wind_speed_probability,
              DirPower, n_turbines, potential_nodes, minimum_distance ):

    """
    Build a Constrained Quadratic Model (CQM) with the specified wind farm optimization objective.

    Parameters:
    - wind_tower_labels (list): List of labels for wind tower locations in polar coords

    Returns:
    - ConstrainedQuadraticModel: The constructed CQM model.
    """

    interactions = {}

    for node1 in wind_tower_labels:
        for node2 in wind_tower_labels:
            # Skip if it's the same node
            if node1 == node2:
                continue

            # Extract angle and radius from node1
            if '.' in node1:  # Check if it's in the format '0.0_433'
                angle_i, radius_i = map(float, node1.split('_'))
            else:  # Assume it's in the format '0_433'
                angle_i, radius_i = map(float, node1.split('_'))

            # Extract angle and radius from node2
            if '.' in node2:  # Check if it's in the format '0.0_433'
                angle_j, radius_j = map(float, node2.split('_'))
            else:  # Assume it's in the format '0_433'
                angle_j, radius_j = map(float, node2.split('_'))

            # Convert coordinates to downwind/crosswind coordinates
            x1, y1 = radius_i * np.cos(np.radians(angle_i)), radius_i * np.sin(np.radians(angle_i))
            x2, y2 = radius_j * np.cos(np.radians(angle_j)), radius_j * np.sin(np.radians(angle_j))
            turb_xc = np.asarray([x1, x2])
            turb_yc = np.asarray([y1, y2])
            coordinate = np.dtype([('x', 'f8'), ('y', 'f8')])
            turb_coords = np.recarray(turb_xc.shape, coordinate)
            turb_coords.x, turb_coords.y = turb_xc, turb_yc

            power_values = []
            for wind_speed_bin in wind_speed_bins:
                # Assuming wind_speed_bin is the actual wind speed at this point
                power_value = DirPower(
                    turb_coords,
                    wind_speed_bin,
                    rated_wind_speed,
                    diameter,
                    cut_in_wind_speed,
                    cut_out_wind_speed,
                    rated_wind_speed,
                    wind_turbine_power,
                )
                power_values.append(power_value)

            # Store the results in the interactions dictionary
            interactions[(node1, node2)] = copy.deepcopy(power_values)


    refactored_interactions = {}

    for bin_index, wind_speed_bin in enumerate(wind_speed_bins):
        bin_key = bin_index  # Using the bin index as the new key
        bin_values = {}

        for pair, actions in interactions.items():
            pair_key = pair  # Keeping the pair as a key in the nested dictionary
            action_value = actions[bin_index]  # Extracting the action for the current bin
            probability = wind_speed_probability[bin_index]  # Extracting the corresponding probability

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

    #  Constraint : Choose exactly n_turbines
    constraint1 = QuadraticModel()
    constraint1.add_linear_from(
                    (
                    (wind_tower_labels[i], 1) for i in range(len(wind_tower_labels)) 
                    ), 
                    default_vartype='BINARY'
                    )
    cqm.add_constraint(constraint1, sense="==", rhs=n_turbines, label='n_turbines', copy=False)

    #  Create a dictionary to associate each tuple of polar coordinates with a BinaryQuadraticModel
    wind_towers = {}
    for angle, radius in potential_nodes:
        label = f'{angle:.1f}_{radius:.0f}'  # Adjust the label format
        key = (round(angle, 1), round(radius, 0))  # Round the values when creating the key
        wind_towers[key] = BinaryQuadraticModel({label: 1.0}, {}, 0.0, 'BINARY')

    pair_list = list(combinations(wind_towers.keys(), 2))

    # Calculate the minimum distance constraint
    for i, (a, b) in enumerate(pair_list):
        angle1, radius1 = a
        angle2, radius2 = b

        # Calculate the distance between two points in polar coordinates
        dist = math.sqrt((radius1**2 + radius2**2) - 2 * radius1 * radius2 * math.cos(math.radians(angle1 - angle2)))

        # Enforce the minimum distance constraint 
        if dist <= minimum_distance:
            cqm.add_constraint(wind_towers[a] + wind_towers[b] <= 1, label=f"min_distance_{a}_{b}")


    return cqm
