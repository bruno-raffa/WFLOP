#  Wind Farm Layout Optimization Problem


The Wind Farm Layout Optimization Problem (WFLOP) is a well-known problem in the field of renewable energy, particularly in the design and planning of wind farms. The goal of WFLOP is to determine the optimal placement of wind turbines within a designated area to maximize the overall energy output of the wind farm. This involves considering various factors and constraints to find an arrangement that balances efficiency, cost, and environmental impact.

Key factors and considerations in WFLOP typically include:

- Wind Turbine Spacing: Determining the optimal distance between individual wind turbines to avoid interference and maximize energy capture. Turbine wakes can reduce the - efficiency of downstream turbines.

- Wind Direction and Speed: Taking into account prevailing wind patterns and variations in wind speed across the designated area to position turbines for optimal energy production.

- Terrain and Topography: Considering the geographical features of the site, such as hills, valleys, and obstacles, which can affect wind patterns and turbine performance.

- Land Use and Environmental Constraints: Factoring in environmental regulations, wildlife concerns, and other land-use constraints to ensure that the wind farm's impact is minimized.

- Costs and Infrastructure: Considering the costs associated with turbine installation, grid connection, and other infrastructure requirements.

The IEA37 (https://github.com/byuflowlab/iea37-wflo-casestudies) submits a case-study aimed at optimize the optimal placement for different scenarios.


The approach used in this notebook involves the use of a Constrained Quadratic Model (CQM) from Dwave (https://cloud.dwavesys.com/). The solution is computed using the LeapHybridCQMSampler.


### Objective
**Highest Power Ouput**. Maximize total output of the Wind Farm in terms of AEP (Annual Energy Production) for all different wind directions

### Constraints

1. All turbines will have to be confined to the assigned area: a circular-shaped area incase of sub-problem CS1.
2. No turbines can be positioned at a distance lesser than two rotor diameters from any other turbines.


### User Inputs
- Nr of Turbines and rotor radius. 
- Radius of the areas


### Output
AEP of the optimized layout and a visual diagram with their coordinates (both polar and cartesian)

## References

Baker, Nicholas & Stanley, Andrew & Thomas, Jared & Ning, Andrew & Dykes, Katherine. (2019). Best Practices for Wake Model and Optimization Algorithm Selection in Wind Farm Layout Optimization. 10.2514/6.2019-0540. 