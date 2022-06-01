# Ethyl-Acetate

This project's goal was to perform a techno-economical analysis of an Ethyl Acetate chemical plant.

The relevant scripts in this project are as follows:

MATLAB Files
1. EtaNPV_Reactor_Optimization184B.m
• This code iterates through temperatures and pressures of interest, calling EtAcetateReactor_PBR_Ideal_184B m. It ignores the separation system as that is only costed at discrete points. It then performs an economic analysis, costing all equipment, chemicals, and energy required to generate a series of cashflows which are then used to calculate the NPV of the project.
2. Discrete_EtaNPV_Reactor_Optimization184B.m
• This code operates on the same chasis that the former code does, except that it only calculates economic parameters at discrete points. These points line up for specific conversions at 240C and 2 atm only. These discrete points allow the code to include data about the separations system from Aspen plus models in order to compare the cost of separating reactor effluents at different conversions.
3. EtAcetateReactor_PBR_Ideal_184B.m
• This function calculates the plugged flow reactor conditions when called using an ODE solver, outputting molar flows at various volumes in the PFR.
4. Sensitivity_HYSYS_EthylAcetate.m
• This function uses the base case parameters from the HYSYS flowsheet to calculate the NPV of the project as a function of various parameters such as chemical costs and tax rates. This function can be used for sensitivity analysis and is called by MonteCarlo_TriangleDist_EthylAcetate.m.
5. MonteCarlo_TriangleDist_EthylAcetate.m
• This code performs a Monte Carlo simulation on the NPV of the project by assigning triangle distributions to each of the variables of interest, randomizing them N times, and calling Sensitivity_HYSYS_EthylAcetate.m to calculate the NPV. Then a probability distribution function and cumulative distribution function are created and plotted.
