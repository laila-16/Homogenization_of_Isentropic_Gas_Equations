# Homogenization of Isentropic Gas Equations

In this repository, you can find all the codes used in the paper 

> Homogenized Equations for Isentropic Gas in a Pipe with Periodically Varying Cross-Section

by Laila Busaleh and David I. Ketcheson.

The directories contain the following:

    - `Homogenized_Equations_Derivation_Mathematica`: A Mathematica notebook in which we derive the homogenized equations
    - `Pseudospectral_Code_Homogenized_System`: A pseudospectral solver for the homogenized equations
    - `RiemannSolver_IsentropicGas_Equations`: Code to set up and solve the variable-coefficient problem with PyClaw.

Use the following commands to run the simulations for the finite volume solution:

1. In the terminal we use this command to compile the Riemann Solver:
   ```
    cd RiemannSolver_IsentropicGas_Equations
    python3 -m numpy.f2py -c rp1_IsentropicGas_1d_capacity_function.f90 -m RS_variable_coeff
   ```

2. To generate the plots
   
   For piece-wise constant case:
   ```
   mpirun -np 32 python3 IsentropicGas_1d_pwc.py htmlplot=1
   ```
   For sinusoidal case:
   ```
   mpirun -np 32 python3 IsentropicGas_1d_sin.py htmlplot=1
   ```


To reproduce the homogenized system solution you can run the code included in Pseudospectral_Code_Homogenized_System .
