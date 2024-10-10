# Homogenization of Isentropic Gas Equations
In this repository, you will find all the codes used in the paper Homogenized Equations for Isentropic Gas in a Pipe with Periodically Varying Cross-Section. Follow these steps to run the simulations for the finite volume solution:
1. From the terminal, navigate to the directory where the code is located.

2. Use this command to compile the Riemann Solver:
   ```
    python3 -m numpy.f2py -c rp1_IsentropicGas_1d_capacity_function.f90 -m RS_variable_coeff
   ```

3. To generate the plots
   
   For piece-wise constant case:
   ```
   mpirun -np 32 python3 IsentropicGas_1d_pwc.py htmlplot=1
   ```
   For sinusoidal case:
   ```
   mpirun -np 32 python3 IsentropicGas_1d_sin.py htmlplot=1
   ```
