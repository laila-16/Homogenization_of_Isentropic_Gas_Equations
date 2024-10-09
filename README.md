# Homogenization of Isentropic Gas Equations
In this repository, you can find all the codes used in the paper Homogenized Equations for Isentropic Gas in a Pipe with Periodically Varying Cross-Section.
The used commands to run the simulations of finite volume solution:
1. On the terminal we use this command to compile the Reimann Solver:
   ```
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
