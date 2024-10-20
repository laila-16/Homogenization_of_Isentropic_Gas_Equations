#!/usr/bin/env python
# encoding: utf-8

r"""
One-dimensional isothermalgas
=============================

Solve the isothermalgas equations with varying-cross section:

.. math::
    a\rho_t + (q)_x & = 0 \\
    q_t + (q^2 /(a\rho) +a p(\rho))_x & = a_x p(\rho).

Here p is the pressure where p=kappa \rho^gamma, q is the flux, a is the cross section area,
and :math:`\rho` is the density. 

The initial condition is a Gaussian and the boundary conditions are periodic.
The final solution is identical to the initial data because both waves have
crossed the domain exactly once.
"""
from __future__ import absolute_import
import numpy as np
from clawpack import riemann
import RS_variable_coeff

def setup(use_petsc=True, kernel_language='Fortran', solver_type='classic',
          outdir='./_output', ptwise=False, weno_order=5, order=2,
          time_integrator='SSP104', disable_output=False, output_style=1,
          L=1600, mx=1000000, bc='wall', tmax=900.0, num_output_times=450, CFL=0.5):

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if kernel_language == 'Fortran':
            riemann_solver = RS_variable_coeff

    if solver_type == 'classic':
        solver = pyclaw.ClawSolver1D(riemann_solver)
        solver.limiters = pyclaw.limiters.tvd.MC
        solver.order = order
        solver.fwave=True
        solver.step_source=classic_source_step

    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver1D(riemann_solver)
        solver.weno_order = weno_order
        solver.time_integrator = time_integrator
        if time_integrator == 'SSPLMMk3':
            solver.lmm_steps = 4
        solver.dq_src = sharpclaw_source_step
    else:
        raise Exception('Unrecognized value of solver_type.')

    solver.kernel_language = kernel_language


    x = pyclaw.Dimension(0., L/2., mx, name='x')
    domain = pyclaw.Domain(x)
    num_eqn = 2
    state = pyclaw.State(domain, num_eqn,num_aux=1)

    #Auxiliary vector will contain a(x)
    if bc=='periodic':
        solver.aux_bc_lower[0] = pyclaw.BC.periodic
        solver.aux_bc_upper[0] = pyclaw.BC.periodic

        solver.bc_lower[0] = pyclaw.BC.periodic
        solver.bc_upper[0] = pyclaw.BC.periodic
        
    elif bc=='wall':
        solver.aux_bc_lower[0] = pyclaw.BC.wall
        solver.aux_bc_upper[0] = pyclaw.BC.wall

        solver.bc_lower[0] = pyclaw.BC.wall
        solver.bc_upper[0] = pyclaw.BC.wall
        
    else: 
        solver.aux_bc_lower[0] = pyclaw.BC.extrap
        solver.aux_bc_upper[0] = pyclaw.BC.extrap

        solver.bc_lower[0] = pyclaw.BC.extrap
        solver.bc_upper[0] = pyclaw.BC.extrap

    solver.num_waves = 2
    solver.num_eqn = 2
    
    solver.cfl_desired = CFL
    solver.cfl_max = 1.0

    
    xc = domain.grid.x.centers
    state.index_capa = 0
    init(state, xc)

    solver.max_steps=5000000

    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state, domain)
    claw.keep_copy = True
    claw.solver = solver
    claw.outdir = outdir
    claw.output_style = output_style
    claw.write_aux_init = True
    claw.tfinal = tmax
    claw.num_output_times = num_output_times
    claw.setplot = setplot

    return claw

def a_function(x):
    return (3/5) + (2/5) *np.sin(2.*np.pi*x)

def a_function_prime(x):
    return (2/5)*2*np.pi*np.cos(2.*np.pi*x)

def init(state, xc):
    A= 1/12
    rho0 = 0.3 + A*np.exp(-1*(xc/5)**2)
    q0 = np.zeros_like(rho0)
    state.q[0, :] = rho0.copy()
    state.q[1, :] = q0.copy()
    #Initialize capacity array
    a_values = []
    state.aux[0,:] = a_function(xc)

def classic_source_step(solver,state,dt):
    q = state.q
    xc = state.grid.x.centers
    
    def S(u1,u2,x):
        return  0.*u1, 0.*u2

    #Get variables
    u1 = np.copy(q[0,:])
    u2 = np.copy(q[1,:])

    #RK4 time integration of source term
    k11,k12 = S(u1,u2,xc)
    k21,k22 = S(u1+0.5*dt*k11,u2+0.5*dt*k12,xc)
    k31,k32 = S(u1+0.5*dt*k21,u2+0.5*dt*k22,xc)
    k41,k42 = S(u1+dt*k31,u2+dt*k32,xc)
    
    q[0,:] = u1+(1./6.)*dt*(k11+2*k21+2*k31+k41)
    q[1,:] = u2+(1./6.)*dt*(k12+2*k22+2*k32+k42)

def sharpclaw_source_step(solver,state,dt):
    q = state.q
    xc = state.grid.x.centers

    #Get variables
    u1 = np.copy(q[0,:])
    u2 = np.copy(q[1,:])

    dq = np.empty(q.shape)
    dq[0,:] = 0.*u1
    dq[1,:] = 0.*u2
    return dq


def setplot(plotdata):
    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """
    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for density
    plotfigure = plotdata.new_plotfigure(name='Solution', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    #plotaxes.xlimits = [0,4.0]
    #plotaxes.ylimits = [0.28, 0.4]
    plotaxes.title = 'Rho'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth': 1, 'markersize': 5}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    #plotaxes.xlimits = [-100.0,100.0]
    #plotaxes.ylimits = [-1., 1.]
    plotaxes.title = 'q'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 1
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth': 1, 'markersize': 5}
    
    return plotdata


def run_and_plot(**kwargs):
    claw = setup(kwargs)
    claw.run()
    from clawpack.petclaw import plot
    #plot.interactive_plot(setplot=setplot) use this command for interactive plots
    plot.html_plot(setplot=setplot) #to create HTML plots

if __name__ == "__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup, setplot)
