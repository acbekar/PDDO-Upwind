import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size' : 14})

"""
Inviscid Burgers - PDDO - Explicit - Upwind
Ali Can Bekar
2021
"""

"""
1D PD Derivative Stencil Generator for Given Horizon Radius.
Generates the Derivative Operators for u>0 and u<0. Since the
Discretization is Uniform and the Boundary Conditions are Periodic,
Created Stencil can be used for all Material Points without Modification.
"""
##################################################################################
##################################################################################
# G Function Generator for 1D
def InvGFunc(xi, delta, amat, order):
    GFS = np.zeros(order + 1)
    Avec = np. zeros(order + 1)
    for i in range(order + 1):
        Avec[i] = pow(xi, i)
    w = np.exp(-4.0 * (np.abs(xi) / delta) ** 2)
    for i in range(order + 1):
        C = amat[:, i] * Avec
        C *= w
        GFS[i] = np.sum(C)
    return GFS        
# A Matrix Generator
def GenAmat(xi, delta, order):
    Avec = np. zeros(order + 1)
    for i in range(order + 1):
        Avec[i] = pow(xi, i)
    w = np.exp(-4.0 * (np.abs(xi) / delta) ** 2)
    Amat = np.outer(Avec, Avec) * w
    return Amat
# b Matrix Generator
def Genbmat(order):
    bmat = np.eye(order + 1)
    for i in range(1, order + 1):
        bmat[i, i] = bmat[i - 1, i - 1] * i
    return bmat
# Derivative Operator Using PDDO
def GenDVec(dx, delta, order):
    nump = round(delta) + 1
    # Positive and Negative Flux Stencils
    # are Created Sequentially. Discretization
    # is Uniform therefore No Family Search is
    # Needed.
    ptsp = np.arange(0, nump)
    ptsn = np.arange(-nump + 1, 1)
    pts = [ptsn, ptsp]
    DVecs = [[], []]
    # Dvecs Include Derivative Stencils up to 1st order.(0th and 1st).
    # First Order Derivative is not Needed for this Problem.
    for i in range(order + 1):
        DVecs[0].append(np.zeros([round(delta)+1]))
        DVecs[1].append(np.zeros([round(delta)+1]))
    for k in range(2):
        bmat = Genbmat(order)
        Amat = np.zeros([order + 1, order + 1])
        for num in pts[k]:
            xi = num * dx
            Amat += GenAmat(xi, delta * dx, order)
        Amat *= dx
        amat = np.linalg.solve(Amat, bmat)
        for kk in range(len(ptsp)):
            xi = pts[k][kk] * dx
            for derord in range(order + 1):
                # Multiplying with dx Beforehand for Integration.
                DVecs[k][derord][kk] = InvGFunc(xi, delta * dx, amat, order).T[derord] * dx
    return DVecs
##################################################################################
##################################################################################


# Main Loop for Explicit Integration
##################################################################################
##################################################################################
def EeInt(u, dt, dx, delta, numtstp, numpt, DerVec):
    # Copying Field Variables for Integration at
    # Every Timestep.
    u_prev = u[:]
    nump = round(delta) + 1
    ptsp = np.arange(0, nump)
    ptsn = np.arange(-nump + 1, 1)
    pts = []
    for i in range(numtstp):
        for j in range(numpt):
            # Checking the Flux Direction.
            if u[j] > 0.0:
                # Mod is used to Enforce Periodicity.
                pts = np.mod(j + ptsn, numpt)
                # Explicit Integration.
                u[j] -= (dt / 2.0) * np.dot(u_prev[pts] * u_prev[pts], DerVec[0][1])
            if u[j] < 0.0:
                pts = np.mod(j + ptsp, numpt)
                u[j] -= (dt / 2.0) * np.dot(u_prev[pts] * u_prev[pts], DerVec[1][1])
        # Copying Field Variables for the Next Timestep.        
        u_prev = u[:]
    return u
##################################################################################
##################################################################################

# Input Values are the same as the Inviscid Burgers' Equation
# Given in the Paper.
numpt = 201
xmin = 0.0
xmax = 2.0
dx = (xmax - xmin) / (numpt - 1)
delta = 2.015
order = 1
tmax = 0.6
dt = 0.001
numtstp = round(tmax / dt)
x = np.linspace(xmin, xmax, numpt)
# Initial Condition.
u = np.sin(np.pi * x)
DV = GenDVec(dx, delta, order)
u = EeInt(u, dt, dx, delta, numtstp, numpt, DV)

# Plotting The Result.
fig = plt.figure(1, figsize=(5, 5))
plt.plot(x, u, marker='s', markersize=4, mfc='white',
        linestyle='-', alpha=0.8, label='$t=0.6$')
plt.xlabel('x(m)')
plt.ylabel('u')
plt.legend(loc='upper right')
plt.show()
