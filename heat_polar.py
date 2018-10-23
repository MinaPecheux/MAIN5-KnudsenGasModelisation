# [MAIN5] Polytech Sorbonne: GASPART Project
# Autumn 2018
# ===========================================
# I. Bentoumi, A. Khizar, V. Nicol, M. Pecheux, S. Sleiman
#
# Discrete approximation of the heat equation
# -------------------------------------------
import sys
import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.pylab import subplot2grid


# DESCRIPTION
# -----------
# u = approximation:
# 2-D array: (nb_frames, 1 (for central point) + (P-1)*Q points)

# explicit scheme : ( u^n+1_i,j - u^n_i,j )/dt
# centered : ( u^n_i+1,j - 2*u^n_i,j + u^n_i-1,j )/dx^2

# This script plots 2 animated graphs:
#  - a top-view of the particle with polar-discretized space points
# (i.e.: concentric rings of points!)
#  - a side-view of one particular particle radius: we only consider
# the angle theta = 0 and we look at the behavior of each u_(i,0) on this line
# through time


# VARIABLES
# ---------
# problem variables
R = 1                           # radius of the particle
l = 1                           # heat diffusion coefficient
max_time      = 10              # number of frames to simulate
Thot, Tborder = 1., 0.          # temperatures at the hotspot
                                # .. and on the borders of the particle
                                # .. (for Dirichlet boundary condition)
TP = 0.1                        # external temperature (for phantom points)
border_condition = 'Neumann'    # type of boundary condition:
                                # .. 'Dirichlet' or 'Neumann'
k = 1.

# discretization variables
P       = 10 # radius: nb of discretization points
Q       = 30 # angle: nb of discretization points
dr      = R / P
dtheta  = 2*np.pi / Q
dt      = 0.5 / (l/(dr**2) + l/(dtheta**2)) # to verify the CFL condition

# visualization variables
EXPORT    = False
ANIM_RATE = 400
tx, ty    = 7.*R/10., 8.*R/10.

# FUNCTIONS
# ---------
# util function to easily map between (i,j) coordinates and array index
idx     = lambda p, q: 1 + (p-1)*Q + ((q-1) % Q) if p > 0 else 0 # (skip first point in middle)

# function to compute the interactions with the environment
def g():
    return 1.

# function to compute the value of the particle's border points
# (Dirichlet or Neumann boundary conditions)
def border_func(x):
    if border_condition == 'Dirichlet': return Tborder
    elif border_condition == 'Neumann':
        return TP - g()*dr        # with arbitrary phantom point
        # return g()*dr - x       # with interior point
        # return (TP - x)/(2*dr)  # centered?
    else:
        print('[Error] Unknown border condition.')
        sys.exit(1)

# function to initialize the particle
# (for now, sets a simple Gaussian)
def init():
    u = np.zeros((max_time, 1 + (P-1)*Q))
    # init point with Gaussian values
    c = 1.
    u[0,0] = c
    for p in range(1, P):
        for q in range(1, Q+1):
            # u[0,idx(p, q)] = c * np.exp(-(p*dr*np.cos(q*dtheta))**2)
            u[0,idx(p, q)] = c * np.exp(-k*(p*dr)**2)

    return u

# function to update the approximation array
def update(u, t):
    next_u = np.zeros_like(u)

    # middle points
    for p in range(1, P - 1):
        for q in range(1, Q+1):
            A = (u[idx(p+1,q)] - 2*u[idx(p,q)] + u[idx(p-1,q)]) / (dr**2) \
                + (u[idx(p+1,q)] - u[idx(p-1,q)]) / (2*dr**2*p) \
                + (u[idx(p,q+1)] - 2*u[idx(p,q)] + u[idx(p,q-1)]) / ((p*dr)**2 * dtheta**2)
            next_u[idx(p, q)] = u[idx(p, q)] + l * dt * A
    # special treatment: center point
    # take mean of first circle
    # next_u[0] = u[0]
    # next_u[0] = np.mean(next_u[idx(1,1):idx(1,Q)])
    #next_u[0] = 2.*(u[idx(1,1)] - 2*u[0] + u[idx(1,Q//2)])/(dr**2)
    integrand = lambda x, t, p: np.exp(-(p - x)**2 / (4*t*dt) - k*x**2)
    next_u[0] = integrate.quad(
        integrand, -np.inf, np.inf, args=(t, 0)
    )[0]
    next_u[0] *= 1 / np.sqrt(4*np.pi*t*dt)

    # border points
    for q in range(1, Q+1):
        next_u[idx(P-1, q)] = border_func(next_u[idx(P-2, q)])

    return next_u

# function to update the exact solution u
def update_exact(u, t):
    next_u_exact = np.zeros_like(u)

    integrand = lambda x, t, p: np.exp(-(p - x)**2 / (4*t*dt) - k*x**2)
    for p in range(0, P):
        next_u_exact[idx(p, 0)] = integrate.quad(
            integrand, -np.inf, np.inf, args=(t, p)
        )[0]
        next_u_exact[idx(p, 0)] *= 1 / np.sqrt(4*np.pi*t*dt)

    # # middle points
    # for p in range(0, P):
    #     # approximation of the exact solution with a quadrature method
    #     next_u_exact[idx(p, 0)] = 1
    #     for yi in np.arange(-1000,1000, 0.1):
    #         next_u_exact[idx(p, 0)] += np.exp( - (p - yi)**2 / (4*t*dt) - yi**2 ) * 0.1
    #     next_u_exact[idx(p, 0)] *= 1 / np.sqrt(4*np.pi*t*dt)
        
    return next_u_exact


# function to compute and output the simulation
def simulate(t, u, u_exact, text, p01, p02, p03):
    # compute new simulation iteration
    if t > 0:
        u[t] = update(u[t-1], t)
        u_exact[t] = update_exact(u_exact[t-1], t)

    # plot new image
    data = np.zeros((len(u[t]), 3))
    data[0, -1] = u[t, 0]
    xx, yy = [0], [u[t, 0]]
    xx_exact, yy_exact = [0], [u_exact[t, 0]]
    i = 0
    for p in range(1, P):
        for q in range(1, Q+1):
            data[idx(p, q),0] = (p * dr) * np.cos(q * dtheta)
            data[idx(p, q),1] = (p * dr) * np.sin(q * dtheta)
            data[idx(p, q),2] = u[t, idx(p, q)]
            i += 1
        xx.append(p)
        yy.append(u[t, idx(p, 0)])
        
        xx_exact.append(p)
        yy_exact.append(u_exact[t, idx(p, 0)])

    # update figures
    p01.set_offsets(data[:,:2])
    p01.set_array(data[:,2])
    p01.set_clim(np.min(data[:,2]), np.max(data[:,2]))
    p02.set_data((xx, yy))
    p03.set_data((xx_exact, yy_exact))
    text.set_text('Iter. %d' % (t+1))

    return u, u_exact, text, p01, p02, p03,

# MAIN ROUTINE
# ------------
# init grid
u = init()
u_exact = init()

# prepare figure
fig = plt.figure(figsize=(14, 6))
ax01 = subplot2grid((1, 2), (0, 0))
ax02 = subplot2grid((1, 2), (0, 1))
ax01.set_xlim((-R, R))
ax01.set_ylim((-R, R))
ax02.set_xlim((0, P-1))
ax02.set_ylim((-np.max(u), np.max(u)))

# set plots
p011 = ax01.text(tx, ty, 'Iter. 0', bbox=dict(facecolor='white'), label='p011')
p012 = ax01.scatter([], [], c=[], cmap='coolwarm', label='p012', vmin=0., vmax=Thot)
p021, = ax02.plot([], [], 'b--', label='p021')
p022, = ax02.plot([], [], 'r-', label='p022')
ax02.legend(['Numerical scheme', 'Exact solution'])

# create final animation
im_ani = animation.FuncAnimation(
    fig, simulate, fargs=(u, u_exact, p011, p012, p021, p022), frames=max_time,
    interval=ANIM_RATE, blit=False
)
# output anim (to file, or directly)
if EXPORT: im_ani.save('anim.gif', writer='imagemagick')
else: plt.show()
