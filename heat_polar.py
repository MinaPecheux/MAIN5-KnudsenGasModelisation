# [MAIN5] Polytech Sorbonne: GASPART Project
# Autumn 2018
# ===========================================
# I. Bentoumi, A. Khizar, V. Nicol, M. Pecheux, S. Sleiman
#
# Discrete approximation of the heat equation
# -------------------------------------------
import sys
import numpy as np
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
max_time      = 50              # number of frames to simulate
Thot, Tborder = 1., 0.          # temperatures at the hotspot
                                # .. and on the borders of the particle
                                # .. (for Dirichlet boundary condition)
TP = 0.1
border_condition = 'Neumann'    # type of boundary condition:
                                # .. 'Dirichlet' or 'Neumann'

# discretization variables
P       = 10 # radius: nb of discretization points
Q       = 30 # angle: nb of discretization points
dr      = R / P
dtheta  = 2*np.pi / Q
dt      = 0.5 / (l/(dr**2) + l/(dtheta**2)) # to verify the CFL condition

# visualization variables
EXPORT    = False
ANIM_RATE = 200
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
    elif border_condition == 'Neumann': return TP - g()*dr
    else:
        print('[Error] Unknown border condition.')
        sys.exit(1)

# function to initialize the particle
# (for now, sets an arbitrary hotspot in the middle of the particle,
# zero everywhere else)
def init():
    u = np.zeros((max_time, 1 + (P-1)*Q))

    # set max temp point in the middle
    u[0,0] = Thot
    # set border points
    for q in range(1, Q+1):
        u[0, idx(P-1, q)] = border_func(u[0, idx(P-2, q)])

    return u

# function to update the approximation array
def update(u):
    next_u = np.zeros_like(u)

    # middle points
    for p in range(1, P - 1):
        for q in range(1, Q+1):
            A = (u[idx(p+1,q)] - 2*u[idx(p,q)] + u[idx(p-1,q)]) / (dr**2) \
                + (u[idx(p+1,q)] - u[idx(p-1,q)]) / (2*dr**2*p) \
                + (u[idx(p,q+1)] - 2*u[idx(p,q)] + u[idx(p,q-1)]) / ((p*dr)**2 * dtheta**2)
            next_u[idx(p, q)] = u[idx(p, q)] + l * dt * A
    # special treatment: middle point
    # take mean of first circle
    next_u[0] = np.mean(next_u[idx(1,1):idx(1,Q)])

    # border points
    for q in range(1, Q+1):
        next_u[idx(P-1, q)] = border_func(next_u[idx(P-2, q)])

    return next_u

# function to compute and output the simulation
def simulate(t, u, text, p01, p02):
    # compute new simulation iteration
    if t > 0:
        u[t] = update(u[t-1])

    # plot new image
    data = np.zeros((len(u[t]), 3))
    data[0, -1] = u[t, 0]
    xx, yy = [0], [u[t, 0]]
    i = 0
    for p in range(1, P):
        for q in range(1, Q+1):
            data[idx(p, q),0] = (p * dr) * np.cos(q * dtheta)
            data[idx(p, q),1] = (p * dr) * np.sin(q * dtheta)
            data[idx(p, q),2] = u[t, idx(p, q)]
            i += 1
        xx.append(p)
        yy.append(u[t, idx(p, 0)])

    # update figures
    p01.set_offsets(data[:,:2])
    p01.set_array(data[:,2])
    p01.set_clim(np.min(data[:,2]), np.max(data[:,2]))
    p02.set_data((xx, yy))
    text.set_text('Iter. %d' % (t+1))

    return u, text, p01, p02,

# MAIN ROUTINE
# ------------
# init grid
u = init()

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
p02, = ax02.plot([], [], label='p02')

# create final animation
im_ani = animation.FuncAnimation(
    fig, simulate, fargs=(u, p011, p012, p02), frames=max_time,
    interval=ANIM_RATE, blit=False
)
# output anim (to file, or directly)
if EXPORT: im_ani.save('anim.gif', writer='imagemagick')
else: plt.show()
