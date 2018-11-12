# [MAIN5] Polytech Sorbonne: GASPART Project
# Autumn 2018
# ===========================================
# I. Bentoumi, A. Khizar, V. Nicol, M. Pecheux, S. Sleiman
#
# Discrete approximation of the heat equation
# -------------------------------------------
import sys
import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.pylab import subplot2grid


# DESCRIPTION
# -----------
# u = approximation:
# 2-D array: (nb_frames, 1 (for central point) + (P-1)*Q points)

# implicit scheme
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
k = 1                           # heat diffusion coefficient
max_time      = 30              # number of frames to simulate
Thot, Tborder = 1., 0.          # temperatures at the hotspot
                                # .. and on the borders of the particle
                                # .. (for Dirichlet boundary condition)
TP = 0.1                        # external temperature (for phantom points)
border_condition = 'Dirichlet'    # type of boundary condition:
                                # .. 'Dirichlet' or 'Neumann'
m = 0.
sigma = 0.1

# discretization variables
P       = 20 # radius: nb of discretization points
Q       = 30 # angle: nb of discretization points
dr      = R / P
dtheta  = 2*np.pi / Q
dt      = 0.5 / (k/(dr**2) + k/(dtheta**2)) # to verify the CFL condition

# visualization variables
EXPORT    = False
ANIM_RATE = 300
tx, ty    = 7.*R/10., 8.*R/10.

# FUNCTIONS
# ---------
# util function to easily map between (i,j) coordinates and array index
idx = lambda p, q: 1 + (p-1)*Q + ((q-1) % Q) if p > 0 else 0 # (skip first point in middle)

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
        
def M():
    t    = k*dt
    Ai   = lambda i: 1 + 2*t/(dr**2) + 2*t/((i*dr)**2 * (dtheta**2))
    Bi   = lambda i: -t/((i*dr)**2 * (dtheta**2))
    Ci_p = lambda i: -t/(dr**2) * (1 + 1/(2*i))
    Ci_m = lambda i: -t/(dr**2) * (1 - 1/(2*i))
    
    M = np.zeros((1 + (P-1)*Q, 1 + (P-1)*Q)) # square matrix
    # central point: centered scheme
    M[0,0]      = 1. + 4*t/(dr**2)
    M[0,1]      = -2*t/(dr**2)
    M[0,Q//2+1] = -2*t/(dr**2)
    # specific treatment for first ring (p = 1)
    for q in range(1, Q+1):
        M[q,0]   = Ci_m(1)
        M[q,q]   = Ai(1)
        if q == 1: M[q,Q]   = Bi(1)   # modulo on q indices
        else:      M[q,q-1] = Bi(1)
        if q == Q: M[q,1]   = Bi(1)
        else:      M[q,q+1] = Bi(1)
        M[q,q+Q] = Ci_p(1)
    # middle points
    for p in range(2, P-1):
        for tmp_q in range(1, Q+1):
            q = tmp_q + (p-1)*Q           # get real row index
            M[q,q-Q] = Ci_m(p)
            M[q,q]   = Ai(p)
            if tmp_q == 1: M[q,q+Q-1] = Bi(p)   # modulo on q indices
            else:          M[q,q-1]   = Bi(p)
            if tmp_q == Q: M[q,1+(p-1)*Q] = Bi(p)
            else:          M[q,q+1]       = Bi(p)
            M[q,q+Q] = Ci_p(p)
    # specific treatment for last ring (p = P-1)
    for tmp_q in range(1, Q+1):
        q = tmp_q + (P-2)*Q           # get real row index
        M[q,q-Q] = Ci_m(P-1)
        M[q,q]   = Ai(P-1)
        if tmp_q == 1: M[q,q+Q-1] = Bi(P-1)   # modulo on q indices
        else:          M[q,q-1]   = Bi(P-1)
        if tmp_q == Q: M[q,1+(P-2)*Q] = Bi(P-1)
        else:          M[q,q+1]       = Bi(P-1)
                
    M = sparse.csr_matrix(M)
    return M

# function to initialize the particle
# (for now, sets a simple Gaussian)
def init():
    u = np.zeros((max_time, 1 + (P-1)*Q))
    # init point with Gaussian values
    c = 1./(np.sqrt(2*np.pi)*sigma)
    u[0,0] = c
    for p in range(1, P):
        for q in range(1, Q+1):
            u[0,idx(p, q)] = c * np.exp(-((p*dr)-m)**2/(2*sigma**2))

    return u

# function to update the approximation array
def update(u, M, t):
    next_u = np.zeros_like(u)
    
    # first step: solve linear system for middle points
    next_u[:] = sparse.linalg.spsolve(M, u)

    # special treatment: center point
    # next_u[0] = u[0] # no change
    # next_u[0] = np.mean(next_u[idx(1,1):idx(1,Q)]) # take mean of first circle
    # next_u[0] = (u[idx(1,1)] - u[0] + u[idx(1,Q//2+1)])/(dr**2)
    # next_u[0] = 1. / np.sqrt(2*np.pi*(sigma**2 + 2*k*t*dt)) # exact solution

    # border points
    for q in range(1, Q+1):
        next_u[idx(P-1, q)] = border_func(next_u[idx(P-2, q)])

    return next_u

# function to update the exact solution u
def update_exact(u, t):
    next_u_exact = np.zeros_like(u)

    c = 1. / np.sqrt(2*np.pi*(sigma**2 + 2*k*t*dt))
    for p in range(0, P):
        next_u_exact[idx(p, 0)] = np.exp(-(p*dr-m)**2 / (2*(sigma**2 + 2*k*t*dt)))
        next_u_exact[idx(p, 0)] *= c
        
    return next_u_exact


# function to compute and output the simulation
def simulate(t, u, M, u_exact, text, p01, p02, p03):
    # compute new simulation iteration
    if t > 0:
        u[t]       = update(u[t-1], M, t)
        u_exact[t] = update_exact(u_exact[t-1], t)

    # plot new image
    data               = np.zeros((len(u[t]), 3))
    data[0, -1]        = u[t, 0]
    xx, yy             = [0], [u[t, 0]]
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
u       = init()
u_exact = init()
M       = M()

# prepare figure
fig  = plt.figure(figsize=(14, 6))
ax01 = subplot2grid((1, 2), (0, 0))
ax02 = subplot2grid((1, 2), (0, 1))
ax01.set_xlim((-R, R))
ax01.set_ylim((-R, R))
ax02.set_xlim((0, P-1))
ax02.set_ylim((-np.max(u), np.max(u)))

# set plots
p011  = ax01.text(tx, ty, 'Iter. 0', bbox=dict(facecolor='white'), label='p011')
p012  = ax01.scatter([], [], c=[], cmap='coolwarm', label='p012', vmin=0., vmax=Thot)
p021, = ax02.plot([], [], 'b--', label='p021')
p022, = ax02.plot([], [], 'r-', label='p022')
ax02.legend(['Numerical scheme', 'Exact solution'])

# create final animation
im_ani = animation.FuncAnimation(
    fig, simulate, fargs=(u, M, u_exact, p011, p012, p021, p022), frames=max_time,
    interval=ANIM_RATE, blit=False
)
# output anim (to file, or directly)
if EXPORT: im_ani.save('anim.gif', writer='imagemagick')
else: plt.show()
