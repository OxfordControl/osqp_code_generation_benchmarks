from __future__ import print_function
import numpy as np
import scipy.sparse as spa
import osqp
from builtins import range


# Model of the system
dyn_A = spa.csc_matrix([[1., 0.02, 0., 0.], [0., 1., 0., 0.],
                        [0., 0., 1.0049, 0.02], [0., 0., 0.4913, 1.0049]])
dyn_B = spa.csc_matrix([[0.0002], [0.02], [-0.0005], [-0.0501]])
[nx, nu] = dyn_B.shape

# Objective function
obj_Q = spa.diags([1., 0.3, 0.3, 0.1])
obj_R = spa.diags([0.1])
# QN = dare(A,B,Q,R)
obj_QN = spa.csc_matrix([[69.4161, 39.9917, 111.6822, 22.3337],
                         [39.9917, 43.8824, 131.5453, 26.2255],
                         [111.6822, 131.5453, 618.2991, 121.9162],
                         [22.3337, 26.2255, 121.9162, 24.6494]])

# Input constraints
constr_D = spa.diags([1.0])
constr_dl = -5.0
constr_du = 5.0

# State constraints
constr_E = spa.eye(nx)
constr_el = np.array([-0.5, -1.0, -0.2, -0.5])
constr_eu = np.array([0.5, 1.0, 0.2, 0.5])

# Terminal state constraints
constr_F = constr_E
constr_fl = constr_el
constr_fu = constr_eu

# Prediction horizon
N = 10

'''
MPC Matrices
'''
cond_A = spa.csc_matrix(((N+1)*nx, nx))
cond_A[:nx, :] = spa.eye(nx)
for i in range(N):
    cond_A[(i+1)*nx:(i+2)*nx, :] = cond_A[i*nx:(i+1)*nx, :].dot(dyn_A)

B_tmp = dyn_B
cond_B = spa.csc_matrix(((N+1)*nx, N*nu))
for i in range(N):
    cond_B += spa.kron(spa.vstack([spa.csc_matrix((1, N)),
                       spa.eye(N, k=-i)]), B_tmp)
    B_tmp = dyn_A.dot(B_tmp)

cond_R = spa.kron(spa.eye(N), obj_R)
cond_Q = spa.block_diag([spa.kron(spa.eye(N), obj_Q), obj_QN])

u_low = np.tile(constr_dl, N)
u_upp = np.tile(constr_du, N)
x_low = np.hstack([np.tile(constr_el, N), constr_fl])
x_upp = np.hstack([np.tile(constr_eu, N), constr_fu])

# Pass the data to OSQP
x = np.array([0., 0., 0.15, 0.])  # initial state: [p, p_dot, theta, theta_dot]
P = cond_B.T.dot(cond_Q).dot(cond_B) + cond_R
q = cond_B.T.dot(cond_Q).dot(cond_A).dot(x)
A = spa.vstack([spa.eye(N*nu), cond_B])
l = np.hstack([u_low, x_low - cond_A.dot(x)])
u = np.hstack([u_upp, x_upp - cond_A.dot(x)])

m = osqp.OSQP()
m.setup(P, q, A, l, u, eps_rel=1e-3, eps_abs=1e-3,
        rho=1e-1, sigma=1e-5, alpha=1.95, max_iter=3000)

# Generate the code
m.codegen("code", project_type="Makefile", embedded=1,
          python_ext_name='emosqp', force_rewrite=True)

'''
Apply MPC in closed loop
'''

import emosqp

# Apply MPC to the system
sim_steps = 20

for i in range(sim_steps):

    # Solve
    sol = emosqp.solve()
    u = sol[0][(N+1)*nx:]
    status_val = sol[2]
    numofiter = sol[3]
    runtime = 1000*sol[4]

    print("iter = %4d,    runtime = %7.4f ms" % (numofiter, runtime))

    # Compute details
    if status_val != 1:
        print('Problem unsolved')
        break

    # Apply first control input to the plant
    x = dyn_A.dot(x) + dyn_B.dot(u)

    # Update initial state
    q_new = cond_B.T.dot(cond_Q).dot(cond_A).dot(x)
    l_new = np.hstack([u_low, x_low - cond_A.dot(x)])
    u_new = np.hstack([u_upp, x_upp - cond_A.dot(x)])
    emosqp.update_lin_cost(q_new)
    emosqp.update_bounds(l_new, u_new)
