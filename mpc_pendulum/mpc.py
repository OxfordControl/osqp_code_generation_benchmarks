from __future__ import print_function
from __future__ import division
import numpy as np
import scipy.sparse as spa
from builtins import range
import os

# Import subprocess to run matlab script
from subprocess import call

# Import scipy io to write/read mat file
import scipy.io as io

# For importing python modules from string
import importlib

# Import osqp
import osqp

# Import qpoases
# import qpoases as qpoases

# Plotting
import matplotlib.pylab as plt
plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=15)   # fontsize of the tick labels
plt.rc('ytick', labelsize=15)   # fontsize of the tick labels
plt.rc('legend', fontsize=15)   # legend fontsize
plt.rc('text', usetex=True)     # use latex
plt.rc('font', family='serif')

colors = {'b': '#1f77b4',
          'g': '#2ca02c',
          'o': '#ff7f0e',
          'r': '#d62728'}

'''
SYSTEM DESCRIPTION
'''
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


# Iterations
# import tqdm


class QPmatrices(object):
    """
    QP problem matrices

    q_vecs is the matrix containing different linear costs
    """
    def __init__(self, P, q,
                 A, l, u, N, lx=None, ux=None):
        self.P = P
        self.q = q
        self.A = A
        self.l = l
        self.u = u
        self.N = N
        self.lx = lx
        self.ux = ux


class Statistics(object):
    """
    Solve statistics
    """
    def __init__(self, x):
        self.x = x
        self.avg = np.mean(x)
        self.median = np.median(x)
        self.max = np.max(x)
        self.min = np.min(x)


def b(x, n, N):
    b = np.zeros((N+1)*n)
    b[:n] = -x
    return b


def gen_qp_matrices(N, x0, version):
    """
    Generate QP matrices for portfolio optimization problem
    """

    # Construct the MPC matrices
    if version == 'sparse':
        # Objective
        Hx = spa.kron(spa.eye(N), obj_Q)
        Hu = spa.kron(spa.eye(N), obj_R)
        H = spa.block_diag([Hx, obj_QN, Hu])

        # Dynamics
        Mx = spa.kron(spa.eye(N+1), -spa.eye(nx)) + \
            spa.kron(spa.eye(N+1, k=-1), dyn_A)
        Mu = spa.kron(spa.vstack([spa.csc_matrix((1, N)), spa.eye(N)]), dyn_B)
        M = spa.hstack([Mx, Mu])

        # Constraints
        Gx = spa.kron(spa.eye(N), constr_E)
        Gu = spa.kron(spa.eye(N), constr_D)
        G = spa.block_diag([Gx, constr_F, Gu])
        gl = np.hstack([np.tile(constr_el, N), constr_fl,
                        np.tile(constr_dl, N)])
        gu = np.hstack([np.tile(constr_eu, N), constr_fu,
                        np.tile(constr_du, N)])

        # QP data
        P = H
        A = spa.vstack([M, G])
        q = np.zeros((N+1)*nx + N*nu)
        l = np.hstack([b(x0, nx, N), gl])
        u = np.hstack([b(x0, nx, N), gu])

        # Store data
        qp_matrices = QPmatrices(P, q, A, l, u, N)

    # Return QP matrices
    return qp_matrices


def solve_loop(qp_matrices, x0, nsim, solver='emosqp'):
    """
    Solve portfolio optimization loop for all gammas
    """

    # Shorter name for qp_matrices
    qp = qp_matrices

    print('\nSolving MPC pendulum problem loop for N = %d and solver %s' %
          (qp.N, solver))

    # Initialize time vector
    time = np.zeros(nsim)

    # Initialize number of iterations vector
    niter = np.zeros(nsim)

    if solver == 'emosqp':
        # Pass the data to OSQP
        m = osqp.OSQP()
        m.setup(qp.P, qp.q, qp.A, qp.l, qp.u, eps_rel=1e-2, eps_abs=1e-2,
                rho=1e-1, sigma=1e-5, alpha=1.6, verbose=False, max_iter=3000)

        # Get extension name
        module_name = 'emosqpn%s' % str(qp.N)

        # Generate the code
        m.codegen("code", python_ext_name=module_name, force_rewrite=True)

        # Import module
        emosqp = importlib.import_module(module_name)

        for i in range(nsim):

            # Solve
            x, y, status, niter[i], time[i] = emosqp.solve()

            # Check if status correct
            if status != 1:
                raise ValueError('OSQP did not solve the problem!')

            # Apply first control input to the plant
            u = x[-qp.N*nu:-(qp.N-1)*nu]
            x0 = dyn_A.dot(x0) + dyn_B.dot(u)

            # DEBUG
            # solve with gurobi
            # import mathprogbasepy as mpbpy
            # prob = mpbpy.QuadprogProblem(qp.P, q, qp.A, qp.l, qp.u)
            # res = prob.solve(solver=mpbpy.GUROBI)

    # elif solver == 'qpoases':
    #
    #     n_dim = qp.P.shape[0]
    #     m_dim = qp.A.shape[0]
    #
    #     # Initialize qpoases and set options
    #     qpoases_m = qpoases.PyQProblem(n_dim, m_dim)
    #     options = qpoases.PyOptions()
    #     options.printLevel = qpoases.PyPrintLevel.NONE
    #     qpoases_m.setOptions(options)
    #
    #     for i in range(n_prob):
    #
    #         # Get linera cost as contiguous array
    #         q = np.ascontiguousarray(qp.q_vecs[:, i])
    #
    #         # Reset cpu time
    #         qpoases_cpu_time = np.array([20.])
    #
    #         # Reset number of of working set recalculations
    #         nWSR = np.array([1000])
    #
    #         if i == 0:
    #             res_qpoases = qpoases_m.init(qp.P, q, qp.A, qp.lx, qp.ux,
    #                                          qp.l, qp.u,
    #                                          nWSR, qpoases_cpu_time)
    #         else:
    #             # Solve new hot started problem
    #             res_qpoases = qpoases_m.hotstart(q, qp.lx, qp.ux,
    #                                              qp.l, qp.u, nWSR,
    #                                              qpoases_cpu_time)
    #
    #         # # DEBUG Solve with gurobi
    #         # qpoases solution
    #         # sol_qpoases = np.zeros(qp.n)
    #         # qpoases_m.getPrimalSolution(sol_qpoases)
    #         # import mathprogbasepy as mpbpy
    #         # Agrb = spa.vstack((spa.csc_matrix(qp.A),
    #         #                    spa.eye(qp.n))).tocsc()
    #         # lgrb = np.append(qp.l, qp.lx)
    #         # ugrb = np.append(qp.u, qp.ux)
    #         # prob = mpbpy.QuadprogProblem(spa.csc_matrix(qp.P), q,
    #         #                              Agrb, lgrb, ugrb)
    #         # res = prob.solve(solver=mpbpy.GUROBI, verbose=True)
    #         # print("Norm difference x qpoases - GUROBI = %.4f" %
    #         #       np.linalg.norm(sol_qpoases - res.x))
    #         # print("Norm difference objval qpoases - GUROBI = %.4f" %
    #         #       abs(qpoases_m.getObjVal() - res.obj_val))
    #
    #         if res_qpoases != 0:
    #             raise ValueError('qpoases did not solve the problem!')
    #
    #         # Save time
    #         time[i] = qpoases_cpu_time[0]
    #
    #         # Save number of iterations
    #         niter[i] = nWSR[0]

    else:
        raise ValueError('Solver not understood')

    # Return statistics
    return Statistics(time), Statistics(niter)


'''
Solve problems
'''
# Generate gamma parameters and cost vectors
nsim = 100

# Prediction horizon
N_vec = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
x_init = np.array([0., 0., 0.15, 0.])

# Define statistics for osqp and qpoases
osqp_timing = []
osqp_iter = []
# qpoases_timing = []
# qpoases_iter = []


for i in range(len(N_vec)):

    # # Generate QP dense matrices
    # qp_matrices_dense = gen_qp_matrices(k_vec[i], n_vec[i],
    #                                     gammas, 'dense')
    #
    # # Solving loop with qpoases
    # timing, niter = solve_loop(qp_matrices_dense, 'qpoases')
    # qpoases_timing.append(timing)
    # qpoases_iter.append(niter)

    # Generate QP sparsematrices
    qp_matrices_sparse = gen_qp_matrices(N_vec[i], x_init, 'sparse')

    # Solve loop with emosqp
    timing, niter = solve_loop(qp_matrices_sparse, x_init, nsim, 'emosqp')
    osqp_timing.append(timing)
    osqp_iter.append(niter)


# Plot timings
osqp_avg = np.array([x.avg for x in osqp_timing])
# qpoases_avg = np.array([x.avg for x in qpoases_timing])

plt.figure()
ax = plt.gca()
plt.semilogy(N_vec, osqp_avg, color=colors['b'], label='OSQP')
# plt.semilogy(n_vec, qpoases_avg, color=colors['o'], label='qpOASES')
plt.legend()
plt.grid()
ax.set_xlabel(r'Prediction horizon $N$')
ax.set_ylabel(r'Time [s]')
plt.tight_layout()
plt.show(block=False)
plt.savefig('results.pdf')
