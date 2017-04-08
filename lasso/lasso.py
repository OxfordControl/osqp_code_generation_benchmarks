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
import qpoases
import pickle

# Plotting
import matplotlib.pylab as plt
plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=15)   # fontsize of the tick labels
plt.rc('ytick', labelsize=15)   # fontsize of the tick labels
plt.rc('legend', fontsize=15)   # legend fontsize
plt.rc('text', usetex=True)     # use latex
plt.rc('font', family='serif')

colors = { 'b': '#1f77b4',
           'g': '#2ca02c',
           'o': '#ff7f0e',
           'r': '#d62728'}


class QPmatrices(object):
    """
    QP problem matrices

    q_vecs is the matrix containing different linear costs
    """
    def __init__(self, P, q_vecs, A, l, u, n, m):
        self.P = P
        self.q_vecs = q_vecs
        self.A = A
        self.l = l
        self.u = u
        self.n = n
        self.m = m


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


def gen_qp_matrices(m, n, gammas):
    """
    Generate QP matrices for lasso problem
    """
    # Reset random seed for repetibility
    np.random.seed(1)

    # Problem parameters
    dens_lvl = 0.4

    # Generate data
    Ad = spa.random(m, n, density=dens_lvl, format='csc')
    x_true = np.multiply((np.random.rand(n) > 0.5).astype(float),
                         np.random.randn(n)) / np.sqrt(n)
    bd = Ad.dot(x_true) + .5*np.random.randn(m)

    #       minimize	y.T * y + gamma * np.ones(n).T * t
    #       subject to  y = Ax - b
    #                   -t <= x <= t
    P = spa.block_diag((spa.csc_matrix((n, n)), spa.eye(m),
                        spa.csc_matrix((n, n))), format='csc')
    # q = np.append(np.zeros(m + n), gamma*np.ones(n))
    In = spa.eye(n)
    Onm = spa.csc_matrix((n, m))
    A = spa.vstack([spa.hstack([Ad, -spa.eye(m), Onm.T]),
                    spa.hstack([In, Onm, In]),
                    spa.hstack([-In, Onm, In])]).tocsc()
    l = np.hstack([bd, np.zeros(2*n)])
    u = np.hstack([bd, np.inf * np.ones(2*n)])

    # Create linear cost vectors
    q_vecs = np.empty((2*n + m, 0))
    for gamma in gammas:
        q_vecs = np.column_stack(
            (q_vecs, np.append(np.zeros(n+m), gamma*np.ones(n))))
    qp_matrices = QPmatrices(P, q_vecs, A, l, u, n, m)

    # # Save matrices for CVXGEN
    # if n <= 50:
    #     cvxgen_dir = 'cvxgen/n%d/' % n
    #     if not os.path.exists(cvxgen_dir):
    #         os.mkdir(cvxgen_dir)
    #     io.savemat(cvxgen_dir + 'datafilen%d.mat' % n,
    #                {'gammas': gammas,
    #                 'Ad': Ad,
    #                 'bd': bd})
    #
    # # Save matrices for FiOrdOs
    # fiordos_dir = 'fiordos/n%d/' % n
    # if not os.path.exists(fiordos_dir):
    #     os.mkdir(fiordos_dir)
    # io.savemat(fiordos_dir + 'datafilen%d.mat' % n,
    #            {'gammas': gammas,
    #             'Ad': Ad,
    #             'bd': bd})

    # Return QP matrices
    return qp_matrices


def solve_loop(qp_matrices, solver='emosqp'):
    """
    Solve portfolio optimization loop for all gammas
    """

    # Shorter name for qp_matrices
    qp = qp_matrices

    print('\nSolving lasso problem loop for n = %d and solver %s' %
          (qp.n, solver))

    # Get number of problems to solve
    n_prob = qp.q_vecs.shape[1]

    # Initialize time vector
    time = np.zeros(n_prob)

    # Initialize number of iterations vector
    niter = np.zeros(n_prob)

    if solver == 'emosqp':
        # Pass the data to OSQP
        m = osqp.OSQP()
        m.setup(qp.P, qp.q_vecs[:, 0], qp.A, qp.l, qp.u,
                rho=0.01, auto_rho=False, verbose=False)

        # Get extension name
        module_name = 'emosqpn%s' % str(qp.n)

        # Generate the code
        m.codegen("code", python_ext_name=module_name,
                  force_rewrite=True)

        # Import module
        emosqp = importlib.import_module(module_name)

        for i in range(n_prob):
            q = qp.q_vecs[:, i]

            # Update linear cost
            emosqp.update_lin_cost(q)

            # Solve
            x, y, status, niter[i], time[i] = emosqp.solve()

            # Check if status correct
            if status != 1:
                print('OSQP did not solve the problem!')
                import ipdb; ipdb.set_trace()
                # raise ValueError('OSQP did not solve the problem!')

            # DEBUG
            # solve with gurobi
            # import mathprogbasepy as mpbpy
            # prob = mpbpy.QuadprogProblem(qp.P, q, qp.A, qp.l, qp.u)
            # res = prob.solve(solver=mpbpy.GUROBI)

    elif solver == 'qpoases':

        n_dim = qp.P.shape[0]
        m_dim = qp.A.shape[0]

        # Initialize qpoases and set options
        qpoases_m = qpoases.PyQProblem(n_dim, m_dim)
        options = qpoases.PyOptions()
        options.printLevel = qpoases.PyPrintLevel.NONE
        qpoases_m.setOptions(options)


        # Setup matrix P and A
        P = np.ascontiguousarray(qp.P.todense())
        A = np.ascontiguousarray(qp.A.todense())

        for i in range(n_prob):

            # Get linera cost as contiguous array
            q = np.ascontiguousarray(qp.q_vecs[:, i])

            # Reset cpu time
            qpoases_cpu_time = np.array([10.])

            # Reset number of of working set recalculations
            nWSR = np.array([10000])

            if i == 0:
                res_qpoases = qpoases_m.init(P, q, A,
                                             None, None,
                                             qp.l, qp.u,
                                             nWSR, qpoases_cpu_time)
            else:
                # Solve new hot started problem
                res_qpoases = qpoases_m.hotstart(q, None, None,
                                                 qp.l, qp.u, nWSR,
                                                 qpoases_cpu_time)

            # if res_qpoases != 0:
            #     import ipdb; ipdb.set_trace()
            #     raise ValueError('qpoases did not solve the problem!')

            # Save time
            time[i] = qpoases_cpu_time[0]

            # Save number of iterations
            niter[i] = nWSR[0]

    elif solver == 'gurobi':

        n_dim = qp.P.shape[0]
        m_dim = qp.A.shape[0]

        for i in range(n_prob):

            # Get linera cost as contiguous array
            q = np.ascontiguousarray(qp.q_vecs[:, i])

            # solve with gurobi
            import mathprogbasepy as mpbpy
            prob = mpbpy.QuadprogProblem(qp.P, q, qp.A, qp.l, qp.u)
            res = prob.solve(solver=mpbpy.GUROBI, verbose=False)

            # Save time
            time[i] = res.cputime

            # Save number of iterations
            niter[i] = res.total_iter

    else:
        raise ValueError('Solver not understood')

    # Return statistics
    return Statistics(time), Statistics(niter)


'''
Problem parameters
'''
# Generate gamma parameters and cost vectors
n_gamma = 21
gammas = np.logspace(-2, 2, n_gamma)


# Variables
# n_vec = np.array([20, 30, 50])
# n_vec = np.array([100, 200, 300, 400, 500])
n_vec = np.array([10, 20, 30, 50, 80, 100, 150, 200, 250, 300, 350, 400])

# Measurements
m_vec = (10 * n_vec).astype(int)

# Setup if solve with gurobi/qpoases or not
solve_qpoases = False
solve_gurobi = False

# Define statistics for osqp and qpoases
osqp_timing = []
osqp_iter = []
if solve_gurobi:
    gurobi_iter = []
    gurobi_timing = []
if solve_qpoases:
    qpoases_timing = []
    qpoases_iter = []

'''
Solve problems
'''
for i in range(len(n_vec)):

    # Generate QP sparsematrices
    qp_matrices = gen_qp_matrices(m_vec[i], n_vec[i], gammas)

    # Solve loop with emosqp
    timing, niter = solve_loop(qp_matrices, 'emosqp')
    osqp_timing.append(timing)
    osqp_iter.append(niter)


    if solve_qpoases:
        # Solving loop with qpoases
        timing, niter = solve_loop(qp_matrices, 'qpoases')
        qpoases_timing.append(timing)
        qpoases_iter.append(niter)

    if solve_gurobi:
        # Solve loop with gurobi
        timing, niter = solve_loop(qp_matrices, 'gurobi')
        gurobi_timing.append(timing)
        gurobi_iter.append(niter)


if solve_qpoases:
    # Dump qpoases timings
    with open('qpoases/timing.pickle', 'wb') as f:
                pickle.dump({'timing': qpoases_timing,
                             'iter': qpoases_iter}, f)

if solve_gurobi:
    # Dump gurobi timings
    with open('gurobi/timing.pickle', 'wb') as f:
                pickle.dump({'timing': gurobi_timing,
                             'iter': gurobi_iter}, f)


if not solve_qpoases:
    # Load qpoases timings
    with open('qpoases/timing.pickle', 'rb') as f:
        qpoases_structs = pickle.load(f)
    qpoases_timing = qpoases_structs['timing']
    qpoases_iter = qpoases_structs['iter']

if not solve_gurobi:
    # Load gurobi timings
    with open('gurobi/timing.pickle', 'rb') as f:
        gurobi_structs = pickle.load(f)
    gurobi_timing = gurobi_structs['timing']
    gurobi_iter = gurobi_structs['iter']


# Plot timings
osqp_avg = np.array([x.avg for x in osqp_timing])
qpoases_avg = np.array([x.avg for x in qpoases_timing if x.avg < 9])
gurobi_avg = np.array([x.avg for x in gurobi_timing])


plt.figure()
ax = plt.gca()
plt.semilogy(n_vec, osqp_avg, color='C0', label='OSQP')
plt.semilogy(n_vec[:len(qpoases_avg)], qpoases_avg, color='C1', label='qpOASES')
plt.semilogy(n_vec, gurobi_avg, color='C5', label='GUROBI')
plt.legend()
plt.grid()
ax.set_xlabel(r'Number of parameters $n$')
ax.set_ylabel(r'Time [s]')
plt.tight_layout()
plt.show(block=False)
plt.savefig('results.pdf')
