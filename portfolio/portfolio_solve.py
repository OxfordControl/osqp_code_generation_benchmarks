from __future__ import print_function
from __future__ import division
import numpy as np
import scipy.sparse as spa
from builtins import range
import os
import pandas as pd

# Import scipy io to write/read mat file
import scipy.io as io

# Import subprocess to run matlab script
from subprocess import call
from platform import system

# For importing python modules from string
import importlib


class QPmatrices(object):
    """
    QP problem matrices

    q_vecs is the matrix containing different linear costs
    """
    def __init__(self, P, q_vecs,
                 A, l, u, n, k, lx=None, ux=None):
        self.P = P
        self.q_vecs = q_vecs
        self.A = A
        self.l = l
        self.u = u
        self.n = n
        self.k = k
        self.lx = lx
        self.ux = ux


def gen_qp_matrices(k, n, gammas):
    """
    Generate QP matrices for portfolio optimization problem
    """
    # Reset random seed for repetibility
    np.random.seed(1)

    # Problem parameters

    dens_lvl = 0.5

    # Generate data
    F = spa.random(n, k, density=dens_lvl, format='csc')
    D = spa.diags(np.random.rand(n) * np.sqrt(k), format='csc')
    mu = np.random.randn(n)

    # Construct the problem
    #       minimize	x' D x + y' I y - (1/gamma) * mu' x
    #       subject to  1' x = 1
    #                   F' x = y
    #                   0 <= x <= 1
    P = spa.block_diag((2*D, 2*spa.eye(k)), format='csc')

    A = spa.vstack([
            spa.hstack([spa.csc_matrix(np.ones((1, n))),
                       spa.csc_matrix((1, k))]),
            spa.hstack([F.T, -spa.eye(k)])
        ]).tocsc()
    l = np.hstack([1., np.zeros(k)])   # Linear constraints
    u = np.hstack([1., np.zeros(k)])
    lx = np.zeros(n)   # Bounds
    ux = np.ones(n)

    # Create linear cost vectors
    q_vecs = np.empty((k + n, 0))
    for gamma in gammas:
        q_vecs = np.column_stack((q_vecs,
                                  np.append(-mu / gamma, np.zeros(k))))
    qp_matrices = QPmatrices(P, q_vecs, A, l, u, n, k, lx, ux)

    # Save matrices for CVXGEN (n<= 120)
    if n <= 120:
        io.savemat(os.path.join('cvxgen', 'n%d' % n, 'datafilen%d.mat' % n),
                   {'gammas': gammas,
                    'F': F,
                    'D': D,
                    'mu': mu})

    # Save matrices for FiOrdOs
    io.savemat(os.path.join('fiordos', 'n%d' % n, 'datafilen%d.mat' % n),
               {'gammas': gammas,
                'F': F,
                'D': D,
                'mu': mu})

    # Return QP matrices
    return qp_matrices


def solve_loop(qp_matrices, solver='emosqp'):
    """
    Solve portfolio optimization loop for all gammas
    """

    # Shorter name for qp_matrices
    qp = qp_matrices

    print('\nSolving portfolio problem loop for n = %d and solver %s' %
          (qp.n, solver))

    # Get number of problems to solve
    n_prob = qp.q_vecs.shape[1]

    # number of assets and factors
    n = qp.n
    k = qp.k

    # Results list
    results = []

    if solver == 'emosqp':
        # Construct qp matrices
        Aosqp = spa.vstack([qp.A,
                            spa.hstack([spa.eye(n), spa.csc_matrix((n, k))])
                            ]).tocsc()
        losqp = np.append(qp.l, qp.lx)
        uosqp = np.append(qp.u, qp.ux)

        # Pass the data to OSQP
        m = osqp.OSQP()
        m.setup(qp.P, qp.q_vecs[:, 0], Aosqp, losqp, uosqp,
                rho=10, verbose=False)

        # Get extension name
        module_name = 'emosqpn%s' % str(qp.n)

        # Generate the code
        m.codegen("code", python_ext_name=module_name, force_rewrite=True)

        # Import module
        emosqp = importlib.import_module(module_name)

        for i in range(n_prob):
            q = qp.q_vecs[:, i]

            # Update linear cost
            emosqp.update_lin_cost(q)

            # Solve
            x, y, status, niter, time = emosqp.solve()

            # Check if status correct
            if status != 1:
                import ipdb; ipdb.set_trace()
                raise ValueError('OSQP did not solve the problem!')

            # Solution statistics
            solution_dict = {'solver': [solver],
                             'runtime': [time],
                             'iter': [niter],
                             'n': [qp.n]}
            results.append(pd.DataFrame(solution_dict))

    elif solver == 'qpoases':
        '''
        Sparse problem formulation qpOASES
        '''
        n_dim = qp.P.shape[0]  # Number of variables
        m_dim = qp.A.shape[0]  # Number of constraints without bounds

        # Initialize qpoases and set options
        qpoases_m = qpoases.PyQProblem(n_dim, m_dim)
        options = qpoases.PyOptions()
        options.printLevel = qpoases.PyPrintLevel.NONE
        qpoases_m.setOptions(options)

        # Construct bounds for qpoases
        lx = np.append(qp.lx, -np.inf * np.ones(k))
        ux = np.append(qp.ux, np.inf * np.ones(k))

        # Setup matrix P and A
        P = np.ascontiguousarray(qp.P.todense())
        A = np.ascontiguousarray(qp.A.todense())

        for i in range(n_prob):

            # Get linera cost as contiguous array
            q = np.ascontiguousarray(qp.q_vecs[:, i])

            # Reset cpu time
            qpoases_cpu_time = np.array([20.])

            # Reset number of of working set recalculations
            nWSR = np.array([1000])

            if i == 0:
                res_qpoases = qpoases_m.init(P, q, A,
                                             np.ascontiguousarray(lx),
                                             np.ascontiguousarray(ux),
                                             np.ascontiguousarray(qp.l),
                                             np.ascontiguousarray(qp.u),
                                             nWSR, qpoases_cpu_time)
            else:
                # Solve new hot started problem
                res_qpoases = qpoases_m.hotstart(q,
                                                 np.ascontiguousarray(lx),
                                                 np.ascontiguousarray(ux),
                                                 np.ascontiguousarray(qp.l),
                                                 np.ascontiguousarray(qp.u),
                                                 nWSR,
                                                 qpoases_cpu_time)

            if res_qpoases != 0:
                raise ValueError('qpoases did not solve the problem!')

            # Solution statistics
            solution_dict = {'solver': [solver],
                             'runtime': [qpoases_cpu_time[0]],
                             'iter': [nWSR[0]],
                             'n': [qp.n]}
            results.append(pd.DataFrame(solution_dict))

    elif solver == 'gurobi':

        # Construct qp matrices
        Agurobi = spa.vstack((qp.A,
                              spa.hstack((spa.eye(n), spa.csc_matrix((n, k)))
                                         ))).tocsc()
        lgurobi = np.append(qp.l, qp.lx)
        ugurobi = np.append(qp.u, qp.ux)

        for i in range(n_prob):

            # Get linera cost as contiguous array
            q = np.ascontiguousarray(qp.q_vecs[:, i])

            # solve with gurobi
            prob = mpbpy.QuadprogProblem(qp.P, q, Agurobi, lgurobi, ugurobi)
            res = prob.solve(solver=mpbpy.GUROBI, verbose=False)

            # Solution statistics
            solution_dict = {'solver': [solver],
                             'runtime': [res.cputime],
                             'iter': [res.total_iter],
                             'n': [qp.n]}
            results.append(pd.DataFrame(solution_dict))

    else:
        raise ValueError('Solver not understood')

    return pd.concat(results)


'''
Solve problems
'''
# Generate gamma parameters and cost vectors
n_gamma = 11
gammas = np.logspace(2, -2, n_gamma)

# Assets
n_vec = np.array([50, 80, 100, 120, 150, 200, 250, 300, 400, 500])

# Factors
k_vec = (n_vec / 10).astype(int)

# Setup if solve with gurobi/qpoases or not
solve_osqp = True
solve_gurobi = True
solve_qpoases = True
solve_cvxgen = True
solve_fiordos = True

# Define statistics for osqp, gurobi and qpoases
if solve_osqp:
    import osqp
    osqp_stats = []
    problem_stats = []
if solve_gurobi:
    import mathprogbasepy as mpbpy
    gurobi_stats = []
if solve_qpoases:
    import qpoases
    qpoases_stats = []

# Size of the exe file generated by OSQP
if solve_osqp:
    if system() == 'Windows':
        cmdsep = '&'
        makefile = '"MinGW Makefiles"'
        example_fullname = 'example.exe'
    else:
        cmdsep = ';'
        makefile = '"Unix Makefiles"'
        example_fullname = 'example'


'''
Solve problems
'''
for i in range(len(n_vec)):

    # Generate QP sparse matrices
    qp_matrices = gen_qp_matrices(k_vec[i], n_vec[i], gammas)

    if solve_osqp:
        # Solving loop with emosqp
        stats = solve_loop(qp_matrices, 'emosqp')
        osqp_stats.append(stats)

        # Get size of the generated exe file in KB
        call('cd code %s ' % (cmdsep) +
             'mkdir build %s ' % (cmdsep) +
             'cd build %s ' % (cmdsep) +
             'cmake -G %s .. %s ' % (makefile, cmdsep) +
             ' cmake --build .',
             shell=True)
        example_path = os.path.join('code', 'build', 'out', example_fullname)
        example_size = os.path.getsize(example_path) / 1024.

        # Problem statistics
        N = qp_matrices.P.nnz + qp_matrices.A.nnz
        problem_dict = {'n': [qp_matrices.n],
                        'k': [qp_matrices.k],
                        'N': [N],
                        'filesize': example_size}
        problem_stats.append(pd.DataFrame(problem_dict))

    if solve_qpoases:
        # Solving loop with qpoases
        stats = solve_loop(qp_matrices, 'qpoases')
        qpoases_stats.append(stats)

    if solve_gurobi:
        # Solve loop with gurobi
        stats = solve_loop(qp_matrices, 'gurobi')
        gurobi_stats.append(stats)

if solve_cvxgen:
    call('matlab -nodesktop -nodisplay -nosplash -r ' +
         '"cd cvxgen; run run_all; exit;"')

if solve_fiordos:
    call('matlab -nodesktop -nodisplay -nosplash -r ' +
         '"cd fiordos; run run_all; exit;"')

'''
Store results in CSV files
'''
if solve_osqp:
    # Combine OSQP stats and store them in a CSV file
    df = pd.concat(osqp_stats)
    df.to_csv('osqp_stats.csv', index=False)

    # Combine problem stats and store them in a CSV file
    df = pd.concat(problem_stats)
    df.to_csv('problem_stats.csv', index=False)

if solve_gurobi:
    # Combine GUROBI stats and store them in a CSV file
    df = pd.concat(gurobi_stats)
    df.to_csv('gurobi_stats.csv', index=False)

if solve_qpoases:
    # Combine QPOASES stats and store them in a CSV file
    df = pd.concat(qpoases_stats)
    df.to_csv('qpoases_stats.csv', index=False)
