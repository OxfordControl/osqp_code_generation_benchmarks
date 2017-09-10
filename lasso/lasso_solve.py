from __future__ import print_function
from __future__ import division
import numpy as np
import scipy.sparse as spa
from builtins import range
import os
import pandas as pd

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
    def __init__(self, P, q_vecs, A, l, u, n, m):
        self.P = P
        self.q_vecs = q_vecs
        self.A = A
        self.l = l
        self.u = u
        self.n = n
        self.m = m


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

    # Results list
    results = []

    if solver == 'emosqp':
        # Pass the data to OSQP
        m = osqp.OSQP()
        m.setup(qp.P, qp.q_vecs[:, 0], qp.A, qp.l, qp.u,
                rho=10., verbose=False)

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
                print('OSQP did not solve the problem!')
                import ipdb
                ipdb.set_trace()
                raise ValueError('OSQP did not solve the problem!')

            # Solution statistics
            solution_dict = {'solver': [solver],
                             'runtime': [time],
                             'iter': [niter],
                             'n': [qp.n]}
            results.append(pd.DataFrame(solution_dict))

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
                res_qpoases = qpoases_m.init(P, q, A, None, None, qp.l, qp.u,
                                             nWSR, qpoases_cpu_time)
            else:
                # Solve new hot started problem
                res_qpoases = qpoases_m.hotstart(q, None, None, qp.l, qp.u,
                                                 nWSR, qpoases_cpu_time)

            # Solution statistics
            solution_dict = {'solver': [solver],
                             'runtime': [qpoases_cpu_time[0]],
                             'iter': [nWSR[0]],
                             'n': [qp.n]}
            results.append(pd.DataFrame(solution_dict))

    elif solver == 'gurobi':

        n_dim = qp.P.shape[0]
        m_dim = qp.A.shape[0]

        for i in range(n_prob):

            # Get linera cost as contiguous array
            q = np.ascontiguousarray(qp.q_vecs[:, i])

            # solve with gurobi
            prob = mpbpy.QuadprogProblem(qp.P, q, qp.A, qp.l, qp.u)
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
Problem parameters
'''
# Generate gamma parameters and cost vectors
n_gamma = 21
gammas = np.logspace(2, -2, n_gamma)

# Number of parameters
n_vec = np.array([10, 20, 30, 50, 80, 100, 150, 200, 250, 300, 350, 400])

# Measurements
m_vec = (10 * n_vec).astype(int)

# Setup if solve with gurobi/qpoases or not
solve_osqp = True
solve_gurobi = True
solve_qpoases = True

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
    qp_matrices = gen_qp_matrices(m_vec[i], n_vec[i], gammas)

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
        example_size = int(round(os.path.getsize(example_path) / 1024.))

        # Problem statistics
        N = qp_matrices.P.nnz + qp_matrices.A.nnz
        problem_dict = {'n': [qp_matrices.n],
                        'm': [qp_matrices.m],
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
