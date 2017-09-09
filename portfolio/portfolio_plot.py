from __future__ import print_function
from __future__ import division
import pandas as pd

# Import scipy io to write/read mat file
import scipy.io as io
import os

# Plotting
import matplotlib.pylab as plt
plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=15)   # fontsize of the tick labels
plt.rc('ytick', labelsize=15)   # fontsize of the tick labels
plt.rc('legend', fontsize=15)   # legend fontsize
plt.rc('text', usetex=True)     # use latex
plt.rc('font', family='serif')

# Import data from mat files to pandas dataframe
cvxgen_results = io.loadmat(os.path.join('cvxgen',
                                         'cvxgen_stats.mat'))['cvxgen_stats']
cvxgen_stats = pd.DataFrame(cvxgen_results, columns=['n', 'cvxgen_runtime'])
cvxgen_stats['n'] = cvxgen_stats['n'].astype(int)

fiordos_results = io.loadmat(os.path.join('fiordos',
                             'fiordos_stats.mat'))['fiordos_stats']
fiordos_stats = pd.DataFrame(fiordos_results, columns=['n', 'fiordos_runtime'])
fiordos_stats['n'] = fiordos_stats['n'].astype(int)

# Import data from CSV files
osqp_stats = pd.read_csv('osqp_stats.csv')
gurobi_stats = pd.read_csv('gurobi_stats.csv')
qpoases_stats = pd.read_csv('qpoases_stats.csv')

# Get mean number of iterations and runtime for each number of parameters
osqp_mean = osqp_stats.groupby('n').mean()
gurobi_mean = gurobi_stats.groupby('n').mean()
qpoases_mean = qpoases_stats.groupby('n').mean()
cvxgen_mean = cvxgen_stats.groupby('n').mean()
fiordos_mean = fiordos_stats.groupby('n').mean()

# Save mean stats in a CSV file
osqp_mean.columns = ['osqp_iter', 'osqp_runtime']
gurobi_mean.columns = ['gurobi_iter', 'gurobi_runtime']
qpoases_mean.columns = ['qpoases_iter', 'qpoases_runtime']
mean_stats = pd.concat([osqp_mean, gurobi_mean, qpoases_mean,
                        cvxgen_mean, fiordos_mean], axis=1)
mean_stats.to_csv('portfolio_mean_stats.csv')

# Plot mean runtime
plt.figure()
ax = plt.gca()
plt.semilogy(mean_stats.index.values,
             mean_stats['osqp_runtime'].get_values(),
             color='C0', label='OSQP')
plt.semilogy(mean_stats.index.values,
             mean_stats['qpoases_runtime'].get_values(),
             color='C2', label='qpOASES')
plt.semilogy(mean_stats.index.values,
             mean_stats['cvxgen_runtime'].get_values(),
             color='C3', label='CVXGEN')
plt.semilogy(mean_stats.index.values,
             mean_stats['fiordos_runtime'].get_values(),
             color='C4', label='FiOrdOs')
plt.semilogy(mean_stats.index.values,
             mean_stats['gurobi_runtime'].get_values(),
             color='C5', label='GUROBI')
plt.legend()
plt.grid()
ax.set_xlabel(r'Number of assets $n$')
ax.set_ylabel(r'Time [s]')
plt.tight_layout()
plt.show(block=False)
plt.savefig('results.pdf')
