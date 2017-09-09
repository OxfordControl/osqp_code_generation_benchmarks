from __future__ import print_function
from __future__ import division
import pandas as pd

# Plotting settings
import matplotlib.pylab as plt
plt.rc('axes', labelsize=20)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=15)   # fontsize of the tick labels
plt.rc('ytick', labelsize=15)   # fontsize of the tick labels
plt.rc('legend', fontsize=15)   # legend fontsize
plt.rc('text', usetex=True)     # use latex
plt.rc('font', family='serif')

# Import data from CSV files
osqp_stats = pd.read_csv('osqp_stats.csv')
gurobi_stats = pd.read_csv('gurobi_stats.csv')
qpoases_stats = pd.read_csv('qpoases_stats.csv')

# Get mean number of iterations and runtime for each number of parameters
osqp_mean = osqp_stats.groupby('n').mean()
gurobi_mean = gurobi_stats.groupby('n').mean()
qpoases_mean = qpoases_stats.groupby('n').mean()

# Save mean stats in a CSV file
osqp_mean.columns = ['osqp_iter', 'osqp_runtime']
gurobi_mean.columns = ['gurobi_iter', 'gurobi_runtime']
qpoases_mean.columns = ['qpoases_iter', 'qpoases_runtime']
mean_stats = pd.concat([osqp_mean, gurobi_mean, qpoases_mean], axis=1)
mean_stats.to_csv('lasso_mean_stats.csv')

# Plot mean runtime
plt.figure()
ax = plt.gca()
plt.semilogy(mean_stats.index.values,
             mean_stats['osqp_runtime'].get_values(),
             color='C0', label='OSQP')
plt.semilogy(mean_stats.index.values,
             mean_stats['qpoases_runtime'].get_values(),
             color='C1', label='qpOASES')
plt.semilogy(mean_stats.index.values,
             mean_stats['gurobi_runtime'].get_values(),
             color='C5', label='GUROBI')
plt.legend()
plt.grid()
ax.set_xlabel(r'Number of parameters $n$')
ax.set_ylabel(r'Time [s]')
plt.tight_layout()
plt.show(block=False)
plt.savefig('results.pdf')
