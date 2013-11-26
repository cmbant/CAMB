# sample script to compare two different sets of cross-correlation spectra from scalar_covariance_output_file output
# e.g. use PyDev in Eclipse to easily edit and result run python scripts

from cambPlots import *
from pylab import *


root = r'../test_lensing_scalCovCls.dat'
root2 = r'../test_lensing2_scalCovCls.dat'

# Plot just one output
plot_array([root], y='log', x='log', half=True, diff=False)

# Fractional difference between two outputs
# plot_array([root, root2], y=None, x='log', half=True, diff=True)

# Fractional difference between two outputs, restricted to correlations between sources 4 and 5 (CAMB sources, first two window functions)
# plot_array([root, root2], indices=[3, 4], y=None, x='log', half=True, diff=True)

# indices[0,2] would compare TT and CMB lensing potentials
# plot_array([root, root2], indices=[0, 2], y=None, x='log', half=True, diff=True)


show()


