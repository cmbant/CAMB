# sample script to compare two standard _lensedCl results
# ncl=1 does just TT, ncl=3,4 does more
# see camPlots.py for function declaration

from cambPlots import *
from pylab import *


root = r'z:\test3_lensedCls.dat_1_3'
root2 = r'z:\test_zerosource_lensedCls.dat_1_3'


plot_compare([root, root2], ncl=2, x='log', y='log')


show()

