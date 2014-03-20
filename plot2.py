# -*- encoding: utf-8 -*-

# -------------------------------------------------------------------------------
# Copyright (c) 2014 Vincent Gauthier Telecom SudParis.
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
# -------------------------------------------------------------------------------

__author__ = """\n""".join(['Vincent Gauthier <vgauthier@luxbulb.org>'])

import pickle as p
import pylab as plt
import numpy as np

###############################################################################
#
# Begining of global definition
#
simulation_end_time=10.0
tau = 1.0/30
filelist=['Results/mult/1', 'Results/mult/2', 'Results/mult/3', 'Results/mult/4',
  'Results/mult/5', 'Results/mult/6', 'Results/mult/7', 'Results/mult/8', 'Results/mult/9']
###############################################################################
#
# End of global definition
#

def matplotlib_setup(figsize_x=10, figsize_y=6):
    import matplotlib as mpl
    mpl.rcParams['font.size'] = 9.0
    mpl.rcParams['font.weight'] = 'bold'
    mpl.rcParams['xtick.labelsize'] = 10
    mpl.rcParams['ytick.labelsize'] = 10
    mpl.rcParams['axes.labelsize'] = 'large'
    mpl.rcParams['axes.labelweight'] = 'bold'
    mpl.rcParams['axes.linewidth'] = 0.75
    mpl.rcParams['lines.linewidth'] = 2
    mpl.rcParams['lines.markersize'] = 8
    mpl.rcParams['legend.numpoints'] = 1
    # figure size in inch
    mpl.rcParams['figure.figsize'] = figsize_x, figsize_y
    # figure dots per inch
    mpl.rcParams['figure.dpi'] = 300

def plot(I):
  x = np.arange(0,simulation_end_time,tau)
  plt.loglog(x, I, 'r', alpha=0.8)
  plt.xlim((10**(-1), simulation_end_time))
  plt.xlabel ('Time in days ')
  plt.ylabel ('Infectious')
  plt.savefig('diffusion1.svg')
  plt.savefig('diffusion1.pdf')

def load_files(directory):
  with open(directory+'/I.p', 'rb') as fp:
      I = p.load(fp)
  return I

if __name__ == '__main__':
  matplotlib_setup()
  plt.figure()
  for filename in filelist:
    I = load_files(filename)
    plot(I)