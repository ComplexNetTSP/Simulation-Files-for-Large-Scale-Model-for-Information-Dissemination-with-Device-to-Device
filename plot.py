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
simulation_end_time = 30.0
tau = 1.0/60
directory = 'Results/mult/50'
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

def plot(S, I, R, A):
  dim = A.shape[0]
  matplotlib_setup()
  x = np.arange(0, simulation_end_time, tau)
  plt.figure()
  for i in xrange(dim):
      plt.loglog(x, A[i, :], alpha=0.15)
      plt.ylabel('Population')
      plt.xlabel('Time in days')
      plt.xlim((10**(-1), 40))
  plt.savefig('diffusion.svg')
  plt.savefig('diffusion.pdf')

  plt.figure()
  # plt.subplot(311)
  # plt.plot(S, 'g')
  # plt.xlabel ('Time (years)')
  # plt.ylabel ('Susceptible')
  # plt.subplot(312)
  plt.loglog(x, I, 'r', alpha=0.8)
  plt.xlim((10**(-1), 40))
  plt.xlabel('Time in days ')
  plt.ylabel('Infectious')
  plt.savefig('diffusion1.svg')
  plt.savefig('diffusion1.pdf')

  #plt.subplot(313)
  # plt.plot(R, 'k')
  # plt.xlabel ('Time (years)')
  # plt.ylabel ('Recovered')
  # plt.show()

def load_files(directory):
  with open(directory+'/S.p', 'rb') as fp:
      S = p.load(fp)
  with open(directory+'/I.p', 'rb') as fp:
      I = p.load(fp)
  with open(directory+'/R.p', 'rb') as fp:
      R = p.load(fp)
  with open(directory+'/A.p', 'rb') as fp:
      A = p.load(fp)

  return S, I, R, A

if __name__ == '__main__':
  S, I, R, A = load_files(directory)
  plot(S, I, R, A)
