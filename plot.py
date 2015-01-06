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
import os
import argparse

###############################################################################
#
# Begining of global definition
#
simulation_end_time = 30.0
tau = 1.0/5
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

def plot(I, A, output_dir):
  dim = A.shape[0]
  matplotlib_setup()
  x = np.arange(0, simulation_end_time, tau)
  plt.figure()
  for i in xrange(dim):
      plt.loglog(A[i, :], alpha=0.15)
      plt.ylabel('Population')
      plt.xlabel('Time in days')
      #plt.xlim((10**(-1), 40))
  plt.savefig(output_dir + '/diffusion.svg')
  plt.savefig(output_dir + '/diffusion.pdf')

  plt.figure()
  plt.loglog(I, 'r', alpha=0.8)
  #plt.xlim((10**(-1), 40))
  plt.xlabel('Time in days ')
  plt.ylabel('Infectious')
  plt.savefig(output_dir + '/diffusion1.svg')
  plt.savefig(output_dir + '/diffusion1.pdf')

def load_files(directory):
  with open(directory+'/I.p', 'rb') as fp:
      I = p.load(fp)
  with open(directory+'/A.p', 'rb') as fp:
      A = p.load(fp)
  return I, A

if __name__ == '__main__':
  #
  # Parse argument
  #
  parser = argparse.ArgumentParser(description='Process SIR simulation with latent states.')
  parser.add_argument('--output', help='output directory', required=True)
  parser.add_argument('--input', help='input directory', required=True)
  args = parser.parse_args()
  argsdict = vars(args)
  if args.output and args.input:
    output_dir = argsdict['output']
    input_dir = argsdict['input']
    # Remove the last backslash of the string if exist
    if output_dir.endswith('\\'):
      output_dir = output_dir[:-1]
    # Remove the last backslash of the string if exist
    if input_dir.endswith('\\'):
      input_dir = input_dir[:-1]

    # if output dir doesn' extist create it
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)

  I, A = load_files(input_dir)
  plot(I, A, output_dir)
