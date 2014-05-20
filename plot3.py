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

def plot(I_1, A_1, I_2, A_2, out_dir):
  dim = A_1.shape[0]
  matplotlib_setup()
  plt.figure()
  plt.loglog(I_1, 'r', alpha=0.8)
  plt.loglog(I_2, 'r', alpha=0.8)
  plt.xlabel('Time in days ')
  plt.ylabel('Infectious')
  plt.savefig(out_dir + '/diffusion1.svg')
  plt.savefig(out_dir + '/diffusion1.pdf')

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
  #
  # Parse argument
  #
  parser = argparse.ArgumentParser(description='Process SIR simulation with latent states.')
  parser.add_argument('--output', help='output directory', required=True)
  parser.add_argument('--latent', help='directory of the latent simualtion results', required=True)
  parser.add_argument('--nolatent', help='directory of the nolatent simualtion results', required=True)
  args = parser.parse_args()
  argsdict = vars(args)
  if args.latent and  args.nolatent and args.output:
    output_dir = argsdict['output']
    latent = argsdict['latent']
    nolatent = argsdict['nolatent']
    # Remove the last backslash of teh sting if exist
    if output_dir.endswith('\\'):
      output_dir = output_dir[:-1]
    # Remove the last backslash of teh sting if exist
    if latent.endswith('\\'):
      latent = latent[:-1]
    # Remove the last backslash of teh sting if exist
    if nolatent.endswith('\\'):
      nolatent = nolatent[:-1]
    # if output dir doesn' extist create it
    if not os.path.exists(output_dir):
      os.makedirs(output_dir)
    S_L, I_L, R_L, A_L = load_files(latent)
    S_NL, I_NL, R_NL, A_NL = load_files(nolatent)
    plot(I_L, A_L, I_NL, A_NL, output_dir)
