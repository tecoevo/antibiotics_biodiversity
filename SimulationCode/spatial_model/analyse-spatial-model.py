import this
import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import seaborn as sns
import ternary
import math
import pickle

##################################################### CONFIG PARAMETERS ##########################################################

# all parameters exactly correspond to those in eco-dyn.py, see there for more information. 

no_runs = 3 

update_iterations = 10

N = 3
M = 1

parameter_values = np.arange(0.7,18.1,0.6) 

end_stats = np.zeros([len(parameter_values),len(parameter_values),3])

community_identifier = 'PSD-1-antibiotic'

##################################################### FUNCTION DEFINITIONS ##########################################################

def color_point(x, y, z, scale):
      w = 255
      x_color = x * w / float(scale)
      y_color = y * w / float(scale)
      z_color = z * w / float(scale)
      r = z_color / w
      g = y_color / w
      b = x_color / w
      return (r, g, b, 1.)

def generate_heatmap_data(scale=90):
      from ternary.helpers import simplex_iterator
      d = dict()
      for (i, j, k) in simplex_iterator(scale):
          d[(i, j, k)] = color_point(i, j, k, 90)
      return d

def tuple2rgb(tup):
      rgb = []
      for number in tup:
           number = 256*float(format(number, '.3f'))
           rgb.append(round(number))
      return rgb

def UniformSampleSimplex():
     """Returns a uniform random point on the 2-simplex i.e., the triangle."""
     random = [0.0]
     random.append(np.random.uniform(0,1))
     random.append(np.random.uniform(0,1))
     random.append(1.0)
     random = sorted(random)
     simplex_number = [random[1]-random[0], random[2]-random[1], random[3]-random[2]]
     return simplex_number

##################################################### ANALYSIS ##########################################################

parameter_values = [3.1,18.0]
for kd in range(len(parameter_values)):
      for kp in range(len(parameter_values)):
           jobname = community_identifier + '_KD-' + str(format(parameter_values[kd], '.1f')) + '_KP-' + str(format(parameter_values[kp], '.1f')) 
           all_run_stats = []
           for r in range(1, no_runs):
                filename = "results/simulation_statistics/" + community_identifier  + '/stats/' + jobname + 'run' + str(r) + "_stats"
                with open(filename, 'rb') as fil:
                     this_run_stats = pickle.load(fil)
                all_run_stats.append(this_run_stats)

           # now we calculate the average trajectory 
           avg_trajectory = np.zeros([update_iterations,N]) # update_iterations rows, N columns

           for t in range(update_iterations): 
                # for each time point t, find the average abundance of each strain at t
                for strain_idx in range(N):
                     avg_trajectory[t][strain_idx] = np.mean([traj[t][strain_idx] for traj in all_run_stats])

           # long run behaviour
           end_slice = avg_trajectory[update_iterations-5:update_iterations]
           end_abundance = end_slice.mean(axis=0)
           end_stats[len(parameter_values)-1-kd][kp][:] = end_abundance
           

# plotting behaviour in the K_P-K_D parameter space

plt.imshow(end_stats[:,5:], extent=[3.1, 18.1, 0.7, 18.1])
plt.xlabel("K_P")
plt.ylabel("K_D")
plt.suptitle("Region of coexistence in the K_P-K_D parameter space")
plt.savefig("../Manuscript/Working_draft/Figures/region-of-coexistence.pdf")