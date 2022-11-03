# This script analyses the results of runs of eco-dyn.py - it can be used to do 3 things:
# 1. Plot abundance trajectories for individual KP, KD values - each run + mean
# 2. Plot a heatmap showing the fate of the ecological dynamics across the KP-KD parameter space
# 3. Plot a measure of variance between independent runs across the KP-KD parameter space

import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import seaborn as sns
import ternary
import math
import pickle
from matplotlib.colors import ListedColormap
import pickle
import os

##################################################### CONFIG PARAMETERS ##########################################################

# all parameters exactly correspond to those in eco-dyn.py, see there for more information. 

no_runs = 15

update_iterations = 6000

N = 3
M = 1

parameter_values = np.arange(0.7,18.5,0.6) 

end_stats = np.zeros([len(parameter_values),len(parameter_values),3])
traj_std = np.zeros([len(parameter_values), len(parameter_values)])

community_identifier = 'PSD-1-antibiotic'
strain_labels = ['P', 'S', 'D']
grid_colors = [[0.0,0.0,0.0], [0.172549, 0.482353, 0.713725], [0.843137, 0.29, 0.29], [0.670588, 0.85098, 1], [0.992157, 0.682353, 0.380392]] 

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


## 1. plotting trajectories for a particular parameter combination:
K_D = 15.7
K_P = 16.3
generations_cutoff = 4000 # when, if at all, do you want to cut off the trajectories to enhance the visualisation of the dynamics. 
                          # If this cutoff is not necessary, set generations_cutoff to the same number as update_iterations.

jobname = community_identifier + '_KD-' + str(K_D) + '_KP-' + str(K_P) # job identifier that will be used to store figures/movies related to this job

all_run_stats = []
for r in range(1, no_runs+1):
       filename = "../SimulationCode/spatial_model/results/simulation_statistics/" + community_identifier  + '/stats/' + jobname + 'run' + str(r) + "_stats"
       with open(filename, 'rb') as fil:
             this_run_stats = pickle.load(fil)
       all_run_stats.append(this_run_stats[:generations_cutoff])

# now we calculate the average trajectory 
avg_trajectory = np.zeros([generations_cutoff,N]) # update_iterations rows, N columns
for t in range(generations_cutoff): 
       # for each time point t, find the average abundance of each strain at t
       for strain_idx in range(N):
             avg_trajectory[t][strain_idx] = np.mean([traj[t][strain_idx] for traj in all_run_stats])

fig, ax = plt.subplots(1, 1)
ax.set_ylim([0,1])
legend_labels = strain_labels

# plot the average trajectory

for run in range(1,no_runs+1):
       run_stats = all_run_stats[run-1]
       for strain_idx in range(N):
         strain_abundance = [timestep[strain_idx] for timestep in run_stats]
         ax.plot(range(generations_cutoff), strain_abundance, 
                       color = grid_colors[strain_idx+1],
                       linewidth=0.3, alpha=0.5)  

for strain_idx in range(N):
     strain_abundance = [timestep[strain_idx] for timestep in avg_trajectory]
     ax.plot(range(generations_cutoff), strain_abundance, 
                   label = legend_labels[strain_idx], 
                   color = grid_colors[strain_idx+1],linewidth=3.0, alpha=1.0)

plt.ylim([-0.05, 1.05])
plt.savefig(jobname + '.pdf', format='pdf')
plt.show()
# plt.close(fig)

# 2. Analysing data for all KP-KD values. Given time-abundance data for each run of all parameter values, this block will find the long-term fate 
# of the ecological dynamics i.e., the time-average of the last 300 timesteps of the run-averaged trajectory and store it as a matrix with shape 
# identical to the parameter space. 

# If you already have such a list that you want to analyse, you can comment this block and only keep the "opening.." data block below

for kd in range(len(parameter_values)):
      for kp in range(len(parameter_values)):
           jobname = community_identifier + '_KD-' + str(format(parameter_values[kd], '.1f')) + '_KP-' + str(format(parameter_values[kp], '.1f')) 
           all_run_stats = []
           for r in range(1, no_runs+1):
                filename = "../SimulationCode/spatial_model/results/simulation_statistics/" + community_identifier  + '/stats/' + jobname + 'run' + str(r) + "_stats"
                with open(filename, 'rb') as fil:
                     this_run_stats = pickle.load(fil)
                all_run_stats.append(this_run_stats)
           # now we calculate the average trajectory 
           avg_trajectory = np.zeros([update_iterations,N]) # update_iterations rows, N columns
          
           for t in range(update_iterations): 
                # for each time point t, find the average abundance of each strain at t
                for strain_idx in range(N):
                     avg_trajectory[t][strain_idx] = np.mean([traj[t][strain_idx] for traj in all_run_stats])

          #  # calculating standard deviation in trajectories. do this only if there is no file named trajectory_deviation_data in this directory
          #  # or if you have generated new trajectories. 
           overall_std = []
           for strain_idx in range(N):
                # finding standard deviation of each strain trajectory
                std_t = list(range(update_iterations))
                for t in range(update_iterations):
                     # for each time point t, find the standard deviation of runs about the mean value at time t 
                     std_t[t] = np.std([traj[t][strain_idx] for traj in all_run_stats])
                strain_std = sum(std_t)
                overall_std.append(std_t)
           traj_std[len(parameter_values)-1-kd][kp] = np.mean(overall_std)

           # long run behaviour
           end_slice = avg_trajectory[update_iterations-300:update_iterations]
           end_abundance = end_slice.mean(axis=0)
           end_stats[len(parameter_values)-1-kd][kp][:] = np.array(end_abundance)
      print(kd, 'done')

for i in range(len(parameter_values)):
      for j in range(len(parameter_values)):
           [r,g,b] = end_stats[i][j]
           end_stats[i][j][:] = r*np.array([0.172549, 0.482353, 0.713725]) + g*np.array([0.843137, 0.29, 0.29]) + b*np.array([0.670588, 0.85098, 1])

end_stats = end_stats[:,5:]
# The pruning above is done since we don't want K_P<3 since that corresponds to a standard deviation of less than one lattice site (see main text). 

# # Saving the data
# with open("trajectory_deviation_data",'wb') as fil:
#       pickle.dump(traj_std[:,5:],fil)

# with open("region-of-coexistence_data", 'wb') as fil2:
#       pickle.dump(end_stats[:,5:],fil2)

# # 3. Opening the data if you already have it
# with open("trajectory_deviation_data",'rb') as fil:
#       traj_std = pickle.load(fil)

# with open("region-of-coexistence_data",'rb') as fil2:
#       end_stats = pickle.load(fil2)

num_ticks = 10
# the index of the position of yticks
yticks = np.linspace(0, len(parameter_values) - 1, num_ticks, dtype=np.int)
xticks = np.linspace(0, len(parameter_values[5:]) - 1, num_ticks, dtype=np.int)
# the content of labels of these yticks
yticklabels = [format(parameter_values[idx], '.1f') for idx in yticks]
yticklabels.reverse()
xticklabels = [format(parameter_values[4+idx], '.1f') for idx in xticks]

# Plotting the variation (traj_std) between independent runs of the simulation for all parameter combinations.
s = sns.heatmap(traj_std, yticklabels=yticklabels, xticklabels=xticklabels)
s.set(xlabel='Diffusivity of antibiotic (K_P)', ylabel='Diffusivity of degrader (K_D)', rasterized=True)
s.set_yticks(yticks)
s.set_xticks(xticks)
plt.savefig("trajectory_deviation.pdf", format='pdf')
plt.show()

# Plotting the fate (end_stats) for all parameter combinations
plt.imshow(end_stats, extent = [parameter_values[5], parameter_values[-1], parameter_values[0], parameter_values[-1]], aspect='equal')
plt.xlabel("Diffusivity of antibiotic (K_P)")
plt.ylabel("Diffusivity of degrader (K_D)")
plt.savefig("region-of-coexistence.pdf", format='pdf')
plt.show()
