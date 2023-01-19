import numpy as np 
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import scipy as sp
from matplotlib import colors
import itertools as it
import os
from scipy.ndimage import gaussian_filter
from copy import deepcopy
from matplotlib.colors import ListedColormap
from math import pi, exp, log 
import multiprocessing as mp
import random
import argparse
import pickle

if __name__ == "__main__": 
     #this check makes sure that arguments are asked for only when this file is the main program($ python3 ...) and not when it is imported in another piece of code
     parser = argparse.ArgumentParser(description="Individual-based spatial model of antibitic-mediated interactions.  ")

     parser.add_argument("-kd", "--kd", dest="K_D", help="K_D, the diffusivity of degrader molecules.  ", required=True)

     parser.add_argument("-kp", "--kp", dest="K_P", help="K_P, the diffusitivity of antibiotic molecules.  ", required=True)
     
     parser.add_argument("-ss", "--ss", dest="save_snapshots", type=bool, help="Binary, do you want to save snapshots after each \
         lattice update? Requires a good amount of memory space.  Default False. DO NOT pass an argument on the command line if you want to keep it False.", default=False)

     parser.add_argument("-sm", "--sm", dest="save_movie", type=bool, help="Binary, do you want to save snapshots after \
         each lattice update? Requires a good amount of memory space.  Default False. DO NOT pass an argument on the command line if you want to keep it False.", default=False)

     parser.add_argument("-v", "--v", dest="verbose", type=bool, help="Binary, prints a line stating the number of generations \
         completed after every lattice update. Note that if you are running this simulation on a cluster, the output files will \
         likely take up more space than expected if this flag is set to true. Default False. DO NOT pass an argument on the command line if you want to keep it False. ", default=False)


     args = parser.parse_args() 
     K_D = float(args.K_D)
     K_P = float(args.K_P)
     save_snapshots = args.save_snapshots
     save_movie = args.save_movie
     verbose = args.verbose

############################################### STRAIN CLASS DEFINITION ##########################################################

class Strain(object):
     """Strain class."""

     def __init__(self,
                  phenotype):
         """Initialize strain instance.

         Args:
             phenotype: row vector of length M describing phenotype with respect to each antibiotic.
                        each element of the vector is taken from {1,2,3,4}, standing for producer,
                        sensitive, degrader, intrinsically resistant respectively.
         """
         self.phenotype = phenotype

         # calculate growth rates for the strain by taking into account the metabolic cost of
         # antibiotic-related phenotypes.

         rate = g
         for antibiotic in range(M):
             phenotype_wrt_antibiotic_i = self.phenotype[antibiotic]
             if phenotype_wrt_antibiotic_i == 1:   # if the strain is a producer of antibiotic 1
                 rate = rate - c_prod[antibiotic]
             elif phenotype_wrt_antibiotic_i == 3: # if the strain is a degrader of antibiotic 1
                 rate = rate - c_deg[antibiotic]
             elif phenotype_wrt_antibiotic_i == 4: # if the strain is intrinsically resistant to antibiotic 1
                 rate = rate - c_res[antibiotic]
                 
         self.growth_rate = rate

##################################################### CONFIG PARAMETERS ##########################################################

# This file contains the parameters necessary for running the agent-based model in lattice_model.py

# number of runs of the simulation in parallel

no_runs = 3

# set the jobname and the names of the strains
 
community_identifier = 'PSD-1-antibiotic'
jobname = community_identifier + '_KD-' + str(K_D) + '_KP-' + str(K_P) # job identifier that will be used to store figures/movies related to this job

# setting up directories to store results
try: # error thrown if directory already exists
     os.mkdir("results/simulation_snapshots/" + community_identifier)
except: 
     pass

try: # error thrown if directory already exists
     os.mkdir("results/simulation_statistics/" + community_identifier)
except: 
     pass

try: 
     os.mkdir("results/simulation_snapshots/" + community_identifier + '/' + jobname)
except: 
     pass

try: 
     os.mkdir("results/simulation_statistics/" + community_identifier + '/plots')
     os.mkdir("results/simulation_statistics/" + community_identifier + '/stats')
except: 
     pass


# Create the strains using the Strain class. Note that the number of strains created here must be the same as the value of N, and 
# the number of columns in the row vector provided as argument provided to Strain() must be equal to M. The values of N and M can 
# be found in config_parameters.py

N = 3 # number of species
M = 1 # number of antibiotics
g = 0.7 # base fitness of each strain
d = 0.3 # base death rate of each strain

# metabolic costs of:
c_res = np.array([1/20.0]) # intrinsic resistance. uniformly chosen from [0.01, 0.06]. 
c_deg = 2.1*c_res  # resistance via degradation
c_prod = 3*c_res # production

strain_1 = Strain([1])
strain_2 = Strain([2])
strain_3 = Strain([3])

strain_labels = ['P', 'S', 'D']

strains = [strain_1,strain_2,strain_3]

initial_occupation_fraction = 0.3 # Fraction of grid cells that we want to be occupied when initializing the simulation
                                  # Probability of any grid cell being occupied by a certain strain is given by initial_occupation_fraction/N.

update_iterations = 10 # number of generations for which the grid must be updated. 
# Note that 1 generation is considered to be the total time taken to update the state of each individual currently on the lattice. 
# An update step is defined to be the process of updating the state of a single individual.

grid_colors = [[0.0,0.0,0.0], [0.172549, 0.482353, 0.713725], [0.843137, 0.29, 0.29], [0.670588, 0.85098, 1], [0.992157, 0.682353, 0.380392]] 
# used to colour the grid cells. First element must always be "black" (i.e., empty), 
# length of the array must be N+1 (one colour for each strain + 1 for empty cell)

cmap = ListedColormap(grid_colors) # color map to colour the grid cells according to the identity of their 
                                                          # occupants (or lack thereof). 

size = 200 # length (and breadth) of the lattice on which we are performing simulations

# K_D = 2.1   # degrading area, enter float value
# K_P = 1.5*K_D     # killing area

u_p = 10 # amount of antibiotic being produced at every update step

u_d = 10 # amount of degrader being produced at every update step

def G_2d(x, y, sigma):
     """ This function returns the value of the gaussian density at (x, y), with mean = 0 and sigma taken as input. """
     return exp( -(x**2 + y**2) / (2* (sigma**2)) ) / (2*pi*sigma**2)

c1 = 0.99
cs = 0.9

cstar = log((1-cs)/(1-d)) / log((1-c1)/(1-d))

threshold = (u_p*(cstar*G_2d(0,1,K_P/3) - G_2d(0,K_P/3, K_P/3))) / (cstar - 1) # see the main text to understand how this parameter has been set. 
susceptibility = ((cstar-1)*(log((1-c1)/(1-d)))) / (u_p*(G_2d(0,1,K_P/3)-G_2d(0,K_P/3,K_P/3))) # see the main text to understand how this parameter has been set.


############################################### FUNCTIONS AND MODEL DEFINITIONS #####################################################

# This file defines the classes of objects that we will use to build an agent-based model.

def create_movie(figuredir, movie_name, iterations, savedir, run):
    """Compile a sequence of pictures into a movie and save it.  \n 

    figuredir: directory from which the images are being taken. All images 
    must be named "generation_i.jpg" for some number i. 
    movie_name: name given to the resulting movie created from the images. 
    iterations: how many iterations has the lattice been updated for? 
    savedir: directory in which the movie must be saved. """

    images = ['run_' + str(run) + '_generation_'+str(generation)+'.jpg' for generation in range(iterations)]

    frame = cv2.imread(os.path.join(figuredir, images[0]))
  
    # setting the frame width, height width 
    # the width, height of first image
    height, width, layers = frame.shape  
  
    video = cv2.VideoWriter(filename=savedir + '/' + movie_name, fourcc=cv2.VideoWriter_fourcc(*"mp4v"), fps=100, frameSize=(width, height)) 
  
    # Appending the images to the video one by one
    for image in images: 
        video.write(cv2.imread(os.path.join(figuredir, image))) 
      
    # Deallocating memories taken for window creation
    cv2.destroyAllWindows() 
    video.release()  # releasing the video generated

def antibiotic_degrader_interaction(A, D, size):
     """A: matrix of antibiotic concentrations at each grid cell.
     D: matrix of degrader concentrations at each grid cell.
     size: side length of the square matrices.

     Implement the inactivation of the antibiotic by its degrader molecule, if present. This amounts to 
     subtracting the matrix D from the matrix A only if antibiotic concentration is non-zero at that point,
     and changing the degrader concentration to reflect how much has been used up in the reaction. 
     Returns modified versions of the same pair of matrices as a list.   

     Note: The two matrices must be of identical shape. """

     A_new, D_new = np.zeros([size, size]), np.zeros([size,size])
     for row in range(size):
         for col in range(size):
             A_new[row][col] = np.max([0, A[row][col] - D[row][col]])
             D_new[row][col] = np.max([0, D[row][col] - A[row][col]])
                          
     return [A_new, D_new]  
            
def dose_response(strain, A, M, row, col, tau, k):
     """Decides the death probability of a focal individual depending on the grid-concentration of all antibiotics that it is sensitive to.
     strain: instance of the Strain class, describing the phenotype of the focal individual.
     A: dict of antibiotic concentrations on the lattice, and in particular, on the grid cell that the focal individual is currently present.
     row, col: position of the individual on the lattice.
     tau: minimum concentration of antibiotic required for inhibition.
     k: rate of convergence of death probability to 1. """

     # the base value of death rate is the intrinsic death probability, which we set here to be d (defined above)
     phen = strain.phenotype
     ant_conc = 0
     for i in range(M):
         A_i_conc = A['antibiotic'+str(i)][row][col]
         if phen[i]==2:
             ant_conc += A_i_conc
      
     if ant_conc < tau:
        #  print('below threshold')
         return d
     if ant_conc >= tau:
        #  print('above threshold')
         return d + (1-d)*(1-(exp(k*tau)*exp(-k*ant_conc))) # one could calculate exp(k*tau - k*conc) instead of multiplying the two exponential terms 
                                                # together but that will lead to a subtraction term, causing loss of significance when  
                                                # ant_conc ~ threshold.

def DispersalNbhd(lattice, row, col, size):
     """row, col: coordinates of a given grid cell \n
     Returns the indices of the empty grid cells in the Moore neighbourhood of this grid cell."""

     indices = [[row % size, (col+1) % size],
             [row % size, (col-1) % size],
             [(row+1) % size, col % size],
             [(row-1) % size, col % size],
             [(row-1) % size, (col-1) % size],
             [(row+1) % size, (col+1) % size]]

     if np.prod([lattice[idxx][idxy] for [idxx,idxy] in indices]) != 0:
         return -1 # no empty cells in the Moore neighbourhood
     else:
         return [idx for idx in indices if lattice[idx[0]][idx[1]]==0]

class Lattice_Simulation():
     """Lattice on which the simulations take place. """

     def __init__(self,
                 N,
                 M,
                 size,
                 iterations,
                 strains,
                 cmap,
                 jobname):
         """Initialize lattice instance.

         Args:
             N: number of strains.
             M: number of antibiotics.
             size: positive integer. length(and breadth) of the square lattice being created.
             iterations: positive integer. number of iterations for which the simulation is performed.
             strains: 1xN array, with each element being an instance of the Strain class defined above.
             cmap: ListedColormap object containing colours to be assigned to the N strains. No strain is always assigned black
             jobname: Identifier for this job. Will be used to name the folder that stores snapshots as well. 
         """
         self.N = N
         self.M = M
         self.size = size
         self.iterations = iterations
         self.strains = strains
         self.cmap = cmap
         self.jobname = jobname

         self.antibiotic_concentrations = dict()
         self.degrader_concentrations = dict()
         self.lattice = np.zeros([self.size, self.size])
         self.abundance_statistics = []

     def initialize_lattice(self):
         """Sporulation of individuals at random points on the lattice. This is the (random) initial condition for the ecological dynamics.
         Also construct the matrices recording antibiotic and degrader concentrations at each lattice point. """

         # initializing matrices for antibiotic and degrader concentration
         for i in range(self.M):
             akey = 'antibiotic' + str(i)
             self.antibiotic_concentrations[akey] = np.zeros([self.size, self.size])
 
             dkey = 'degrader' + str(i)
             self.degrader_concentrations[dkey] = np.zeros([self.size, self.size])
 
         # populate the lattice with individuals
         for row in range(self.size):
             for column in range(self.size):
                 # generate a random number to decide which (if any) individual occupies lattice point (row,column)
                 occupation_decider = random.uniform(0, 1)
                 if occupation_decider > 1 - initial_occupation_fraction:
                     continue  # empty - self.lattice[i][j] remains at 0
                 for strain_idx in range(1, N+1):
                     # individuals are represented by numbers, with each number corresponding uniquely to a strain
                     if (strain_idx-1)*(initial_occupation_fraction/N) <= occupation_decider < strain_idx*(initial_occupation_fraction/N):
                         self.lattice[row][column] = strain_idx
                         break
  
     def update_lattice(self):
         """Updating the state of each grid cell on the lattice. """
         current_state = deepcopy(self.lattice)
         updated_state = deepcopy(self.lattice)

         # First, we wipe out all the small molecules from the previous timestep. This is done to ensure that there is no unrealistic 
         # buildup around producers and degraders. This can also be done by setting up a mechanism by which these small molecules decay. 
         for i in range(self.M):
             self.antibiotic_concentrations['antibiotic' + str(i)] = np.zeros([self.size, self.size])

             self.degrader_concentrations['degrader' + str(i)] = np.zeros([self.size, self.size])

         # Now, given the positions of all current individuals on the lattice (current_state), we let the microbes all secrete their 
         # respective small molecules. Here we assume that within each timestep, secretion takes place on a small enough timescale compared 
         # to diffusion that we can separate these processes - first secretion, then diffusion.

         for row in range(self.size):
             for col in range(self.size):
                 cell_occupant = int(current_state[row][col])

                 # If the cell is empty, then there is nothing to be updated
                 if cell_occupant == 0:
                     continue
                 # If it is not empty, then we record the strain that is present on this grid cell in `individual_type`
                 else:
                     individual_strain = self.strains[cell_occupant-1] # strains number 1 to N, but Python indexing begins from 0(0 in the lattice means empty)
                     individual_phenotype = individual_strain.phenotype # phenotype of this strain
     
                     # secretion of antibiotics and degrader molecules: The individual produces the antibiotic and degrader   
                     # molecules - of volume u_p or u_d resp. - on the grid cell that the  individual occupies
     
                     for i in range(self.M): 
                         phen_i = individual_phenotype[i]
                         if phen_i == 1: # add u_p units of this antibiotic to the grid cell where this individual is present
                             antibiotic_lattice = self.antibiotic_concentrations['antibiotic' + str(i)]
                             antibiotic_lattice[row][col] += u_p
                         elif phen_i == 3:  # add u_d units of this degrader molecule to the grid cell where this individual is present
                             degrader_lattice = self.degrader_concentrations['degrader' + str(i)]
                             degrader_lattice[row][col] += u_d

         # We diffuse each small molecule independently using a Gaussian filter. Here we assume that the diffusion of one molecule 
         # does not depend on the diffusion of another, and that diffusion is not impeded by the density of microbes. 
         # Post diffusion, we allow the antibiotics and degraders to react with each other. 
         
         for i in range(self.M):
             # these molecules must now diffuse to nearby locations. 
             # diffusivities sigma are divided by 3 to coincide with the meaning of K_P and K_D in the well-mixed model 
             # i.e., the region around which the antibiotic is not effective. We assume that the antibiotic concentration 
             # is not sufficiently high to cause inhibition 3 standard deviations away from the centre.
               
             A_conc_new = gaussian_filter(self.antibiotic_concentrations['antibiotic' + str(i)], sigma=K_P/3.0, mode='wrap')
             D_conc_new = gaussian_filter(self.degrader_concentrations['degrader' + str(i)], sigma=K_D/3.0, mode='wrap')
 
             # the degrader molecules (if present) must react with the antibiotics and inhibit them
             [A_conc_new, D_conc_new] = antibiotic_degrader_interaction(A_conc_new, D_conc_new, self.size)
 
             self.antibiotic_concentrations['antibiotic' + str(i)] = A_conc_new
             self.degrader_concentrations['degrader' + str(i)] = D_conc_new
         
         # Now we iterate over the grid cells uniformly randomly to decide the fate of each individual in the presence
         # of these small molecules. We pick a random grid cell size^2 times, which implies that each grid cell is 
         # updated on average once. This is the implementation of the overlapping generations assumption. 

         randrows = list(range(self.size)) 
         randcols = list(range(self.size))
         random.shuffle(randrows)
         random.shuffle(randcols) 
         
         for row in randrows:
             for col in randcols: 
                 cell_occupant = int(current_state[row][col])
             
                 # If the cell is empty, then there is nothing to be updated
                 if cell_occupant == 0:
                     continue
                 # If it is not empty, then we record the strain that is present on this grid cell in `individual_type`
                 else: 
                     individual_strain = self.strains[cell_occupant-1] # strains number 1 to N, but Python indexing begins from 0(0 in the lattice means empty)
                     individual_phenotype = individual_strain.phenotype # phenotype of this strain

                     # death and birth are implemented below (in that order). 

                     if random.uniform(0,1) < dose_response(individual_strain, self.antibiotic_concentrations, 
                                                                            self.M, row, col, tau=threshold, k=susceptibility):
                         updated_state[row][col] = 0 # death
                 
                     # if the focal individual dies in the death step, then it obviously cannot subsequently give birth. 
                     if updated_state[row][col] == 0:
                         continue
                     else:
                         birth_prob = individual_strain.growth_rate
                         D = DispersalNbhd(updated_state, row, col, self.size) # D is the dispersal neighbourhood
     
                         if D == -1: # if there are no empty cells in the Moore neighbourhood
                             continue
                         else:
                             if random.uniform(0,1) < birth_prob:
                                 dispersal_idx = D[random.randrange(0, len(D))]
                                 updated_state[dispersal_idx[0]][dispersal_idx[1]] = current_state[row][col] #birth
         
         self.lattice = updated_state
         return self
 
     def update_statistics(self):
         """Finding abundances of each strain after a given generation. """
 
         all_gridcells = list(it.chain.from_iterable(self.lattice))
         # total population size is number of cells that are occupied i.e., number of cells not empty
         #  total_popn_size = self.size**2 - all_gridcells.count(0)
         current_abundance = [0.0]*N

         for strain_idx in range(1, N+1):
             # count number of cells which have individuals of a given strain
             current_abundance[strain_idx - 1] = all_gridcells.count(strain_idx)/self.size**2
         self.abundance_statistics.append(current_abundance)

         return self

     def save_lattice_snapshot(self, run, generation):
         """Save an image of the lattice into a separate folder after each generation so that it is 
         easy to check progress while the code is running. """
        
         data = self.lattice
         # create figures of the same size - will be helpful when creating movie
         fig, axs = plt.subplots(1, 1,
                             constrained_layout=True, squeeze=False)
         # plot the colormesh using the array self.lattice
         for [ax, cmap] in zip(axs.flat, [self.cmap]):
             psm = ax.pcolormesh(data, cmap=cmap, rasterized=True, vmin=0, vmax=4)
         plt.axis('off')
         plt.savefig('results/simulation_snapshots/' + community_identifier  + '/' + jobname + '/' + 'run_' + str(run) + '_generation_' + str(generation) + '.jpg')
         plt.close(fig)

         return self


###################################################### RUN MODEL ##################################################################

# This is the main section that combines the classes defined and parameters set to perform the simulations. 


def run_model(R):
     # Runs one simulation of the lattice model and saves the snapshots, movie, and returns an array with the abundances of the N strains at 
     # each time point i.e., an array of shape [update_iteration,N]
     
     # call the Lattice_Simulation class and construct the grid.
     grid = Lattice_Simulation(N=N,
                 M=M,
                 size=size,
                 iterations=update_iterations,
                 strains=strains,
                 cmap=cmap,
                 jobname=jobname)

     # place spores on the grid

     grid.initialize_lattice() 
     grid.save_lattice_snapshot(R,0)

     # update the lattice for update_iteration many generations

     if save_snapshots==True and verbose==True:
         for generation in range(1,update_iterations+1):
             print('generation number ', generation, ' started out of ', update_iterations)
             grid.update_lattice()
             grid.update_statistics()
             grid.save_lattice_snapshot(R,generation)

     if save_snapshots==True and verbose==False:
         for generation in range(1,update_iterations+1):
             grid.update_lattice()
             grid.update_statistics()
             grid.save_lattice_snapshot(R,generation)

     if save_snapshots==False and verbose==True:
         for generation in range(1,update_iterations+1):
             print('generation number ', generation, ' started out of ', update_iterations)
             grid.update_lattice()
             grid.update_statistics()

     if save_snapshots==False and verbose==False:
         for generation in range(1,update_iterations+1):
             grid.update_lattice()
             grid.update_statistics()
  
     # compile snapshots into a movie only if the save_snapshots flag is set to True (default False)
     if save_movie==True and R==0:
         create_movie(figuredir = 'results/simulation_snapshots/' + community_identifier  + '/' + jobname, 
                                movie_name = jobname + '_run' + str(R) + '.mp4', 
                                iterations=update_iterations, 
                                savedir = 'results/simulation_movies/',
                                run=R) 
     
     return [R,grid.abundance_statistics]


def record_result(result):
     # required for the "callback" arg in the multiprocessing.apply_async function to record the result given by each CPU 
     global abundance_stats_over_runs
     abundance_stats_over_runs.append(result)

# Now we run no_runs independent simulations of the lattice model on a parallelized for loop
# using the multiprocessing module

abundance_stats_over_runs = []
pool = mp.Pool(no_runs) # number of CPUs

for run in range(1,no_runs+1):
     #  print('run number ',run, 'started out of ', no_runs)
     pool.apply_async(run_model, args=[run], callback=record_result)
     
pool.close() # close the processes
pool.join()  # finalise the processes by recording results

# plot abundance statistics

fig, ax = plt.subplots(1, 1)
ax.set_ylim([0,1])
legend_labels = strain_labels

# first we find the "average" trajectory for each of the strains

avg_trajectory = np.zeros([update_iterations,N]) # update_iterations rows, N columns

for t in range(update_iterations): 
     # for each time point t, find the average abundance of each strain at t
     for strain_idx in range(N):
         avg_trajectory[t][strain_idx] = np.mean([traj[1][t][strain_idx] for traj in abundance_stats_over_runs])
       
# plot the average trajectory

for strain_idx in range(N):
     strain_abundance = [timestep[strain_idx] for timestep in avg_trajectory]
     ax.plot(range(update_iterations), strain_abundance, 
                   label = legend_labels[strain_idx], 
                   color = grid_colors[strain_idx+1],
                   linewidth=3.3)

# now we plot the statistics for all runs

for run in range(1,no_runs+1):
     run_stats = [run_result[1] for run_result in abundance_stats_over_runs if run_result[0]==run][0]
     for strain_idx in range(N):
         strain_abundance = [timestep[strain_idx] for timestep in run_stats]
         ax.plot(range(update_iterations), strain_abundance, 
                       color = grid_colors[strain_idx+1],
                       linewidth=0.3)  

ax.legend()
plt.suptitle("Abundance statistics of the strains over time")
plt.savefig('results/simulation_statistics/' + community_identifier  + '/plots/' + jobname + '.jpg')
plt.close(fig)

# record average trajectory of abundances. This will take the average abundances of each strain over all the runs and pickle the list. 
# Each jobname makes one file named *_endstats". This is not a text file and cannot be directly human-read. To get it back, follow the syntax 

for r in range(1,no_runs+1):
     with open("results/simulation_statistics/" + community_identifier  + '/stats/' + jobname + 'run' + str(r) + "_stats", 'wb') as filname:
         pickle.dump([stats for [R, stats] in abundance_stats_over_runs if R==r][0], filname)     
 
# Note: This is not a text file and cannot be directly human-read. To get it back and see it, follow the syntax: 
# with open ('outfile', 'rb') as filname:
#      itemlist = pickle.load(filname)
# print(itemlist)
