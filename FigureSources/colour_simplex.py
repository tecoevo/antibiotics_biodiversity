import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import seaborn as sns
import ternary
import math

def color_point(x, y, z, scale):
      w = 255
      x_color = x * w / float(scale)
      y_color = y * w / float(scale)
      z_color = z * w / float(scale)
      r = z_color / w
      g = y_color / w
      b = x_color / w
      return [r, g, b]

def generate_heatmap_data(scale=90):
      from ternary.helpers import simplex_iterator
      d = dict()
      for (i, j, k) in simplex_iterator(scale):
          [r,g,b] = color_point(i, j, k, 90) 
          d[(i, j, k)] = r*np.array([0.172549, 0.482353, 0.713725]) + g*np.array([0.843137, 0.29, 0.29]) + b*np.array([0.670588, 0.85098, 1])
      return d


#legend 
scale = 90
data = generate_heatmap_data()
figure,tax = ternary.figure(scale=scale)
tax.heatmap(data, style="hexagonal", use_rgba=True, colorbar=False)
tax.left_corner_label("P", fontsize=12)
tax.top_corner_label("S", fontsize=12)
tax.right_corner_label("D", fontsize=12)

# Remove default Matplotlib Axes
tax.clear_matplotlib_ticks()
tax.get_axes().axis('off')
tax.boundary()
plt.savefig("../Manuscript/Working_draft/Figures/colour_simplex.pdf")
plt.show()