import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from typing import List
import numpy as np




def generate_cmap(color_boundary: List = None, colors: List = None) :
    
    if color_boundary is None :
        boundaries = [0, -np.log10(0.05), 10]
    else :
        boundaries = color_boundary

    if colors is None :
        colors = ["#97C3FB", "#EEEEE1","#FF7777"]
    else :
        colors = colors
    
    norm=plt.Normalize(min(boundaries),max(boundaries))
    tuples = list(zip(map(norm,boundaries), colors))
    cmap = LinearSegmentedColormap.from_list("", tuples)

    return cmap