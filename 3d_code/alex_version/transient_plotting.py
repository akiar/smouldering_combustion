import numpy as np 
from astropy.table import Table
from matplotlib import pyplot as plt
'''Global variables'''
# Set number of time steps
time_steps = 60
#
# Set number of control volumes 
x_cv = 21 - 2 + 1
y_cv = 61 - 2 + 1
z_cv = 21 - 2 + 1 
colours = ["b", "k", "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b",
           "k", "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b", "k",
           "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b", "k", "r", 
           "b", "k", "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b",
           "k", "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b", "k",
           "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b", "k", "r", 
           "y", "g", "m", "b", "k", "r", "y", "g", "m"]
markers = ["o", ".", "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o",
           ".", "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o", ".",
           "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o", ".", "*",
           "o", ".", "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o",
           ".", "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o", ".",
           "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o", ".", "*",
           ">", "^", "s", "o", ".", "*", ">", "^", "s"]
'''========================================================================='''

def y_plot(): 
    '''main routine for plotting transient results from alex_version.prj'''
    # Set up figure
    fig = plt.figure(figsize=(8,15))
    ax = fig.add_subplot(111)

    # Loop over all time steps. Set number in "Global variables" 
    for i in range(1,time_steps+1,2):
        # Load py_i.txt file
        transient = load_data(i)

        # assign variables
        yp = transient["yp"]
        T = transient["tf"] # "tf" or "ts"

        # Make plot
        if i <= 5: 
            ax.scatter(T, yp, color = colours[i], marker=markers[i-1],
                       label=i)
        else:
            ax.plot(T, yp, color = colours[i-1], label=i)
    # End of loop 

    # format plot
    ax.legend(loc='upper right', fontsize=11, scatterpoints=1)
    ax.set_ylabel("y [m]", fontweight='bold', fontsize=14)
    ax.set_xlabel("TF [K]", fontweight='bold', fontsize=14)
    return

def load_data(time_step):
    '''subroutine for loading data files. 
       Modify thermocouple location you are plotting by changing the 
           center_data criteria'''
    data = Table.read('py_'+str(time_step)+'.txt',
                      format='ascii.commented_header', guess=False)
    # Select only data through the center of the column
    center_data = data[data["i"]==(x_cv/2)]

    return(center_data)