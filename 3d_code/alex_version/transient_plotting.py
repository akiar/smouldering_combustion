import numpy as np 
from astropy.table import Table
from matplotlib import pyplot as plt
'''Global variables'''
# Set number of time steps
time_steps = 240
start_step = 1
num_steps_t = 5
num_steps_y = 10
delta_t = 50.0

# Set number of control volumes
temp = "ts"
x_cv = 21 - 2 + 1
y_cv = 81 - 2 + 1
z_cv = 21 - 2 + 1 
colours = ["b", "k", "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b",
           "k", "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b", "k",
           "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b", "k", "r", 
           "b", "k", "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b",
           "k", "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b", "k",
           "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b", "k", "r", 
           "y", "g", "m", "b", "k", "r", "y", "g", "m","b", "k", "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b",
           "k", "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b", "k",
           "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b", "k", "r", 
           "b", "k", "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b",
           "k", "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b", "k",
           "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b", "k", "r", 
           "y", "g", "m", "b", "k", "r", "y", "g", "m","b", "k", "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b",
           "k", "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b", "k",
           "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b", "k", "r", 
           "b", "k", "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b",
           "k", "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b", "k",
           "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b", "k", "r", 
           "y", "g", "m", "b", "k", "r", "y", "g", "m","b", "k", "r", "y", "g", "m", "b", "k", "r", "y", "g", "m", "b",
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
           ">", "^", "s", "o", ".", "*", ">", "^", "s","o", ".", "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o",
           ".", "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o", ".",
           "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o", ".", "*",
           "o", ".", "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o",
           ".", "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o", ".",
           "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o", ".", "*",
           ">", "^", "s", "o", ".", "*", ">", "^", "s","o", ".", "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o",
           ".", "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o", ".",
           "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o", ".", "*",
           "o", ".", "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o",
           ".", "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o", ".",
           "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o", ".", "*",
           ">", "^", "s", "o", ".", "*", ">", "^", "s","o", ".", "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o",
           ".", "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o", ".",
           "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o", ".", "*",
           "o", ".", "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o",
           ".", "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o", ".",
           "*", ">", "^", "s", "o", ".", "*", ">", "^", "s", "o", ".", "*",
           ">", "^", "s", "o", ".", "*", ">", "^", "s"]
'''========================================================================='''

def plot_all():
    print "y/T"
    y_plot()
    print "t/T"
    time_plot()
    return

def y_plot(): 
    '''main routine for plotting transient results from alex_version.prj'''
    # Set up figure
    fig = plt.figure(figsize=(15,8))
    ax = fig.add_subplot(111)

    # Loop over all time steps. Set number in "Global variables" 
    for i in range(start_step,time_steps+1,num_steps_y):
        # Load py_i.txt file
        transient = load_data(i)

        # assign variables
        yp = transient["yp"]
        T = transient[temp] # "tf" or "ts"

        # Make plot
        if i <= 5: 
            ax.scatter(yp, T, color = colours[i], marker=markers[i-1],
                       label=i)
        else:
            ax.plot(yp, T, color = colours[i-1], label=i)
    # End of loop 

    # format plot
#    ax.legend(loc='upper right', fontsize=11, scatterpoints=1)
    ax.set_ylabel("TF [K]", fontweight='bold', fontsize=14)
    ax.set_xlabel("y [m]", fontweight='bold', fontsize=14)
    return

def time_plot():
    '''main routine for plotting transient results from alex_version.prj'''
    fig = plt.figure(figsize=(15,8))
    ax = fig.add_subplot(111)

    for i in range(start_step,time_steps+1):   # set desired looping
        transient = load_data(i)
        time = i*delta_t
        T = transient[temp]
        for j in range(0,y_cv,num_steps_t):
            ax.plot(time, T[j], color = colours[j], marker=markers[j],
                       label=j)
#    ax.legend(loc='upper right', fontsize=11, scatterpoints=1)
    ax.set_ylabel("TF [K]", fontweight='bold', fontsize=14)
    ax.set_xlabel("Time [s]", fontweight='bold', fontsize=14)        
    return


def load_data(time_step):
    '''subroutine for loading data files. 
       Modify thermocouple location you are plotting by changing the 
           center_data criteria
       plot    1: y/T 2: t/T'''
    data = Table.read('py_'+str(time_step)+'.txt',
                      format='ascii.commented_header', guess=False)
    # Select only data through the center of the column
    center_data = data[data["i"]==(x_cv/2)]
    return(center_data)