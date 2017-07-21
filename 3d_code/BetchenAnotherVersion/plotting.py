'''Plot from CVELO.dat file from BetchenAnotherVersion'''
import numpy as np
import os
import pylab as pylab
from astropy.table import Table, join
from matplotlib import pyplot as plt

'''Global variables'''
directory_path = 'C:\\Users\\Alex\\Documents\\GitHub\\smouldering_combustion\\3d_code\\BetchenAnotherVersion\\'
U_infinity = 5.6e-5 #5.6e-5 6.3167
H = 0.02    #0.02 0.01
rho = 1000.0

'''Main plotting function'''
def porous(file_name, scheme):
    data = Table.read(directory_path+file_name,
                      format='ascii.commented_header', guess=False)

    # Declare variables
    pressure = data['P']
    u_velocity = data['U']
    xp = data['XP']
    print '-------------------------------'
    print 'Left Boundary U: ', u_velocity[0]/U_infinity
    print 'Left Boundary P: ', pressure[0]/rho/U_infinity**2
    print 'Max: ', xp[pressure/rho/U_infinity**2==max(pressure/rho/U_infinity**2)]/H
    print '-------------------------------'

    '''Porous Plug'''
    # Velocity plot
    figure_name = "velocity_plot_pp_{}.png".format(scheme)
    vel_plot = plt.figure(figsize=(10,10))
    ax = vel_plot.add_subplot(111)
    ax.plot(xp/H, u_velocity/U_infinity)
    ax.set_xlabel('x/H', fontweight='bold', fontsize=14)
    ax.set_ylabel('u/U', fontweight='bold', fontsize=14)
#    plt.ylim(1.2,1.6)
#    plt.xlim(0,8)
    pylab.savefig(os.path.join(directory_path, figure_name))
    plt.show()
    plt.close()

    # Pressure plot
    figure_name = "pressure_plot_pp_{}.png".format(scheme)
    press_plot = plt.figure(figsize=(10,10))
    ax2 = press_plot.add_subplot(111)
    ax2.plot(xp/H, pressure/rho/U_infinity**2)
    ax2.set_xlabel('x/H', fontweight='bold', fontsize=14)
    ax2.set_ylabel('P/(rho*u^2)', fontweight='bold', fontsize=14)
    plt.ylim(0,350)
    plt.xlim(0,8)
    pylab.savefig(os.path.join(directory_path, figure_name))
 #   plt.show()
 #   plt.close()

    return

def beavers(file_name, scheme):
    data = Table.read(directory_path+file_name,
                      format='ascii.commented_header', guess=False)

    # Declare variables
    u_velocity = data['U']
    yp = data['YP']
    print '-------------------------------'
    print 'Left Boundary U: ', u_velocity[0]/U_infinity
    print '-------------------------------'

    '''Beavers Johnson Flow'''
    # Velocity plot
    figure_name = "velocity_plot_beaver_{}.png".format(scheme)
    vel_plot = plt.figure(figsize=(10,10))
    ax = vel_plot.add_subplot(111)
    ax.scatter(u_velocity/U_infinity, yp/H)
    ax.set_xlabel('u/U', fontweight='bold', fontsize=14)
    ax.set_ylabel('y/H', fontweight='bold', fontsize=14)
    pylab.savefig(os.path.join(directory_path, figure_name))
    plt.show()
    plt.close()
#    plt.ylim(1.2,1.6)
#    plt.xlim(0,2.6)

    return