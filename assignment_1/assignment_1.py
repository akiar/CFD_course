"""Plot output for assignment 1"""
import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt

def results_1d():
    # Load outpy.txt, formatted output file from fortran assignment_1.prj
    output = Table.read("outpy.txt", format='ascii.commented_header',
                        guess=False)

    # Plot temperature as a function of node value along the bar
    temp_fig = plt.figure(figsize=(15,10))
    ax = temp_fig.add_subplot(111)
    ax.plot(output['I'], output['T'])
    ax.scatter(output['I'], output['T'])

    # Format figure
    ax.set_xlabel("Node Number")
    ax.set_ylabel("Node Temperature")
    ax.xaxis.set_major_locator(plt.MultipleLocator(1.0))
    ax.yaxis.set_major_locator(plt.MultipleLocator(10.0))
    plt.xlim(0, max(output['I']+1))

    return()