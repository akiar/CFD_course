"""Plot output for assignment 1"""
import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
import os
import pylab

def results_1d(question_num):
    # iteration variables
    di = "001"
    t_o = "100"
    '''---------------------------------------------------------------------'''

    # Load outpy.txt, formatted output file from fortran assignment_1.prj
    output = Table.read("outpy.txt", format='ascii.commented_header',
                        guess=False)
    path = "C:\\Users\\Alex\\Documents\\GitHub\\CFD_course\\assignment_1\\{}".format(question_num)
    file_name = "{}_CV-{}_CVL-{}_To-{}".format(question_num,
                                               len(output['I'])-2, di, t_o)

    print "Number of CVs: ", len(output['I']) - 2   # -2 for IB-1 and IE+1
    print "Max temp: ", max(output['T']), output['I'][output['T']==max(output['T'])]

    # Plot temperature as a function of node value along the bar
    temp_fig = plt.figure(figsize=(15,10))
    ax = temp_fig.add_subplot(111)
    ax.plot(output['I'], output['T'])
    ax.scatter(output['I'], output['T'])

    # Format figure
    ax.set_xlabel("Node Number")
    ax.set_ylabel("Node Temperature")
    ax.xaxis.set_major_locator(plt.MultipleLocator(1.0))
    plt.ylim(min(output['T'])-0.01*min(output['T']),
             max(output['T'])+0.01*max(output['T']))
    plt.xlim(0, max(output['I']+1))

    print "temperature:"
    print output['T']

    if not os.path.exists(path):
        print "New_path"
        os.makedirs(path)
    pylab.savefig(os.path.join(path, file_name))
    plt.show()
    # plt.close()

    return()