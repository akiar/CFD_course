"""Plot output for assignment 1"""
import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
import os
import pylab
import shutil

def results_1d(question_num):
    # iteration variables
    di = "00025"
    t_o = "100"
    lin = "1"
    '''---------------------------------------------------------------------'''

    # Load outpy.txt, formatted output file from fortran assignment_1.prj
    output = Table.read("outpy.txt", format='ascii.commented_header',
                        guess=False)
    path = "C:\\Users\\Alex\\Documents\\GitHub\\CFD_course\\assignment_1\\{}".format(question_num)
    file_name = "{}_CV-{}_CVL-{}_To-{}".format("Node",
                                               len(output['I'])-2, di, t_o)
    data_name = "output_linearization_{}_CV_{}".format(lin,
                                                       len(output['I'])-2)
    shutil.copy2('outpy.txt', path +"\\"+data_name+'_output.txt')

    print "Number of CVs: ", len(output['I']) - 2   # -2 for IB-1 and IE+1
    print "max temp: ", max(output['T']), output['XP'][output['T']==max(output['T'])]

    # Plot temperature as a function of node value along the bar
    node_fig = plt.figure(figsize=(10,10))
    ax = node_fig.add_subplot(111)
    ax.plot(output['I'], output['T'])
    ax.scatter(output['I'], output['T'])

    # Format figure
    ax.set_xlabel("Node Position", fontweight='bold', fontsize=14)
    ax.set_ylabel("Node Temperature", fontweight='bold', fontsize=14)
    #ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    plt.ylim(min(output['T'])-0.01*min(output['T']),
             max(output['T'])+0.01*max(output['T']))
   # plt.xlim(0, max(output['XP']))

    print "temperature:"
    print output['T']

    if not os.path.exists(path):
        print "New_path"
        os.makedirs(path)
    pylab.savefig(os.path.join(path, file_name))
    plt.show()
    
    # Plot temperature as a function of x position along the bar
    x_fig = plt.figure(figsize=(10,10))
    ax1 = x_fig.add_subplot(111)
    ax1.plot(output['XP'], output['T'])
    ax1.scatter(output['XP'], output['T'])

    # Format figure
    ax1.set_xlabel("X Position (m)", fontweight='bold', fontsize=14)
    ax1.set_ylabel("Temperature (K)", fontweight='bold', fontsize=14)
    #ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    plt.ylim(min(output['T'])-0.01*min(output['T']),
             max(output['T'])+0.01*max(output['T']))
    # plt.xlim(0, max(output['XP']))

    file_name = "{}_CV-{}_CVL-{}_To-{}".format("Xpos",len(output['I'])-2,
                                               di, t_o)
    if not os.path.exists(path):
        print "New_path"
        os.makedirs(path)
    pylab.savefig(os.path.join(path, file_name))
    plt.show()
    
    return()