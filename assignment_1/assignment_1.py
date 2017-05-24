"""Plot output for assignment 1"""
import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
import os
import pylab
import shutil

symbol = ['*', '+', '^', 's', 'o', '>', '<', '.']
colours = ['k', 'b', 'y', 'g', 'm', 'c', 'r', 'k']

def results_1d(question_num):
    # iteration variables
    di = "00025"
    t_o = "100"
    lin = "3"
    '''---------------------------------------------------------------------'''

    # Load outpy.txt, formatted output file from fortran assignment_1.prj
    output = Table.read("outpy.txt", format='ascii.commented_header',
                        guess=False)
    path = "C:\\Users\\Alex\\Documents\\GitHub\\CFD_course\\assignment_1\\{}".format(question_num)
    file_name = "{}_CV-{}_CVL-{}_To-{}".format("Node",
                                               len(output['I'])-2, di, t_o)
    data_name = "output_linearization_{}_CV_{}".format(lin,
                                                       len(output['I'])-2)
    if not os.path.exists(path):
        print "New_path"
        os.makedirs(path)
    shutil.copy2('outpy.txt', path +"\\"+data_name+'_output.txt')

    print "Number of CVs: ", len(output['I']) - 2   # -2 for IB-1 and IE+1
    print "min temp: ", min(output['T']), output['XP'][output['T']==min(output['T'])]

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

def composite_plots(problem):
    '''Make composites of each mesh with analytic solution'''
    lin = "3"
    path = "C:\\Users\\Alex\\Documents\\GitHub\\CFD_course\\assignment_1\\{}\\".format(problem)

    composite = plt.figure(figsize=(12,12))
    ax = composite.add_subplot(111)
    num_cvs = ["1", "2", "4", "8", "16", "32", "64"] #, "65"]
    converged = "4"

    for i in range(0, len(num_cvs)):
        results = Table.read(path + "output_linearization_"+lin+"_CV_"+num_cvs[i]+"_output.txt",
                            format='ascii.commented_header', guess=False)
        if i==0:
            legend = num_cvs[i]+" CV"
        else:
            legend = num_cvs[i]+" CVs"

        if num_cvs[i] == converged:
            ax.plot(results['XP'], results['T'], color=colours[i],
                    label=legend+", Converged")
        else:
            ax.scatter(results['XP'], results['T'], color=colours[i],
                       marker=symbol[i], s=700/((i+1)*5), label=legend)

    ax.legend(loc='lower left', fontsize=14, scatterpoints=1)
    ax.set_xlabel('x position [m]', fontweight='bold', fontsize=14)
    ax.set_ylabel('T(x) [K]', fontweight='bold', fontsize=14)

    file_name = "composite_{}_lin{}".format(problem, lin)
    pylab.savefig(os.path.join(path, file_name))
    plt.show()
    plt.close()

    if problem != "problem_4":  
        analytic = Table.read(path + "analytic.txt",
                              format="ascii.commented_header", guess=False)
        analytic_plot = plt.figure(figsize=(12,12))
        ax1 = analytic_plot.add_subplot(111)
        converged_data = Table.read(path + "output_linearization_"+lin+"_CV_"+converged+"_output.txt",
                                    format='ascii.commented_header', guess=False)

        ax1.scatter(converged_data["XP"], converged_data["T"], marker='o',
                    color='b', s=50, label=converged+" CVs, Converged Solution")
        ax1.plot(analytic["XP"], analytic["T"], color='k',
                 label="Analytic Solution")
        ax1.legend(loc='lower left', fontsize=14, scatterpoints=1)
        ax1.set_xlabel('x position [m]', fontweight='bold', fontsize=14)
        ax1.set_ylabel('T(x) [K]', fontweight='bold', fontsize=14)
        file_name = "analytic_{}_lin{}".format(problem, lin)
        pylab.savefig(os.path.join(path, file_name))
        plt.show()
        plt.close()
        
    return ()
