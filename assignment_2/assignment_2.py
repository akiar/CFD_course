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
    di = "025"
    t_o = "100"
    lin = "3"
    con_vol = 40
    '''---------------------------------------------------------------------'''

    # Load outpy.txt, formatted output file from fortran assignment_2.prj
    output = Table.read("outpy.txt", format='ascii.commented_header',
                        guess=False)
    path = "C:\\Users\\Alex\\Documents\\GitHub\\CFD_course\\assignment_2\\{}".format(question_num)
    file_name = "{}_TS-{}_CVL-{}_To-{}".format("Xpos",
                                               (len(output['I'])-2)/con_vol - 1,
                                               di, t_o)
    data_name = "output_linearization_{}_TS_{}".format(lin,
                                                       (len(output['I'])-2)/con_vol - 1)
    if not os.path.exists(path):
        print "New_path"
        os.makedirs(path)
    shutil.copy2('outpy.txt', path +"\\"+data_name+'_output.txt')

    print "Number of TS: ", (len(output['I']) - 2)/con_vol-2  # -2 for IB-1 and IE+1
    print "max temp: ", max(output['T']), output['XP'][output['T']==max(output['T'])]
    print "temperature:"
    print output['T']
    
    # Plot temperature as a function of x position along the bar
    x_fig = plt.figure(figsize=(10,10))
    ax1 = x_fig.add_subplot(111)
    #ax1.plot(output['XP'], output['T'])
    ax1.scatter(output['XP'], output['T'])

    # Format figure
    ax1.set_xlabel("x* [x/L]", fontweight='bold', fontsize=14)
    ax1.set_ylabel("T (x,t) [K]", fontweight='bold', fontsize=14)
    #ax.xaxis.set_major_locator(plt.MultipleLocator(1))

    plt.ylim(min(output['T'])-0.01*min(output['T']),
             max(output['T'])+0.01*max(output['T']))
    # plt.xlim(0, max(output['XP']))

    if not os.path.exists(path):
        print "New_path"
        os.makedirs(path)
    pylab.savefig(os.path.join(path, file_name))
    plt.show()
    
    print "Center temp: ", output['T'][output['I']==1]
    
    return

def composite_plots():
    '''Make composites of each mesh with analytic solution'''
    path = "C:\\Users\\Alex\\Documents\\GitHub\\CFD_course\\assignment_2\\results\\"

    composite = plt.figure(figsize=(12,12))
    ax = composite.add_subplot(111)
    num_ts = ["2", "4", "8", "16", "32"] 
    analytic_start = Table.read(path + "analytic_start.txt",
                                format="ascii.commented_header", guess=False)
    analytic_end = Table.read(path + "analytic_end.txt",
                              format="ascii.commented_header", guess=False)
    ax.plot(analytic_start["XP"], analytic_start["T"], color='k',
            label="Initialized Analytic")
    ax.plot(analytic_end["XP"], analytic_end["T"], color='r',
            label="End Analytic t*")

    for i in range(0, len(num_ts)):
        results = Table.read(path + num_ts[i] + "_ts_final.txt",
                             format='ascii.commented_header', guess=False)
        legend = num_ts[i]+" Time steps"
        ax.scatter(results['XP'], results['T'], color=colours[i],
                   marker=symbol[i], s=700/((i+1)*5), label=legend)

    ax.legend(loc='upper right', fontsize=12, scatterpoints=1)
    ax.set_xlabel('x*', fontweight='bold', fontsize=14)
    ax.set_ylabel('T(x, t) [C]', fontweight='bold', fontsize=14)

    file_name = "composite_plot"
    pylab.savefig(os.path.join(path, file_name))
    plt.show()
    plt.close()
        
    return ()
