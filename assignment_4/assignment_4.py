"""Plot output for assignment 3"""
import numpy as np
from astropy.table import Table
from matplotlib import pyplot as plt
import os
import pylab
import shutil

symbol = ['*', '+', '^', 's', 'o', '>', '<', '.']
colours = ['k', 'r', 'y', 'g', 'm', 'c', 'r', 'k']

def results_1d(question_num, adv):
    # iteration variables
#    adv = "QUICK"
    lin = "3"
    '''---------------------------------------------------------------------'''

    # Load outpy.txt, formatted output file from fortran assignment_4.prj
    output = Table.read("outpym.txt", format='ascii.commented_header',
                        guess=False)
    con_vol = len(output['I'])-2
    path = "C:\\Users\\Alex\\Documents\\GitHub\\CFD_course\\assignment_4\\{}\\".format(question_num)
    #if an=="rev":
    #    analytic = Table.read(path+"analytic_rev.txt", format='ascii.commented_header',
    #                          guess=False)
    #else:
    #    analytic = Table.read(path+"analytic.txt", format='ascii.commented_header',
    #                          guess=False)
    file_name = "{}_CVL-{}_Advection-{}".format("PDist", con_vol, adv)
    data_name = "output_advection_{}_CV_{}".format(adv, con_vol)

    if not os.path.exists(path):
        print "New_paths"
        os.makedirs(path)
    shutil.copy2('outpym.txt', path +"\\"+data_name+'_output.txt')

    print "------- ", adv, " -------"
    print "Number of CVs: ", con_vol  # -2 for IB-1 and IE+1
    print "pressure, velocity"
    print output['P'], output['U']
    print "max velocity: ", max(output['U'])
    
    # Plot temperature as a function of x position along the bar
    x_fig = plt.figure(figsize=(10,10))
    ax1 = x_fig.add_subplot(111)
    #ax1.plot(output['XP'], output['T'])
    ax1.scatter(output['XP'], output['P'])
    #ax1.plot(analytic['XP'], analytic['T'])

    # Format figure
    ax1.set_xlabel("x [m]", fontweight='bold', fontsize=14)
    ax1.set_ylabel("P (x) [Pa]", fontweight='bold', fontsize=14)
    #ax.xaxis.set_major_locator(plt.MultipleLocator(1))

    plt.ylim(min(output["P"]),#-0.01*min(output['T']),
             max(output["P"]))
    # plt.xlim(0, max(output['XP']))

    if not os.path.exists(path):
        print "New_path"
        os.makedirs(path)
    pylab.savefig(os.path.join(path, file_name))
    plt.show()
    
    u_fig = plt.figure(figsize=(10,10))
    ax2 = u_fig.add_subplot(111)
    ax2.scatter(output['XP'], output['U'])
    # Format figure
    ax2.set_xlabel("x [m]", fontweight='bold', fontsize=14)
    ax2.set_ylabel("U (x) [m/s]", fontweight='bold', fontsize=14)
    plt.ylim(min(output['U']),
             max(output['U']))

    file_name = "{}_CVL-{}_Advection-{}".format("UDist", con_vol, adv)
    if not os.path.exists(path):
        print "New_path"
        os.makedirs(path)
    pylab.savefig(os.path.join(path, file_name))
    plt.show()
    
    return

def composite_plots(problem, cv):
    '''Make composites of each mesh with analytic solution'''
    path = "C:\\Users\\Alex\\Documents\\GitHub\\CFD_course\\assignment_3\\{}\\".format(problem)

    composite = plt.figure(figsize=(12,12))
    ax = composite.add_subplot(111)
    cvs = cv
    scheme = ["UDS", "CDS", "QUICK"] 
    analytic = Table.read(path + "analytic_rev.txt",
                                format="ascii.commented_header", guess=False)

    ax.plot(analytic["XP"], analytic["T"], color='b',
            label="Analytic")

    for i in range(0, len(scheme)):
        results = Table.read(path + "output_advection_"+scheme[i]+"_CV_"+cvs+"_output.txt",
                             format='ascii.commented_header', guess=False)
        legend = scheme[i]+" Scheme"
        if scheme[i]=="CDS":
            order=3
        else:
            order=i
        ax.scatter(results['XP'], results['T'], color=colours[i],
                   marker=symbol[i], s=700/((i+1)*5), label=legend, zorder=order)

    ax.legend(loc='lower left', fontsize=12, scatterpoints=1)
    ax.set_xlabel('x [m]', fontweight='bold', fontsize=14)
    ax.set_ylabel('T(x) [K]', fontweight='bold', fontsize=14)

    file_name = "composite_plot_{}_{}_rev".format(problem, cvs)
    pylab.savefig(os.path.join(path, file_name))
    plt.show()
    plt.close()

    return ()
    
def final_temp():
    path = "C:\\Users\\Alex\\Documents\\GitHub\\CFD_course\\assignment_3\\results\\"

    final_plt = plt.figure(figsize=(12,12))
    ax = final_plt.add_subplot(111)
    num_ts = ["2", "4", "8", "16", "32"]
    imp_final = [19.22635, 14.9879, 12.59924, 11.32753, 10.67393]
    cn_final = [8.0186, 9.52, 9.8809, 9.9725, 9.9938]
    analytic_end = [10, 10, 10, 10, 10]

    ax.scatter(num_ts, imp_final, color='r', marker='*', s=100,
               label="Final Implicit Center Temperatures")
    ax.scatter(num_ts, cn_final, color='b', marker='*', s=100,
               label="Final Crank-Nicolson Center Temperatures")
    ax.plot(num_ts, analytic_end, color='k',
            label="Final Analytic Center Temperatures")
    
    ax.legend(loc='lower right', fontsize=12, scatterpoints=1)
    ax.set_xlabel('Number of Time-steps', fontweight='bold', fontsize=14)
    ax.set_ylabel('Final Center Temperature', fontweight='bold', fontsize=14)
    
    file_name = "end_temp_plot"
    pylab.savefig(os.path.join(path, file_name))
    plt.show()
    plt.close()
    return
