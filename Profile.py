# 1) Call profile_dftb with the input molecule
# 2) Manually edit the dftb_in.hsd file for the input specifications
# 3) Set Input method properly because it is used to construct the folder name where results get stored
# 4) Set iteration count to the no. of times you want to simulate and profile
# 5) Result directory is the one where the result folders would be saved


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import gc
import os

Input_method = 'QR'  # The eigensolver routine present in the input file
Iteration_count = 10  # Number of times for profiling the source code w.r.t. the same input file
Result_directory = 'gpu_utilization_testing/'
Input_list = '2ndInput.txt'  # The file containing the list of molecules for which you want to run your simulation


def profile_dftb(molname):
    Input = molname  # The name of the input molecule whose corresponding xyz file is present in the directory

    # Converting the input xyz file to gen format which will be placed in the dftb_in.hsd file.
    read_file = Input + '.xyz'
    write_file = 'geo.gen'
    fp_read = open(read_file, "r")
    fp_write = open(write_file, "w")
    count = 1
    linecount = 1
    for ln in fp_read:
        ln = ln.replace('\n', '')
        words = ln.split(' ')
        words = list(filter(None, words))
        # print(words)
        if linecount == 1:
            fp_write.write(ln + ' C\n')
            fp_write.write('C H\n')
        if linecount > 2:
            if len(words) > 2:
                if words[0] == 'C':
                    fp_write.write('\t' + str(count) + '\t1\t' + str(words[1]) + '\t' + str(words[2]) + '\t'
                                   + str(words[3]) + '\n')
                    count = count + 1
                else:
                    fp_write.write('\t' + str(count) + '\t2\t' + str(words[1]) + '\t' + str(words[2]) + '\t'
                                   + str(words[3]) + '\n')
                    count = count + 1
            del words[:]
        linecount = linecount + 1
    fp_read.close()
    fp_write.close()

    # Executing the source code
    for i in range(0, Iteration_count):
        Output_directory = Result_directory + Input + '_' + Input_method + '_' + str(i+1)
        command = 'mkdir ' + Output_directory
        os.system(command)
        run_dftb = 'dftb+/dftb+ dftb_in.hsd'
        print(str(i+1) + ') ' + run_dftb)
        os.system(run_dftb)

    # Profiling the source code
        profile_dftb = 'gprof dftb+/dftb+ > ProfilingReport_' + Input + '_' + str(i+1)
        print(str(i + 1) + ') ' + profile_dftb)
        os.system(profile_dftb)

    # Moving the output files to the output directory
        command = 'mv ProfilingReport_' + Input + '_' + str(i + 1) + ' ' + Output_directory
        print(command)
        os.system(command)
        command = 'cp band.out ' + Output_directory
        print(command)
        os.system(command)
        command = 'cp detailed.out ' + Output_directory
        print(command)
        os.system(command)
        command = 'cp dftb_in.hsd ' + Output_directory
        print(command)
        os.system(command)
        command = 'cp ' + Input + '.xyz ' + Output_directory
        print(command)
        os.system(command)
        command = 'cp dftb_pin.hsd ' + Output_directory
        print(command)
        os.system(command)
        command = 'cp geom.out.xyz ' + Output_directory
        print(command)
        os.system(command)
        command = 'cp geo.gen ' + Output_directory
        print(command)
        os.system(command)
        command = 'cp gmon.out ' + Output_directory
        print(command)
        os.system(command)
        command = 'cp md.out ' + Output_directory
        print(command)
        os.system(command)


def line_plot(x_axis, y_axis, method):
    x_axis = np.array(x_axis)
    # print(x_axis)
    x_axis = list(map(int, x_axis))
    # print(x_axis)
    y_axis = np.array(y_axis)
    # print(y_axis)
    y_axis = list(map(float, y_axis))
    # print(y_axis)
    y_axis = [x for _, x in sorted(zip(x_axis, y_axis))]
    # print(y_axis)
    x_axis.sort()
    # print(x_axis)
    if method == 'QR':
        plt.plot(x_axis, y_axis, color='blue')
    elif method == 'DC':
        plt.plot(x_axis, y_axis, color='red')
    elif method == 'RR':
        plt.plot(x_axis, y_axis, color='green')
    else:
        plt.plot(x_axis, y_axis, color='orange')


def draw_graph():
    dir_list = next(os.walk(Result_directory))[1]
    # print(dir_list)
    # print(len(dir_list))
    molecule = []
    method = []
    graph_filename = Result_directory + '/scalability_plot.png'
    for i in dir_list:
        words = i.split('_')
        molecule.append(words[0])
        method.append(words[1])
    molecule = list(set(molecule))
    method = list(set(method))
    for i in range(0, Iteration_count):
        for j in method:
            x_axis = []
            y_axis = []
            for k in molecule:
                folder_name = Result_directory + k + '_' + j + '_' + str(i+1)
                with open(folder_name + '/geom.out.xyz') as f:
                    first_line = f.readline()  # No. of atoms
                first_line = first_line.replace('\n', '')
                first_line = first_line.replace(' ', '')
                x_axis.append(first_line)
                f = open(folder_name + '/ProfilingReport_' + k + '_' + str(i+1))
                lines = f.readlines()
                f.close()
                lines[5] = lines[5].replace('\n', '')
                sent = lines[5].split(' ')
                sent = list(filter(None, sent))
                if len(sent) > 5:
                    y_axis.append(sent[2])  # Time
                else:
                    y_axis.append('0.00')
            # print(x_axis)
            # print(y_axis)
            line_plot(x_axis, y_axis, j)
    plt.grid()
    blue_patch = mpatches.Patch(color='blue', label='QR')
    red_patch = mpatches.Patch(color='red', label='DC')
    green_patch = mpatches.Patch(color='green', label='RR')
    orange_patch = mpatches.Patch(color='orange', label='MAGMA')
    plt.figlegend(handles=[blue_patch, red_patch, green_patch, orange_patch], loc='upper right', fontsize=9)
    plt.title('Scalability plot (scc=on, v19.1)', fontsize=20)
    plt.xlabel('#Atoms', fontsize=16)
    plt.ylabel('Time (in seconds)', fontsize=16)
    # plt.show()
    plt.savefig(graph_filename, bbox_inches='tight', dpi=100)
    plt.close()
    gc.collect()


def draw_average_graph():
    dir_list = next(os.walk(Result_directory))[1]
    # print(dir_list)
    # print(len(dir_list))
    molecule = []
    method = []
    graph_filename = Result_directory + '/scalability_plot_avg.png'
    for i in dir_list:
        words = i.split('_')
        molecule.append(words[0])
        method.append(words[1])
    molecule = list(set(molecule))
    method = list(set(method))
    for i in method:
        x_axis = []
        y_axis = []
        for j in molecule:
            count = 0
            for k in range(0, Iteration_count):
                folder_name = Result_directory + j + '_' + i + '_' + str(k+1)
                with open(folder_name + '/geom.out.xyz') as f:
                    first_line = f.readline()  # No. of atoms
                first_line = first_line.replace('\n', '')
                first_line = first_line.replace(' ', '')
                f = open(folder_name + '/ProfilingReport_' + j + '_' + str(k+1))
                lines = f.readlines()
                f.close()
                lines[5] = lines[5].replace('\n', '')
                sent = lines[5].split(' ')
                sent = list(filter(None, sent))
                if len(sent) > 5:
                    count = count + float(sent[2])  # Sum up time to take average
                else:
                    count = count + 0
            x_axis.append(first_line)
            y_axis.append(count/Iteration_count)  # Average time
        line_plot(x_axis, y_axis, i)
    plt.grid()
    blue_patch = mpatches.Patch(color='blue', label='QR')
    red_patch = mpatches.Patch(color='red', label='DC')
    green_patch = mpatches.Patch(color='green', label='RR')
    orange_patch = mpatches.Patch(color='orange', label='MAGMA')
    plt.figlegend(handles=[blue_patch, red_patch, green_patch, orange_patch], loc='upper right', fontsize=9)
    plt.title('Scalability plot (scc=on, v19.1) [Avg]', fontsize=20)
    plt.xlabel('#Atoms', fontsize=16)
    plt.ylabel('Time (in seconds)', fontsize=16)
    # plt.show()
    plt.savefig(graph_filename, bbox_inches='tight', dpi=100)
    plt.close()
    gc.collect()


def main():
    fp_read = open(Input_list, "r")
    for ln in fp_read:
        ln = ln.replace('\n', '')
        profile_dftb(ln)  # Profile dftb+
    fp_read.close()
    # profile_dftb('Propane')
    # draw_graph()  # Plot scalability graph
    # draw_average_graph()  # Plot scalability graph by taking the average of the 10 cases


main()






