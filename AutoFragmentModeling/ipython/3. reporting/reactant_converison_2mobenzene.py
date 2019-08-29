#~/usr/bin/env python
#-*- coding: utf-8 -*-
import os
import numpy as np
import matplotlib.pyplot as plt

# set global settings
def init_plotting():
    plt.rcParams['figure.figsize'] = (4, 3)
    plt.rcParams['font.size'] = 8
    #plt.rcParams['font.family'] = 'Helvetica'
    plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = 1.5*plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']
    #plt.rcParams['savefig.dpi'] = 2*plt.rcParams['savefig.dpi']
    plt.rcParams['xtick.major.size'] = 3
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['xtick.minor.width'] = 1
    plt.rcParams['ytick.major.size'] = 3
    plt.rcParams['ytick.minor.size'] = 3
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['ytick.minor.width'] = 1
    plt.rcParams['legend.frameon'] = True
    plt.rcParams['legend.loc'] = 'best'
    plt.rcParams['axes.linewidth'] = 1

    plt.gca().spines['right'].set_color('none')
    plt.gca().spines['top'].set_color('none')
    plt.gca().xaxis.set_ticks_position('bottom')
    plt.gca().yaxis.set_ticks_position('left')

def load_comparison_data(detailed_model, frag_model1, frag_model2=None):

    v0_csv = os.path.join('../', 'data', 'pdd_chemistry', 
                          'detailed', detailed_model,
                          'results', 'reactant_conv.csv')

#    v0_csv = os.path.join('../', 'data', 'pdd_chemistry', 
#                          'two-sided_newcut1',
#                          'results', 'reactant_conv.csv')

    v0_data = []
    with open(v0_csv, 'r') as read_in:
        for line in read_in:
            tokens = line.split(' ')
            entries = [float(token) for token in tokens]
            v0_data.append(entries)

    assert len(v0_data) == 2

    v1_csv = os.path.join('../', 'data', '2mobenzene', 
                          frag_model1,
                          'results', 'reactant_conv.csv')
    v1_data = []
    with open(v1_csv, 'r') as read_in:
        for line in read_in:
            tokens = line.split(' ')
            entries = [float(token) for token in tokens]
            v1_data.append(entries)

    assert len(v1_data) == 2

    if frag_model2:
        v2_csv = os.path.join('../', 'data', 'pdd_chemistry', 
                          frag_model2,
                          'results', 'reactant_conv.csv')
        v2_data = []
        with open(v2_csv, 'r') as read_in:
            for line in read_in:
                tokens = line.split(' ')
                entries = [float(token) for token in tokens]
                v2_data.append(entries)

        assert len(v2_data) == 2

        return np.array(v0_data), np.array(v1_data), np.array(v2_data)

    else:
        return np.array(v0_data), np.array(v1_data), None

def plot_comparison(v0_data, v1_data, v2_data=None,
                    detailed_model=None,
                    frag_model1=None,
                    frag_model2=None,
                    xlabel='', 
                    ylabel='',
                    figure_name='', 
                    xlim=10, ylim=1.0):

    init_plotting()

    plt.figure()
    plt.plot(v0_data[0]/3600.0, v0_data[1], label='Detailed: {0}'.format(detailed_model))
    plt.plot(v1_data[0]/3600,v1_data[1], label='Fragment: {0}'.format(frag_model1))
    if v2_data:
        plt.plot(v2_data[0]/3600,v2_data[1], label='Fragment: {0}'.format(frag_model2))

    plt.gca().set_xscale('log')


    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim((.1, xlim))
    plt.ylim((0, ylim))

    plt.gca().legend(scatterpoints=1)

    plt.tight_layout()
    plt.savefig(figure_name)


detailed_model = 'pdd_2014_pruning4_s4_a3ene_c11'
frag_model1 = 'one-sided'
frag_model2 = None

if frag_model2:
    figure_name = 'reactant_conversion_{0}_vs_{1}'.format(frag_model1, frag_model2)
else:
    figure_name = 'reactant_conversion_{0}'.format(frag_model1)

# plot reactant conversion
xlabel = 'Time / hr'
ylabel = 'Conversion'
detailed, frag1, frag2 = load_comparison_data(detailed_model, frag_model1, frag_model2)
plot_comparison(detailed, frag1, frag2,
                detailed_model,
                frag_model1, 
                frag_model2,
                xlabel, ylabel, 
                '{0}.pdf'.format(figure_name),
                xlim=14)
