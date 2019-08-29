#~/usr/bin/env python
#-*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import os
import numpy as np

# set global settings
def init_plotting():
    plt.rcParams['figure.figsize'] = (4, 3)
    plt.rcParams['font.size'] = 8
    #plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['axes.titlesize'] = plt.rcParams['font.size']
    plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
    plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
    plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']
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

def plot(detailed_mech_mwd_path, frag_mech_mwd_path, mw_cuts, model):

    y1 = get_mw_distri_data(detailed_mech_mwd_path, mw_cuts)
    y2 = get_mw_distri_data(frag_mech_mwd_path, mw_cuts)

    x1 = np.arange(len(y1))*2+0.8
    x2 = np.arange(len(y1))*2
    
    init_plotting()

    plt.figure()

    fig, ax = plt.subplots()
    plt.bar(x1, y1, label='Detailed Mechansim', color='#2c7fb8')
    plt.bar(x2, y2, label='Fragment Mechanism', color='#66c2a4')

    ax.set_xticks(x1)
    
    # xtick_labels = ['<C5', 'C5-C15', 'C15-C25', 'C25-C35', '>C35']
    xtick_labels = []
    for i, mw_cut in enumerate(mw_cuts):
        if i == 0:
            xtick_labels.append('<{0}'.format(mw_cut))
        else:
            xtick_labels.append('{0}-{1}'.format(mw_cuts[i-1], mw_cut))

        if i == len(mw_cuts) - 1:
            xtick_labels.append('>{0}'.format(mw_cut))
    
    ax.set_xticklabels(xtick_labels)
    plt.xlabel('Molecular Weight (g/mol)')
    plt.ylabel('Mole Fraction')

    plt.gca().legend()

    plt.tight_layout()
    plt.savefig('mwd_comparison_{0}.pdf'.format(model))


def get_mw_distri_data(data_filepath, mw_cuts):

    full_path = os.path.join(data_filepath)

    from numpy import genfromtxt
    model_data = genfromtxt(full_path, delimiter=' ')

    mws, molfracs = model_data[0], model_data[1]

    # initialize
    agg_mw_distri = [0.0] * (len(mw_cuts)+1)

    for i, mw in enumerate(mws):
        mw_cut_idx = find_which_mw_cut(mw, mw_cuts)
        agg_mw_distri[mw_cut_idx] += molfracs[i]

    return agg_mw_distri

def find_which_mw_cut(mw, mw_cuts):

    for i, mw_cut in enumerate(mw_cuts):
        if mw <= mw_cut:
            return i
    
    return i + 1


model = 'two-sided_newcut1'
# mw_cuts = [70, 210, 350, 490]
mw_cuts = [30, 150, 350, 550]
detailed_mech_mwd_path = os.path.join('../', 'data', 'pdd_chemistry', 
                                      'detailed', 'pdd_2014_pruning4_s4_a3ene_c11',
                                      'results', 'mwd.csv')
frag_mech_mwd_path = os.path.join('../', 'data', 'pdd_chemistry', 
                                  model,
                                  'results', 'mwd_14hr.csv')

plot(detailed_mech_mwd_path, frag_mech_mwd_path, mw_cuts, model)
