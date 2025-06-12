#!/Users/sergiortizropero/miniconda3/envs/md_traj/bin/python
# Sergi Ortiz @ UAB
# 8 April 2025
# Perform the routine MD analysis of one (or a set of trajectories). Note that
# multiple replicas are not supported (i.e. the analysis is independent to 
# other replicas). If multiple replicas, use the tools available in 
# `analysis`.ipynb (easy analysis)

#===================#
#       IMPORTS     #
#===================#

import os
from numpy import average, std
import numpy as np
from argparse import ArgumentParser

# md analysis
from EMDA import EMDA, __version__
print("EMDA version is:", __version__)

# plotting
import matplotlib.pyplot as plt
plt.style.use('default')
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
rc('lines', lw=1, color='b')
rc('legend', loc='best')
plt.rcParams["legend.fancybox"] = False
plt.rcParams["legend.edgecolor"] = 'black'
plt.rcParams['legend.borderpad'] = 0.25
plt.rcParams['legend.fontsize'] = 11
plt.rcParams.update({'pgf.preamble': r'\usepackage{amsmath}'})


#=====================#
#       FUNCTIONS     #
#=====================#

def parser():
    parser = ArgumentParser(description="Script to adapt AMBER topology and coordinates to be used in ChemShell.")

    parser.add_argument('-i', '--input_',
                        help='MD source directory (i.e. the one that contains the prmtop, inpcrd, pdb, preprod and prod tree)',
                        required=True,
                        type=str
                        )
    parser.add_argument('-load', '--load',
                        help="Wether to load the .pkl file from previous analysis (preprocessed EMDA info). Default False. Do use --load if all the MD analysis have been preprocessed (pkl created).",
                        action='store_true',
                        default=False,
                        required=False,
                        )
    parser.add_argument('-p', '--protein',
                        help='Name of the protein. Default hLOX15',
                        required=False,
                        default='hLOX15', 
                        type=str
                        )
    parser.add_argument('-m', '--multiple',
                        help='Perform the analysis on multiple MD directories. The input -i directory must contain the MD directories to be analyzed.',
                        action='store_true',
                        default=False,
                        required=False,
                        )
    
    args = parser.parse_args()

    return args


# CUSTOM PLOTTING FUNCTIONS

def analyze_MD(md_base_dir, parent_dir):
    def plot_two_metrics_hist(replicate_data, save_name, metric_name=r'', convergence_line_last=5, x_time=True, bin_size=0.1, display_mean=True, skip=None, replica_names=[], legend='', xlims=[]):

        replicate_data = np.array(replicate_data)

        # all have the same length, traj length
        frames = np.arange(1, len(replicate_data[0])+1, 1)
        
        if x_time:
            # each frame corresponds to 10 ps
            # 10 ns, produces 1000 frames, i.e. 0.01 ns / frame
            times = frames * 0.01 

        bins = np.arange(np.min(replicate_data), np.max(replicate_data), bin_size)

        # generate time figure
        fig, axes = plt.subplots(ncols=1, nrows=1, figsize=(3.5, 3.5))
        _axes = [axes]

        # color lists
        colors = ['darkgreen', 'deeppink', 'mediumblue', 'orange']
        colors = ['mediumblue', 'deeppink']

        # histogram
        _axes[0].tick_params(bottom=True, top=True, left=True, right=False,direction="in")
        _axes[0].tick_params(labelbottom=True, labeltop=False, labelleft=True,labelright=False, direction="in")

        # for each replicate
        for j in range(len(replicate_data)):
            # filling
            _axes[0].hist(replicate_data[j], bins=bins, density=True, edgecolor=colors[j], color='w', histtype='stepfilled', alpha=1)

        # for each replicate
        for j in range(len(replicate_data)):

            # contour
            if len(replica_names) != 0:
                if j==0:
                    _axes[0].hist(replicate_data[j], bins=bins, density=True, color=colors[j], alpha=.1, label=r'$\textnormal{'+f'{replica_names[j]}'+r' }$')  
                else: 
                    _axes[0].hist(replicate_data[j], bins=bins, density=True, color=colors[j], alpha=.02, label=r'$\textnormal{'+f'{replica_names[j]}'+r' }$')  
            else:
                _axes[0].hist(replicate_data[j], bins=bins, density=True, color=colors[j], alpha=.6, label=r'$\textnormal{replicate '+f'{j+1}'+r' }$')    
            if display_mean:
                if len(replica_names) != 0:
                    _axes[0].axvline(x=np.average(replicate_data[j]), color=colors[j], alpha=0.8, ls='dotted', label=r'$\textnormal{average r'+f'{replica_names[j]}'+'}$') 
                else: 
                    _axes[0].axvline(x=np.average(replicate_data[j]), color=colors[j], alpha=0.8, ls='dotted', label=r'$\textnormal{average r'+f'{j+1}'+'}$') 

        # labels
        _axes[0].set(xlabel=r'$\textnormal{'+f'{metric_name}'+'}$', ylabel=r'$\textnormal{density (adim.)}$')
        _axes[0].legend(loc='best', prop={'size': 10}, title=r'$\textsf{'+f'{legend}'+'}$', title_fontsize=11)
        if len(xlims) != 0:
            _axes[0].set_xlim(left=xlims[0], right=xlims[1])

        plt.tight_layout()
        plt.savefig(os.path.join(md_analysis_base_dir, parent_dir, save_name), dpi=400)
        #plt.show(fig)

    def plot_replicates_time(replicate_data, save_name, metric_name=r'', convergence_line_last=None, x_time=True, skip=None):
        '''
        Single plot with a metric wrt time for several replicates. MAX 4 REPLICATES.
        '''

        # all have the same length, traj length
        frames = np.arange(1, len(replicate_data[0])+1, 1)
        
        if x_time:
            # each frame corresponds to 10 ps
            # 10 ns, produces 1000 frames, i.e. 0.01 ns / frame
            times = frames * 0.01 

        # generate time figure
        fig, axes = plt.subplots(ncols=1, nrows=1, figsize=(6, 5))
        _axes = [axes]
        #_axes = axes[0,0], axes[0, 1], axes[1, 0], axes[1, 1], axes[2, 0], axes[2,1], axes[3,0], axes[3,1]

        # color lists
        colors = ['darkgreen', 'deeppink', 'mediumblue', 'orange']

        for i in range(len(_axes)):

            _axes[i].tick_params(bottom=True, top=True, left=True, right=False,direction="in")
            _axes[i].tick_params(labelbottom=True, labeltop=False, labelleft=True,labelright=False, direction="in")

            # for each replicate
            for j in range(len(replicate_data)):
                if x_time:
                    if skip is not None:
                        _axes[i].plot(times[::skip], replicate_data[j, ::skip], color=colors[j], alpha=0.9, label=r'$\textnormal{replicate '+f'{j+1}'+r' }$')
                    else: 
                        _axes[i].plot(times, replicate_data[j], color=colors[j], alpha=0.9, label=r'$\textnormal{replicate '+f'{j+1}'+r' }$')
                else:
                    if skip is not None:
                        _axes[i].plot(frames[::skip], replicate_data[j, ::skip], color=colors[j], alpha=0.9, label=r'$\textnormal{replicate '+f'{j+1}'+r' }$')
                    else: 
                        _axes[i].plot(frames, replicate_data[j], color=colors[j], alpha=0.9, label=r'$\textnormal{replicate '+f'{j+1}'+r' }$')
                if convergence_line_last is not None:
                    _axes[i].axhline(y=np.average(replicate_data[j][-convergence_line_last:]), color=colors[j], alpha=0.8, ls='dotted', label=r'$\textnormal{converged r'+f'{j+1}'+'}$') 
                
            # labels
            if x_time:
                _axes[i].set(xlabel=r'$\textnormal{time (ns)}$', ylabel=r'$\textnormal{'+f'{metric_name}'+'}$')
            else:
                _axes[i].set(xlabel=r'$\textnormal{frame (adim.)}$', ylabel=r'$\textnormal{'+f'{metric_name}'+'}$')
            _axes[i].legend(loc='best', prop={'size': 10})

        plt.tight_layout()
        plt.savefig(os.path.join(md_analysis_base_dir, parent_dir, save_name), dpi=400)
        plt.close(fig)


    def plot_replicates_histogram(replicate_data, save_name, metric_name='', bin_size=0.1, display_mean=True):
        '''
        Single histogram with a metric (statistic) for several replicates. MAX 4 REPLICATES.
        '''

        replicate_data = np.array(replicate_data)
        bins = np.arange(np.min(replicate_data), np.max(replicate_data), bin_size)

        # generate time figure
        fig, axes = plt.subplots(ncols=1, nrows=1, figsize=(6, 5))
        _axes = [axes]
        #_axes = axes[0,0], axes[0, 1], axes[1, 0], axes[1, 1], axes[2, 0], axes[2,1], axes[3,0], axes[3,1]

        # color lists
        colors = ['darkgreen', 'deeppink', 'mediumblue', 'orange']

        for i in range(len(_axes)):

            _axes[i].tick_params(bottom=True, top=True, left=True, right=False,direction="in")
            _axes[i].tick_params(labelbottom=True, labeltop=False, labelleft=True,labelright=False, direction="in")

            # for each replicate
            for j in range(len(replicate_data)):
            
                # filling
                _axes[i].hist(replicate_data[j], bins=bins, density=True, edgecolor=colors[j], color='w', histtype='stepfilled')

            # for each replicate
            for j in range(len(replicate_data)):

                # contour
                _axes[i].hist(replicate_data[j], bins=bins, density=True, color=colors[j], alpha=.08, label=r'$\textnormal{replicate '+f'{j+1}'+r' }$')    
                if display_mean:
                    _axes[i].axvline(x=np.average(replicate_data[j]), color=colors[j], alpha=0.8, ls='dotted', label=r'$\textnormal{average r'+f'{j+1}'+r'}$') 
                

            # labels
            _axes[i].set(xlabel=r'$\textnormal{'+f'{metric_name}'+'}$', ylabel=r'$\textnormal{density (adim.)}$')
            _axes[i].legend(loc='best', prop={'size': 10})

        plt.tight_layout()
        plt.savefig(os.path.join(md_analysis_base_dir, parent_dir, save_name), dpi=400)
        plt.close(fig)


    def plot_replicates_both(replicate_data, save_name, metric_name=r'', convergence_line_last=5, x_time=True, bin_size=0.1, display_mean=True, skip=None, replica_names=[]):
        '''
        Single plot with a metric wrt time for several replicates. MAX 4 REPLICATES.
        '''
        '''
        Single histogram with a metric (statistic) for several replicates. MAX 4 REPLICATES.
        '''


        replicate_data = np.array(replicate_data)

        # all have the same length, traj length
        frames = np.arange(1, len(replicate_data[0])+1, 1)
        
        if x_time:
            # each frame corresponds to 10 ps
            # 10 ns, produces 1000 frames, i.e. 0.01 ns / frame
            times = frames * 0.01 

        bins = np.arange(np.min(replicate_data), np.max(replicate_data), bin_size)

        # generate time figure
        fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(7, 3.5))
        _axes = axes

        # color lists
        colors = ['darkgreen', 'deeppink', 'mediumblue', 'orange']

        
        # wrt time
        _axes[0].tick_params(bottom=True, top=True, left=True, right=False,direction="in")
        _axes[0].tick_params(labelbottom=True, labeltop=False, labelleft=True,labelright=False, direction="in")

        # for each replicate
        for j in range(len(replicate_data)):
            if x_time:
                if skip is not None:
                    if len(replica_names)!= 0:
                        _axes[0].plot(times[::skip], replicate_data[j, ::skip], color=colors[j], alpha=0.9, label=r'$\textnormal{'+f'{replica_names[j]}'+r'}$')
                    else:
                        _axes[0].plot(times[::skip], replicate_data[j, ::skip], color=colors[j], alpha=0.9, label=r'$\textnormal{replicate '+f'{j+1}'+r'}$')
                else:
                    if len(replica_names)!= 0:
                        _axes[0].plot(times, replicate_data[j], color=colors[j], alpha=0.9, label=r'$\textnormal{'+f'{replica_names[j]}'+r'}$')  
                    else:
                        _axes[0].plot(times, replicate_data[j], color=colors[j], alpha=0.9, label=r'$\textnormal{replicate '+f'{j+1}'+r' }$')
            else:
                if skip is not None:
                    _axes[0].plot(frames[::skip], replicate_data[j, ::skip], color=colors[j], alpha=0.9, label=r'$\textnormal{replicate '+f'{j+1}'+r' }$')
                else: 
                    _axes[0].plot(frames, replicate_data[j], color=colors[j], alpha=0.9, label=r'$\textnormal{replicate '+f'{j+1}'+r' }$')
            
            if convergence_line_last is not None:
                _axes[0].axhline(y=np.average(replicate_data[j][-convergence_line_last:]), color=colors[j], alpha=0.8, ls='dotted', label=r'$\textnormal{converged r'+f'{j+1}'+'}$') 
            
        if x_time:
            _axes[0].set(xlabel=r'$\textnormal{time (ns)}$', ylabel=r'$\textnormal{'+f'{metric_name}'+'}$')
            #_axes[0].set_ylim(0, 10)
        else:
            _axes[0].set(xlabel=r'$\textnormal{frame (adim.)}$', ylabel=r'$\textnormal{'+f'{metric_name}'+'}$')
        _axes[1].legend(loc='best', prop={'size': 10}, title=r'$\textsf{'+f'Simulations'+'}$', title_fontsize=11)



        # histogram
        _axes[1].tick_params(bottom=True, top=True, left=True, right=False,direction="in")
        _axes[1].tick_params(labelbottom=True, labeltop=False, labelleft=True,labelright=False, direction="in")

        # for each replicate
        for j in range(len(replicate_data)):
            print(len(replicate_data))
            print(f'algo {j} histogram explode')
            print(replicate_data[j])
        
            # filling
            _axes[1].hist(replicate_data[j], bins=bins, density=True, edgecolor=colors[j], color='w', histtype='stepfilled')

        # for each replicate
        for j in range(len(replicate_data)):

            # contour
            if len(replica_names) != 0:
                _axes[1].hist(replicate_data[j], bins=bins, density=True, color=colors[j], alpha=.08, label=r'$\textnormal{'+f'{replica_names[j]}'+r' }$')  
            else:
                _axes[1].hist(replicate_data[j], bins=bins, density=True, color=colors[j], alpha=.08, label=r'$\textnormal{replicate '+f'{j+1}'+r' }$')    
            if display_mean:
                if len(replica_names) != 0:
                    _axes[1].axvline(x=np.average(replicate_data[j]), color=colors[j], alpha=0.8, ls='dotted', label=r'$\textnormal{average r'+f'{replica_names[j]}'+'}$') 
                else: 
                    _axes[1].axvline(x=np.average(replicate_data[j]), color=colors[j], alpha=0.8, ls='dotted', label=r'$\textnormal{average r'+f'{j+1}'+'}$') 

        # labels
        _axes[1].set(xlabel=r'$\textnormal{'+f'{metric_name}'+'}$', ylabel=r'$\textnormal{density (adim.)}$')
        _axes[1].legend(loc='best', prop={'size': 10}, title=r'$\textsf{'+f'Simulation'+'}$', title_fontsize=11)

    
        plt.tight_layout()
        plt.savefig(os.path.join(md_analysis_base_dir, parent_dir, save_name), dpi=400)

        
    def plot_replicates_time_separate(replicate_data, save_name, metric_name=r'', convergence_line_last=5, x_time=True):
        '''
        Single plot with a metric wrt time for several replicates. MAX 4 REPLICATES.
        '''

        # all have the same length, traj length
        frames = np.arange(1, len(replicate_data[0])+1, 1)
        if x_time:
            # each frame corresponds to 10 ps
            # 10 ns, produces 1000 frames, i.e. 0.01 ns / frame
            times = frames * 0.01 


        # generate time figure
        fig, axes = plt.subplots(ncols=len(replicate_data), nrows=1, figsize=(4.5*len(replicate_data), 4))
        _axes = axes

        # color lists
        colors = ['darkgreen', 'deeppink', 'mediumblue', 'orange']

        for i in range(len(_axes)):

            _axes[i].tick_params(bottom=True, top=True, left=True, right=False,direction="in")
            _axes[i].tick_params(labelbottom=True, labeltop=False, labelleft=True,labelright=False, direction="in")

            if x_time:
                _axes[i].plot(times, replicate_data[i], color=colors[i], alpha=0.9, label=r'$\textnormal{replicate '+f'{i+1}'+r' }$')
            else:
                _axes[i].plot(frames, replicate_data[i], color=colors[i], alpha=0.9, label=r'$\textnormal{replicate '+f'{i+1}'+r' }$')
            
            if convergence_line_last is not None:
                _axes[i].axhline(y=np.average(replicate_data[i][-convergence_line_last:]), color=colors[i], alpha=0.8, ls='dotted', label=r'$\textnormal{converged r'+f'{i+1}'+'}$') 
            
            # labels
            if i==0:
                _axes[i].set(ylabel=r'$\textnormal{'+f'{metric_name}'+'}$')
            if x_time:
                _axes[i].set(xlabel=r'$\textnormal{time (ns)}$')
            else:
                _axes[i].set(xlabel=r'$\textnormal{frame (adim.)}$')
            _axes[i].legend(loc='best', prop={'size': 10})

        plt.tight_layout()
        plt.savefig(os.path.join(md_analysis_base_dir, parent_dir, save_name), dpi=400)
        plt.close(fig)


    def plot_replicates_histogram_separate(replicate_data, save_name, metric_name=r'', bin_size=0.2, display_mean=True):
        '''
        Single plot with a metric wrt time for several replicates. MAX 4 REPLICATES.
        '''

        replicate_data = np.array(replicate_data)
        bins = np.arange(np.min(replicate_data), np.max(replicate_data), bin_size)

        # generate time figure
        # generate time figure
        fig, axes = plt.subplots(ncols=len(replicate_data), nrows=1, figsize=(4.5*len(replicate_data), 4))
        _axes = axes

        # color lists
        colors = ['darkgreen', 'deeppink', 'mediumblue', 'orange']

        for i in range(len(_axes)):

            _axes[i].tick_params(bottom=True, top=True, left=True, right=False,direction="in")
            _axes[i].tick_params(labelbottom=True, labeltop=False, labelleft=True,labelright=False, direction="in")
    
            _axes[i].hist(replicate_data[i], bins=bins, density=True, edgecolor=colors[i], color='w', histtype='stepfilled')
            _axes[i].hist(replicate_data[i], bins=bins, density=True, color=colors[i], alpha=.08, label=r'$\textnormal{replicate '+f'{i+1}'+r' }$')    
            if display_mean:
                _axes[i].axvline(x=np.average(replicate_data[i]), color=colors[i], alpha=0.8, ls='dotted', label=r'$\textnormal{average r'+f'{i+1}'+'}$') 
            
            # labels
            if i==0:
                _axes[i].set(ylabel=r'$\textnormal{density (adim.)}$')
            _axes[i].set(xlabel=r'$\textnormal{'+f'{metric_name}'+'}$')
            _axes[i].legend(loc='best', prop={'size': 10})

        plt.tight_layout()
        plt.savefig(os.path.join(md_analysis_base_dir, parent_dir, save_name), dpi=400)
        plt.close(fig)




    def plot_contacts_histogram(dict_list, residue_names):
        '''
        plot bar plot for each label, with replicas of the contacts of a given group of atoms 
        '''

        # for each contact
        for i  in range(len(dict_list)):
            
            # define the figure (subfigures = num of replicas)

            #print(len(dict_list))
            fig, axes = plt.subplots(ncols=len(dict_list[i]), nrows=1, figsize=(4.5*len(dict_list[i]), 4))
            if len(dict_list[i]) == 1:
                _axes = [axes]
            else:
                _axes = axes

            # for each replica
            for j in range(len(dict_list[i])):
                data_dict = dict_list[i][j]
                filtered_dict = data_dict.copy()

                for key in data_dict.keys():
                    if key not in residue_names[i]:
                        #print(f'{key} not in {residue_names[i]}')
                        filtered_dict.pop(key, None)

                for res in residue_names[i]:
                    if res not in data_dict.keys():
                        filtered_dict[res] = 0.0

                residue_list = [r'$\textnormal{'+f'{x}'+r'}$' for x in list(filtered_dict.keys())]
                value_list =  list(filtered_dict.values())
                #print(residue_list)
                #print(value_list)


                _axes[j].tick_params(bottom=True, top=True, left=True, right=False,direction="in")
                _axes[j].tick_params(labelbottom=True, labeltop=False, labelleft=True,labelright=False, direction="in")
                _axes[j].tick_params(axis='x', labelrotation=60) 
                
                _axes[j].bar(residue_list, value_list, width=4.5/(len(residue_names[i])+ 0.1*len(residue_names[i])), facecolor='darkgreen', alpha=0.1)
                _axes[j].bar(residue_list, value_list, width=4.5/(len(residue_names[i])+ 0.1*len(residue_names[i])), linewidth=1, edgecolor='darkgreen', fill=False)

                # labels
                if j==0:
                    _axes[j].set(ylabel=r'$\textnormal{'+f'contact frequency (percent)'+'}$')
                _axes[j].set(xlabel=r'$\textnormal{residue ID}$')

                #print(data_dict)

        plt.tight_layout()
        plt.savefig(os.path.join(md_analysis_base_dir, parent_dir, f'contacts_hist {label}'), dpi=400)





    
    def plot_contacts_evolution(contact_dict_list, residue_names, traj_list):
        '''
        PLOT SELECTED CONTACTS (A PRIORI, RESXXX) EVOLUTION. ONE REPLICA SUPPORTED!

        contact_dict_list is a list of contact dictionaries, with the frequency of contact for each mini-traj file. 


        dictionary containing

        '''

        # small traj files, for each of these plot the contacts
        x = np.arange(1, len(traj_list)+1) * 10 

        fig, axes = plt.subplots(ncols=len(contact_dict_list), nrows=1, figsize=(4*len(contact_dict_list), 4))
        if len(contact_dict_list) == 1:
            _axes = [axes]
        else:
            _axes = axes


        # color lists
        colors = ['darkviolet', 'crimson', 'deeppink', 'darkgreen', 'mediumseagreen', 'royalblue', 'darkblue']

        # marker list
        marker = ['o', '^', 'd', 's', 'h',  'x', 'o', 'o', 'o']

        k=0
        for key, value in contact_dict_list.items():
            
            data = []

            for res_name in residue_names[k]:
                temp = [res_name]
                for i in range(len(value)):
                    try:
                        temp.append(value[i][res_name])
                    except: 
                        temp.append(0.0)


                data.append(temp)
            #print('data', data)

            _axes[k].tick_params(bottom=True, top=True, left=True, right=False,direction="in")
            _axes[k].tick_params(labelbottom=True, labeltop=False, labelleft=True,labelright=False, direction="in")

            for i in range(len(data)):
                _label = data[i][0]
                #print('printing label', _label)
                _values = data[i][1:]   
                #print('printing values', _values)
                
                _axes[k].plot(x, _values, label=r'$\textnormal{'+f'{_label}'+'}$', marker=marker[i], color=colors[i], alpha=0.7)
                
            # labels
            _axes[k].set(xlabel=r'$\textnormal{time (ns)}$', ylabel=r'$\textnormal{'+f'{key} contact frequency (percent)'+'}$')
            _axes[k].legend(loc='best', prop={'size': 10}, title=r'$\textsf{'+f'Contacts'+'}$', title_fontsize=11)

            k+=1
        
        plt.tight_layout()
        plt.savefig(os.path.join(md_analysis_base_dir, parent_dir, 'contacts.png'), dpi=400)




    lig_ID_list = ['AA', '5-HpETE', '5-HETE']
    if lig_ID not in lig_ID_list:
        raise ValueError(f'Introduced lig_ID {lig_ID} is not among the supported ligands: {lig_ID_list}')
    else:
        if lig_ID == 'AA':
            lig_name = 'AA'

        elif lig_ID == '5-HpETE':
            lig_name = '5-HpETE'

        elif lig_ID == '5-HETE':
            lig_name = '5-HETE'

        import re
        md_number = re.findall(r'\d+', parent_dir)[-2]
        replica_names = [f'{lig_name} {md_number}']


    md_analysis_base_dir = '/Volumes/white_HDD/md_analysis/'
    if not os.path.exists(os.path.join(md_analysis_base_dir, parent_dir)):
        os.makedirs(os.path.join(md_analysis_base_dir, parent_dir))

    # md source directory (all data, preprod, prod, prmtop, pdb, incrd)
    #src_dir_1 = os.listdir(os.path.join(md_base_dir, parent_dir))

    # topology
    topologies_1 = [i for i in os.listdir(os.path.join(md_base_dir, parent_dir)) if i.endswith('.prmtop')]
    top = os.path.join(md_base_dir, parent_dir, topologies_1[0])
    if len(topologies_1) > 1:
        raise ValueError('bruh. More than one topology files? Check md source directory.')
    print(f'loading prmtop    \t{top}')

    # production trajectory
    prod_dir_1 = os.path.join(md_base_dir, parent_dir, 'prod/out')                    # locate the production traj
    prod_files = [int(x.strip('.nc').strip('prod_')) for x in os.listdir(prod_dir_1) if '.nc' in x and ('.' not in x[0])]
    traj_1 = [os.path.join(prod_dir_1, f'prod_{i}.nc')  for i in prod_files]      # assume constant length 250 ns, 25 
    for i in traj_1:
        print(f'loading trajectory\t{i}')

    # load base MD (r=1)
    emda = EMDA()
    emda.load_variant(
        parameters=top,
        trajectory=traj_1,
        variant_name=protein_ID
    )

    total_time_simulated = 10 * len(traj_1)
    print(f'total time simulated: {total_time_simulated} ns')



    print('performing atom selections')
    if lig_ID == lig_ID_list[0]:

        #===============#
        #       AA      #
        #===============#

        print(f'selecting atoms for {lig_ID_list[0]}')

        # LIGAND
        emda.select('LIG', '665', sel_type='res_num')

        emda.select('C10', 'name C10 and resid 665')
        emda.select('H10S', 'name H10A and resid 665')
        emda.select('H10R', 'name H10B and resid 665')

        emda.select('C13', 'name C13 and resid 665')
        emda.select('H13S', 'name H13A and resid 665')
        emda.select('H13R', 'name H13B and resid 665')

        emda.select('Met', 'name C20 and resid 665')

        # planarity dihedrals
        # C10 planarity
        emda.select('p10_A', 'name C12 and resid 665')      # A i B formen el cis double bond
        emda.select('p10_B', 'name C11 and resid 665')
        emda.select('p10_C', 'name C10 and resid 665')      # C conté H10A i H10B
        emda.select('p10_D', 'name C9 and resid 665')       # D és el carboni més proxim al COO

        # C13 planarity
        emda.select('p13_A', 'name C11 and resid 665')      # A i B formen el cis double bond
        emda.select('p13_B', 'name C12 and resid 665')
        emda.select('p13_C', 'name C13 and resid 665')      # C conté H13A i H13B
        emda.select('p13_D', 'name C14 and resid 665')      # D és el carboni més llunyà al COO

        # functional groups
        #                   C (C1) O (O1) O (O2)
        emda.select('COO', [10554, 10555, 10556], sel_type='at_num')
        emda.select('COO_C', 10554, sel_type='at_num')    
        emda.select('COO_O1', 10555, sel_type='at_num')    
        emda.select('COO_O2', 10556, sel_type='at_num')  

        # PROTEIN hLOX15
        # OH1 OH cofactor :664
        #                   O
        emda.select('OH1', 10552, sel_type='at_num') 

        # special residues (sidechain analysis)
        emda.select('ILE417_CD1', 6708, sel_type='at_num')      # cavity depth ILE417

        emda.select('GLU175_OE1', 2827, sel_type='at_num')      # GLU175 O -- ARG402 N, H bond tracking.

        emda.select('ARG402_NE', 6479, sel_type='at_num')       # ARG402 N epsilon -- COO, H bond tracking.
        emda.select('ARG402_NH1', 6482, sel_type='at_num')      # ARG402 NH -- COO and GLU175, H bond tracking.
        emda.select('ARG402_NH2', 6485, sel_type='at_num')      # ARG402 NH -- COO and GLU175, H bond tracking.    

        # backbone
        emda.select('backbone', 'backbone')

        # domains
        #emda.select('PLAT', 'resid 1-110')


        print(emda.selections)


    elif lig_ID == lig_ID_list[1]:

        #====================#
        #       5-HpETE      #
        #====================#
        
        print(f'selecting atoms for {lig_ID_list[1]}')

        # LIGAND
        emda.select('LIG', '665', sel_type='res_num')

        emda.select('C10', 'name C11 and resid 665')
        emda.select('H10S', 'name H19 and resid 665')
        emda.select('H10R', 'name H18 and resid 665')

        emda.select('C13', 'name C8 and resid 665')
        emda.select('H13S', 'name H15 and resid 665')
        emda.select('H13R', 'name H14 and resid 665')

        emda.select('Met', 'name C1 and resid 665')

        # planarity dihedrals
        # C10 planarity
        emda.select('p10_A', 'name C9 and resid 665')       # A i B formen el cis double bond
        emda.select('p10_B', 'name C10 and resid 665')    
        emda.select('p10_C', 'name C11 and resid 665')      # C conté H10A i H10B
        emda.select('p10_D', 'name C12 and resid 665')      # D és el carboni més proxim al COO

        # C13 planarity
        emda.select('p13_A', 'name C10 and resid 665')      # A i B formen el cis double bond
        emda.select('p13_B', 'name C9 and resid 665')
        emda.select('p13_C', 'name C8 and resid 665')       # C conté H13A i H13B
        emda.select('p13_D', 'name C7 and resid 665')       # D és el carboni més llunyà al COO


        # functional groups 
        #                  C (C20) O (O1) O (O2)
        emda.select('COO', [10604, 10605, 10606], sel_type='at_num')
        emda.select('COO_C', 10604, sel_type='at_num')    
        emda.select('COO_O1', 10605, sel_type='at_num')    
        emda.select('COO_O2', 10606, sel_type='at_num')    

        #                   O (O3) O (O4) H (H23)
        emda.select('OOH', [10607, 10608, 10590], sel_type='at_num')


        # PROTEIN hLOX15
        # OH1 OH cofactor :664
        #                   O
        emda.select('OH1', 10552, sel_type='at_num') 

        # special residues (sidechain analysis)
        emda.select('ILE417_CD1', 6708, sel_type='at_num')      # cavity depth ILE417

        emda.select('GLU175_OE1', 2827, sel_type='at_num')      # GLU175 O -- ARG402 N, H bond tracking.

        emda.select('ARG402_NE', 6479, sel_type='at_num')       # ARG402 N epsilon -- COO, H bond tracking.
        emda.select('ARG402_NH1', 6482, sel_type='at_num')      # ARG402 NH -- COO and GLU175, H bond tracking.
        emda.select('ARG402_NH2', 6485, sel_type='at_num')      # ARG402 NH -- COO and GLU175, H bond tracking.    


        # backbone
        emda.select('backbone', 'backbone')

        # domains
        #emda.select('PLAT', 'resid 1-110')


        print(emda.selections)


    elif lig_ID == lig_ID_list[2]:

        #===================#
        #       5-HETE      #
        #===================#

        print(f'selecting atoms for {lig_ID_list[2]}')

        # ligand atoms
        emda.select('LIG', '665', sel_type='res_num')

        emda.select('C10', 'name C11 and resid 665')
        emda.select('H10S', 'name H19 and resid 665')
        emda.select('H10R', 'name H18 and resid 665')

        emda.select('C13', 'name C8 and resid 665')
        emda.select('H13S', 'name H15 and resid 665')
        emda.select('H13R', 'name H14 and resid 665')

        emda.select('Met', 'name C1 and resid 665')

        # planarity dihedrals
        # C10 planarity
        emda.select('p10_A', 'name C9 and resid 665')       # A i B formen el cis double bond
        emda.select('p10_B', 'name C10 and resid 665')    
        emda.select('p10_C', 'name C11 and resid 665')      # C conté H10A i H10B
        emda.select('p10_D', 'name C12 and resid 665')      # D és el carboni més proxim al COO

        # C13 planarity
        emda.select('p13_A', 'name C10 and resid 665')      # A i B formen el cis double bond
        emda.select('p13_B', 'name C9 and resid 665')
        emda.select('p13_C', 'name C8 and resid 665')       # C conté H13A i H13B
        emda.select('p13_D', 'name C7 and resid 665')       # D és el carboni més llunyà al COO


        # functional groups 
        #                  C (C20) O (O1) O (O2)
        emda.select('COO', [10603, 10604, 10605], sel_type='at_num')
        emda.select('COO_C', 10603, sel_type='at_num')    
        emda.select('COO_O1', 10604, sel_type='at_num')    
        emda.select('COO_O2', 10605, sel_type='at_num')    

        #                  O (O3) H (H31)
        emda.select('OH', [10606, 10607], sel_type='at_num')

        # PROTEIN hLOX15
        # OH1 OH cofactor :664
        #                   O
        emda.select('OH1', 10552, sel_type='at_num') 

        # special residues (sidechain analysis)
        emda.select('ILE417_CD1', 6708, sel_type='at_num')      # cavity depth ILE417
        
        emda.select('GLU175_OE1', 2827, sel_type='at_num')      # GLU175 O -- ARG402 N, H bond tracking.

        emda.select('ARG402_NE', 6479, sel_type='at_num')       # ARG402 N epsilon -- COO, H bond tracking.
        emda.select('ARG402_NH1', 6482, sel_type='at_num')      # ARG402 NH -- COO and GLU175, H bond tracking.
        emda.select('ARG402_NH2', 6485, sel_type='at_num')      # ARG402 NH -- COO and GLU175, H bond tracking.    

        # backbone
        emda.select('backbone', 'backbone')

        # domains
        #emda.select('PLAT', 'resid 1-110')


        #print(emda.selections)


    else:
        raise ValueError('Something went wrong!')


    # hydrogen abstraction -- OH1 cofactor distances
    emda.add_distance('distance H10S -- O OH1', 'OH1', 'H10S')
    emda.add_distance('distance H10R -- O OH1', 'OH1', 'H10R')
    emda.add_distance('distance C10 -- O OH1', 'OH1', 'C10')

    emda.add_distance('distance H13S -- O OH1', 'OH1', 'H13S')
    emda.add_distance('distance H13R -- O OH1', 'OH1', 'H13R')
    emda.add_distance('distance C13 -- O OH1', 'OH1', 'C13')

    distance_H_set = ['distance H10S -- O OH1', 'distance H10R -- O OH1', 'distance C10 -- O OH1', 'distance H13S -- O OH1', 'distance H13R -- O OH1', 'distance C13 -- O OH1']

    # ARG402 H bonds
    emda.add_distance('distance GLU175 OE1 -- ARG402 NH1', 'GLU175_OE1', 'ARG402_NH1')      # GLU175 -- ARG402 (NH1 and NH2)
    emda.add_distance('distance GLU175 OE1 -- ARG402 NH2', 'GLU175_OE1', 'ARG402_NH2')     # GLU175 -- ARG402

    distance_GLU_ARG = ['distance GLU175 OE1 -- ARG402 NH1', 'distance GLU175 OE1 -- ARG402 NH2'] 

    emda.add_distance('distance COO O1 -- ARG402 NE', 'COO_O1', 'ARG402_NE')     # COO (both oxygens) -- ARG402 (NE, NH1, NH2)
    emda.add_distance('distance COO O1 -- ARG402 NH1', 'COO_O1', 'ARG402_NH1')
    emda.add_distance('distance COO O1 -- ARG402 NH2', 'COO_O1', 'ARG402_NH2')
    emda.add_distance('distance COO O2 -- ARG402 NE', 'COO_O2', 'ARG402_NE')
    emda.add_distance('distance COO O2 -- ARG402 NH1', 'COO_O2', 'ARG402_NH1')
    emda.add_distance('distance COO O2 -- ARG402 NH2', 'COO_O2', 'ARG402_NH2')

    distance_COO_set = ['distance COO O1 -- ARG402 NE', 'distance COO O1 -- ARG402 NH1', 'distance COO O1 -- ARG402 NH2', 'distance COO O2 -- ARG402 NE', 'distance COO O2 -- ARG402 NH1', 'distance COO O2 -- ARG402 NH2']

    # cavity depth 
    emda.add_distance('distance Methyl -- ILE417 CD1', 'Met', 'ILE417_CD1')

    distance_Met_set = ['distance Methyl -- ILE417 CD1']

    # dihedrals
    emda.add_dihedral('dihedral C12--C11--C10--C9 (C10 planarity)', 'p10_A', 'p10_B', 'p10_C', 'p10_D', domain=360)     # around C10
    emda.add_dihedral('dihedral C11--C12--C13--C14 (C13 planarity)', 'p13_A', 'p13_B', 'p13_C', 'p13_D', domain=360)    # around C13

    dihedral_set = ['dihedral C12--C11--C10--C9 (C10 planarity)', 'dihedral C11--C12--C13--C14 (C13 planarity)']

    # rmsd
    emda.add_RMSD(f'{lig_ID} RMSD', 'LIG')
    #emda.add_RMSD('RMSD-alpha', 'alpha')    # TODO
    emda.add_RMSD('backbone RMSD', 'backbone')

    rmsd_set = [f'{lig_ID} RMSD', 'backbone RMSD']

    # contacts
    emda.add_contacts('COO contacts', 'COO', sel_env=4, include_WAT=True)
    if lig_ID == lig_ID_list[0]:
        contacts_set = ['COO contacts']
    elif lig_ID == lig_ID_list[1]:
        emda.add_contacts('OOH contacts', 'OOH', sel_env=4, include_WAT=True)
        contacts_set = ['COO contacts', 'OOH contacts']
    elif lig_ID == lig_ID_list[2]:
        emda.add_contacts('OH contacts', 'OH', sel_env=4, include_WAT=True)
        contacts_set = ['COO contacts', 'OH contacts']

    print('selected measurements:\n')
    print(emda.measures)


    if not load_results:
        # compute the metrics
        emda.run()

        # save the analysis
        emda.save(os.path.join(md_analysis_base_dir, parent_dir, f'{protein_ID}_{lig_ID}_analysis'))    

    else: 
        # load precomputed analysis 
        emda.load(os.path.join(md_analysis_base_dir, parent_dir, f'{protein_ID}_{lig_ID}_analysis.pkl'))

    print(f'loaded systems:\t{emda.universe}')
    print(f'\ntotal measures conducted:\n{emda.measures}')


    print(f"\nprotein \t{protein_ID}\n")
    for label in distance_H_set:

        for variant, variant_data in emda.measures[label].result.items():
            

            print(f'\n\033[1mMETRIC: \t{label}\n')
            replicate_data = list(variant_data.values())

            for replica, values in variant_data.items():
                if values:  
                    # if values are present
                    if label == 'distance H10S -- O OH1':
                        distance_H10S_value = values
                    if label == 'distance H13S -- O OH1':
                        distance_H13S_value = values

                    print(f"\033[0mreplica \t{replica}")
                    print(f"avg     \t{average(values):.3f}")
                    print(f"std     \t{std(values):.3f}")
                    print(f"max     \t{max(values):.3f}")
                    print(f"min     \t{min(values):.3f}")

                else:
                    print(f"! replica {replica} has no values !")
                    raise ValueError('Some replicas dont have any values! Abort.')
                
            
            # plot wrt time
            # display all available replicas 
            metric_name = f'{label} (Å)'
            #plot_replicates_time(replicate_data, save_name=f'{label}_plot.png', metric_name=metric_name, convergence_line_last=100)
            #plot_replicates_histogram(replicate_data, save_name=f'{label}_hist.png', metric_name=metric_name, display_mean=True)
            plot_replicates_both(replicate_data, save_name=f'{label}_both.png', metric_name=metric_name, display_mean=False, skip=None, replica_names=replica_names, convergence_line_last=None)


            print(f"\nprotein \t{protein_ID}\n")

    plot_two_metrics_hist(
            [distance_H10S_value, distance_H13S_value], 
            save_name='H10_H13_hist.png', 
            metric_name=r'H$_x$ distance (Å)', x_time=True, bin_size=0.1, skip=None, replica_names=[r'\textnormal{H$_{10S}$}', r'\textnormal{H$_{13S}$}'], display_mean=False, legend='Hydrogen', xlims=[2,7])
    for label in distance_GLU_ARG:

        for variant, variant_data in emda.measures[label].result.items():
            

            print(f'\n\033[1mMETRIC: \t{label}\n')
            replicate_data = list(variant_data.values())

            for replica, values in variant_data.items():
                if values:  
                    # if values are present

                    print(f"\033[0mreplica \t{replica}")
                    print(f"avg     \t{average(values):.3f}")
                    print(f"std     \t{std(values):.3f}")
                    print(f"max     \t{max(values):.3f}")
                    print(f"min     \t{min(values):.3f}")

                else:
                    print(f"! replica {replica} has no values !")
                    raise ValueError('Some replicas dont have any values! Abort.')
                
            
            # plot wrt time
            # display all available replicas 
            metric_name = f'{label} (Å)'
            #plot_replicates_time(replicate_data, save_name=f'{label}_plot.png', metric_name=metric_name, convergence_line_last=100)
            #plot_replicates_histogram(replicate_data, save_name=f'{label}_hist.png', metric_name=metric_name, display_mean=True)
            plot_replicates_both(replicate_data, save_name=f'{label}_both.png', metric_name=metric_name, display_mean=False, skip=None, replica_names=replica_names, convergence_line_last=None)



            print(f"\nprotein \t{protein_ID}\n")
    for label in distance_COO_set:

        for variant, variant_data in emda.measures[label].result.items():
            

            print(f'\n\033[1mMETRIC: \t{label}\n')
            replicate_data = list(variant_data.values())

            for replica, values in variant_data.items():
                if values:  
                    # if values are present

                    print(f"\033[0mreplica \t{replica}")
                    print(f"avg     \t{average(values):.3f}")
                    print(f"std     \t{std(values):.3f}")
                    print(f"max     \t{max(values):.3f}")
                    print(f"min     \t{min(values):.3f}")

                else:
                    print(f"! replica {replica} has no values !")
                    raise ValueError('Some replicas dont have any values! Abort.')
                
            
            # plot wrt time
            # display all available replicas 
            metric_name = f'{label} (Å)'
            #plot_replicates_time(replicate_data, save_name=f'{label}_plot.png', metric_name=metric_name, convergence_line_last=100)
            #plot_replicates_histogram(replicate_data, save_name=f'{label}_hist.png', metric_name=metric_name, display_mean=True)
            plot_replicates_both(replicate_data, save_name=f'{label}_both.png', metric_name=metric_name, display_mean=False, skip=None, replica_names=replica_names, convergence_line_last=None)




            print(f"\nprotein \t{protein_ID}\n")
    for label in distance_Met_set:

        for variant, variant_data in emda.measures[label].result.items():
            

            print(f'\n\033[1mMETRIC: \t{label}\n')
            replicate_data = list(variant_data.values())

            for replica, values in variant_data.items():
                if values:  
                    # if values are present

                    print(f"\033[0mreplica \t{replica}")
                    print(f"avg     \t{average(values):.3f}")
                    print(f"std     \t{std(values):.3f}")
                    print(f"max     \t{max(values):.3f}")
                    print(f"min     \t{min(values):.3f}")

                else:
                    print(f"! replica {replica} has no values !")
                    raise ValueError('Some replicas dont have any values! Abort.')
                
            
            # plot wrt time
            # display all available replicas 
            metric_name = f'{label} (Å)'
            #plot_replicates_time(replicate_data, save_name=f'{label}_plot.png', metric_name=metric_name, convergence_line_last=100)
            #plot_replicates_histogram(replicate_data, save_name=f'{label}_hist.png', metric_name=metric_name, display_mean=True)
            plot_replicates_both(replicate_data, save_name=f'{label}_both.png', metric_name=metric_name, display_mean=False, skip=None, replica_names=replica_names, convergence_line_last=None)



            print(f"\nprotein \t{protein_ID}\n")
    for label in dihedral_set:

        for variant, variant_data in emda.measures[label].result.items():
            

            print(f'\n\033[1mMETRIC: \t{label}\n')
            replicate_data = list(variant_data.values())

            for replica, values in variant_data.items():
                if values:  
                    # if values are present
                    if label == 'dihedral C12--C11--C10--C9 (C10 planarity)':
                        dih_C10_value = values
                    if label == 'dihedral C11--C12--C13--C14 (C13 planarity)':
                        dih_C13_value = values

                    print(f"\033[0mreplica \t{replica}")
                    print(f"avg     \t{average(values):.3f}")
                    print(f"std     \t{std(values):.3f}")
                    print(f"max     \t{max(values):.3f}")
                    print(f"min     \t{min(values):.3f}")

                else:
                    print(f"! replica {replica} has no values !")
                    raise ValueError('Some replicas dont have any values! Abort.')
                
            
            # plot wrt time
            # display all available replicas 
            metric_name = f'{label} (Å)'
            #plot_replicates_time(replicate_data, save_name=f'{label}_plot.png', metric_name=metric_name, convergence_line_last=100)
            #plot_replicates_histogram(replicate_data, save_name=f'{label}_hist.png', metric_name=metric_name, display_mean=True)
            plot_replicates_both(replicate_data, save_name=f'{label}_both.png', metric_name=metric_name, display_mean=False, skip=None, replica_names=replica_names, convergence_line_last=None, bin_size=5)


            print(f"\nprotein \t{protein_ID}\n")
        
    plot_two_metrics_hist(
        [dih_C10_value, dih_C13_value], 
        save_name='C10_C13_DIH_hist.png', 
        metric_name=r'C$_{x-1}$--C$_{x+2}$ dihedral ($^\circ$, deg)', x_time=True, bin_size=5, skip=None, replica_names=[r'\textnormal{C$_{10}$}', r'\textnormal{C$_{13}$}'], display_mean=False, legend='Planarity', xlims=[0, 360])


    for label in rmsd_set:

        for variant, variant_data in emda.measures[label].result.items():
            

            print(f'\n\033[1mMETRIC: \t{label}\n')
            replicate_data = list(variant_data.values())

            for replica, values in variant_data.items():
                if values:  
                    # if values are present

                    print(f"\033[0mreplica \t{replica}")
                    print(f"avg     \t{average(values):.3f}")
                    print(f"std     \t{std(values):.3f}")
                    print(f"max     \t{max(values):.3f}")
                    print(f"min     \t{min(values):.3f}")

                else:
                    print(f"! replica {replica} has no values !")
                    raise ValueError('Some replicas dont have any values! Abort.')
                
            
            # plot wrt time
            # display all available replicas 
            metric_name = f'{label} (Å)'
            #plot_replicates_time(replicate_data, save_name=f'{label}_plot.png', metric_name=metric_name, convergence_line_last=100)
            #plot_replicates_histogram(replicate_data, save_name=f'{label}_hist.png', metric_name=metric_name, display_mean=True)
            plot_replicates_both(replicate_data, save_name=f'{label}_both.png', metric_name=metric_name, display_mean=False, skip=None, replica_names=replica_names, convergence_line_last=None)



    # 2 tabs
    contacts_emda = EMDA()

    # for each complete tajectory
    for t in range(len(traj_1)):

        contacts_emda.load_variant(
            parameters=top,
            trajectory=traj_1[t],
            variant_name=f'partial_prod_{t+1}'
        )



    # contacts

    if lig_ID == lig_ID_list[0]:

        contacts_emda.select('COO', [10554, 10555, 10556], sel_type='at_num')
        contacts_emda.add_contacts('COO contacts', 'COO', sel_env=4, include_WAT=True)
        contacts_set = ['COO contacts']

    elif lig_ID == lig_ID_list[1]:

        contacts_emda.select('COO', [10604, 10605, 10606], sel_type='at_num')
        contacts_emda.select('OOH', [10607, 10608, 10590], sel_type='at_num')
        contacts_emda.add_contacts('OOH contacts', 'OOH', sel_env=4, include_WAT=True)
        contacts_emda.add_contacts('COO contacts', 'COO', sel_env=4, include_WAT=True)
        contacts_set = ['COO contacts', 'OOH contacts']

    elif lig_ID == lig_ID_list[2]:

        contacts_emda.select('COO', [10603, 10604, 10605], sel_type='at_num')
        contacts_emda.select('OH', [10606, 10607], sel_type='at_num')
        contacts_emda.add_contacts('OH contacts', 'OH', sel_env=4, include_WAT=True)
        contacts_emda.add_contacts('COO contacts', 'COO', sel_env=4, include_WAT=True)
        contacts_set = ['COO contacts', 'OH contacts']

    if load_results:
        contacts_emda.load(os.path.join(md_analysis_base_dir, parent_dir, f'{protein_ID}_{lig_ID}_analysis_contacts.pkl'))
    else:
        contacts_emda.run()
        contacts_emda.save(os.path.join(md_analysis_base_dir, parent_dir, f'{protein_ID}_{lig_ID}_analysis_contacts'))

    print(contacts_emda.universe)




    # contacts data (global info)
    print(f'total number of frames:\t{total_time_simulated * 100}')


    contact_set_freq = [f'{x}_freq' for x in contacts_set]
    contact_set_freq_percent = [f'{x}_freq_percent' for x in contacts_set]

    print(contact_set_freq_percent)
    contact_set_amount = [f'{x}_freq_amount' for x in contacts_set]
    for contact in contacts_set:
        emda.analyse_contacts_frequency(f'{contact}_freq', contact)
        emda.analyse_contacts_frequency(f'{contact}_freq_percent', contact, percentage=True)
        emda.analyse_contacts_amount(f'{contact}_amount', contact)

    emda.measures




    # DISPLAY GLOBAL INFO
    display_water = False


    print(f"\nprotein \t{protein_ID}")
    print('global contacts analysis')

    # label, replica, per cada replica un dictionary
    dict_list = []


    for label in contact_set_freq_percent:

        label_dict_list = []

        for variant, variant_data in emda.analyses[label].result.items():
            
            print(f'\n\033[1mMETRIC: \t{label}\n\033[0m')
            replicate_data = list(variant_data.values())
            
            for replica, values in variant_data.items():
            
                if values:  
                    print(f"\033[0mreplica \t{replica}")
                    sorted_dict = dict(sorted(values.items(), key=lambda item: item[1], reverse=True))
                    new_dict = sorted_dict.copy()

                    for residue, frequency in sorted_dict.items():
                        
                        if display_water:
                            print(f'  {residue}\t{frequency}')
                        else: 
                            if 'WAT' not in residue:
                                print(f'  {residue}\t{frequency}')
                            else: 
                                new_dict.pop(residue, None)
                else:
                    print(f"! replica {replica} has no values !")
                    raise ValueError('Some replicas dont have any values! Abort.')

                print(new_dict)
                label_dict_list.append(new_dict)
            
            dict_list.append(label_dict_list)
                    
    residue_names = [['ARG402', 'GLN595', 'LEU407', 'GLY406', 'ARG598', 'GLN175'], ['ARG402', 'LEU407', 'GLY406', 'GLN595', 'ARG598', 'ILE399']]
    plot_contacts_histogram(dict_list, residue_names)






    # plotting contacts for every file 





        # ultra shitty code
    # contacts evolution

    algo_dict = {}
    for label in contacts_set:

        # shitty code    
        dict_list = []
        contacts_emda.analyse_contacts_frequency(f'{label}_freq_percent', label, percentage=True)

        for variant, variant_data in contacts_emda.analyses[f'{label}_freq_percent'].result.items():
            
            #print(f'\n\033[1mMETRIC: \t{label}\n\033[0m')
            replicate_data = list(variant_data.values())
            

            for replica, values in variant_data.items():
            
                if values:  
                    #print(f"\033[0mreplica \t{replica}")
                    sorted_dict = dict(sorted(values.items(), key=lambda item: item[1], reverse=True))
                    new_dict = sorted_dict.copy()

                    for residue, frequency in sorted_dict.items():
                        
                        if display_water:
                            #print(f'  {residue}\t{frequency}')
                            pass
                        else: 
                            if 'WAT' not in residue:
                                #print(f'  {residue}\t{frequency}')
                                pass
                            else: 
                                new_dict.pop(residue, None)

                    dict_list.append(new_dict)

        algo_dict[label] = dict_list
                


    plot_contacts_evolution(algo_dict, residue_names, traj_list=traj_1)

    # plot contacts OH and OOH with 407, 406, ILE399, 






    emda.analyse_probability_density('pdf H13S', ['distance H13S -- O OH1', 'dihedral C11--C12--C13--C14 (C13 planarity)'], num_of_points=100, get_basins=True)
    emda.analyses['pdf H13S']
    emda.plot_probability_densities(
        'pdf H13S', 
        width_per_replica=4.5, 
        plot_measures=False,  
        plot_minima=False, 
        height_per_variant=4, 
        show_contour_lines=False, 
        levels_fill=100, 
        color_map='BuPu_r', 
        max_level=10, 
        set_names_in_axis=False, 
        xname='distance H13S -- O OH1', 
        yname='dihedral C11--C12--C13--C14 (C13 planarity)', 
        out_name=os.path.join(md_analysis_base_dir, parent_dir, 'pdf_H13S.png'), 
        xlims=[2,7],
        ylims=[0, 360],
        md_label=replica_names[0],
    )





    emda.analyse_probability_density('pdf H13R', ['distance H13R -- O OH1', 'dihedral C11--C12--C13--C14 (C13 planarity)'], num_of_points=100, get_basins=True)
    emda.analyses['pdf H13R']
    emda.plot_probability_densities(
        'pdf H13R', 
        width_per_replica=4.5, 
        plot_measures=False,  
        plot_minima=False, 
        height_per_variant=4, 
        show_contour_lines=False, 
        levels_fill=100, 
        color_map='BuPu_r', 
        max_level=10, 
        set_names_in_axis=False, 
        xname='distance H13R -- O OH1', 
        yname='dihedral C11--C12--C13--C14 (C13 planarity)', 
        out_name=os.path.join(md_analysis_base_dir, parent_dir, 'pdf_H13R.png'),
        xlims=[2,7],
        ylims=[0, 360],
        md_label=replica_names[0],
    )


    emda.analyse_probability_density('pdf H10S', ['distance H10S -- O OH1', 'dihedral C12--C11--C10--C9 (C10 planarity)'], num_of_points=100, get_basins=True)
    emda.analyses['pdf H10S']
    emda.plot_probability_densities(
        'pdf H10S', 
        width_per_replica=4.5, 
        plot_measures=False,  
        plot_minima=False, 
        height_per_variant=4, 
        show_contour_lines=False, 
        levels_fill=100, 
        color_map='BuPu_r', 
        max_level=10, 
        set_names_in_axis=False, 
        xname='distance H10S -- O OH1', 
        yname='dihedral C12--C11--C10--C9 (C10 planarity)', 
        out_name=os.path.join(md_analysis_base_dir, parent_dir, 'pdf_H10S.png'), 
        xlims=[2,7],
        ylims=[0, 360],
        md_label=replica_names[0],
    )



    emda.analyse_probability_density('pdf H10R', ['distance H10R -- O OH1', 'dihedral C12--C11--C10--C9 (C10 planarity)'], num_of_points=100, get_basins=True)
    emda.analyses['pdf H10R']
    emda.plot_probability_densities(
        'pdf H10R', 
        width_per_replica=4.5, 
        plot_measures=False,  
        plot_minima=False, 
        height_per_variant=4, 
        show_contour_lines=False, 
        levels_fill=100, 
        color_map='BuPu_r', 
        max_level=10, 
        set_names_in_axis=False, 
        xname='distance H10R -- O OH1', 
        yname='dihedral C12--C11--C10--C9 (C10 planarity)', 
        out_name=os.path.join(md_analysis_base_dir, parent_dir, 'pdf_H10R.png'),
        xlims=[2,7],
        ylims=[0, 360],
        md_label=replica_names[0],
    )



# CLI parsing
args = parser()
input_ = args.input_
load = args.load
protein = args.protein
multiple = args.multiple




if not multiple:
    # var definition
    md_base_dir = os.path.split(input_)[0]                  # base MD directory (source dir)
    parent_dir = os.path.split(input_)[1]                   # parent name of the dir that contains the md run of a given replicate
    protein_ID = protein                                    # protein name
    load_results = load                                     # if results have been calculated, load them instead of performing the calculations again

    if 'AA' in parent_dir:
            lig_ID = 'AA'
    elif '5-HETE' in parent_dir:
        lig_ID = '5-HETE'
    elif '5-HpETE' in parent_dir:
        lig_ID = '5-HpETE'
    else:
        raise ValueError(f'ligand coult not be detected in simulation {parent_dir}')

    analyze_MD(md_base_dir, parent_dir)

    #try:
    #    analyze_MD(md_base_dir, parent_dir)
    #except:
    #    print(f'Simulation analysis failed for {parent_dir}')

else:
    input_MDs = [os.path.join(input_, x) for x in os.listdir(input_) if ('precat' not in x) and ('.' not in x[0]) ]

    for input_MD in input_MDs:
        md_base_dir = os.path.split(input_MD)[0]                  # base MD directory (source dir)
        parent_dir = os.path.split(input_MD)[1]                   # parent name of the dir that contains the md run of a given replicate
        protein_ID = protein                                    # protein name
        load_results = load                                     # if results have been calculated, load them instead of performing the calculations again
        
        if 'AA' in parent_dir:
            lig_ID = 'AA'
        elif '5-HETE' in parent_dir:
            lig_ID = '5-HETE'
        elif '5-HpETE' in parent_dir:
            lig_ID = '5-HpETE'
        else:
            raise ValueError(f'ligand coult not be detected in simulation {parent_dir}')

        analyze_MD(md_base_dir, parent_dir)


        '''try:
            
        except:
            print(f'\n\n\nSimulation failed for {parent_dir}\n\n\n')'''
        
# absolute garbage of code lmao