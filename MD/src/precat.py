#!/Users/sergiortizropero/miniconda3/envs/md_traj/bin/python
# Sergi Ortiz @ UAB
# 8 April 2025
# Perform the routine MD analysis of one (or a set of trajectories). Note that
# multiple replicas are not supported (i.e. the analysis is independent to 
# other replicas). If multiple replicas, use the tools available in 
# `analysis`.ipynb (easy analysis)

# TODO add support for replicates

#===================#
#       IMPORTS     #
#===================#

# general 
import os
import sys
from numpy import average, std
import numpy as np
from argparse import ArgumentParser

from EMDA import EMDA, __version__
from EMDA import frame_extractor
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
                        help='Directory containing either a single or multiple MDs, with the tree as generated with create_md.sh',
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
    parser.add_argument('-s', '--selection',
                        help='Hydrogen abstraction selection. Either 10 or 13 supported.',
                        required=True,
                        type=int
                        )
    parser.add_argument('-dist', '--distance',
                        help='H -- O OH1 distance for the filtering. Default is 3 A', 
                        required=False,
                        default=3,
                        type=float
                        )
    parser.add_argument('-rc', '--rxcrd',
                        help='RC filter. Default is RC r > -2.5', 
                        required=False,
                        default=-2.5,
                        type=float
                        )
    parser.add_argument('-dih', '--dihedral',
                        help='C_sel degree planarity criterium. Default is 5 to mean 180 +- 5.', 
                        required=False,
                        default=5,
                        type=float
                        )
    
    args = parser.parse_args()

    return args

# input el dir gran 
# guarda precat els pkl i el txt
# define a function from here

def compute_precat(md_base_dir, parent_dir):
        
    lig_ID_list = ['AA', '5-HpETE', '5-HETE']
    if lig_ID not in lig_ID_list:
        raise ValueError(f'Introduced lig_ID {lig_ID} is not among the supported ligands: {lig_ID_list}')

    # topology
    topologies_1 = [i for i in os.listdir(os.path.join(md_base_dir, parent_dir)) if i.endswith('.prmtop')]
    top = os.path.join(md_base_dir, parent_dir, topologies_1[0])
    if len(topologies_1) > 1:
        raise ValueError('bruh. More than one topology files? Check md source directory.')
    print(f'\nloading prmtop\t{top}')

    # production trajectory
    prod_dir_1 = os.path.join(md_base_dir, parent_dir, 'prod/out')                    # locate the production traj
    prod_files = [int(x.strip('.nc').strip('prod_')) for x in os.listdir(prod_dir_1) if '.nc' in x]
    traj_1 = [os.path.join(prod_dir_1, f'prod_{i}.nc')  for i in prod_files]      # assume constant length 250 ns, 25 

    # load base MD (r=1)
    emda = EMDA()
    emda.load_variant(
        parameters=top,
        trajectory=traj_1,
        variant_name=protein_ID
    )

    total_time_simulated = 10 * len(traj_1)
    print(f'total time simulated: {total_time_simulated} ns')


    if lig_ID == lig_ID_list[0]:

        #===============#
        #       AA      #
        #===============#

        print(f'selecting atoms for {lig_ID_list[0]}')

        # LIGAND
        emda.select('LIG', '665', sel_type='res_num')

        emda.select('H10A', 'name H10A and resid 665')
        emda.select('H10B', 'name H10B and resid 665')

        emda.select('H13A', 'name H13A and resid 665')
        emda.select('H13B', 'name H13B and resid 665')

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

        # PROTEIN hLOX15
        # OH1 OH cofactor :664
        emda.select('OH1', 10552, sel_type='at_num') 


    elif lig_ID == lig_ID_list[1]:

        #====================#
        #       5-HpETE      #
        #====================#
        
        print(f'selecting atoms for {lig_ID_list[1]}')

        # LIGAND
        emda.select('LIG', '665', sel_type='res_num')
        
        emda.select('H10A', 'name H19 and resid 665')
        emda.select('H10B', 'name H18 and resid 665')

        emda.select('H13A', 'name H15 and resid 665')
        emda.select('H13B', 'name H14 and resid 665')

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

        # PROTEIN hLOX15
        # OH1 OH cofactor :664
        #                   O
        emda.select('OH1', 10552, sel_type='at_num') 


    elif lig_ID == lig_ID_list[2]:

        #===================#
        #       5-HETE      #
        #===================#

        print(f'selecting atoms for {lig_ID_list[2]}')

        # ligand atoms
        emda.select('LIG', '665', sel_type='res_num')

        emda.select('H10A', 'name H19 and resid 665')
        emda.select('H10B', 'name H18 and resid 665')

        emda.select('H13A', 'name H15 and resid 665')
        emda.select('H13B', 'name H14 and resid 665')

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

        # PROTEIN hLOX15
        # OH1 OH cofactor :664
        #                   O
        emda.select('OH1', 10552, sel_type='at_num') 

    else:
        raise ValueError('Something went wrong!')

    emda.add_distance('distance H10A -- O OH1', 'OH1', 'H10A')
    emda.add_distance('distance H10B -- O OH1', 'OH1', 'H10B') 

    emda.add_distance('distance H13A -- O OH1', 'OH1', 'H13A')
    emda.add_distance('distance H13B -- O OH1', 'OH1', 'H13B')

    emda.add_distance('distance H10A -- C', 'p10_C', 'H10A')
    emda.add_distance('distance H13A -- C', 'p13_C', 'H13A')  

    emda.add_dihedral('dihedral C12--C11--C10--C9 (planarity around C10)', 'p10_A', 'p10_B', 'p10_C', 'p10_D', domain=360)     # around C10
    emda.add_dihedral('dihedral C11--C12--C13--C14 (planarity around C13)', 'p13_A', 'p13_B', 'p13_C', 'p13_D', domain=360)    # around C13


    if not load_results:
        emda.run()
        
        # save the analysis
        emda.save(os.path.join(frames_dir, f'{parent_dir}_frames'))    

    else: 
        # load precomputed analysis 
        emda.load(os.path.join(frames_dir, f'{parent_dir}_frames.pkl'))

    # perform the analyses
    if selection == 10:
        print('analyzing C10 dist and dihedral')
        emda.analyse_value('dihedral_criterion', 'dihedral C12--C11--C10--C9 (planarity around C10)', mode='thres', val1=180-dih_criterium, val2=180+dih_criterium)
        emda.analyse_value('distance_criterion', 'distance H10A -- O OH1', mode='thres', val1=dist_criterium)
        # obtain the results
        criterion_1 = list(emda.analyses['dihedral_criterion'].result.values())[0]['R1']
        criterion_2 = list(emda.analyses['distance_criterion'].result.values())[0]['R1']

        print('\n\n\n COMPUTING REACTION COORDINATE \n\n\n')
        new_r_CH = list(emda.measures['distance H10A -- C'].result.values())[0]['R1']
        new_r_OH = list(emda.measures['distance H10A -- O OH1'].result.values())[0]['R1']
        rc_array = np.array(new_r_CH) - np.array(new_r_OH)

    elif selection == 13:
        print('analyzing C13 dist and dihedral')
        
        emda.analyse_value('dihedral_criterion', 'dihedral C11--C12--C13--C14 (planarity around C13)', mode='thres', val1=180-dih_criterium, val2=180+dih_criterium)
        emda.analyse_value('distance_criterion', 'distance H13A -- O OH1', mode='thres', val1=dist_criterium)
        print(emda.measures['dihedral C11--C12--C13--C14 (planarity around C13)'].result.values())

        # obtain the results
        criterion_1 = list(emda.analyses['dihedral_criterion'].result.values())[0]['R1']
        criterion_2 = list(emda.analyses['distance_criterion'].result.values())[0]['R1']

        print('\n\n\n COMPUTING REACTION COORDINATE \n\n\n')
        new_r_CH = list(emda.measures['distance H13A -- C'].result.values())[0]['R1']
        new_r_OH = list(emda.measures['distance H13A -- O OH1'].result.values())[0]['R1']
        rc_array = np.array(new_r_CH) - np.array(new_r_OH)


    else:
        raise ValueError('ep! Analysis step exploded!')

    
    # AND-type restriction
    criteria = list(np.array([criterion_1, criterion_2]).min(axis=0))


    # new 2 May 2025
    # select 2000 equidistant frames and check how many of these satisfy the conditions
    # reaction coordinate is provided in the txt file
    #N_new = int(len(traj_1)*2000/25)
    N_new = 2500
    #N_new=20
    if len(rc_array) % 2 == 0:
        trimmed_rc_array = rc_array[::int(len(rc_array)/N_new)]
        trimmed_new_r_OH = np.array(new_r_OH[::int(len(rc_array)/N_new)])
        trimmed_new_r_CH = np.array(new_r_CH[::int(len(rc_array)/N_new)])
        trimmed_criterion_1 = criterion_1[::int(len(rc_array)/N_new)]

    else: 
        trimmed_rc_array = rc_array[::int(len(rc_array-1)/N_new)]
        trimmed_new_r_OH = np.array(new_r_OH[::int(len(rc_array-1)/N_new)])
        trimmed_new_r_CH = np.array(new_r_CH[::int(len(rc_array-1)/N_new)])
        trimmed_criterion_1 = criterion_1[::int(len(rc_array-1)/N_new)]

    # new criteria of r > -2.5 then precatalytic
    precat_rc = rc_criterium
    trimmed_rc_array_copy = np.copy(trimmed_rc_array)
    filtered_rc = trimmed_rc_array_copy[trimmed_rc_array > precat_rc]
    n_new = len(filtered_rc)
    N_new = len(trimmed_rc_array)

    # new criteria: r>-2.5 and a dihedral between ...
    precat_rc = rc_criterium
    trimmed_rc_array_copy = np.copy(trimmed_rc_array)
    filtered_rc_C = trimmed_rc_array_copy[(trimmed_rc_array > precat_rc) & (trimmed_criterion_1) ]
    n_new_C = len(filtered_rc_C)

    # if no dihedral restriction is used, simply use -dih 180. 
    with open(os.path.join(frames_dir, f'precat_result_H{selection}_{dist_criterium:.1f}_{dih_criterium:.1f}_{rc_criterium:.1f}.txt'), 'a') as f:
        MD_name = f'{parent_dir}'.strip(protein_ID).split('_')[1:6]
        N = len(traj_1)*1000
        
        n = int(criteria.count(True))
        f.write(f'{MD_name[0]}'.ljust(10) +
                f'{MD_name[1]}'.ljust(10) +
                f'{MD_name[2]}'.ljust(10) +
                f'{MD_name[3]}'.ljust(10) +
                f'{MD_name[4]}'.ljust(15) +
                f'{(criterion_2.count(True)*100/(1000 * len(traj_1))):.2f}'.ljust(10) +
                f'{(criterion_1.count(True)*100/(1000 * len(traj_1))):.2f}'.ljust(10) +
                f'{N}'.ljust(8) +
                f'{n}'.ljust(8) +
                f'{(criteria.count(True)*100/(1000 * len(traj_1))):.2f}'.ljust(12) + 
                f'{-8.3145 * 300 * np.log(n/N) * (1/4.184) * (1e-3):.1f}'.ljust(12) +
                f'{N_new}'.ljust(8) +
                f'{n_new}'.ljust(8) +
                f'{-8.3145 * 300 * np.log(n_new/N_new) * (1/4.184) * (1e-3):.1f}'.ljust(12) + 
                f'{N_new}'.ljust(8) +
                f'{n_new_C}'.ljust(8) +
                f'{-8.3145 * 300 * np.log(n_new_C/N_new) * (1/4.184) * (1e-3):.1f}\n'

        )
        #f.write('MD run'.ljust(40)+'Criterion 1 (%)'.ljust(20)+'Criterion 2 (%)'.ljust(20)+'Frames'.ljust(10)+'Precat'.ljust(10)+'Precat (%)\n')

    print(parent_dir)
    #print(f'Frames that satisfy the dihedral criterion:\t{(criterion_1.count(True)*100/(1000 * len(traj_1))):.2f}\t%')
    #print(criterion_1.count(True))
    #print(f'Frames that satisfy the distance criterion:\t{(criterion_2.count(True)*100/(1000 * len(traj_1))):.2f}\t%')
    #print(criterion_2.count(True))
    #print(f'Frames that satisfy both criteria:\t\t{(criteria.count(True)*100/(1000 * len(traj_1))):.2f}\t%')
    #print(criteria.count(True))


    # TODO crear una carpeta i enganxar tots els frames precatalitics per cada simulació 
    '''
    # generate new analysis attribute with the combined results (EMDA is lacking lmao)
    emda.analyses['combined_criteria'] = emda.analyses['distance_criterion']
    emda.analyses['combined_criteria'].result = {protein_ID : {'R1' : criteria}}

    # TODO
    # randomly select a set of frames that fulfill the criteria and extract them

    # create directory to store the frames
    if not os.path.isdir('frames'):
        os.makedirs('frames')

    # extract and save the frames
    emda.export_frames_by_analysis(variant=protein_ID, replica='R1', analysis_name='combined_criteria', folder='./frames', out_name='*')
    #print(list(emda.analyses['combined_criteria'].result.values())[0]['R1'])


    #total_length = 1000 * len(traj_1)
    #print(f'total length:\t{total_length}')
    #selected_frames = frame_extractor.frame_selector(criteria, bins=10, frames_per_bin=10)
    #frame_extractor.trajectory_frame_extractor(emda.universe[protein_ID]['R1'], selected_frames, folder='./') 
    '''




# CLI parsing
args = parser()
input_ = args.input_
load = args.load
protein = args.protein
selection = args.selection
dist_criterium = args.distance
dih_criterium = args.dihedral
rc_criterium = args.rxcrd

# var definition


input_MDs = [os.path.join(input_, x) for x in os.listdir(input_) if ('precat' not in x) and ('.' not in x[0])]
frames_dir = os.path.join(input_, 'precat')


for i in range(len(input_MDs)):


    md_base_dir = os.path.split(input_MDs[i])[0]                  # base MD directory (source dir)
    parent_dir = os.path.split(input_MDs[i])[1]                   # parent name of the dir that contains the md run of a given replicate
    protein_ID = protein                                    # protein name
    load_results = load                                     # if results have been calculated, load them instead of performing the calculations again
    
    if 'AA' in parent_dir:
        lig_ID = 'AA'
    elif '5-HETE' in parent_dir:
        lig_ID = '5-HETE'
    elif '5-HpETE' in parent_dir:
        lig_ID = '5-HpETE'
    else:
        print(input_MDs[i])
        raise ValueError(f'ligand coult not be detected in simulation {parent_dir}')


    if i == 0:

        if not os.path.exists(frames_dir):
            os.makedirs(frames_dir)

        with open(os.path.join(frames_dir, f'precat_result_H{selection}_{dist_criterium:.1f}_{dih_criterium:.1f}_{rc_criterium:.1f}.txt'), 'w') as f:
            pass

        with open(os.path.join(frames_dir, f'precat_result_H{selection}_{dist_criterium:.1f}_{dih_criterium:.1f}_{rc_criterium:.1f}.txt'), 'a') as f:
            f.write(f'# Precatalytic structure computation from MD trajectory\n')
            f.write(f'# Sergi Ortiz @ UAB\n')
            f.write(f'# {protein_ID} protein\n')
            f.write(f'# Reaction coordinate r = r_C{selection}-H{selection} - r_H{selection}-O\n')  
            f.write(r'# \Delta E (kcal/mol) refers to the precatalytic structure formation probability contirbution to \Delta E^\ddag,' + '\n') 
            f.write(r'# this contibution depends heavily on how a \'precatalytic\' structure is defined (distance and dihedral filters)!' + '\n') 
            f.write(f'# Energy contribution in kcal/mol\n\n')
            f.write(f'# Three definitions of precatalytic structures:\n')
            f.write(f'# A1. Criterium 1: Distance  H{selection} < {dist_criterium:.1f} A\n')
            f.write(f'# A2. Criterium 2: Planarity C{selection} = 180 +- {dih_criterium:.1f} deg\n')
            f.write(f'# B. N_spaced equispaced frames, with the criterion r > {rc_criterium}.\n') 
            f.write(f'# C1. N_spaced equispaced frames. Criterium 1: RC r > {rc_criterium}\n')
            f.write(f'# C2. N_spaced equispaced frames. Criterium 2: Planarity C{selection} = 180 +- {dih_criterium:.1f} deg\n\n')
            f.write('\n') 
            f.write('MD run info\n')
            f.write('Ligand'.ljust(10) +
                    'Docking'.ljust(10) +
                    'Radius'.ljust(10) +
                    'Type'.ljust(10) +
                    'Receptor'.ljust(15) +
                    'A1 (%)'.ljust(10) + 
                    'A2 (%)'.ljust(10) +
                    'A_N'.ljust(8)+
                    'A_n'.ljust(8) +
                    'Precat (%)'.ljust(12) + 
                    r'A \Delta E'.ljust(12) +
                    'B_N'.ljust(8) + 
                    'B_n'.ljust(8) +
                    r'B \Delta E'.ljust(12) +
                    'C_N'.ljust(8) + 
                    'C_n'.ljust(8) +
                    r'C \Delta E' + '\n'
            )
    

    compute_precat(md_base_dir, parent_dir)
    
    '''
    try: 
        
    except:
        print(f'Something went wrong with MD run {parent_dir}')   
    '''

