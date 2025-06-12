# __Human ALOX15 positional regioselectivity: Scripts__
This repository contains all the scripts used in the work "Resolving inflammation: Theoretical study of the positional regioselectivity of 5S-HETE and 5S-HpETE in human ALOX15", as the Chemistry BSc thesis at Universitat Autònoma de Barcelona. 

The protocol involves molecular docking, MD and QM/MM and the scripts are organized to provide an overview on the setup as follows.

## __Molecular Docking__
Two subdirectories, `analysis` and `final_poses`. The former contains the analysis scripts to easily visualize the docking poses in Chimera (`pose_visualization.ipynb`), while the latter contains proper analysis based on structural filtering of poses (`pose_filtering.ipynb`). This filtering works by first preprocessing the poses and extracting relevant metrics (done in Chimera) and are then filtered with a dedicated notebook. 

The resulting best poses are stored in `final_poses` directory. 

## __Parametrization (Non-standard residues)__
Three subdirectories containing the protein models used (`protein_models`), iron coordination sphere (obtained from [MolBioMed GitHub](https://github.com/MolBioMedUAB/AMBER_parameters)) and ligands' `.frcmod`. For the latter, two scripts are provided to build the RESP charges job and to extract the charges (`resp_job_gen.sh` and `resp_build.sh`, respectively).     

## __Molecular Dynamics__
All the source code related to MD is stored in `src` subdirectory. The following subdirectories contain examples for each MD step:

- **`1_preprocessing`**: The `.pdb` preprocessing is performed via `AMBER` with a `preprocess.sh` script. 

- **`2_setup`**: The setup consists in building the topology via `tleap` and the actual `pmemd` job for the preprod and production steps. The former can be performed with `tleap` and the provided reference `tleap.in` files for each ligand. The latter is generated from a `create_md.sh` file that builds the directory tree from the topology (`prmtop` and `inpcrd` files).   

- **`3_trajectories`**: Contains the trajectory files (not included in the repo). 

- **`4_analysis`**: Contains `analysis.ipynb` and `analysis.py` scripts that perform the analysis of a single run or a set of runs, respectively, fetching the MD trajectory from `3_trajectories` and storing it in `4_analysis`. 

- **`5_precat_analysis`**: Performs a precatalytic analysis on the trajectories stored in `3_trajectories`. This analysis is not required for QM/MM calculations but provides insight on the sampling and precatalytic structure ratios. 


## __QM/MM__
Finally, all the source code for QM/MM calculations is stored in `src` subdirectory. The QM/MM setup is highly automated (further code polishing is needed). It extracts all the precatalytic frames and does the preprocessing to perform an optimization, scan and TS jobs for each frame. 

The resulting scan data can be used to compute barrier heights and plot the results via `scan_analysis.ipynb`. 
