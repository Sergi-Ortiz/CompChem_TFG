#!/bin/bash
#SBATCH -J opt_<FRAME>
#SBATCH -e %x.%j.err
#SBATCH -o %x.%j.log
#SBATCH --mail-user=sergi.ortizr@autonoma.cat
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH -p normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --time=30-00:00

# module load
ml Turbomole
ml Chemshell

# node computing options
export OMP_NUM_THREADS=8    #`echo $SLURM_JOB_NODELIST | wc -w`
export PARNODES=8           #`echo $SLURM_JOB_NODELIST | wc -w`
export PARA_ARCH=SMP
ulimit -s unlimited


# source dir    $SLURM_SUBMIT_DIR
# working dir   $TMP_DIR

# file management (input)
echo "<FRAME>    setting up input files"

cp $SLURM_SUBMIT_DIR/<FRAME>* $TMP_DIR              # prmtop, inpcrd, pdb, act_list
cp $SLURM_SUBMIT_DIR/qm_list $TMP_DIR           # qm_list
cp $SLURM_SUBMIT_DIR/opt.chm $TMP_DIR           # ChemShell job
ls -l $TMP_DIR > $SLURM_SUBMIT_DIR/input_files.txt


# optimization
echo "<FRAME>    starting QM/MM optimization @ `date`"
cd $TMP_DIR/
srun chemsh opt.chm > $SLURM_SUBMIT_DIR/opt.out
echo "<FRAME>    finished QM/MM optimization @ `date`"


# file management (output)
echo "<FRAME>    retrieving input files"
mkdir $SLURM_SUBMIT_DIR/output
cp $TMP_DIR/*.pdb $SLURM_SUBMIT_DIR/output
cp $TMP_DIR/*.c   $SLURM_SUBMIT_DIR/output

ls -l $TMP_DIR > $SLURM_SUBMIT_DIR/output/output_files.txt

echo ""