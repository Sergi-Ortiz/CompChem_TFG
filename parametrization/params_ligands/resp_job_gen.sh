#!/opt/homebrew/bin/bash
## !/bin/bash
# Created by Sergi Ortiz Ropero @ 15/03/2025
# RESP CHARGE NON-STANDARD RESIDUE PARAMETRIZATION
# Inputs the ligand (mol2 or pdb) with Hs
# Generates the Gaussian job needed

# parse flags
input=''
input_format=''
output_dir=''
charge=0
multiplicity=1
processors=16
jobdir='/home/sortiz'            # job directory at picard
verbose=0

print_usage() {
    printf "\e[1mUsage:\e[0m\n"
    printf '%b\n' "-i\tINPUT FILE\tLigand input structure. See -f for the input format"
    printf '%b\n' "-f\tINPUT FORMAT\tLigand input structure format. Supported: pdb or mol2"
    printf '%b\n' "-o\tOUTPUT DIR\tWorking directory for the parametrization protocol. All files will be derivatives of INPUT FILE."
    printf '%b\n' "-c\tCHARGE\t\tGlobal charge of the ligand. Default 0."
    printf '%b\n' "-m\tMULTIPLICITY\tSpin multiplicity of the ligand, computed as 2S+1. Default 1 (singlet)."
    printf '%b\n' "-p\tPROCESSORS\tNumber of processors to use for the Gaussian Job. Default 16."
    printf '%b\n' "-p\tPROCESSORS\tHome directory in picard /home/<username>"
    printf "\n"
    printf '%b\n' "\e[1mObs.\e[0m \e[0m This script will generate a subdir where the ligand structure is located with the Gaussian job and a log file. Use that dir for everything RESP related.\e[0m"
    printf "\n"
    printf '%b\n' "\e[1mExample.\e[0m bash resp_job_gen.sh -i <ligand>.pdb -f pdb -c 0 -m 1 -p 16"
    printf "\n"
}

while getopts 'i:f:c:m:p:j:h' flag; do
  case "${flag}" in
    i) input="${OPTARG}" ;;
    f) input_format="${OPTARG}" ;;
    c) charge="${OPTARG}" ;;
    m) multiplicity="${OPTARG}" ;;
    p) processors="${OPTARG}" ;;
    j) jobdir="${OPTARG}" ;;
    h)
       print_usage
       exit 0 ;;
    *)
       print_usage
       exit 0 ;;
  esac
done


script_name=$(basename "$0")
printf "Generating a RESP Gaussian16 job with $script_name\n"
printf '%(%d-%m-%Y %H:%M:%S)T\n' -1

#==========================#
#       PREPROCESSING      #
#==========================#

# check input exist
if [ "$input_format" != 'pdb' ] && [ "$input_format" != 'mol2' ]; then
    echo "Input format $input_format is not accepted. Try pdb or mol2."
    echo "Script terminated with an error."
    exit 1
fi

# obtain ligand base name
LIG_NAME=$(basename "$input" ".$input_format")

# check input file
if [ ! -f "$input" ]; then
    echo "Input file $input not found at $PWD base dir"
    echo "Script terminated with an error."
    exit 1
fi

# check output directory
#TODO

# create subdirs
SUBDIR="$PWD/"$LIG_NAME"_resp"
echo "creating $SUBDIR"
mkdir $SUBDIR
#touch $SUBDIR/resp_$LIG_NAME.txt
#LOG_FILE="$SUBDIR/resp_$LIG_NAME.txt"


# TODO customize name output files
#if [ ! -f "$output" ]; then
#    echo "No output name given"
#    # assign output names
#    OUTPUT=""$LIG_NAME"_resp.com"
#    GESP_OUTPUT=""$LIG_NAME"_resp.gesp"
#fi

# output names
OUTPUT="${LIG_NAME}_resp.com"
GESP_OUTPUT="${LIG_NAME}_resp.gesp"

#========================#
#       SCRIPT CORE      #
#========================#

# generate Gaussian base output
antechamber -i $input -fi $input_format -o $OUTPUT -fo gcrt -gv 1 -ge $GESP_OUTPUT  -m $multiplicity -nc $charge -pf 1 -gn "%nprocshared=$processors" -gm "%mem=2GB"
mv $OUTPUT $SUBDIR


# modify Gaussian job according to flags

#NEW_LINE=" $charge   $multiplicity"
#G_CHG_MULT=$(sed '8q;d' $SUBDIR/$OUTPUT)
#echo "current chg and mult: $G_CHG_MULT"
#sed -i '' "8s,.*,$NEW_LINE," "$SUBDIR/$OUTPUT"
#G_CHG_MULT=$(sed '8q;d' $SUBDIR/$OUTPUT)
#echo "new chg and mult:    $G_CHG_MULT"

# modify .com input message
sed -i '' "8s,.*,Amber RESP parametrization of ligand $LIG_NAME," "$SUBDIR/$OUTPUT"

# add number of processors and memory
#sed -i '' '2i\'$'\n''%nprocshared=16'$'\n' "$SUBDIR/$OUTPUT"
#sed -i '' '3i\'$'\n''%mem=2GB'$'\n' "$SUBDIR/$OUTPUT"

# change .chk output name
sed -i '' "3s,molecule,"$jobdir"/${LIG_NAME}_resp/${LIG_NAME}_resp.chk," "$SUBDIR/$OUTPUT"

# print final gaussian job
printf "\nGaussian Job:\n\n"
cat "$SUBDIR/$OUTPUT"
echo "Finished generating Gaussian job :)"

# create the slurm script to batch the job for Picard
G_JOB_NAME="job_$LIG_NAME.slm"
echo "generating $G_JOB_NAME"
cat > "$SUBDIR/$G_JOB_NAME" <<EOF
#!/bin/bash
#SBATCH -J ${LIG_NAME}_resp
#SBATCH -e /home/sortiz/${LIG_NAME}_resp/${LIG_NAME}_resp.%j.err
#SBATCH -o /home/sortiz/${LIG_NAME}_resp/${LIG_NAME}_resp.%j.out
#SBATCH -p normal
#SBATCH -n 1
#SBATCH -c $processors
#SBATCH -t 2-02:00

# script generated automatically by resp_job_gen.sh
# Amber RESP charge parametrization
# Gaussian HF/6-31(d) job 

# gaussian settings
ml Gaussian/16.C.02-AVX2
#ml Gaussian/16.B.01-LEGACY 

export GAUSS_SCRDIR=\$SLURM_SUBMIT_DIR
printf "gaussian dir:\t\$GAUSS_SCRDIR\n"

#  Modify the input and output files!
INPUT_FILES=$OUTPUT      
OUTPUT_FILES=${LIG_NAME}_resp.log

printf "Started: \`date\`\n"

cp \$INPUT_FILES \$TMP_DIR
cd \$TMP_DIR

srun g16 < \$SLURM_SUBMIT_DIR/\$INPUT_FILES > \$SLURM_SUBMIT_DIR/\$OUTPUT_FILES

printf "\n\n\`ls -l\`\n\n"
mv *  \$SLURM_SUBMIT_DIR/

printf "Finished: \`date\`\n"
EOF
chmod +x "$G_JOB_NAME"
printf "\nJob script:\n\n"
cat "$SUBDIR/$G_JOB_NAME"


# instruct the user to copy the directory to Picard
printf "\n\n"
printf "Access Picard and create /home/sortiz/"$LIG_NAME"_resp/\n\n"
printf "Use\n\n\tscp -r $SUBDIR sortiz@158.109.174.75:/home/sortiz/\n\nto copy the files to Picard"
printf "\n\n"
printf "Use\n\n\tscp -r sortiz@158.109.174.75:/home/sortiz/"$LIG_NAME"_resp/ $PWD \n\nto copy the results from Picard once the job is finished"
printf "\n\n"

echo $OUTPUT_CHK

exit 0