#!/opt/homebrew/bin/bash
## !/bin/bash
# Created by Sergi Ortiz Ropero @ 22/03/2025
# AMBER PROTEIN-LIGAND PDB PREPROCESSING
# generates a '_4amber' preprocessed pdb ready to be used in `tleap``


# parse flags
input=''        # is the path to that file
no_ligand=1     # no ligand is false (we have ligand). Either 0 or 1
ligand_name=''  # 'AA', 'HETE', 'HpETE'

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
    printf '%b\n' "\e[1mObs.\e[0m \e[0m Do execute this script with in the same directory as the input file, or use bash ../preprocess.sh ... as command.\e[0m"
    printf "\n"
    printf '%b\n' "\e[1mExample.\e[0m bash preprocess.sh -i <complex_name>.pdb -l 0 [only if no ligand is present] -n [AA, HETE, HpETE]"
    printf "\n"
}

while getopts 'i:l:n:h' flag; do
  case "${flag}" in
    i) input="${OPTARG}" ;;
    l) no_ligand="${OPTARG}" ;;
    n) ligand_name="${OPTARG}" ;;
    h)
       print_usage
       exit 0 ;;
    *)
       print_usage
       exit 0 ;;
  esac
done


script_name=$(basename "$0")
printf "Preprocessing ${input} with $script_name\n"
printf '%(%d-%m-%Y %H:%M:%S)T\n' -1

# get absolute path to file
ABS_PATH="`pwd`/${input}"
echo "$ABS_PATH"

# get protein ligand complex name
COMPLEX_NAME=$(basename ${ABS_PATH} ".pdb")
echo "$COMPLEX_NAME"

# get directory where it is
DIR_PATH=$(dirname ${ABS_PATH})
echo "$DIR_PATH"

# PREPROCESSING
cd ${DIR_PATH} 
echo "`ls -l`"

# TODO remove all redundant data from the pdb (pdb4amber does this automatically)
#grep -e '^ATOM\|^HETATM\|^TER\|^END' ${COMPLEX_NAME}.pdb ${COMPLEX_NAME}_4amber.pdb
#sed -i 's/HOH/WAT/g' $thing_clean

if [ "$no_ligand" -eq 1 ]; then
    # add LIG as ligand RESNAME
    # lines 10592 - 10645, add 'MOL B 665' -- necessary for pdb4amber
    if [ "$ligand_name" = "AA" ]; then
      echo "$ligand_name"
      sed -r '10592,10644s,ACD     1,LIG B 665,' "${COMPLEX_NAME}.pdb" > "${COMPLEX_NAME}_modded.pdb"
    fi
    if [ "$ligand_name" = "HETE" ]; then
      sed -r '10592,10660s,1       1,LIG B 665,' "${COMPLEX_NAME}.pdb" > "${COMPLEX_NAME}_modded.pdb"
    fi
    if [ "$ligand_name" = "HpETE" ]; then
      sed -r '10592,10660s,MOL     1,LIG B 665,' "${COMPLEX_NAME}.pdb" > "${COMPLEX_NAME}_modded.pdb"
    fi
    
    pdb4amber -i "${COMPLEX_NAME}_modded.pdb" -o "${COMPLEX_NAME}_4amber.pdb"
else
    pdb4amber -i "${COMPLEX_NAME}.pdb" -o "${COMPLEX_NAME}_4amber.pdb"
fi

printf "\nPreprocessing finished successfully :)\n\n\n"

printf "Following steps:\n"
printf " - Check the ligand name is MOL in both the .pdb and the .frcmod\n"
printf " - Ligand LIG number is 665\n"
printf " - Check the .frcmod parameters of the ligand and iron coord sphere\n"
printf " - Use tleap.in to produce the coordinates and the topology.\n\n"

exit 0