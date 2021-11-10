# can just "run" this file to set up the calibration environment
# same as useful_shell_definitions_sgra, run with . into terminal

# append this directory to PATH so we can use calibrate.py and such
PATH=$(pwd):$PATH
PATH="/n/sgraraid/filaments/data/TEST4/pacs_calibrate":$PATH

# useful for printing the Planck beta files from a given set of regions
print_beta () {
  file_list=$(ls ${1}Herschel/processed/*/beta-stats.txt)
  for i in $file_list ; do
    echo $i
    echo "-----------"
    cat $i
    echo "-----------"
    echo "-----------"
  done
}
