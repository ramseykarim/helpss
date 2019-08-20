#!/bin/bash
# bash script to run manticore, updated
# author: ramsey
# this needs to be run on sgra (or another manticore-ready computer)

manticore_version="1.4.2"
flux_mod="-plus000045"
perr_mod="-plus5.0pct"
serr_mod="-plus1.5pct"
beta_h="1.80"
beta_c="2.10"
T_hot="16.00"
working_dir="$(pwd)/"
n_param="3"
n_cores="5"
region="Per1"

print_usage_exit() {
    printf "manticore.sh: usage: ./manticore.sh [valid arguments]
    -h help (prints this message and exits)
    -x run (even if no arguments are present) (at least one argument is necessary to run)
    -H halo dust model/law (decimal spectral index, or OH5, DL5, or DL3)
        default -h 1.80
    -C core dust model/law (decimal spectral index, or OH5, DL5, or DL3)
        default -c 2.10
    -T halo temperature in Kelvin (number) (only used if running 3 parameter)
        default -T 16.00
    -d working directory containing input data
        NEEDS to end in '/'
        default -d <current directory> ($(pwd)/)
    -2 run 2-parameter (single-component) Manticore fit (no argument)
        default (no flag) is 3 parameter
    -l filename of log file
        see pattern rules for -o
        not overwritten, just appended to. created if does not exist.
        default -l mant_log.log in working directory (-d)
    -s filename of 2 parameter result (only used if running 3 parameter)
        see pattern rules for -o
        the file will only be read
        default -s full-${manticore_version}-Per1-pow-1000-0.1-1.80.fits, in working directory (-d)
    -t number of CPU threads to use (positive integer)
        default -t 5
    -o output filename
        if ends in '/', assumes this specifies a write directory and uses default filename
        if contains '/' but does not end in '/', assumes this specifies full file name and path to which to write
        if not contains '/', assumes this specifies only the filename and writes to working directory
        if not set, writes to default filename in working directory
"
    exit 1
}

# parse arguments
while getopts 'hxH:C:T:d:2l:s:t:o:' flag ; do
    case "${flag}" in
        h) print_usage_exit ;;
        x) : ;;
        H) beta_h="${OPTARG}" ;;
        C) beta_c="${OPTARG}" ;;
        T) T_hot="${OPTARG}" ;;
        d) working_dir="${OPTARG}"
            if [[ "$working_dir" != *\/ ]] ; then
                printf "$0: invalid directory name -- d\n"
                print_usage_exit
            fi ;;
        2) n_param="2" ;;
        l) log="${OPTARG}" ;;
        s) single="${OPTARG}" ;;
        t) n_cores="${OPTARG}" ;;
        o) out="${OPTARG}" ;;
        *) print_usage_exit ;;
    esac
done

if [[ -z $1 ]] ; then
    printf "$0: need at least one argument (-x to run with all defaults)\n"
    print_usage_exit
fi


# A couple helper functions for parsing arguments
parse_dust() {
    # first/only argument is dust model argument
    if [[ "$1" != "DL3" ]] && [[ "$1" != "DL5" ]] && [[ "$1" != "OH5" ]] ; then
        echo "pow-1000-0.1-${1}"
    else
        echo "$1"
    fi
}

parse_filename() {
    # first arg is input
    # second arg is default filename
    # default directory is working directory
    if [[ -z $1 ]] ; then
        echo "${working_dir}$2"
    elif [[ "$1" == *\/ ]] ; then
        echo "${1}$2"
    elif [[ "$1" != *\/* ]] ; then
        echo "${working_dir}${1}"
    else
        echo "$1"
    fi
}

# Format power laws if not using full dust model
beta_h=$(parse_dust $beta_h)
beta_c=$(parse_dust $beta_c)
# Put together some manticore arguments depending on 2/3 parameter run
if [[ "$n_param" == "3" ]] ; then
    dust="${beta_h},${beta_c}"
    Thstub="-Th${T_hot}"
else
    dust="${beta_h}"
    Thstub=""
fi

# Parse log argument
log=$(parse_filename $log mant_log.log)
# Parse single component argument
single=$(parse_filename $single full-${manticore_version}-Per1-pow-1000-0.1-1.80.fits)
# Parse output argument
out=$(parse_filename $out full-${manticore_version}-${region}-${dust}${Thstub}.fits)

# Print out a buch of status text
echo "============================="
echo "MANTICORE via bash wrapper"
echo "Working in ${working_dir}"
if [[ "$n_param" == "3" ]] ; then
    printf "Running 3 parameter (two component) fit\n"
else
    printf "Running 2 parameter (single component) fit\n"
fi
echo "Writing to ${log}"
echo "Halo dust model: ${beta_h}"
if [[ "$n_param" == "3" ]] ; then
    echo "Core dust model: ${beta_c}"
    echo "Halo temperature for 3-param fit: ${T_hot}"
fi
echo "Starting now!" $(date)

# Write a bunch of stuff to log
echo "=============================" >> $log
echo "MANTICORE via bash wrapper" >> $log
echo "Started at $(date)" >> $log
echo "=============================" >> $log
if [[ "$n_param" == "3" ]] ; then
  echo "Running 3 parameter (two component) fit" >> $log
  echo "T_hot: ${T_hot}" >> $log
else
  echo "Running 2 parameter (single component) fit" >> $log
fi
echo "dust: ${dust}" >> $log
echo "working from directory:" >> $log
echo "${working_dir}" >> $log
echo "using ${n_cores} cores on $(hostname) as $(whoami)" >> $log



# Really put together the command line arguemnts
manticore="/sgraraid/filaments/manticore-${manticore_version}/manticore"
rccstub="-remapped-conv"
fitsstub=".fits"
dust=$(echo $dust | sed 's/-/:/g')

img="image${rccstub}"
err="error${rccstub}"

b1="${working_dir}PACS160um-"
b2="${working_dir}SPIRE250um-"
b3="${working_dir}SPIRE350um-"
b4="${working_dir}SPIRE500um-"
flux_args="${b1}${img}${flux_mod}${fitsstub} ${b2}${img}${fitsstub} ${b3}${img}${fitsstub} ${b4}${img}${fitsstub}"
err_args="${b1}${err}${perr_mod}${fitsstub}" # need to prefix with "-e "
for b in $b2 $b3 $b4 ; do
  err_args="${err_args} -e ${b}${err}${serr_mod}${fitsstub}"
done

if [[ "$n_param" == "3" ]] ; then
  call_command="${manticore} -3 -T ${T_hot} -s ${single}  -c 0.0 -g 100,100"
else
  call_command="${manticore} -${n_param} -g 100"
fi
call_command="${call_command} -v -a -t ${n_cores} -D ${dust} -l ${log} -o ${out} -e ${err_args} ${flux_args}"

echo "*****************************"
echo "CALLING MANTICORE AS:"
echo "${call_command}"

$call_command

echo "=============================" >> $log
echo "Finished at $(date)" >> $log
echo "=============================" >> $log
echo "" >> $log
