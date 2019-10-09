#!/bin/bash
# bash script to run manticore, updated
# author: Ramsey Karim
# this needs to be run on sgra (or another manticore-ready computer)

this_script_name="manticore.sh"

# Defaults
manticore_version="1.5"
manticore="/sgraraid/filaments/manticore-${manticore_version}/manticore"
flux_mod="45"
beta_h="1.80"
beta_c="2.10"
T_hot="16.00"
working_dir="$(pwd)/"
n_param="3"
n_cores="5"
region="unassigned_name"

# Check that manticore executable exists
if [[ ! -f $manticore ]] ; then
    printf "$0: manticore executable not found at\n  $manticore\n"
    if [[ ! "$(hostname | cut -d'.' -f1 )" == "sgra" ]] ; then
        printf "You are using the machine $(hostname); try sgra\n"
    fi
    exit 1
fi

print_usage_exit() {
    printf "${this_script_name}: usage: ./${this_script_name} [valid arguments]
    -h help (prints this message and exits)
        use -hh to get manticore help (-hhh for verbose)
    -x run (even if no arguments are present) (at least one argument is necessary to run)
    -n name of region (string)
        default -n unassigned_name
    -H halo dust model/law (decimal spectral index, or OH5, DL5, or DL3)
        if 2 parameter fit, this is the only dust model
        default -h 1.80
    -C core dust model/law (decimal spectral index, or OH5, DL5, or DL3)
       	only used if running 3 parameter
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
    -f PACS flux offset (number)
        either integer or float, shouldn't be more than 6 digits
        the PACS image with this offset must already exist (and be named properly)
        this script will do the formatting (0 padding, etc)
        default -f 45 (good for Per1 region)
    -G grid sample
        See manticore documentation
    -R random sample (integer number)
        See manticore documentation
    -S systematic uncertainty sample
        See manticore documentation
"
    exit 0
}

# Useful regular expression
re_float='^[0-9]{1,5}\.[0-9]$'
re_int='^[0-9]+$'
# Help request tracker
hflag=""

# parse arguments
while getopts 'hxn:H:C:T:d:2l:s:t:o:f:GR:S' flag ; do
    case "${flag}" in
        h) hflag="${hflag}h" ;;
        x) : ;;
        n) region="${OPTARG}" ;;
        H) beta_h="${OPTARG}" ;;
        C) beta_c="${OPTARG}" ;;
        T) T_hot="${OPTARG}" ;;
        d) working_dir="${OPTARG}"
            if [[ "$working_dir" != *\/ ]] ; then
                printf "$0: invalid directory name -- d (-h for instructions)\n"
                exit 1
            fi ;;
        2) n_param="2" ;;
        l) log="${OPTARG}" ;;
        s) single="${OPTARG}" ;;
        t) n_cores="${OPTARG}" ;;
        o) out="${OPTARG}" ;;
        f) flux_mod="${OPTARG}" ;;
        G) GRSflag="G${GRSflag}" ;;
        R) GRSflag="${GRSflag}R ${OPTARG}"
            if [[ ! $OPTARG =~ $re_int ]] ; then
                printf "$0: -R argument should be an integer\n"
                exit 1
            fi ;;
        S) GRSflag="S${GRSflag}" ;;
        *) print_usage_exit ;;
    esac
done

# Handle help request
if [[ $hflag =~ ^h$ ]] ; then
  print_usage_exit
elif [[ $hflag =~ ^h{2}$ ]] ; then
  $manticore -h ; exit 0
elif [[ $hflag =~ ^h{3,}$ ]] ; then
  $manticore -hv ; exit 0
fi

# Exit out if there are no arguments
if [[ -z $1 ]] ; then
    printf "${this_script_name}: need at least one argument (-h for instructions, -x to run with all defaults)\n"
    exit 1
fi


# A couple helper functions for parsing arguments
parse_dust() {
    # first/only argument is dust model argument
    if [[ "$1" == "DL3" ]] || [[ "$1" == "DL5" ]] || [[ "$1" == "OH5" ]] ; then
        echo "$1"
    else
        # assume it is a decimal, like "1.80"
        # manticore should complain if invalid.
        echo "pow-1000-0.1-${1}"
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

# Format flux offset as present in PACS image filename
# RE are defined above the argument parsing case statement
if [[ $flux_mod =~ $re_float ]] ; then
    flux_mod=$(printf "plus%06.1f"  $flux_mod)
elif [[ $flux_mod =~ $re_int ]] ; then
    flux_mod=$(printf "plus%06d"  $flux_mod)
else
    printf "Invalid PACS flux offset: %s (-h for instructions)\n" $flux_mod
    exit 1
fi
flux_mod="-${flux_mod}"


# Parse log argument
log=$(parse_filename $log mant_log.log)
# Parse single component argument
single=$(parse_filename $single full-${manticore_version}-Per1-pow-1000-0.1-1.80.fits)
# Parse output argument
out=$(parse_filename $out full-${manticore_version}-${region}-${dust}${Thstub}.fits)

# Print out a buch of status text
echo "============================="
echo "MANTICORE via bash wrapper"
echo "Working in ${working_dir} on region ${region}"
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
err_args="${b1}${err}${fitsstub}" # need to prefix with "-e "
for b in $b2 $b3 $b4 ; do
  err_args="${err_args} -e ${b}${err}${fitsstub}"
done

if [[ "$n_param" == "3" ]] ; then
  call_command="${manticore} -3 -T ${T_hot} -s ${single}  -c 0.0 -g 100,100"
else
  call_command="${manticore} -${n_param} -g 100"
fi
if [[ ! -z $GRSflag ]] ; then GRSflag="-${GRSflag}" ; fi
call_command="${call_command} -v -a -t ${n_cores} -D ${dust} -l ${log} -o ${out} ${GRSflag} -e ${err_args} ${flux_args}"

echo "*****************************"
echo "CALLING MANTICORE AS:"
echo "${call_command}"

$call_command

echo "=============================" >> $log
echo "Finished at $(date)" >> $log
echo "=============================" >> $log
echo "" >> $log
