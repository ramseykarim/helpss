# bash script to run manticore, updated
# author: ramsey
# this needs to be run on sgra (or another manticore-ready computer)

flux_mod="-plus000045"
perr_mod="-plus5.0pct"
serr_mod="-plus1.5pct"
beta_h="1.80"
beta_c="2.10"
T_hot="16.00"
working_dir="$(pwd)/"
n_param="3"
log="mant_log.log"
n_cores="5"
region="Per1"
single="full-1.4.1-Per1-pow-1000-0.1-1.80.fits"

print_usage() {
    printf "Usage:\n -h for beta_hot\n -c for beta_cold\n -T for hot temperature\n -d for directory\n"
}

# parse arguments
while getopts 'h:c:T:d:2l:s:' flag; do
    case "${flag}" in
        h) beta_h="${OPTARG}" ;;
        c) beta_c="${OPTARG}" ;;
        T) T_hot="${OPTARG}" ;;
        d) working_dir="${OPTARG}" ;;
        2) n_param="2" ;;
        l) log="${OPTARG}" ;;
        s) single="${OPTARG}" ;;
        *) print_usage
            exit 1 ;;
    esac
done
log="${working_dir}${log}"

# format power laws if not using full dust model
if [[ $beta_h != "DL3" ]] && [[ $beta_h != "OH5" ]] ; then
    beta_h="pow-1000-0.1-${beta_h}"
fi
if [[ $beta_c != "DL3" ]] && [[ $beta_c != "OH5" ]] ; then
    beta_c="pow-1000-0.1-${beta_c}"
fi

# put together some manticore arguments
if [[ "$n_param" == "3" ]] ; then
    dust="${beta_h},${beta_c}"
    Thstub="Th${T_hot}"
else
    dust="${beta_h}"
    Thstub=""
fi

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
out="${working_dir}full-1.4.1-${region}-${dust}${Thstub}.fits"
manticore="/sgraraid/filaments/manticore-1.4.1/manticore"
rccstub="-remapped-conv"
fitsstub=".fits"

img="image${rccstub}"
err="error${rccstub}"

b1="${dirstub}PACS160um-"
b2="${dirstub}SPIRE250um-"
b3="${dirstub}SPIRE350um-"
b4="${dirstub}SPIRE500um-"
flux_args="${b1}${img}${flux_mod}${fitsstub} ${b2}${img}${fitsstub} ${b3}${img}${fitsstub} ${b4}${img}${fitsstub}"
err_args="${b1}${err}${perr_mod}${fitsstub}" # need to prefix with "-e "
for b in b2 b3 b4 ; do
  err_args="${err_args} -e ${b}${err}${serr_mod}${fitsstub}"
done

if [[ "${n_param}" == "3" ]] ; then
  call_command="${manticore} -3 -T ${T_hot} -s ${working_dir}${single}  -c 0.0"
else
  call_command="${manticore} -${n_param}"
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
