# bash script to run manticore, updated
# author: ramsey
# this needs to be run on sgra (or another manticore-ready computer)

flux_mod="plus045"
beta_h="1.80"
beta_c="2.10"
T_hot="16.00"
working_dir="$(pwd)/"
three_param=true # two parameter if false. not set up for 1 parameter yet
log="${working_dir}mant_log.log"

print_usage() {
    printf "Usage:\n -h for beta_hot\n -c for beta_cold\n -T for hot temperature\n -d for directory\n"
}

# parse arguments
while getopts 'h:c:T:d:2l:' flag; do
    case "${flag}" in
        h) beta_h="${OPTARG}" ;;
        c) beta_c="${OPTARG}" ;;
        T) T_hot="${OPTARG}" ;;
        d) working_dir="${OPTARG}" ;;
        2) three_param=false ;;
        l) log="${OPTARG}" ;;
        *) print_usage
            exit 1 ;;
    esac
done

# format power laws if not using full dust model
if [[ $beta_h != "DL3" ]] && [[ $beta_h != "OH5" ]] ; then
    beta_h="pow-1000-0.1-${beta_h}"
fi
if [[ $beta_c != "DL3" ]] && [[ $beta_c != "OH5" ]] ; then
    beta_c="pow-1000-0.1-${beta_c}"
fi

echo "Halo dust model: ${beta_h}"
if [[ "$three_param" = true ]] ; then
    echo "Core dust model: ${beta_c}"
    echo "Halo temperature for 3-param fit: ${T_hot}"
fi
echo "Working in ${working_dir}"
if [[ "$three_param" = true ]] ; then
    printf "Running 3 parameter (two component) fit\n"
else
    printf "Running 2 parameter (single component) fit\n"
fi
echo "Writing to ${log}"

# put together manticore arguments
if [[ "$three_param" = true ]] ; then
    dust="${beta_h},${beta_c}"
else
    dust="${beta_h}"
fi




