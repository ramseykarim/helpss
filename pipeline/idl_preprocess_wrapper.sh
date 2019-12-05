#!/bin/bash
# bash wrapper for Tracy's IDL preprocessing for HELPSS Herschel data
# Author: Ramsey Karim


# The idea here is to run Tracy's preprocessing code from command line using variables

# This wrapper should be run with certain command line arguments:
# 1: directory containing the "level2_5/" directory (e.g. /the/directory/)
# 2: directory containing the manticore prep IDL procedures Tracy wrote
# 3: name of the object (needs a name, like "Perseus1" or "NGC1333")
# 4: directory to which to write all these images

this_script_name="idl_preprocess_wrapper.sh"

# command line arguments are set up
# default values for parameters are clearly enumerated below
default_obs_dir="$(pwd)/"
default_mprep_dir="/n/sgraraid/filaments/manticore-prep/"
default_object="unassigned_name"
default_working_dir="$(pwd)/"

# set to defaults; will change if argument specifies such
obs_directory=$default_obs_dir
Pobs_directory=""
Sobs_directory=""
mprep_directory=$default_mprep_dir
object_name=$default_object
working_dir=$default_working_dir
pacs70=false

print_usage_exit() {
    printf "${this_script_name} usage: ./${this_script_name} [valid arguments]
    -h help (prints this message and exits)
    -x run (even if no arguments are present) (at least one argument is necessary to run)
    -d data directory containing the archival Herschel data
        this directory MUST be the \"observation number\" directory
        this directory MUST contain the \"level2_5\" OR \"level2\" directory
        relative or absolute path are both ok
        default -d <current directory> ($(pwd)/)
    -P data directory containing the archival Herschel PACS data
        overrides -d
        see \"-d\" documentation above
    -S data directory containing the archival Herschel SPIRE data
        overrides -d
        see \"-d\" documentation above
    -i IDL preprocessing scripts directory
        this directory must contain MakeHerschelImages.pro, RemapImages.pro, and ConvolveHerschelImages.pro
        if you are running this on sgra, you may specify this without the '/n' prefix
        default -i ${default_mprep_dir}
    -n object name to assign in case the name is missing
        default -n ${default_object}
    -o output directory to which to write the processed FITS files
        the directory MUST already exist
        default -o <current directory> ($(pwd)/)
    Note: try to avoid relative paths. This script will cd into the working directory.
"
    exit 1
}

complain_directory() {
    # first/only arg is offending directory
    printf "$BASH_SOURCE: DIRECTORY ${1} DOES NOT EXIST
"
    exit 1
}

sanitize_directory() {
    # first/only argument is directory name
    # ensures the name ends in '/'
    # ensures this directory exists, and if not, prints usage and exits
    directory=$1
    if [[ "$directory" != *\/ ]] ; then
        directory="${directory}/"
    fi
    # Resolve relative paths to absolute ones (more or less)
    if [[ "$directory" != \/* ]] ; then
        directory="$(pwd)/${directory}"
    fi
    echo $directory
    if [[ ! -d $directory ]] ; then
        exit 1
    fi
}

# parse arguments
while getopts 'hxd:S:P:i:n:o:7' flag ; do
    case "${flag}" in
        h) print_usage_exit ;;
        x) : ;;
        d) obs_directory="$(sanitize_directory ${OPTARG})"
            if [[ $? -eq 1 ]] ; then complain_directory "${obs_directory}" ; fi ;;
        P) Pobs_directory="$(sanitize_directory ${OPTARG})"
            if [[ $? -eq 1 ]] ; then complain_directory "${Pobs_directory}" ; fi ;;
        S) Sobs_directory="$(sanitize_directory ${OPTARG})"
            if [[ $? -eq 1 ]] ; then complain_directory "${Sobs_directory}" ; fi ;;
        i) mprep_directory="$(sanitize_directory ${OPTARG})"
            if [[ $? -eq 1 ]] ; then complain_directory "${mprep_directory}" ; fi ;;
        n) object_name="${OPTARG}" ;;
        o) working_dir="$(sanitize_directory ${OPTARG})"
            if [[ $? -eq 1 ]] ; then complain_directory "${working_dir}" ; fi ;;
        7) pacs70=true ;;
        *) print_usage_exit ;;
    esac
done

if [[ -z $1 ]] ; then
    printf "${this_script_name}: need at least one argument (-x to run with all defaults)\n"
    print_usage_exit
fi


# Get PACS+SPIRE observation directories and figure out level 2/2_5

# If Pobs_directory or Sobs_directory are not set, set them
if [[ -z "$Pobs_directory" ]] ; then Pobs_directory="$obs_directory" ; fi
if [[ -z "$Sobs_directory" ]] ; then Sobs_directory="$obs_directory" ; fi
# The obs_directory (something like ../1342190326/) should contain the
#  level2_5/ or level2/ directory
# Check PACS directory
if [[ -d "${Pobs_directory}level2_5/" ]] ; then
  Plvl2or25_directory="${Pobs_directory}level2_5/"
elif [[ -d "${Pobs_directory}level2/" ]] ; then
  Plvl2or25_directory="${Pobs_directory}level2/"
else
  printf "${this_script_name}: PACS directory not valid\n"
  exit 1
fi
# Check SPIRE directory
if [[ -d "${Sobs_directory}level2_5/" ]] ; then
  Slvl2or25_directory="${Sobs_directory}level2_5/"
elif [[ -d "${Sobs_directory}level2/" ]] ; then
  Slvl2or25_directory="${Sobs_directory}level2/"
else
  printf "${this_script_name}: SPIRE directory not valid\n"
  exit 1
fi

# The directory structure of the PACS and SPIRE data is fairly standard
# We can assume the name of these subdirectories and that they each contain 1 file
if [[ "$pacs70" = true ]] ; then
    p70_source=\"$(find "${Plvl2or25_directory}HPPJSMAPB/" -name "*.*")\"
fi
p160_source=\"$(find "${Plvl2or25_directory}HPPJSMAPR/" -name "*.*")\"
s250_source=\"$(find "${Slvl2or25_directory}extdPSW/" -name "*.*")\"
s350_source=\"$(find "${Slvl2or25_directory}extdPMW/" -name "*.*")\"
s500_source=\"$(find "${Slvl2or25_directory}extdPLW/" -name "*.*")\"


# Construct the IDL call (based on the NOTES Tracy made in mprep_directory)

# Make ".run directory" shorthand and create all the import statements
idlrun=".run ${mprep_directory}"
make_herschel_images_import="${idlrun}MakeHerschelImages"
remap_images_import="${idlrun}RemapImages"
convolve_herschel_images_import="${idlrun}ConvolveHerschelImages"

# MakeHerschelImages setup and call
if [[ "$pacs70" = true ]] ; then
make_herschel_images_setup="filearr=strarr(5)
filearr(0)=${p70_source}
filearr(1)=${p160_source}
filearr(2)=${s250_source}
filearr(3)=${s350_source}
filearr(4)=${s500_source}
"
else
make_herschel_images_setup="filearr=strarr(4)
filearr(0)=${p160_source}
filearr(1)=${s250_source}
filearr(2)=${s350_source}
filearr(3)=${s500_source}
"
fi
# Note that this will dump outputs to current working directory
make_herschel_images_cmd="MakeHerschelImages, filearr, object=\"${object_name}\""
# Get filenames for these newly created files (standard filenaming scheme)
img="image"
err="error"
if [[ "$pacs70" = true ]] ; then
    p70="\"./PACS70um-"
fi
p160="\"./PACS160um-"
s250="\"./SPIRE250um-"
s350="\"./SPIRE350um-"
s500="\"./SPIRE500um-"
fits=".fits\""

# RemapImages setup and call
# Reference is SPIRE500 (largest pixels, so least number of pixels)
# Need to remap other 3 images+errors (6 total files) to the reference
if [[ "$pacs70" = true ]] ; then
remap_images_setup="reference=${s500}${img}${fits}
filearr=strarr(8)
filearr(0)=${p70}${img}${fits}
filearr(1)=${p70}${err}${fits}
filearr(2)=${p160}${img}${fits}
filearr(3)=${p160}${err}${fits}
filearr(4)=${s250}${img}${fits}
filearr(5)=${s250}${err}${fits}
filearr(6)=${s350}${img}${fits}
filearr(7)=${s350}${err}${fits}
"
else
remap_images_setup="reference=${s500}${img}${fits}
filearr=strarr(6)
filearr(0)=${p160}${img}${fits}
filearr(1)=${p160}${err}${fits}
filearr(2)=${s250}${img}${fits}
filearr(3)=${s250}${err}${fits}
filearr(4)=${s350}${img}${fits}
filearr(5)=${s350}${err}${fits}
"
fi
remap_images_cmd="RemapImages, reference, filearr"


# ConvolveHerschelImages setup and call
# Convolving to reference wavelength of 500um (worst resolution)
rmp="-remapped"
if [[ "$pacs70" = true ]] ; then
convolve_herschel_images_setup="wavearr=[70, 160, 250, 350, 500]
imarr=[${p70}\", ${p160}\", ${s250}\", ${s350}\", ${s500}\"]+\"${img}${rmp}${fits}
errarr=[${p70}\", ${p160}\", ${s250}\", ${s350}\", ${s500}\"]+\"${err}${rmp}${fits}
refwave=500
"
else
convolve_herschel_images_setup="wavearr=[160, 250, 350, 500]
imarr=[${p160}\", ${s250}\", ${s350}\", ${s500}\"]+\"${img}${rmp}${fits}
errarr=[${p160}\", ${s250}\", ${s350}\", ${s500}\"]+\"${err}${rmp}${fits}
refwave=500
"
fi
convolve_herschel_images_cmd="ConvolveHerschelImages, wavearr, imarr, errarr, refwave=refwave"

# Change directory to working directory so all file reads/writes are in there
cd $working_dir

# Make the IDL call using a "here document", which emulates interactive mode
idl <<EOF
${make_herschel_images_import}
${remap_images_import}
${convolve_herschel_images_import}

${make_herschel_images_setup}
${make_herschel_images_cmd}

${remap_images_setup}
${remap_images_cmd}

${convolve_herschel_images_setup}
${convolve_herschel_images_cmd}
EOF

printf "done with IDL preprocessing; written to ${working_dir}
"
