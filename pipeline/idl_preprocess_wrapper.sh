#!/bin/bash
# bash wrapper for Tracy's IDL preprocessing for HELPSS Herschel data
# Author: Ramsey Karim


# The idea here is to run Tracy's preprocessing code from command line using variables

# This wrapper should be run with certain command line arguments:
# 1: directory containing the "level2_5/" directory (e.g. /the/directory/)
# 2: directory containing the manticore prep IDL procedures Tracy wrote
# 3: name of the object (needs a name, like "Perseus1" or "NGC1333")
# 4: directory to which to write all these images

# COMMAND LINE ARGS ARE NOT YET SETUP, SO USING THESE DEFAULS BELOW
# delet this eventually
default_obs_dir="/n/sgraraid/filaments/data/Per/1342190326/"
default_mprep_dir="/n/sgraraid/filaments/manticore-prep/"
default_object="NGC1333"
default_working_dir="./test_folder/"

# should reference an argument
obs_directory=$default_obs_dir
mprep_directory=$default_mprep_dir
object_name=$default_object
# Change directory to working directory so all file reads/writes are in there
cd $default_working_dir

# The obs_directory (something like ../1342190326/) should contain the level2_5/ directory
lvl25_directory="${default_obs_dir}level2_5/"
# The directory structure of the PACS and SPIRE data is fairly standard
# We can assume the name of these subdirectories and that they each contain 1 file
p160_source=\"$(find "${lvl25_directory}HPPJSMAPR/" -name "*.*")\"
s250_source=\"$(find "${lvl25_directory}extdPSW/" -name "*.*")\"
s350_source=\"$(find "${lvl25_directory}extdPMW/" -name "*.*")\"
s500_source=\"$(find "${lvl25_directory}extdPLW/" -name "*.*")\"

# Construct the IDL call (based on the NOTES Tracy made in mprep_directory)

# Make ".run directory" shorthand and create all the import statements
idlrun=".run ${mprep_directory}"
make_herschel_images_import="${idlrun}MakeHerschelImages"
remap_images_import="${idlrun}RemapImages"
convolve_herschel_images_import="${idlrun}ConvolveHerschelImages"

# MakeHerschelImages setup and call
make_herschel_images_setup="filearr=strarr(4)
filearr(0)=${p160_source}
filearr(1)=${s250_source}
filearr(2)=${s350_source}
filearr(3)=${s500_source}
"
# Note that this will dump outputs to current working directory
make_herschel_images_cmd="MakeHerschelImages, filearr, object=\"${object_name}\""
# Get filenames for these newly created files (standard filenaming scheme)
img="image"
err="error"
p160="\"./PACS160um-"
s250="\"./SPIRE250um-"
s350="\"./SPIRE350um-"
s500="\"./SPIRE500um-"
fits=".fits\""

# RemapImages setup and call
# Reference is SPIRE500 (largest pixels, so least number of pixels)
# Need to remap other 3 images+errors (6 total files) to the reference
remap_images_setup="reference=${s500}${img}${fits}
filearr=strarr(6)
filearr(0)=${p160}${img}${fits}
filearr(1)=${p160}${err}${fits}
filearr(2)=${s250}${img}${fits}
filearr(3)=${s250}${err}${fits}
filearr(4)=${s350}${img}${fits}
filearr(5)=${s350}${err}${fits}
"
remap_images_cmd="RemapImages, reference, filearr"


# ConvolveHerschelImages setup and call
# Convolving to reference wavelength of 500um (worst resolution)
rmp="-remapped"
convolve_herschel_images_setup="wavearr=[160, 250, 350, 500]
imarr=[${p160}\", ${s250}\", ${s350}\", ${s500}\"]+\"${img}${rmp}${fits}
errarr=[${p160}\", ${s250}\", ${s350}\", ${s500}\"]+\"${err}${rmp}${fits}
refwave=500
"
convolve_herschel_images_cmd="ConvolveHerschelImages, wavearr, imarr, errarr, refwave=refwave"

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

printf "DONE; RETURNED TO BASH
"
