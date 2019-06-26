#!/bin/bash

machine="rkarim@jupiter.astro.umd.edu"
directory="/n/sgraraid/filaments/data/TEST4/helpss_scratch_work/"
fig_filename="Figure_X_current.png"
fig_path="/home/rkarim/Downloads/"
cmd_line_args="PACS160um Per1"
local_image_dir="/home/ramsey/Documents/Research/Filaments/images/"
# python="/jupiter/rkarim/anaconda3/bin/python"
# python="/jupiter/rkarim/anaconda3/envs/py36/bin/python"
python="python"

# CALL SIGNATURE:
# $1 is either "fetch" or "push"
# $2 is the filename (in current directory)
# $3 is optional, either norun for just putting the file there,
#  or noplot to avoid scp-ing the generated plot back here

if [[ $1 != "fetch" ]] && [[ $1 != "push" ]] ; then
    echo "what? fetch or push";
    exit 1
fi
filename=$2
if [[ -z ${filename} ]] ; then
    echo "please supply a filename"
    exit 1
fi

if [[ $1 == "fetch" ]] ; then
    echo "getting $filename as new_${filename} (for safety)";
    scp ${machine}:${directory}${filename} ./new_${filename}
elif [[ $1 == "push" ]] ; then
    echo "pushing ${filename}...";
    scp ./${filename} ${machine}:${directory}${filename}
    if [[ $3 != "norun" ]] ; then
	ssh ${machine} "${python} ${directory}${filename}"
	    if [[ $3 != "noplot" ]] ; then
	        scp ${machine}:${fig_path}${fig_filename} ${local_image_dir}
	        eog ${image_dir}${fig_filename};
	    fi
    fi
fi
