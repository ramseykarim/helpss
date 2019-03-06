#!/bin/bash

sgra="rkarim@sgra.astro.umd.edu"
directory="/sgraraid/filaments/data/TEST4/regridding_stuff/"
fig_filename="Figure_X_current.png"
cmd_line_args="PACS160um Per1"
#python="/jupiter/rkarim/anaconda3/bin/python"
python="python"

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
    scp ${sgra}:${directory}${filename} ./new_${filename}
elif [[ $1 == "push" ]] ; then
    echo "pushing ${filename}...";
    scp ./${filename} ${sgra}:${directory}${filename}
    if [[ $3 != "norun" ]] ; then
	ssh ${sgra} "${python} ${directory}${filename}"
	if [[ $3 != "noplot" ]] ; then
	    scp ${sgra}:~/${fig_filename} ./
	    eog ${fig_filename};
	fi
    fi
fi



