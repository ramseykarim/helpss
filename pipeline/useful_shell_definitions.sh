# this file is not meant to be "run" like a script
# just run these into your current terminal (use the dot thing)
# and then use these to create directories
# Created: June 29, 2020
# Author: rkarim

mprep="/sgraraid/filaments/data/TEST4/pacs70_cal_test/manticore-prep/"
ipw="/sgraraid/filaments/data/TEST4/helpss_scratch_work/pipeline/idl_preprocess_wrapper.sh"

idlproc () {
    foldr="processed/$2/"
    if [ ! -d $foldr ] ; then
        mkdir $foldr
    fi
    cmd="$ipw -i $mprep -7 -n $1 -d archived/$2/ -o processed/$2/"
    $cmd
}

mkobjdir () {
    if [ $(hostnaem -s) == "sgra" ] ; then
        homedir="/sgraraid/filaments/$1/"
    else
        homedir="/n/sgraraid/filaments/$1/"
    fi
    otherdir="$2"
    shift 2
    if [ ! -d $homedir ] ; then
        mkdir $homedir
    fi
    if [ $? -ne 0 ] ; then
        exit 1
    fi
    pushd $homedir

    mkdir Herschel
    mkdir Spitzer
    cd Herschel

    mkdir archived
    mkdir processed
    mkdir products
    mkdir results

    while [ ! -z $1 ] ; do
        mv ${otherdir}$1 ./archived/
        shift 1
    done
    popd
}
