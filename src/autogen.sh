#!/bin/bash
###############################################################################
# Simple script to bootstrap and create the configure script, should be run
# when any control files have been changed (e.g. new source files added which
# led to changes to Makefile.am files) including:
#
#    configure.ac
#    Makefile.am
#    Makefile.in
#
############################################################

SCRIPT_NAME=$(basename $0)
INITIAL_DIR=$(pwd)
SCRIPT_DIR=$(dirname $0)
if ! echo $SCRIPT_DIR | egrep -q "(^)/" ; then
   BASE_DIR=$INITIAL_DIR/$SCRIPT_DIR
else
   BASE_DIR=$SCRIPT_DIR
fi


# Cleans the sub-directory out and refetches the empty placeholder
# from the zmap repository.
#
# Args are:    lib_repository_name  lib_local_dir
#
function clean_lib
{
    echo "Starting removing $1 from $2"

    rm -rf ./$2/*

    git checkout ./$2

    echo "Finished removing $1 from $2"

}


# Sets up a sub-directory to contain the supplied git repository code.
#
# Args are:    full_git_repository_uri lib_local_dir git_branch
#
#
function fetch_lib
{
    tmp_dir='autogen_tmp_git_checkout_dir'

    echo "Starting cloning $1 into $2"

    mkdir $tmp_dir

    git clone $3 $1 $tmp_dir || zmap_message_exit "could not clone $1 into $PWD/$tmp_dir"

    cp -rf ./$tmp_dir/* ./$2

    rm -rf ./$tmp_dir

    # SHOULD WE BE DOING THIS ????? NOT TOO SURE....
    # Make sure the placeholder files (.gitignore, README) are their original zmap versions
    #git checkout ./$2/

    echo "Finished cloning $1 into $2"

}



# Should gbtools be installed ? Default is "no"
gbtools_install='no'



while getopts ":g" opt ; do
    case $opt in
	g  ) gbtools_install='yes' ;;
    esac
done





# set up seqtools version number, this is arcane and horrible and all
# autoconf's fault. See http://sources.redhat.com/automake/automake.html#Rebuilding
# and the stuff about AC_INIT
# NOTE that seqtools_version.m4 is referenced from configure.ac
#
version_macro_file='seqtools_version.m4'
rm -f $version_macro_file

SEQTOOLS_VERSION=`git describe --abbrev=1`

echo 'dnl seqtools_version.m4 generated by autogen.sh script. Do not hand edit.'  > $version_macro_file
echo "m4_define([VERSION_NUMBER], [$SEQTOOLS_VERSION])" >> $version_macro_file
echo "SeqTools version is: $SEQTOOLS_VERSION"

# Should gbtools be installed ?
if [[ "$gbtools_install" == "yes" ]] ; then

  clean_lib 'gbtools' ./gbtools

  fetch_lib git.internal.sanger.ac.uk:/repos/git/annotools/gbtools  'gbtools'

fi


# If the gbtools autogen.sh script exists then run that. This is necessary
# for gbtools to create its gbtools_version.m4 file.
#
if [ -e "$BASE_DIR/gbtools/autogen.sh" ] ; then
  echo "found gbtools autogen.sh so running it"
  cur_dir=`pwd`
  cd $BASE_DIR/gbtools
  ./autogen.sh
  cd $cur_dir
fi


# then for now we just run autoconf
autoreconf -fi -v || echo "autoreconf failed...."



