#!/bin/bash

# Copyright [2018-2022] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


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

    # Copy the entire contents of the temp directory to the destination directory. Note that 
    # the dot on the end of the source directory is essential for including hidden files. We
    # need to include the hidden .git directory so that we can determine the correct git version
    # number for the library.
    cp -rf ./$tmp_dir/. ./$2

    rm -rf ./$tmp_dir

    echo "Finished cloning $1 into $2"

}



checkout_only='no'

#git_location='git.internal.sanger.ac.uk:/repos/git/annotools'
git_location='git@github.com:Ensembl'

# gbtools stuff.
gbtools_install='yes'
gbtools_repo='gbtools'


while getopts ":cgr:z" opt ; do
    case $opt in
        c  ) checkout_only='yes' ;;
        g  ) gbtools_install='yes' ;;
        r  ) git_location=$OPTARG ;;
        z  ) gbtools_install='no' ;;
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

  fetch_lib "$git_location/$gbtools_repo" 'gbtools' '-b production'

fi


# Sometimes we want a tree containing any necessary subdirectories (aceconn etc)
# but don't want to run any autoconf stuff.
#
if [ "$checkout_only" = 'yes' ] ; then
  echo "Subdirectories installed, exiting before autoreconf."

  exit 0
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
echo "seqtools running: autoreconf -fi -v"
autoreconf -fi -v || echo "autoreconf failed...."
echo "Done"


exit 0

