#!/bin/bash
#
# runconfig.sh
#
# Runs ./configure to use the required flags for specific 
# machines architectures. Also sets the CFLAGS to compile
# with debugger info.
#
# Takes a single argument which specifies the --prefix; this
# can be omitted to use the default prefix.
#

opsys="`uname -s`"
arch="`uname -m`"
machine="`uname -s`-`uname -m`"
pkgpath=""


ldflags_args="$LDFLAGS"

# try commenting out....
#cppflags_args="-g -Wall $CPPFLAGS"


export CXXFLAGS=''



case $opsys in
    Linux )
        case $arch in
            x86_64 )
                # gb10: This is no longer required but is left here as an example in case we need 
                # to write similar code in future...
                #set pkgpath=/software/acedb/gtk/lib/pkgconfig
                #set ldflags_args="-Xlinker -rpath -Xlinker /software/acedb/gtk/lib"
                ;;
            ia64 )
                ;;
            * )
                ;;
        esac
        ;;

    Darwin )
        ;;
    
    * )
        ;;
esac


run_dir=`dirname $0`
script_name="configure"
script_exe="$run_dir/$script_name"

echo "Running $script_exe PKG_CONFIG_PATH='$pkgpath' LDFLAGS='$ldflags_args' CPPFLAGS='$cppflags_args' $@"

if [[ $#argv < 1 ]] ; then
  $script_exe PKG_CONFIG_PATH="$pkgpath" LDFLAGS="$ldflags_args" CPPFLAGS="$cppflags_args"  USE_GBTOOLS="yes" 
else
  $script_exe PKG_CONFIG_PATH="$pkgpath" LDFLAGS="$ldflags_args" CPPFLAGS="$cppflags_args" "$@"
fi


exit 0
