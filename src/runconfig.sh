#!/bin/csh
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

set opsys = "`uname -s`"
set arch = "`uname -m`"
set machine = "`uname -s`-`uname -m`"
set pkgpath = ""
set ldflags = ""
set cflags = "-g -Wall"
set configure_args = $*

switch ( $opsys )
  case "Linux":
    switch ( $arch )
      case "x86_64":
      case "ia64":
        # gb10: This is no longer required but is left here as an example in case we need 
        # to write similar code in future...
        #set pkgpath=/software/acedb/gtk/lib/pkgconfig
        #set ldflags="-Xlinker -rpath -Xlinker /software/acedb/gtk/lib"
        breaksw
      default:
        breaksw
    endsw

  case "Darwin":
    breaksw

  default:
    breaksw
endsw


set run_dir = `dirname $0`
set script_name = "configure"
set script_exe = "$run_dir/$script_name"

echo "Running $script_exe $configure_args PKG_CONFIG_PATH='$pkgpath' LDFLAGS='$ldflags' CFLAGS='$cflags'"

if ($#argv < 1 ) then
  $script_exe PKG_CONFIG_PATH="$pkgpath" LDFLAGS="$ldflags" CFLAGS="$cflags"
else
  $script_exe $configure_args PKG_CONFIG_PATH="$pkgpath" LDFLAGS="$ldflags" CFLAGS="$cflags"
endif


