#
# mexopts.sh	Shell script for configuring MEX-file creation script,
#               mex.  These options were tested with the specified compiler.
#
# usage:        Do not call this file directly; it is sourced by the
#               mex shell script.  Modify only if you don't like the
#               defaults after running mex.  No spaces are allowed
#               around the '=' in the variable assignment.
#
# Note: For the version of system compiler supported with this release,
#       refer to Technical Note 1601 at:
#       http://www.mathworks.com/support/tech-notes/1600/1601.html
#
#
# SELECTION_TAGs occur in template option files and are used by MATLAB
# tools, such as mex and mbuild, to determine the purpose of the contents
# of an option file. These tags are only interpreted when preceded by '#'
# and followed by ':'.
#
#SELECTION_TAG_MEX_OPT: Template Options file for building MEX-files via the system ANSI compiler
#
# Copyright 1984-2007 The MathWorks, Inc.
# $Revision: 1.78.4.15 $  $Date: 2007/11/12 22:52:41 $
#----------------------------------------------------------------------------
#
    TMW_ROOT="$MATLAB"
    MFLAGS=''
    if [ "$ENTRYPOINT" = "mexLibrary" ]; then
        MLIBS="-L$TMW_ROOT/bin/$Arch -lmx -lmex -lmat -lmwservices -lut"
    else  
        MLIBS="-L$TMW_ROOT/bin/$Arch -lmx -lmex -lmat"
    fi
    case "$Arch" in
        Undetermined)
#----------------------------------------------------------------------------
# Change this line if you need to specify the location of the MATLAB
# root directory.  The script needs to know where to find utility
# routines so that it can determine the architecture; therefore, this
# assignment needs to be done while the architecture is still
# undetermined.
#----------------------------------------------------------------------------
            MATLAB="$MATLAB"
            ;;
        glnx86)
#----------------------------------------------------------------------------
            RPATH="-Wl,-rpath-link,$TMW_ROOT/bin/$Arch"
            # StorageVersion: 1.0
            # CkeyName: GNU C
            # CkeyManufacturer: GNU
            # CkeyLanguage: C
            # CkeyVersion:
            CC='gcc'
            CFLAGS='-ansi -D_GNU_SOURCE'
            CFLAGS="$CFLAGS -fPIC -pthread -m32"
            CFLAGS="$CFLAGS  -fexceptions"
            CFLAGS="$CFLAGS -D_FILE_OFFSET_BITS=64" 
            CLIBS="$RPATH $MLIBS -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
            CLIBS="$CLIBS -lstdc++"
#           
            # C++keyName: GNU C++
            # C++keyManufacturer: GNU
            # C++keyLanguage: C++
            # C++keyVersion: 
            CXX='g++'
            CXXFLAGS='-ansi -D_GNU_SOURCE'
            CXXFLAGS="$CXXFLAGS -D_FILE_OFFSET_BITS=64" 
            CXXFLAGS="$CXXFLAGS -fPIC -pthread"
            CXXLIBS="$RPATH $MLIBS -lm"
            CXXOPTIMFLAGS='-O -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
#
            # FortrankeyName: g95
            # FortrankeyManufacturer: GNU
            # FortrankeyLanguage: Fortran
            # FortrankeyVersion: 
            FC='g95'
            FFLAGS='-fexceptions'
            FFLAGS="$FFLAGS -fPIC"
            FLIBS="$RPATH $MLIBS -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDEXTENSION='.mexglx'
            LDFLAGS="-pthread -shared -m32 -Wl,--version-script,$TMW_ROOT/extern/lib/$Arch/$MAPFILE -Wl,--no-undefined"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        glnxa64)
#----------------------------------------------------------------------------
            RPATH="-Wl,-rpath-link,$TMW_ROOT/bin/$Arch"
            # StorageVersion: 1.0
            # CkeyName: GNU C
            # CkeyManufacturer: GNU
            # CkeyLanguage: C
            # CkeyVersion:
            CC='gcc'
            CFLAGS='-ansi -D_GNU_SOURCE'
            CFLAGS="$CFLAGS  -fexceptions"
            CFLAGS="$CFLAGS -fPIC -fno-omit-frame-pointer -pthread"
            CLIBS="$RPATH $MLIBS -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
            CLIBS="$CLIBS -lstdc++"
#
            # C++keyName: GNU C++
            # C++keyManufacturer: GNU
            # C++keyLanguage: C++
            # C++keyVersion: 
            CXX='g++'
            CXXFLAGS='-ansi -D_GNU_SOURCE'
            CXXFLAGS="$CXXFLAGS -fPIC -fno-omit-frame-pointer -pthread"
            CXXLIBS="$RPATH $MLIBS -lm"
            CXXOPTIMFLAGS='-O -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
            # FortrankeyName: g95
            # FortrankeyManufacturer: GNU
            # FortrankeyLanguage: Fortran
            # FortrankeyVersion: 
#
            FC='g95'
            FFLAGS='-fexceptions'
            FFLAGS="$FFLAGS -fPIC -fno-omit-frame-pointer"
            FLIBS="$RPATH $MLIBS -lm"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$COMPILER"
            LDEXTENSION='.mexa64'
            LDFLAGS="-pthread -shared -Wl,--version-script,$TMW_ROOT/extern/lib/$Arch/$MAPFILE -Wl,--no-undefined"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        sol64)
#----------------------------------------------------------------------------
            # StorageVersion: 1.0
            # CkeyName: Sun Studio
            # CkeyManufacturer: Sun
            # CkeyLanguage: C
            # CkeyVersion:
            CC='cc -xarch=v9a'
            CFLAGS='-dalign -xlibmieee -D__EXTENSIONS__ -D_POSIX_C_SOURCE=199506L -mt'
            CFLAGS="$CFLAGS -KPIC"
            CLIBS="$MLIBS -lm"
            CLIBS="$CLIBS -lc"
            COPTIMFLAGS='-xO3 -xlibmil -DNDEBUG'
            CDEBUGFLAGS='-xs -g'
#           
            # C++keyName: Sun Studio
            # C++keyManufacturer: Sun
            # C++keyLanguage: C++
            # C++keyVersion:
            CXX='CC -xarch=v9a -compat=5'
            CCV=`CC -xarch=v9a -V 2>&1`
            version=`expr "$CCV" : '.*\([0-9][0-9]*\)\.'`
            if [ "$version" = "4" ]; then
                    echo "SC5.0 or later C++ compiler is required"
            fi
            CXXFLAGS='-dalign -xlibmieee -D__EXTENSIONS__ -library=stlport4,Crun'
            CXXFLAGS="$CXXFLAGS -D_POSIX_C_SOURCE=199506L -mt"
            CXXFLAGS="$CXXFLAGS -KPIC -norunpath"
            CXXLIBS="$MLIBS -lm"
            CXXOPTIMFLAGS='-xO3 -xlibmil -DNDEBUG'
            CXXDEBUGFLAGS='-xs -g'
#
            # FortrankeyName: Sun Studio
            # FortrankeyManufacturer: Sun
            # FortrankeyLanguage: Fortran
            # FortrankeyVersion:
            FC='f90 -xarch=v9a'
            FFLAGS='-dalign'
            FFLAGS="$FFLAGS -KPIC -mt"
            FLIBS="$MLIBS -lfui -lfsu -lsunmath -lm -lc"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-xs -g'
#
            LD="$COMPILER"
            LDEXTENSION='.mexs64'
            LDFLAGS="-G -mt -M$TMW_ROOT/extern/lib/$Arch/$MAPFILE"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-xs -g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        mac)
#----------------------------------------------------------------------------
            # StorageVersion: 1.0
            # CkeyName: GNU C
            # CkeyManufacturer: GNU
            # CkeyLanguage: C
            # CkeyVersion:
            CC='gcc-4.0'
            CFLAGS='-fno-common -no-cpp-precomp'
            CFLAGS="$CFLAGS"
            CLIBS="$MLIBS"
            COPTIMFLAGS='-O3 -fno-loop-optimize -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            CLIBS="$CLIBS -lstdc++"
            # C++keyName: GNU C++
            # C++keyManufacturer: GNU
            # C++keyLanguage: C++
            # C++keyVersion: 
            CXX=g++-4.0
            CXXFLAGS='-fno-common -no-cpp-precomp -fexceptions -arch ppc'
            CXXLIBS="$MLIBS -lstdc++"
            CXXOPTIMFLAGS='-O3 -fno-loop-optimize -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
            # FortrankeyName: g95
            # FortrankeyManufacturer: GNU
            # FortrankeyLanguage: Fortran
            # FortrankeyVersion: 
            FC='g95'
            FFLAGS='-fexceptions'
            FC_LIBDIR=`$FC -print-file-name=libf95.a 2>&1 | sed -n '1s/\/*libf95\.a//p'`
            FLIBS="$MLIBS -L$FC_LIBDIR -lf95"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$CC"
            LDEXTENSION='.mexmac'
            LDFLAGS='-Wl,-flat_namespace -undefined suppress'
            LDFLAGS="$LDFLAGS -bundle -Wl,-exported_symbols_list,$TMW_ROOT/extern/lib/$Arch/$MAPFILE"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        maci)
#----------------------------------------------------------------------------
            # StorageVersion: 1.0
            # CkeyName: GNU C
            # CkeyManufacturer: GNU
            # CkeyLanguage: C
            # CkeyVersion:
            CC='gcc-4.0'
            CFLAGS='-fno-common -no-cpp-precomp'
            CFLAGS="$CFLAGS  -fexceptions"
            CLIBS="$MLIBS"
            COPTIMFLAGS='-O3 -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            CLIBS="$CLIBS -lstdc++"
            # C++keyName: GNU C++
            # C++keyManufacturer: GNU
            # C++keyLanguage: C++
            # C++keyVersion: 
            CXX=g++-4.0
            CXXFLAGS='-fno-common -no-cpp-precomp -fexceptions -arch i386'
            CXXLIBS="$MLIBS -lstdc++"
            CXXOPTIMFLAGS='-O3 -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
            # FortrankeyName: g95
            # FortrankeyManufacturer: GNU
            # FortrankeyLanguage: Fortran
            # FortrankeyVersion: 
            FC='g95'
            FFLAGS='-fexceptions'
            FC_LIBDIR=`$FC -print-file-name=libf95.a 2>&1 | sed -n '1s/\/*libf95\.a//p'`
            FLIBS="$MLIBS -L$FC_LIBDIR -lf95"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$CC"
            LDEXTENSION='.mexmaci'
            LDFLAGS='-Wl,-flat_namespace -undefined suppress'
            LDFLAGS="$LDFLAGS -bundle -Wl,-exported_symbols_list,$TMW_ROOT/extern/lib/$Arch/$MAPFILE"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
        maci64)
#----------------------------------------------------------------------------
            # StorageVersion: 1.0
            # CkeyName: GNU C
            # CkeyManufacturer: GNU
            # CkeyLanguage: C
            # CkeyVersion:
            CC='gcc-4.0'
            CFLAGS='-fno-common -no-cpp-precomp -fexceptions -arch x86_64'
            CLIBS="$MLIBS -lstdc++"
            COPTIMFLAGS='-O3 -DNDEBUG'
            CDEBUGFLAGS='-g'
#
            # C++keyName: GNU C++
            # C++keyManufacturer: GNU
            # C++keyLanguage: C++
            # C++keyVersion: 
            CXX=g++-4.0
            CXXFLAGS='-fno-common -no-cpp-precomp -fexceptions -arch x86_64'
            CXXLIBS="$MLIBS -lstdc++"
            CXXOPTIMFLAGS='-O3 -DNDEBUG'
            CXXDEBUGFLAGS='-g'
#
            # FortrankeyName: Intel Fortran
            # FortrankeyManufacturer: Intel
            # FortrankeyLanguage: Fortran
            # FortrankeyVersion: 
            FC='ifort'
            FFLAGS=''
            FC_LIBDIR=''
            FLIBS="$MLIBS"
            FOPTIMFLAGS='-O'
            FDEBUGFLAGS='-g'
#
            LD="$CC"
            LDEXTENSION='.mexmaci64'
            LDFLAGS='-Wl,-twolevel_namespace -undefined error -arch x86_64'
            LDFLAGS="$LDFLAGS -bundle -Wl,-exported_symbols_list,$TMW_ROOT/extern/lib/$Arch/$MAPFILE"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'
#
            POSTLINK_CMDS=':'
#----------------------------------------------------------------------------
            ;;
    esac
#############################################################################
#
# Architecture independent lines:
#
#     Set and uncomment any lines which will apply to all architectures.
#
#----------------------------------------------------------------------------
#           CC="$CC"
#CC="/mathworks/hub/sol2/apps/SUNWspro_studio11_20070319/opt/SUNWspro/bin/cc"
#           CFLAGS="$CFLAGS"
#           COPTIMFLAGS="$COPTIMFLAGS"
#           CDEBUGFLAGS="$CDEBUGFLAGS"
#           CLIBS="$CLIBS"
#
### vanilla netcdf-3 settings, maci
### NetCDF libraries configured with "--disable-shared"
#          NETCDF="/opt/netcdf-3.6.2"
#          CFLAGS="$CFLAGS -g -I${NETCDF}/include  "
#          CLIBS="-L${NETCDF}/lib -lnetcdf $CLIBS "
#
### vanilla netcdf-3 settings, glnxa64
### NetCDF library configured with "--disable-shared", CFLAGS="-fPIC"
#          NETCDF="/local/netcdf-3.6.3"
#          NETCDF="/opt/netcdf-4.1.3-i386"
#          CFLAGS="$CFLAGS -g -I${NETCDF}/include "
#          CLIBS="${NETCDF}/lib/libnetcdf.a ${NETCDF}/lib/libhdf5_hl.a ${NETCDF}/lib/libhdf5.a $CLIBS "
#
### vanilla netcdf-4 settings, maci
### NetCDF, HDF5 libraries configured with "--disable-shared"
#          HDF5="/opt/hdf5-1.8.2"
#          NETCDF="/opt/hdf5-1.8.2"
#          CFLAGS="$CFLAGS -g -I${NETCDF}/include -I${HDF5}/include "
#          CLIBS="${NETCDF}/lib/libnetcdf.a ${HDF5}/lib/libhdf5_hl.a ${HDF5}/lib/libhdf5.a -lz $CLIBS "
#
### Opendap settings, glnxa64
### Opendap libraries compiled with CXXFLAGS="-fPIC" and "--disable-shared"
#            NETCDF="/local/opendap"
#            CFLAGS="$CFLAGS -I${NETCDF}/include/libnc-dap -I${NETCDF}/include/libdap "
#            CLIBS="-ldl  -lssl -lcrypto -lxml2 -lz  -lm -lpthread $CLIBS"
#            CLIBS="-lcom_err -lresolv  $CLIBS"
#            CLIBS="-lcurl -lgssapi_krb5 -lkrb5 -lk5crypto -lkrb5support $CLIBS"
#            CLIBS="-L${NETCDF}/lib -lnc-dap -ldapclient -ldap $CLIBS"
#
# Opendap settings, maci.  
            NETCDF="/opt/nc4-32"
            CFLAGS="$CFLAGS -I${NETCDF}/include  "
            CLIBS="-L${NETCDF}/lib -lnetcdf -lcurl -lhdf5_hl -lhdf5 $CLIBS"
#
#           FC="$FC"
#           FFLAGS="$FFLAGS"
#           FOPTIMFLAGS="$FOPTIMFLAGS"
#           FDEBUGFLAGS="$FDEBUGFLAGS"
#           FLIBS="$FLIBS"
#
#           LD="$LD"
#           LDFLAGS="$LDFLAGS"
#           LDFLAGS="$LDFLAGS -mt"        # sol64
#           LDOPTIMFLAGS="$LDOPTIMFLAGS"
#           LDDEBUGFLAGS="$LDDEBUGFLAGS"
#----------------------------------------------------------------------------
#############################################################################
