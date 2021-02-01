#!/bin/bash

INCS="-I$MRB_INSTALL/sbnana/$SBNANA_VERSION/include/
      -I$MRB_INSTALL/sbnana/$SBNANA_VERSION/include/sbnana/
      -I$MRB_INSTALL/sbnana/$SBNANA_VERSION/include/sbnana/CAFAna
      -I$SRPROXY_INC"

LIBDIR=$MRB_INSTALL/sbnana/$SBNANA_VERSION/${CET_SUBDIR}.${MRB_QUALS/:/.}/lib/
#LIBS="-L$LIBDIR -lCAFAnaCore"

DIRS="Cuts Vars"

find $DIRS -name '*.h' | sed 's:\(.*\):#include "SBNAna/\1":' > classes.h

genreflex classes.h \
    -s classes_def.xml \
    -o SBNAna_dict.cc \
    --noIncludePaths \
    --interpreteronly \
    $INCS

g++ SBNAna_dict.cc `root-config --cflags --libs` $INCS $LIBS --shared -fPIC -o $MRB_BUILDDIR/sbnana/lib/libSBNAna_dict.so

cp $MRB_BUILDDIR/sbnana/lib/libSBNAna_dict.so $LIBDIR
