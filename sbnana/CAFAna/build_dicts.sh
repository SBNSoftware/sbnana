#!/bin/bash

INCS="-I$MRB_INSTALL/sbnana/$SBNCODE_VERSION/include/
      -I$MRB_INSTALL/sbnana/$SBNCODE_VERSION/include/sbnana/
      -I$MRB_INSTALL/sbnana/$SBNCODE_VERSION/include/sbnana/CAFAna
      -I$SRPROXY_INC"

LIBDIR=$MRB_INSTALL/sbnana/$SBNCODE_VERSION/${CET_SUBDIR}.${MRB_QUALS/:/.}/lib/
#LIBS="-L$LIBDIR -lCAFAnaCore"

DIRS="Analysis Core Cuts Experiment Extrap Prediction Systs Vars"

find $DIRS -name '*.h' | grep -v CalcsNuFit_cdr.h | grep -v Prolog.h | grep -v EpilogFwd.h | sed 's:\(.*\):#include "CAFAna/\1":' > classes.h

genreflex classes.h \
    -s classes_def.xml \
    -o CAFAna_dict.cc \
    --noIncludePaths \
    --interpreteronly \
    $INCS

g++ CAFAna_dict.cc `root-config --cflags --libs` $INCS $LIBS --shared -fPIC -o $MRB_BUILDDIR/sbnana/lib/libCAFAna_dict.so

cp $MRB_BUILDDIR/sbnana/lib/libCAFAna_dict.so $LIBDIR
