#!/bin/bash

cp makefile makefile.old
sed 's/-fPIC//g' makefile.old > makefile.tmp
sed 's/CFLAGS = /CFLAGS = -fPIC /g' makefile.tmp > makefile
rm -rf libcuba.a makefile.tmp makefile.old
make
