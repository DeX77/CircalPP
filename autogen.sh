#!/bin/sh
autoheader
aclocal -I m4 --install
libtoolize --automake --force --copy
automake -a -c
autoconf