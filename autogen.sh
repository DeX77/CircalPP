#!/bin/sh
aclocal -I m4 --install
libtoolize --automake --force --copy
automake -a -c
autoconf