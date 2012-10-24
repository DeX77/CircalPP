SMALL TUTORIAL FOR CIRCALPP
=============

[![Build Status](https://secure.travis-ci.org/DeX77/CircalPP.png)](http://travis-ci.org/DeX77/CircalPP)

How to get it
-------

The up to date version of CircalPP can always be downloaded in form of a tar.bz2
archive from
http://www.bioinf.uni-leipzig.de. Additionally it is possible to get the
latest state of the art version via GIT by:
git clone git://github.com/DeX77/CircalPP.git

Dependencies
-------

CircalPP depends on three of the BioPP Libraries, obtainable from
http://biopp.univ-montp2.fr at least version 2.0.3 is needed, because CircalPP makes
use of the introduced namespaces. It is recommended to use a C++ compiler pro-
viding the at least the tr1 lib or newer: CircalPP will make use of the new un-
ordered map container introduced by TR1.

Installation
-------

Uncompressed the tar.bz2 archive via tar -xjvf <tar-ball> Or, if you got it
via git run ./autogen.sh in the source directory, to create the necessary ﬁles.
Afterward the standard Unix ./configure && make && make install is
needed to build CircalPP.

First run
-------

When you run the tool without any options it will present you the following:

> usage: CircalPP [options]  
options:  
-v verbose mode  
-b brute-force mode  
-m build multiple alignment using t-coffee  
-s output results stepwise (for very large files)  
-t <integer> random sequences statistics up to size n  
-d <integer> delta value  
-D input is dna  
-R input is rna  
-G alphabet is build from scoring scheme  
-S <filename> file including the scoring scheme  
-I <filename> read data from fasta file instead of stdin  
-F <filename> write all pairwise alignments as fasta to  
-O <filename> write output to <filename> instead of stdout

From these only the used Scoring scheme and the used alphabet type are an obliga-
tion. It is recommended to change the delta value, because it defaults to 1. Input
ﬁles must be in Fasta format. All ﬁles mentioned in this thesis, and some smaller
test examples are included in the testcases directory.
Extending CircalPP
Extending CircalPP to your own needs is easy: all you have to change is in
CircalPP.cpp in the src directory.

Licence
-------

I thereby place CircalPP under the GNU Public licence version 2. Anyone can
use it under these terms.