################################################################################
#
#   MRC FGU Computational Genomics Group
#
#   $Id: pipeline_snps.py 2870 2010-03-03 10:20:29Z andreas $
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#################################################################################
"""
====================================
Pipeline pipeline_immunogenicity.py
====================================

:Author: Nick Ilott
:Release: $Id$
:Date: |today|
:Tags: Python


Overview
========

A pipeline to predict the immunogenicity i.e. MHC binding potential of a set of 
peptides. It uses tools from the IEDB epitope prediction resource http://www.iedb.org/.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file. 

Input
-----

Input is a fasta file of peptide sequences. The fasta file must end in .fasta or .fasta.gz


Pipeline output
===============

The output is a table of peptides and the associated predicted IC50 values for
MHC binding. The IC50 is the minimum IC50 for the protein i.e. the 15mer that
has the strongest affinity.

Glossary
========

.. glossary::


Code
====

"""
from ruffus import *

import sys
import os
import glob
import gzip
import itertools
import collections
import time
import optparse
import shutil
import sqlite3
from rpy2.robjects import r as R
import CGAT.Experiment as E
import CGAT.IOTools as IOTools
import CGAT.Database as Database
import CGAT.FastaIterator as FastaIterator
import numpy as np

###################################################
###################################################
###################################################
## Pipeline configuration
###################################################

# load options from the config file

import CGATPipelines.Pipeline as P

P.getParameters( 
    [ "pipeline.ini" ] )

PARAMS = P.PARAMS

@follows(mkdir("fasta.dir"))
@split("*.fa.gz", "fasta.dir/*.fasta")
def splitFasta(infile, outfiles):
    '''
    split fasta file into separate files
    '''
    for fasta in FastaIterator.iterate(IOTools.openFile(infile[0])):
        filename = fasta.title.replace(" ", "_").split("_")[0]
        outf = IOTools.openFile(os.path.join(
            "fasta.dir", filename + ".fasta"), "w")
        outf.write(">%s\n%s\n" % (fasta.title, fasta.sequence))
    outf.close()

###################################################
###################################################
###################################################

@follows(mkdir("prediction.dir"))
@transform(splitFasta, regex("(\S+)/(\S+).fasta"), r"prediction.dir/\2.pred")
def predictEpitope(infile, outfile):
    '''
    make a prediction of IC50 of binding to MHC based on 
    IEDB resources
    '''
    allele=PARAMS.get("prediction_allele")
    method=PARAMS.get("prediction_method")
    seq_name=P.snip(os.path.basename(infile), ".fasta")
    scriptsdir=PARAMS.get("scriptsdir")

    # This is very strange but mhc_II_binding.py does not
    # run using P.run(). Not sure why...very problematic for
    # parallelisation
    statement="""mhc_II_binding.py \
                 %(method)s \
                 %(allele)s \
                 %(infile)s \
                 | python %(scriptsdir)s/predictions2minimumic50.py \
                 --seq-name=%(seq_name)s \
                 --log=%(outfile)s.log \
                 > %(outfile)s""" % locals()
    os.system(statement)

###################################################
###################################################
###################################################

@merge(predictEpitope, "prediction.dir/agg_epitopes.tsv.gz")
def aggregatePredictions(infiles, outfile):
    '''
    concatenate the individual predictions
    '''
    infs = " ".join(infiles)
    statement = '''awk 'FNR==1 && NR !=1 { while (/^allele/) getline; } 1 {print}' prediction.dir/*.pred | gzip > %(outfile)s'''
    P.run()

###################################################
###################################################
###################################################

@follows(aggregatePredictions)
def full():
    pass


if __name__== "__main__":
    sys.exit( P.main(sys.argv) )    



