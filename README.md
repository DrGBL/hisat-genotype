# HISAT-genotype-modified

This is NOT the official HISAT-genotype Repository. It is only a slight modification of the HISAT-genotype software to allow for easy use of the latest IMGT-HLA reference.

For a guide to building your own indexes, including how to add DRB3, DRB4, and other genes not on the main GRCh38 scaffold, please see here: https://github.com/DrGBL/hisat_genotype_qgp

## Everything below is for reference, and was copied directely from the HISAT-genotype git

Please see the official website for HISATgenotype:
https://daehwankimlab.github.io/hisat-genotype/

*This directory no longer contains HISAT2. Please see the link below if you are looking for the newest HISAT2 release*
https://daehwankimlab.github.io/hisat2/

## Current Release Version
v1.3.2 - Update to patch critical extract reads error

## Previous Releases
v1.3.1 - Update database locations and handling, and general stability fixes - BROKEN

v1.3.0 - Python 3 version

## Overview
HISAT-genotype is a next-generation genomic analysis software platform capable of assembling and genotyping human genes and genomic regions. Thie software leverages HISAT2s graph FM index and graph alignemnt algorithm to align reads to a specially constructed graph genome. An Expectation-Maximization (EM) algorithm finds the maximum likelihood estimates for each gene allele and a guided de Bruijn graph is used to construct the allele sequences.
