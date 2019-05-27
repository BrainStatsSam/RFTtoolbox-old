# RFTtoolbox (Beta Version): A toolbox designed for generation and analysis of random fields both continuously sampled and on a lattice.
NOTE: This is a BETA version of this toolbox. Many more random field functions will soon be added.
And some existing features will be tidied up. Watch this space! Feel free to use the functions available in your research.
However if you do please cite us.

## Table of contents
* [Introduction](#introduction)
* [Folder Structure](#folderstruct)
    * [Cluster_Size_Inference](#CLInf)
    * [Random_Field_Generation](#RFfunctions)
    * [RFT_functions](#siggen)
    * [SPM_Functions](#power)
* [Set Up](#setup)
    * [Dependencies](#dependencies)

## Introduction <a name="introduction"></a>
The RFTtoolbox currently contains code to generate smooth Gaussian, t and 
F fields (with a zero or non-zero mean, the peaks of which can be specified) 
on a lattice of arbitrary size accounting for the edge effect.

It will soon contain code to perform clusterwise inference and analysis and thresholding 
using LKCS and to generate convolution fields as well as other functionalities.

## Folder Structure <a name="folderstruct"></a>

### Cluster Size Inference <a name="CLInf"></a>

A collection of functions to perform cluster size inference using random 
field theory. Many more functions will be added to this folder.

### Random Field Generation <a name="RFfunctions"></a>

Functions to generate Gaussian random fields.

#### asdf

## Set Up
If you have any difficulties getting this code to run or have any questions
feel free to get in touch with me via samuel.davenport(AT)stats.ox.ac.uk.