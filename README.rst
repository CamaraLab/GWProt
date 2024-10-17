
# GWProt - Gromov-Wasserstein for Protein Structual Alignment <a href='https://github.com/CamaraLab/GWProt'><img src="docs/images/logo.png" align="right" width="24%"/></a>

GWProt is a Python library for structural alignment of proteins using the Gromov-Wasserstein distance.
Currently it requires PDB (Protein Data Bank) files of proteins.

This project is under active development. Documentation and code will be updated continuously


Installation
=======================================

``pip install GWProt@git+https://github.com/CamaraLab/GWProt``

GWProt uses Pymol 3 for visualization and the fasta36 package for sequence alignment.
These must be installed and added to the PATH to be used.
However all other functionality of GWProt is unaffected if they are not installed.


They can be downloaded from:

``https://pymol.org/``

``https://fasta.bioch.virginia.edu/fasta_www2/fasta_down.shtml``

Documentation
=======================================

Full documentation and tutorials can be found in the the docs/build/html directory

Docker
=================================

A minimal Docker image can be built using the Dockerfile.

It can also be pulled from ``https://hub.docker.com/repository/docker/camaralaboratory/gwprot/general``