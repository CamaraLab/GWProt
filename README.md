![License](https://img.shields.io/github/license/CamaraLab/GWProt)

# GWProt - Gromov-Wasserstein for Protein Structural Alignment <a href='https://github.com/CamaraLab/GWProt'><img src="docs/GWProt_logo.png" align="right" width="25%"/></a>

GWProt is a Python library for structural alignment of proteins using the Gromov-Wasserstein distance.

Installation
=======================================

``pip install git+https://github.com/CamaraLab/GWProt``

GWProt uses [Pymol 3](https://pymol.org/) for visualization and the [fasta36](https://fasta.bioch.virginia.edu/fasta_www2/fasta_down.shtml) package for sequence alignment. These must be installed and added to the ``PATH`` environment variable if you wish to use visualization and sequence-alignment features. However, all other GWProt functionality is unaffected if they are not installed.

Documentation
=======================================

Full documentation and tutorials can be found at [gwprot.readthedocs.io](https://gwprot.readthedocs.io/)

Docker
=================================

A minimal Docker image is available at [DockerHub](https://hub.docker.com/repository/docker/camaralaboratory/gwprot/general)
