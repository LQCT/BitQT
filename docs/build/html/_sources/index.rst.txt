.. bitQT documentation master file, created by
   sphinx-quickstart on Mon Feb  8 05:05:34 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to BitQT's documentation!
=================================
BitQT is a Python command-line interface (CLI) conceived to speed up
the Heyer's Quality Threshold (QT) [1]_ clustering of long Molecular Dynamics.
The package implements a heuristic approach to `this exact variant
of QT <https://doi.org/10.1021/acs.jcim.9b00558>`_.

The construction of a binary-encoded RMSD matrix, instead of the classical (half/single/double)-precision float matrix, led to considerable RAM savings compared to the few existing QT implementations. This binary matrix also allows implementing the significant steps as bitwise operations, which are faster than the corresponding set operations when dealing with considerable amounts of data. 

.. toctree::
   :maxdepth: 2
   :hidden:
   
   description
   installation
   technical
   tutorial
   citation
   changelog
   faq
   
   
.. [1] Heyer, L. J.; Kruglyak, S.; Yooseph, S. Exploring Expression Data Identification and Analysis of Coexpressed Genes. Genome Res. 1999, 9 (11), 1106–1115.
