RNA
===

This module calculates the abundance of each of the 4 RNA bases. 
The main input files by user are the GenBank annotation file as well as the transcriptomic data.
The total RNA weight percentage in the cell is provided as a fraction (number between 0 and 1).
The GenBank files and transcriptomic should be compatible together. This means that the Gene IDs provided by the transcriptomic file should be accessible in the GenBank file 
RNA uses 3 different types of elements in the annotation: CDS, tRNA and rRNA. The location and strand of these element should be provided as well as the locus tag. 
Hence the transcriptomic file should include the GenBank locus tags as gene identifiers in the first column and the relative abundances in the second.  

File example
------------

.. image:: /rna_file.png
