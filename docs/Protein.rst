Protein
=======

This module calculates the abundance of each of the 20 amino acids who compones proteins. 
The main input files by user are the **GenBank annotation file** as well as the **proteomic data**.

.. note:: Quantitative proteomic is hard to obtain. Using transcriptomic data assuming a 1:1 RNA abundance to protein may provide a working estimate for the cell's amino acid composition.


The total protein weight percentage in the cell is provided as a fraction (number between 0 and 1).
The GenBank files and proteomic should be compatible together. This means that the Gene IDs provided by the proteomic file should be accessible in the GenBank file 
Protein uses CDS elements only. The **translation** into amino acid should be available as well as the **protein_id**. 
BOFdat supports multiple handle GenBank files as for eukaryotes. The use of **protein_id** in the input file is mandatory. **protein_id** are chosen instead of locus tag 
because they are unique and present in all GenBank files.
Hence the transcriptomic file should include the GenBank locus tags as gene identifiers in the first column and the relative abundances in the second.  


