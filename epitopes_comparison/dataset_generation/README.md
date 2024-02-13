Pipeline for Generating a Dataset to Train Neural Networks for Local Structural Similarity Predictions:

01_afdb_prep.py - Preprocesses the database of UniProt structures predicted by AlphaFold. The steps involved are:

Split the full sequence of the protein into domains predicted with high confidence.
Cluster sequences by similarity using MMSeqs.
Align proteins within each cluster to the structure from the cluster's center using PyMOL.
02_structure_based_alignment.py - Generates a file with the sequence alignment derived from structure-based alignment.

03_calculate_sema_scores.py - Precalculates SEMA scores for the proteins.

04_generate_ds.py - Extracts protein regions from the prepared structure in step 1 that have the highest average SEMA scores and combines them with the structure-based alignment to the full-length reference structure of the cluster center of homologous proteins.