Wechangethenamelater
=====================

  **Wechangethenamelater** is intended hold a makefile which will automatize a sample Bioinformatics research pipeline for recovering shared transcrips. This work
  is based on the shared transcript recovery algorithm described in *"Heuristic pairwise alignment for shared transcript discovery"* published in PubMed.


  PS. As suggested by the name we are hoping to change the name of this repo later.

# Dependencies:
  + [Biopython](https://github.com/biopython/biopython)
    + Quick install:
      ```python
      pip install numpy
      pip install biopython
      ```

## Components:
  + Trimmomatic
  + khmer
    + Diginorm
  + Trinity
  + Sourmash

## Targets:
  + all
  + specs
  + QualityControl
  + Trinity
  + Comparison
  + clean
