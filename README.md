# TOPMed CHIP Calling Pipeline

This repository contains a Snakemake workflow for making clonal hematopoiesis of indeterminate potential (CHIP) calls. 

The location of GRCh38 reference and other input files must be configured in [config.yml](config.yml). Mutect2 reference files can be downloaded from gs://gatk-best-practices/somatic-hg38/. A Singularity/Apptainer container image can be built using the definition file [chip.def](chip.def).

