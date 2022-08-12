# fools
FISH tools (fools) to reproducibly generate FISH probes with PaintShop and oligominer


# Installation

    pip install git+https://github.com/adrtsc/fools.git
    
# Specification of metadata file for PaintShop probes

The metadata file contains all the information about the probe color, primer pair to use and the oligo pool the probes are associated with. It has to be submitted as .csv file with at least the following columns:


| gene_name  | color | primer_pair | pool_name |
| ---------- | ----- | ----------- | --------- |
| HSP60      | red   | 1           | HSP60-TP53 |
| TP53        | far-red | 2        | HSP60-TP53 |


* The "gene_name" column contains the gene symbol of the target gene.
* The "color" column can have the values green, red and far-red. It defines which barcode will be appended to the probe and therefore which color the probe will be in the end.
* The "primer_pair" can have the values 1 or 2. The column defines which primer sequence will be appended. The primer sequences are used for the clean-up protocol of the probes. 
* The "pool_name" column defines the name of the pool associated to the probes of that gene in the IDT order form that is generated by the scripts.

# Specification of metadata file for oligominer probes

We use oligominer mainly to generate intronic FISH probes. Therefore the metadata file needs to include all columns as above, but also chromosome coordinates of interest:


| gene_name  | color | primer_pair | pool_name | start | end |
| ---------- | ----- | ----------- | --------- | ----- | --- |
| HSP60      | red   | 1           | HSP60-TP53 | 11020230 | 11023230 |
| TP53        | far-red | 2        | HSP60-TP53 | 919239 | 921239 |


* The "start" column defines the start coordinate of the region of interest
* The "end" column defines the end coordrinate of the region of interest

The chromosome number and the coding strand of the gene do not have to be included in the metadata file, they will be inferred from the annotation file.
    
# Generating a multi-input fasta file for oligominer:

    generate_fasta --genome_path /path/to/genome.fasta --annotation_path /path/to/annotations.gtf --metadata_path /path/to/metadata.csv

# Generating an IDT order file from paintshop output:

    generate_probes --probe_path /path/to/paintshop/result.txt --metadata_path /path/to/metadata.csv
    
# Generating an IDT order file from oligominer output:

    generate_probes --probe_path /path/to/oligominer/result.csv --metadata_path /path/to/metadata.csv --probe_design oligominer
    
