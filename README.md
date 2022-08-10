# fools
FISH tools (fools) to reproducibly generate FISH probes with PaintShop and oligominer


# Installation

    pip install git+https://github.com/adrtsc/fools.git
    
# Generating a multi-input fasta file for oligominer:

    generate_fasta --genome_path /path/to/genome.fasta --annotation_path /path/to/annotations.gtf --metadata_path /path/to/metadata.csv

# Generating an IDT order file from paintshop output:

    generate_probes --probe_path /path/to/paintshop/result.csv --metadata_path /path/to/metadata.csv
    
# Generating an IDT order file from oligominer output:

    generate_probes --probe_path /path/to/oligominer/result.txt --metadata_path /path/to/metadata.csv --probe_design oligominer
    
