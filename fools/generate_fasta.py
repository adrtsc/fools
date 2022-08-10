from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd
from pathlib import Path
from functools import reduce
import argparse


def gene_info(x):
    # Extract gene names
    g_name = list(filter(lambda x: 'gene_name' in x,
                         x.split("; ")))[0].split(" ")[1].replace('"', '')
    g_type = list(filter(lambda x: 'gene_type' in x,
                         x.split("; ")))[0].split(" ")[1].replace('"', '')
    return g_name, g_type


def get_sequence(x, record_dict):
    sequence = record_dict[x['seqname']][x['start']:x['end']]
    if x['strand'] == '-':
        sequence = sequence.reverse_complement()
    return str(sequence.seq)


def add_sequence(series):
    return reduce(lambda x, y: x + 50*'N' + y, series)


def main():
    # Parse Arguments
    parser = argparse.ArgumentParser(
        description="Script to construct fasta files for oligominer")
    parser.add_argument('--genome_path', type=str,
                        help='Full path to the genome .fasta file')
    parser.add_argument('--annotation_path', type=str,
                        help='Full path to the annotations .gtf file.')
    parser.add_argument('--metadata_path', type=str,
                        help='Full path to the metadata .csv file with the'
                             'following columns:'
                             '"gene_name": the gene symbol of the gene of interest'
                             '"chromosome": the chromosome on which the gene is located'
                             '"start": start coordinate of the sequence to fetch'
                             '"end" end coordinate of the sequence to fetch')


    args = parser.parse_args()

    # read the metadata file:
    df = pd.read_csv(Path(args.metadata_path))

    # read the genome sequences
    print("reading genome")
    record_dict = SeqIO.to_dict(SeqIO.parse(args.genome_path, "fasta"))
    # read the annotations
    print("reading annotations")
    annotations = pd.read_table(args.annotation_path,
                                comment="#", sep="\t",
                                names=['seqname', 'source', 'feature',
                                       'start', 'end', 'score', 'strand',
                                       'frame', 'attribute'])
    # filter genes
    gencode_genes = annotations[
        annotations.feature == "gene"].reset_index().drop('index', axis=1)
    # add attributes to df
    gencode_genes["gene_name"], gencode_genes["gene_type"] = zip(
        *gencode_genes.attribute.apply(lambda x: gene_info(x)))

    # drop non-unique genes:
    gencode_genes = gencode_genes.drop_duplicates(subset='gene_name', keep=False)

    df = pd.merge(df, gencode_genes[['gene_name', "seqname", "strand"]],
                  on="gene_name")


    # look up the sequence at the coordinates
    df['sequence'] = df.apply(lambda x: get_sequence(x, record_dict), axis=1)

    # aggregate in case multiple sequences for one gene were fetched:
    agg_map = {'gene_name': 'first', 'chromosome': 'first', 'start': 'first',
               'end': 'first', 'seqname': 'first', 'strand': 'first',
               'sequence': add_sequence}

    df = df.groupby('gene_name', as_index=False).aggregate(agg_map)

    # oligominer has a bug where it does not look for probes for the last
    # fasta sequence, add random sequence that will be ignored
    df = pd.concat([df, pd.DataFrame({i:'filler' for i in df.columns}, index=[0])])

    # save fasta file containing all the sequences
    records = (SeqRecord(Seq(seq),
                         id=df['gene_name'].iloc[index],
                         name=df['gene_name'].iloc[index],
                         description=f"{df['seqname'].iloc[index]}: {df['start'].iloc[index]}-{df['end'].iloc[index]}")
               for index, seq in enumerate(list(df["sequence"])))

    out_path = Path(args.metadata_path).with_suffix(".fasta")
    with open(out_path, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")


if __name__ == "__main__":
    main()
