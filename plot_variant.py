#!/usr/bin/env python3

import argparse
from alphagenome.data import gene_annotation, genome, track_data, transcript
from alphagenome.models import dna_client, variant_scorers
from alphagenome.visualization import plot_components
import pickle
import pandas as pd
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(description="Plot variant using AlphaGenome")
    parser.add_argument('--variant', required=True)
    parser.add_argument('--ontology', required=True)
    parser.add_argument('--assay', required=True)
    parser.add_argument('--sequence_length', required=True)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--alternate', required=True)
    parser.add_argument('--gtf', required=True)
    parser.add_argument('--plot', required=True)
    args = parser.parse_args()

    # Parse the variant from the input argument
    chrom, pos, ref, alt = args.variant.split(':')
    variant = genome.Variant(
        chromosome=chrom,
        position=int(pos),
        reference_bases=ref,
        alternate_bases=alt,
    )

    # Load the reference and alternate sequences from the files
    with open(args.reference, 'rb') as f:
        reference = pickle.load(f)
    with open(args.alternate, 'rb') as f:
        alternate = pickle.load(f)
    
    # Plot predictions
    # Extract the attribute specified by args.assay from reference and alternate
    ref_data = getattr(reference, args.assay.lower())
    alt_data = getattr(alternate, args.assay.lower())

    # Get sequence_length: ["2KB", "16KB", "100KB", "500KB", "1MB"] { type:"string" }
    sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[
        f'SEQUENCE_LENGTH_{args.sequence_length}'
    ]

    # Create an interval around the variant position
    interval = genome.Interval(
        chromosome=chrom,
        start=max(0, int(pos) - int(sequence_length / 2)),
        end=int(pos) + int(sequence_length / 2),
    )

    # Load gene annotations (from GENCODE).
    gtf = pd.read_feather(args.gtf)

    # Filter to protein-coding genes and highly supported transcripts.
    gtf_transcript = gene_annotation.filter_transcript_support_level(
        gene_annotation.filter_protein_coding(gtf), ['1']
    )

    # Extractor for identifying transcripts in a region.
    transcript_extractor = transcript.TranscriptExtractor(gtf_transcript)

    # Also define an extractor that fetches only the longest transcript per gene.
    gtf_longest_transcript = gene_annotation.filter_to_longest_transcript(
        gtf_transcript
    )
    longest_transcript_extractor = transcript.TranscriptExtractor(
        gtf_longest_transcript
    )

    # Extract the longest transcripts per gene for this interval.
    longest_transcripts = longest_transcript_extractor.extract(interval)

    plot_components.plot(
        [
            plot_components.TranscriptAnnotation(longest_transcripts),
            plot_components.OverlaidTracks(
                tdata={
                    'REF': ref_data,
                    'ALT': alt_data,
                },
                colors={'REF': 'dimgrey', 'ALT': 'red'},
            ),
        ],
        interval=ref_data.interval.resize(2**15),
        # Annotate the location of the variant as a vertical line.
        annotations=[plot_components.VariantAnnotation([variant], alpha=0.8)],
        title=f'{args.assay} ({args.ontology}): {args.variant} ({args.sequence_length})',
    )

    # Save the plot to the specified output file
    plt.savefig(args.plot)
    plt.close()

if __name__ == "__main__":
    main()