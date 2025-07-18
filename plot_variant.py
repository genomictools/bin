#!/usr/bin/env python3

import argparse
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
from alphagenome.visualization import plot_components
import pickle
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(description="Plot variant using AlphaGenome")
    parser.add_argument('--variant', required=True)
    parser.add_argument('--assay', required=True)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--alternate', required=True)
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

    plot_components.plot(
        [
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
    )

    # Save the plot to the specified output file
    plt.savefig(args.plot)
    plt.close()

if __name__ == "__main__":
    main()