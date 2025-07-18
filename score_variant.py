#!/usr/bin/env python3

import argparse
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
import pickle

def main():
    parser = argparse.ArgumentParser(description="Predict variant using AlphaGenome")
    parser.add_argument('--api_key', required=True)
    parser.add_argument('--variant', required=True)
    parser.add_argument('--organism', required=True)
    parser.add_argument('--sequence_length', required=True)
    parser.add_argument('--scores', required=True)
    args = parser.parse_args()

    # Initialize the DNA client with the provided API key
    model = dna_client.create(args.api_key)

    organism_map = {
        'human': dna_client.Organism.HOMO_SAPIENS,
        'mouse': dna_client.Organism.MUS_MUSCULUS,
    }
    organism = organism_map[args.organism]

    # Get sequence_length: ["2KB", "16KB", "100KB", "500KB", "1MB"] { type:"string" }
    sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[
        f'SEQUENCE_LENGTH_{args.sequence_length}'
    ]

    # Parse the variant from the input argument
    chrom, pos, ref, alt = args.variant.split(':')
    variant = genome.Variant(
        chromosome=chrom,
        position=int(pos),
        reference_bases=ref,
        alternate_bases=alt,
    )

    # Create an interval around the variant position
    interval = variant.reference_interval.resize(sequence_length)

    variant_scores = model.score_variant(
        interval=interval,
        variant=variant,
        variant_scorers=list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values()),
    )

    df_scores = variant_scorers.tidy_scores(variant_scores)
    df_scores.to_csv(f'{args.scores}', index=False)

if __name__ == "__main__":
    main()