#!/usr/bin/env python3

import argparse
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
from tqdm import tqdm
import pysam
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Predict variant using AlphaGenome")
    parser.add_argument('--api_key', required=True)
    parser.add_argument('--vcf_file', required=True)
    parser.add_argument('--organism', required=True)
    parser.add_argument('--sequence_length', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    # Load the model.
    dna_model = dna_client.create(args.api_key)

    # Open VCF file using pysam
    vcf_in = pysam.VariantFile(args.vcf_file)

    # Parse organism specification.
    organism_map = {
        'human': dna_client.Organism.HOMO_SAPIENS,
        'mouse': dna_client.Organism.MUS_MUSCULUS,
    }
    organism = organism_map[args.organism.lower()]

    # Parse sequence length
    sequence_length = args.sequence_length
    sequence_length = dna_client.SUPPORTED_SEQUENCE_LENGTHS[
        f'SEQUENCE_LENGTH_{sequence_length}'
    ]

    # Specify which scorers to use to score your variants:
    scorer_selections = {
        'rna_seq': True,
        'cage': True,
        'procap': True,
        'atac': True,
        'dnase': True,
        'chip_histone': True,
        'chip_tf': True,
        'polyadenylation': True,
        'splice_sites': True,
        'splice_site_usage': True,
        'splice_junctions': True,
    }

    all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    selected_scorers = [
        all_scorers[key]
        for key in all_scorers
        if scorer_selections.get(key.lower(), False)
    ]

    # Remove any scorers or output types that are not supported for the chosen organism.
    selected_scorers = [
        scorer for scorer in selected_scorers
        if (
            organism.value in variant_scorers.SUPPORTED_ORGANISMS[scorer.base_variant_scorer]
            and not (
                scorer.requested_output == dna_client.OutputType.PROCAP
                and organism == dna_client.Organism.MUS_MUSCULUS
            )
        )
    ]

    # Collect results for CSV
    results = []

    for record in tqdm(vcf_in.fetch()):
        # Unique variant_id: chrom:pos:ref:alt
        variant_id = f"{record.chrom}:{record.pos}:{record.ref}:{record.alts[0]}"

        variant = genome.Variant(
            chromosome=str(record.chrom),
            position=int(record.pos),
            reference_bases=record.ref,
            alternate_bases=record.alts[0],  # Only first ALT allele
            name=variant_id,
        )
        interval = variant.reference_interval.resize(sequence_length)

        variant_scores = dna_model.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=selected_scorers,
            organism=organism,
        )

        tidy = variant_scorers.tidy_scores([variant_scores])
        results.append(tidy)

    # Concatenate all results and write to CSV
    if results:
        df = pd.concat(results, ignore_index=True)
        df.to_csv(args.output, index=False, sep='\t', na_rep='.')
        print(f"Results saved to {args.output}")
    else:
        print("No variants processed, no output written.")

if __name__ == "__main__":
    main()