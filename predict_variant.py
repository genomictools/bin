#!/usr/bin/env python3

import argparse
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
import pickle

def main():
    parser = argparse.ArgumentParser(description="Predict variant using AlphaGenome")
    parser.add_argument('--api_key', required=True)
    parser.add_argument('--variant', required=True)
    parser.add_argument('--ontology', required=True)
    parser.add_argument('--assay', required=True)
    parser.add_argument('--sequence_length', required=True)
    parser.add_argument('--reference', required=True)
    parser.add_argument('--alternate', required=True)
    parser.add_argument('--status', required=True)
    args = parser.parse_args()

    # Initialize the DNA client with the provided API key
    model = dna_client.create(args.api_key)

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
    interval = genome.Interval(
        chromosome=chrom,
        start=max(0, int(pos) - int(sequence_length / 2)),
        end=int(pos) + int(sequence_length / 2),
    )

    # Set the ontology and assay terms
    ontology_terms = [args.ontology]
    requested_outputs = [getattr(dna_client.OutputType, args.assay)]
    
    # Always create empty output files first and flush to disk
    with open(args.reference, 'wb') as f:
        f.flush()
    with open(args.alternate, 'wb') as f:
        f.flush()
    with open(args.status, 'w') as f:
        f.flush()

    try:
        outputs = model.predict_variant(
            interval=interval,
            variant=variant,
            ontology_terms=ontology_terms,
            requested_outputs=requested_outputs,
        )
        with open(args.reference, 'wb') as f:
            pickle.dump(outputs.reference, f)

        with open(args.alternate, 'wb') as f:
            pickle.dump(outputs.alternate, f)

    except Exception as e:
        import sys
        if hasattr(e, 'details') and 'Unsupported' in str(e.details()):
            msg = f"WARNING: {e.details()}"
            with open(args.status, 'w') as f:
                f.write(msg + '\n')
            sys.exit(0)
        else:
            msg = f"ERROR: {str(e)}"
            with open(args.status, 'w') as f:
                f.write(msg + '\n')
            sys.exit(1)

if __name__ == "__main__":
    main()