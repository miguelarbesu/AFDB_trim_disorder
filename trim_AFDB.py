#!/usr/bin/env python3

from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Select
from tqdm import tqdm
import pandas as pd
import argparse
import warnings
# warnings.simplefilter(action='ignore', category=UserWarning)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_dir",
        type=str,
        required=True,
        help="Path to AlphaFold model folder",
    )
    parser.add_argument(
        "-p",
        "--prediction_file",
        type=str,
        required=True,
        help="File with AlphaFold-disorder predictions for all models in input_dir",
    )
    parser.add_argument(
        "-o",
        "--out_dir",
        type=str,
        required=False,
        help="Output folder for trimmed models. Default is next to input_dir",
    )
    parser.add_argument(
        "-s",
        "--segment_lenght",
        type=int,
        required=False,
        nargs="?",
        default=25,
        help="Consecutive disordered amino acids to consider a segment disordered",
    )
    parser.add_argument(
        "-t",
        "--disorder_threshold",
        type=float,
        required=False,
        nargs="?",
        default=0.581,
        help="Minimum score to consider disorder",
    )
    return parser.parse_args()


def load_predictions(prediction_file):
    predictions = pd.read_csv(prediction_file, sep="\t")
    grouped_predictions = predictions.groupby("name")
    return grouped_predictions


def get_structures(input_dir):
    structures = Path(input_dir).rglob("*.pdb")
    return structures


def get_ordered_positions(disorder_prediction, segment_length, disorder_threshold):
    # forward and reverse rolling window
    wdw_fdw = disorder_prediction["disorder-25"].rolling(segment_length)
    wdw_rev = disorder_prediction["disorder-25"][::-1].rolling(segment_length)
    ordered_positions = []

    for wdw in [wdw_fdw, wdw_rev]:
        mask = (wdw.min() < disorder_threshold).reindex_like(disorder_prediction)
        ordered = disorder_prediction[mask]
        ordered_positions.extend(ordered["pos"].values)
    ordered_positions = set(ordered_positions)

    return ordered_positions


class ResSelect(Select):
    def accept_residue(self, res):
        if res.id[1] in ordered_positions:
            return True
        else:
            return False


def trim_save_pdb(structure, pdb_id, out_dir):
    parser = PDBParser()
    io = PDBIO()
    s = parser.get_structure("tmp", structure)
    io.set_structure(s)
    # Select residues in ordered_list and write to new PDB file.
    # A Select subclass is defined for this, see below
    io.save(
        f"{out_dir}/{pdb_id}_trimmed.pdb", ResSelect(), preserve_atom_numbering=True
    )
    return


if __name__ == "__main__":
    args = parse_args()

    if not args.out_dir:
        args.out_dir = str(Path(args.input_dir)) + "_trimmed"
    Path.mkdir(Path(args.out_dir), exist_ok=True)

    grouped_predictions = load_predictions(args.prediction_file)
    structures = get_structures(args.input_dir)

    for structure in tqdm(structures):
        pdb_id = Path(structure).stem
        prediction = grouped_predictions.get_group(pdb_id)
        # If the whole protein in disordered do nothing (it would print an empty PDB)
        if (prediction["disorder-25"] >= args.disorder_threshold).all():
            pass
        else: 
            ordered_positions = get_ordered_positions(
                prediction, args.segment_lenght, args.disorder_threshold
            )
            trim_save_pdb(structure, pdb_id, args.out_dir)
