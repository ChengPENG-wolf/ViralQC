#!/usr/bin/env python3 -u
# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.

import argparse
import pathlib
from tqdm import tqdm
from torch import nn
import torch
import GPUtil

from esm import Alphabet, FastaBatchedDataset, ProteinBertModel, pretrained, MSATransformer


def create_parser():
    parser = argparse.ArgumentParser(
        description="Extract per-token representations and model outputs for sequences in a FASTA file"  # noqa
    )

    parser.add_argument(
        "model_location",
        type=str,
        help="PyTorch model file OR name of pretrained model to download (see README for models)",
    )
    parser.add_argument(
        "fasta_file",
        type=pathlib.Path,
        help="FASTA file on which to extract representations",
    )
    parser.add_argument(
        "output_dir",
        type=pathlib.Path,
        help="output directory for extracted representations",
    )

    parser.add_argument("--toks_per_batch", type=int, default=4096, help="maximum batch size")
    parser.add_argument(
        "--repr_layers",
        type=int,
        default=[-1],
        nargs="+",
        help="layers indices from which to extract representations (0 to num_layers, inclusive)",
    )
    parser.add_argument(
        "--include",
        type=str,
        nargs="+",
        choices=["mean", "per_tok", "bos", "contacts"],
        help="specify which representations to return",
        required=True,
    )
    parser.add_argument(
        "--truncation_seq_length",
        type=int,
        default=1022,
        help="truncate sequences longer than the given value",
    )

    parser.add_argument("--nogpu", action="store_true", help="Do not use GPU even if available")
    return parser


def run(args):
    model, alphabet = pretrained.load_model_and_alphabet(args.model_location)
    model.eval()
    print(model)

    if isinstance(model, MSATransformer):
        raise ValueError(
            "This script currently does not handle models with MSA input (MSA Transformer)."
        )

    if torch.cuda.is_available() and not args.nogpu:
        model = model.cuda()
        print("Transferred model to GPU")
        if torch.cuda.device_count() > 1:
            print(f'\nRunning on {torch.cuda.device_count()} GPUs.')
            model = nn.DataParallel(model)
    
    dataset = FastaBatchedDataset.from_file(args.fasta_file)
    batches = dataset.get_batch_indices(args.toks_per_batch, extra_toks_per_seq=1)
    data_loader = torch.utils.data.DataLoader(
        dataset, collate_fn=alphabet.get_batch_converter(args.truncation_seq_length), batch_sampler=batches
    )
    print(f"Read {args.fasta_file} with {len(dataset)} sequences")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    return_contacts = "contacts" in args.include

    #assert all(-(model.module.num_layers + 1) <= i <= model.module.num_layers for i in args.repr_layers)
    #repr_layers = [(i + model.module.num_layers + 1) % (model.module.num_layers + 1) for i in args.repr_layers]
    assert all(-(model.num_layers + 1) <= i <= model.num_layers for i in args.repr_layers)
    repr_layers = [(i + model.num_layers + 1) % (model.num_layers + 1) for i in args.repr_layers]

    result = []
    args.output_file = args.output_dir / f"embed.pt"
    args.output_file.parent.mkdir(parents=True, exist_ok=True)
    
    with torch.no_grad():
        pbar = tqdm(data_loader)
        for batch_idx, (labels, strs, toks) in enumerate(pbar):
            #print(f"Processing {batch_idx + 1} of {len(batches)} batches ({toks.size(0)} sequences)")
            #print(torch.cuda.get_device_name(0))
            #GPUtil.showUtilization()
            if torch.cuda.is_available() and not args.nogpu:
                toks = toks.to(device="cuda", non_blocking=True)

            out = model(toks, repr_layers=repr_layers, return_contacts=return_contacts)

            #logits = out["logits"].to(device="cpu")
            representations = {
                layer: t.to(device="cpu") for layer, t in out["representations"].items()
            }
            if return_contacts:
                contacts = out["contacts"].to(device="cpu")

            for i, label in enumerate(labels):
                #args.output_file = args.output_dir / f"{label}.pt"
                #args.output_file.parent.mkdir(parents=True, exist_ok=True)
                #result = {"label": label}
                truncate_len = min(args.truncation_seq_length, len(strs[i]))
                # Call clone on tensors to ensure tensors are not views into a larger representation
                # See https://github.com/pytorch/pytorch/issues/1995
                
                if "mean" in args.include:
                    embed = {
                        layer: t[i, 1 : truncate_len + 1].mean(0).clone()
                        for layer, t in representations.items()
                    }
                    result.append([label, embed[list(embed.keys())[0]]])
                if "per_tok" in args.include:
                    embed = {
                        layer: t[i, 1 : truncate_len + 1].clone()
                        for layer, t in representations.items()
                    }
                    result.append([label, embed])
                    
            if batch_idx % 1 == 0:
                pbar.set_description(f'Processing {batch_idx}')
                #GPUtil.showUtilization()

    torch.save(result, args.output_file)
    print(len(result))


def main():
    parser = create_parser()
    args = parser.parse_args()
    run(args)

if __name__ == "__main__":
    main()
