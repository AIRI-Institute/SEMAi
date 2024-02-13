import os
import pickle

import torch
from torch.utils.data import Dataset
from tqdm import tqdm
from transformers import EsmTokenizer


class SaProtSeqDataset(Dataset):
    def __init__(
        self,
        pairs_dataset,
        epitopes_df,
        antigens_df,
        tokenizer,
        max_target_size=999999,
        max_antigen_size=300,
        max_epitope_size=30,
    ):
        self.epitopes_df = epitopes_df
        self.antigens_df = antigens_df

        # filter out antigen-epitope pairs with large structures
        self.pairs_dataset = []
        for item in tqdm(pairs_dataset):
            epitope_id = item["epitope_id"]
            antigen_id = item["antigen_id"]

            epi_seq = self.epitopes_df.loc[epitope_id]["seq"]
            if len(epi_seq) > max_epitope_size:
                continue

            ag_seq = self.antigens_df.loc[antigen_id]["seq"]
            if len(ag_seq) > max_antigen_size:
                continue

            match_matrix_shape = item["match"].shape
            if match_matrix_shape[0] * match_matrix_shape[1] > max_target_size:
                continue

            self.pairs_dataset.append(item)
        print(
            f"{len(pairs_dataset)} total, {len(self.pairs_dataset)} left after filtering"
        )

        self.tokenizer = tokenizer

    def __getitem__(self, idx):
        item = self.pairs_dataset[idx]
        epitope_id = item["epitope_id"]
        antigen_id = item["antigen_id"]

        epi_saprot_seq = self.epitopes_df.loc[epitope_id]["saprot_seq"]
        ag_saprot_seq = self.antigens_df.loc[antigen_id]["saprot_seq"]

        ag_inputs = self.tokenizer([ag_saprot_seq], return_tensors="pt")
        epi_inputs = self.tokenizer([epi_saprot_seq], return_tensors="pt")

        match_matrix = torch.tensor(item["match"]).to(torch.float32)
        return {
            "antigen_epitope": [ag_inputs, epi_inputs],
            "label": torch.flatten(match_matrix).unsqueeze(0),
        }

    def __len__(self):
        return len(self.pairs_dataset)
