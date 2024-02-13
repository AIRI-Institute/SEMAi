import torch
from torch import nn
from transformers import EsmForMaskedLM


class SaProtEpitopeCompareNN(nn.Module):
    def __init__(
        self,
        protein_out_dim=64,
        target_dim=1,
        fc_layer_num=2,
        no_encoder_grad=False,
        saprot_version="SaProt_35M_AF2",
    ):
        super().__init__()

        self.encoder = EsmForMaskedLM.from_pretrained(f"westlake-repl/{saprot_version}")
        self.saprot_emb_dim = 446  # SaProt model embedding size

        if no_encoder_grad:
            self.encoder.eval()
            for param in self.encoder.parameters():
                param.requires_grad = False
        else:
            self.encoder.train()
            for param in self.encoder.parameters():
                param.requires_grad = True

        self.h1 = nn.Linear(self.saprot_emb_dim, protein_out_dim)
        self.batchnorm = nn.LayerNorm(protein_out_dim, elementwise_affine=True)
        self.activation = nn.ReLU()
        self.sigmoid = nn.Sigmoid()
        self.logit = nn.Linear(protein_out_dim, target_dim)

    def forward(self, antigen_epitope, **kwargs):
        ag_inputs = antigen_epitope[0]
        epi_inputs = antigen_epitope[1]

        ag_inputs = {k: v.to(self.encoder.device) for k, v in ag_inputs.items()}
        ag_outputs = self.encoder(**ag_inputs).logits[:, 1:-1, :]

        epi_inputs = {k: v.to(self.encoder.device) for k, v in epi_inputs.items()}
        epi_outputs = self.encoder(**epi_inputs).logits[:, 1:-1, :]

        z0 = ag_outputs.transpose(1, 2)  # (batch_size, emb_size, seq1_len)
        z1 = epi_outputs.transpose(1, 2)  # (batch_size, emb_size, seq2_len)

        z_mul = z0.unsqueeze(3) * z1.unsqueeze(
            2
        )  # (batch_size, emb_size, seq1_len, seq2_len)
        z_cat = torch.transpose(z_mul, 1, 3)
        z_cat = torch.transpose(
            z_cat, 1, 2
        )  # (batch_size, seq1_len, seq2_len, emb_size)
        res = self.h1(z_cat)
        res = self.batchnorm(res)
        res = self.activation(res)
        res = self.logit(res)
        return res
