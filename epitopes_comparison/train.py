import os
import pickle

import pandas as pd
import torch
from sklearn.metrics import (
    roc_curve,
    auc,
    precision_recall_curve,
    average_precision_score,
)
from torch import nn
from torch.utils.tensorboard import SummaryWriter
from transformers import Trainer, TrainingArguments
from transformers import EsmTokenizer

from saprot_dataset import SaProtSeqDataset
from saprot_epi_compare_model import SaProtEpitopeCompareNN


class MaskedTrainer(Trainer):
    def __init__(self, *args, **kwargs):
        self.my_logger = kwargs.pop("my_logger")
        super(MaskedTrainer, self).__init__(*args, **kwargs)
        self.loss_fn = nn.BCEWithLogitsLoss()
        self.n = {"train": 0, "val": 0}

    def compute_loss(self, model, inputs, return_outputs=False):
        labels = inputs.pop("label")
        output = model(**inputs)
        # output = torch.flatten(output)
        labels = torch.flatten(labels)
        loss = self.loss_fn(torch.flatten(output), labels)

        # Log loss
        if self.model.training:
            mode = "train"
        else:
            mode = "val"
        self.n[mode] += 1
        self.my_logger.add_scalar(f"Loss/{mode}", loss.item(), self.n[mode])
        return (
            (loss, {"res": torch.flatten(output[..., 0]).unsqueeze(0)})
            if return_outputs
            else loss
        )


def compute_metrics(pred):
    labels = pred.label_ids
    preds = pred.predictions

    # remove padding
    non_padding_mask = labels != -100.0
    labels = labels[non_padding_mask]
    preds = preds[non_padding_mask]

    # calculate metrics
    fpr, tpr, threshold = roc_curve(labels, preds)  # predictions, labels)
    true_roc_auc = auc(fpr, tpr)

    precision, recall, thresholds = precision_recall_curve(labels, preds)
    prc_auc = auc(recall, precision)

    ap = average_precision_score(labels, preds)
    return {
        "true_roc_auc": true_roc_auc,
        "prc_auc": prc_auc,
        "ap": ap,
    }


def train_me(
    dataset_train,
    dataset_test,
    cur_model,
    seed=42,
    log_folder="logs",
    run_name="test",
    output_dir="res",
    n_epochs=20,
    eval_steps=1500,
    lr=10 * 1e-05,
):

    if not os.path.exists(log_folder):
        os.mkdir(log_folder)

    log_folder = os.path.join(log_folder, run_name)
    my_logger = SummaryWriter(log_folder)
    if not os.path.exists(log_folder):
        os.mkdir(log_folder)

    training_args = TrainingArguments(
        output_dir=output_dir,
        num_train_epochs=n_epochs,
        per_device_train_batch_size=1,
        per_device_eval_batch_size=1,
        warmup_steps=0,
        learning_rate=lr,
        logging_steps=100,
        logging_strategy="steps",
        logging_dir=log_folder,
        evaluation_strategy="steps",
        eval_steps=eval_steps,  # 3000,
        save_strategy="steps",
        save_steps=eval_steps,  # 3000,
        save_safetensors=False,
        gradient_accumulation_steps=2,
        fp16=False,
        run_name=run_name,
        seed=seed,
        remove_unused_columns=False,
        metric_for_best_model="eval_loss",
        load_best_model_at_end=True,
        label_names=["label"],
    )
    trainer = MaskedTrainer(
        model=cur_model,
        args=training_args,
        train_dataset=dataset_train,
        eval_dataset=dataset_test,
        data_collator=collator_fn,
        my_logger=my_logger,
        compute_metrics=compute_metrics,
    )
    trainer.train()
    return cur_model, trainer


def collator_fn(x):
    """This collator function currently only works with batch size of 1,
    otherwise it will produce an error."""
    if len(x) == 1:
        return x[0]
    for i in x:
        print(i)
    return x


if __name__ == "__main__":
    # Set up run parameters
    device = torch.device("cuda:0")
    seed = 42

    dataset_type = "small"
    if dataset_type == "big":
        n_epochs = 20
        eval_steps = 1500
    else:
        n_epochs = 10
        eval_steps = 3000

    saprot_version = "SaProt_35M_AF2"
    saprot_name = "saprot"
    if "650" in saprot_version:
        saprot_name = "saprot_650"

    run_name = f"{saprot_name}_{dataset_type}_long"
    output_dir = f"checkpoints/{run_name}"
    trained_model_path = os.path.join(output_dir, "best_model.pt")


    # Load datasets
    datasets_path = "small_dataset/"
    # load antigen structural SaProt sequences
    antigen_structures = pd.read_csv(
        os.path.join(datasets_path, f"{dataset_type}_full_structures.csv"),
        index_col=-1,
    )
    # load epitope structural SaProt sequences
    epitope_structures = pd.read_csv(
        os.path.join(datasets_path, f"{dataset_type}_epitopes_structures.csv"),
        index_col=-1,
    )
    # load dataset with antigen and epitope ids and alignment matrix
    with open(os.path.join(datasets_path, f"{dataset_type}_dataset.pkl"), "rb") as f:
        dataset = pickle.load(f)
    # split into train and test
    test_set_size = int(0.1 * len(dataset))
    train_set = dataset[:-test_set_size]
    test_set = dataset[-test_set_size:]
    # create datasets
    tokenizer = EsmTokenizer.from_pretrained("westlake-repl/SaProt_35M_AF2")
    antigen_epitope_trainset = SaProtSeqDataset(
        train_set,
        epitope_structures,
        antigen_structures,
        tokenizer,
        max_target_size=5000,
    )
    antigen_epitope_testset = SaProtSeqDataset(
        test_set,
        epitope_structures,
        antigen_structures,
        tokenizer,
        max_target_size=5000,
    )
    print(f"Loaded {dataset_type} dataset")


    # Run training
    cur_model = SaProtEpitopeCompareNN(saprot_version=saprot_version).to(device)
    cur_model.main_input_name = "antigen_epitope"

    cur_model, trainer = train_me(
        antigen_epitope_trainset,
        antigen_epitope_testset,
        cur_model,
        run_name=run_name,
        output_dir=output_dir,
        n_epochs=n_epochs,
        eval_steps=eval_steps,
    )
    torch.save(trainer.model, trained_model_path)
