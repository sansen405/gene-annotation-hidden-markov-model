from __future__ import annotations

import torch
from torch import nn


class StartCNN(nn.Module):
    """1D CNN for translation-start scoring.

    Mirrors the splice CNN's position-aware design (dilated convs + a few
    positional pooling bins instead of one global value) so the head keeps
    coarse information about where context sits relative to the candidate ATG at
    the window center. The crop is asymmetric and downstream-biased, following
    translation-initiation biology:
      * ~20 bp upstream — the Kozak / initiation context (core -6..+4, with the
        -3 purine and +4 G the strongest positions) plus proximal 5' UTR
        composition that separates a real start from intergenic/UTR sequence.
      * the ATG at the window center.
      * ~60 bp (~20 codons) downstream — coding-frame periodicity at the ORF
        onset, the single most discriminative cue between a true start and a
        downstream internal in-frame Met (which the gene-start penalty cannot
        distinguish).
    The crop (80 bp) is the same order as the splice heads (donor 31 bp,
    acceptor 61 bp) and well inside the 121 bp input window.
    Output shape: (batch, 1) — the start logit.
    """

    POOL_BINS = 8

    def __init__(self, window_size: int = 121, hidden_channels: int = 128) -> None:
        super().__init__()
        self.window_size = window_size
        center = window_size // 2
        # ATG occupies [center, center + 3). Upstream Kozak/5' UTR context plus a
        # wider downstream coding-onset window (downstream-biased).
        self.start_slice = (max(0, center - 20), min(window_size, center + 60))

        def _backbone() -> nn.Sequential:
            return nn.Sequential(
                nn.Conv1d(4, hidden_channels, kernel_size=7, padding=3),
                nn.BatchNorm1d(hidden_channels),
                nn.ReLU(),
                nn.Dropout(p=0.2),
                nn.Conv1d(hidden_channels, hidden_channels, kernel_size=5, padding=2),
                nn.BatchNorm1d(hidden_channels),
                nn.ReLU(),
                # Stacked dilated convolutions grow the receptive field
                # exponentially (d=2: 15bp, d=4: 23bp, d=8: 39bp) so the head can
                # see Kozak context and the first codons without extra parameters.
                nn.Conv1d(hidden_channels, hidden_channels, kernel_size=3, padding=2, dilation=2),
                nn.BatchNorm1d(hidden_channels),
                nn.ReLU(),
                nn.Conv1d(hidden_channels, hidden_channels, kernel_size=3, padding=4, dilation=4),
                nn.BatchNorm1d(hidden_channels),
                nn.ReLU(),
                nn.Conv1d(hidden_channels, hidden_channels, kernel_size=3, padding=8, dilation=8),
                nn.BatchNorm1d(hidden_channels),
                nn.ReLU(),
                nn.AdaptiveMaxPool1d(StartCNN.POOL_BINS),
            )

        self.start_features = _backbone()
        head_in = hidden_channels * StartCNN.POOL_BINS
        self.start_head = nn.Sequential(nn.Dropout(p=0.3), nn.Linear(head_in, 1))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        s_lo, s_hi = self.start_slice
        start_feat = self.start_features(x[:, :, s_lo:s_hi]).flatten(1)
        return self.start_head(start_feat)


def one_hot_encode_windows(windows: list[str]) -> torch.Tensor:
    encoded = torch.zeros((len(windows), 4, len(windows[0])), dtype=torch.float32)
    channel_by_base = {"A": 0, "C": 1, "G": 2, "T": 3}

    for row, window in enumerate(windows):
        for col, base in enumerate(window.upper()):
            channel = channel_by_base.get(base)
            if channel is not None:
                encoded[row, channel, col] = 1.0

    return encoded
