from __future__ import annotations

import torch
from torch import nn


class SpliceCNN(nn.Module):
    """1D CNN for splice-site scoring with separate donor and acceptor backbones.

    Each backbone: Conv7 -> Conv5 -> dilated Conv3 (receptive field ~23 bp),
    with BatchNorm + Dropout throughout and a global max-pool at the end.
    Separate heads produce independent donor and acceptor logits.
    Output shape: (batch, 2) — column 0 donor, column 1 acceptor.
    """

    def __init__(self, window_size: int = 121, hidden_channels: int = 128) -> None:
        super().__init__()
        self.window_size = window_size

        def _backbone() -> nn.Sequential:
            return nn.Sequential(
                nn.Conv1d(4, hidden_channels, kernel_size=7, padding=3),
                nn.BatchNorm1d(hidden_channels),
                nn.ReLU(),
                nn.Dropout(p=0.2),
                nn.Conv1d(hidden_channels, hidden_channels, kernel_size=5, padding=2),
                nn.BatchNorm1d(hidden_channels),
                nn.ReLU(),
                # dilation=2 grows receptive field without adding parameters
                nn.Conv1d(hidden_channels, hidden_channels, kernel_size=3, padding=2, dilation=2),
                nn.BatchNorm1d(hidden_channels),
                nn.ReLU(),
                nn.AdaptiveMaxPool1d(1),
            )

        self.donor_features = _backbone()
        self.acceptor_features = _backbone()
        self.donor_head = nn.Sequential(nn.Dropout(p=0.3), nn.Linear(hidden_channels, 1))
        self.acceptor_head = nn.Sequential(nn.Dropout(p=0.3), nn.Linear(hidden_channels, 1))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        donor_logit = self.donor_head(self.donor_features(x).squeeze(-1))
        acceptor_logit = self.acceptor_head(self.acceptor_features(x).squeeze(-1))
        return torch.cat([donor_logit, acceptor_logit], dim=1)


def one_hot_encode_windows(windows: list[str]) -> torch.Tensor:
    encoded = torch.zeros((len(windows), 4, len(windows[0])), dtype=torch.float32)
    channel_by_base = {"A": 0, "C": 1, "G": 2, "T": 3}

    for row, window in enumerate(windows):
        for col, base in enumerate(window.upper()):
            channel = channel_by_base.get(base)
            if channel is not None:
                encoded[row, channel, col] = 1.0

    return encoded
