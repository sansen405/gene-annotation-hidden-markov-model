from __future__ import annotations

import torch
from torch import nn


class SpliceCNN(nn.Module):
    """Small 1D CNN for donor/acceptor splice-site scoring."""

    def __init__(self, window_size: int = 81, hidden_channels: int = 64) -> None:
        super().__init__()
        self.window_size = window_size

        self.features = nn.Sequential(
            nn.Conv1d(in_channels=4, out_channels=hidden_channels, kernel_size=7, padding=3),
            nn.ReLU(),
            nn.Conv1d(in_channels=hidden_channels, out_channels=hidden_channels, kernel_size=5, padding=2),
            nn.ReLU(),
            nn.AdaptiveMaxPool1d(1),
        )
        self.classifier = nn.Linear(hidden_channels, 2)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        features = self.features(x).squeeze(-1)
        return self.classifier(features)


def one_hot_encode_windows(windows: list[str]) -> torch.Tensor:
    encoded = torch.zeros((len(windows), 4, len(windows[0])), dtype=torch.float32)
    channel_by_base = {"A": 0, "C": 1, "G": 2, "T": 3}

    for row, window in enumerate(windows):
        for col, base in enumerate(window.upper()):
            channel = channel_by_base.get(base)
            if channel is not None:
                encoded[row, channel, col] = 1.0

    return encoded
