from __future__ import annotations

import torch
from torch import nn


class SpliceCNN(nn.Module):
    POOL_BINS = 8

    def __init__(self, window_size: int = 121, hidden_channels: int = 128) -> None:
        super().__init__()
        self.window_size = window_size
        center = window_size // 2
        self.donor_slice = (max(0, center - 10), min(window_size, center + 21))
        self.acceptor_slice = (max(0, center - 50), min(window_size, center + 11))

        def _backbone() -> nn.Sequential:
            return nn.Sequential(
                nn.Conv1d(4, hidden_channels, kernel_size=7, padding=3),
                nn.BatchNorm1d(hidden_channels),
                nn.ReLU(),
                nn.Dropout(p=0.2),
                nn.Conv1d(hidden_channels, hidden_channels, kernel_size=5, padding=2),
                nn.BatchNorm1d(hidden_channels),
                nn.ReLU(),
                nn.Conv1d(hidden_channels, hidden_channels, kernel_size=3, padding=2, dilation=2),
                nn.BatchNorm1d(hidden_channels),
                nn.ReLU(),
                nn.Conv1d(hidden_channels, hidden_channels, kernel_size=3, padding=4, dilation=4),
                nn.BatchNorm1d(hidden_channels),
                nn.ReLU(),
                nn.Conv1d(hidden_channels, hidden_channels, kernel_size=3, padding=8, dilation=8),
                nn.BatchNorm1d(hidden_channels),
                nn.ReLU(),
                nn.AdaptiveMaxPool1d(SpliceCNN.POOL_BINS),
            )

        self.donor_features = _backbone()
        self.acceptor_features = _backbone()
        head_in = hidden_channels * SpliceCNN.POOL_BINS
        self.donor_head = nn.Sequential(nn.Dropout(p=0.3), nn.Linear(head_in, 1))
        self.acceptor_head = nn.Sequential(nn.Dropout(p=0.3), nn.Linear(head_in, 1))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        d_lo, d_hi = self.donor_slice
        a_lo, a_hi = self.acceptor_slice
        donor_feat = self.donor_features(x[:, :, d_lo:d_hi]).flatten(1)
        acceptor_feat = self.acceptor_features(x[:, :, a_lo:a_hi]).flatten(1)
        donor_logit = self.donor_head(donor_feat)
        acceptor_logit = self.acceptor_head(acceptor_feat)
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
