from __future__ import annotations

import torch
from torch import nn


class StartCNN(nn.Module):
    POOL_BINS = 8

    def __init__(self, window_size: int = 121, hidden_channels: int = 128) -> None:
        super().__init__()
        self.window_size = window_size
        center = window_size // 2
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
