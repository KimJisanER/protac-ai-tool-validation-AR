"""
PROTAC-STAN torch.load monkeypatch (C17 대응).

원본 리포(~/PROTAC-STAN)를 수정하지 않고, torch.load의 기본값
weights_only를 False로 강제한다. torch 2.8에서 weights_only=True가
기본값이 되며 data.py:181 / data_loader.py:47,49,52 / inference.py:69의
torch.load(... torch_geometric Data/slices, state_dict ...)가
UnpicklingError로 실패하는 문제를 해결한다.

사용: STAN 모듈 import 전에 이 모듈을 먼저 import.
    import stan_patch  # noqa
"""
import functools
import torch

_orig_torch_load = torch.load


@functools.wraps(_orig_torch_load)
def _patched_load(*args, **kwargs):
    kwargs.setdefault("weights_only", False)
    return _orig_torch_load(*args, **kwargs)


torch.load = _patched_load
print("[stan_patch] torch.load monkeypatched: weights_only=False (default)")
