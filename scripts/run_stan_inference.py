"""
PROTAC-STAN class-1(활성) 확률 덤프 추론 러너.

원본 inference.py는 torch.max(...,dim=1)로 argmax 라벨만 반환 → 랭킹 불가.
이 러너는 동일 모델/로더를 쓰되 F.log_softmax 출력의 class-1 성분을
exp 하여 활성확률을 덤프한다. 원본 리포는 미수정(monkeypatch + 별도 러너).

실행:
  ~/anaconda3/envs/protac/bin/python scripts/run_stan_inference.py \
      --root data/custom --name custom --out outputs/stan_probs.csv

cwd는 ~/PROTAC-STAN 이어야 한다(config.toml, saved_models 상대경로).
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
import toml

# --- torch.load monkeypatch (반드시 STAN 모듈 import 전) ---
_THIS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _THIS)
import stan_patch  # noqa: F401

STAN_DIR = os.path.expanduser("~/PROTAC-STAN")
sys.path.insert(0, STAN_DIR)

from data_loader import PROTACLoader  # noqa: E402
from model import PROTAC_STAN  # noqa: E402


def setup_seed(seed=21332):
    import random
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=str, default="data/custom")
    parser.add_argument("--name", type=str, default="custom")
    parser.add_argument("--out", type=str, default="stan_probs.csv")
    args = parser.parse_args()

    setup_seed()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("device:", device)

    cfg = toml.load(os.path.join(STAN_DIR, "config.toml"))
    model = PROTAC_STAN(cfg["model"])
    ckpt = os.path.join(STAN_DIR, "saved_models", "protac-stan.pt")
    model.load_state_dict(torch.load(ckpt))
    model = model.to(device).eval()
    print("model loaded:", ckpt)

    _, test_loader = PROTACLoader(root=args.root, name=args.name,
                                  batch_size=1, train_ratio=0.0)

    probs, preds = [], []
    smiles_list = []
    with torch.no_grad():
        for data in test_loader:
            protac = data["protac"].to(device)
            e3 = data["e3_ligase"].to(device)
            poi = data["poi"].to(device)
            out, _ = model(protac, e3, poi, mode="eval")  # log_softmax
            p = torch.exp(out)  # [B,2] 확률
            probs.extend(p[:, 1].cpu().numpy().tolist())
            preds.extend(torch.argmax(out, dim=1).cpu().numpy().tolist())
            smiles_list.extend(list(protac.smiles) if isinstance(protac.smiles, (list, tuple)) else [protac.smiles])

    # 입력 csv에서 식별자 회수
    raw = pd.read_csv(os.path.join(args.root, f"{args.name}.csv"))
    res = pd.DataFrame({
        "row": range(len(probs)),
        "Smiles": smiles_list,
        "stan_active_prob": probs,
        "stan_pred": preds,
    })
    # ID 컬럼 있으면 병합
    if "ID" in raw.columns and len(raw) == len(res):
        res.insert(0, "ID", raw["ID"].values)
    out_path = args.out
    res.to_csv(out_path, index=False)
    print(f"saved {len(res)} rows -> {out_path}")
    print(res.to_string(index=False))


if __name__ == "__main__":
    main()
