"""
Phase TS-B — time-split 코호트에 PROTAC-STAN 추론(class-1 활성확률 덤프).
순서 보존(loader train_ratio=0.0 → no shuffle) + 입력 CSV(ik14/label/Uniprot 등)와
행순서 join. 원본 리포 미수정(stan_patch monkeypatch + 별도 러너).

실행 (cwd = ~/PROTAC-STAN):
  ~/anaconda3/envs/protac/bin/python <thispath>/run_stan_ts.py --sub ts_new
"""
import argparse, os, sys
import numpy as np, pandas as pd, torch, toml

_THIS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _THIS)
import stan_patch  # noqa: F401  (torch.load weights_only=False 패치, STAN import 전)

STAN_DIR = os.path.expanduser('~/PROTAC-STAN')
sys.path.insert(0, STAN_DIR)
from data_loader import PROTACLoader  # noqa: E402
from model import PROTAC_STAN  # noqa: E402

def setup_seed(seed=21332):
    import random
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True; torch.backends.cudnn.benchmark = False

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sub', required=True, help='data/<sub>/<sub>.csv')
    ap.add_argument('--batch', type=int, default=8)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    out = a.out or os.path.join('/home/kimjisan95/ar_protac_project/outputs', f'stan_{a.sub}_probs.csv')

    setup_seed()
    dev = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print('device:', dev)
    cfg = toml.load(os.path.join(STAN_DIR,'config.toml'))
    model = PROTAC_STAN(cfg['model'])
    model.load_state_dict(torch.load(os.path.join(STAN_DIR,'saved_models','protac-stan.pt')))
    model = model.to(dev).eval()

    root = os.path.join(STAN_DIR,'data',a.sub)
    _, loader = PROTACLoader(root=root, name=a.sub, batch_size=a.batch, train_ratio=0.0)

    probs, preds = [], []
    with torch.no_grad():
        for data in loader:
            protac = data['protac'].to(dev); e3 = data['e3_ligase'].to(dev); poi = data['poi'].to(dev)
            o, _ = model(protac, e3, poi, mode='eval')   # log_softmax [B,2]
            p = torch.exp(o)
            probs.extend(p[:,1].cpu().numpy().tolist())
            preds.extend(torch.argmax(o, dim=1).cpu().numpy().tolist())

    raw = pd.read_csv(os.path.join(root, f'{a.sub}.csv'))
    assert len(raw) == len(probs), f'행수 불일치 {len(raw)} vs {len(probs)}'
    raw = raw.copy()
    raw['stan_active_prob'] = probs
    raw['stan_pred'] = preds
    keep = ['ik14','Uniprot','E3 ligase Uniprot','label','stan_active_prob','stan_pred',
            'Molecular Weight','XLogP3','Topological Polar Surface Area','Smiles']
    raw[[c for c in keep if c in raw.columns]].to_csv(out, index=False)
    print(f'saved {len(raw)} rows -> {out}')
    print('  양성률(label):', round(raw["label"].mean(),3),
          '| 평균 확률:', round(float(np.mean(probs)),3))

if __name__ == '__main__':
    main()
