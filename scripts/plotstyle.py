"""matplotlib 한글 폰트 설정(Noto Sans CJK KR). 그림 스크립트에서 import."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

_PATH = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
try:
    fm.fontManager.addfont(_PATH)
    _name = fm.FontProperties(fname=_PATH).get_name()
    plt.rcParams["font.family"] = _name
except Exception as e:
    print("[plotstyle] 폰트 로드 실패:", e)
plt.rcParams["axes.unicode_minus"] = False
