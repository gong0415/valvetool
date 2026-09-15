"""驗證 site/index.html 參照的每個檔案都真的存在於 site/ 底下。

在 CI 建置時擋下「頁面上線但模組缺檔」這種只有使用者才會發現的錯誤。
執行： python tools/check_site.py
"""
import json
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
SITE = ROOT / "site"
INDEX = SITE / "index.html"


def main():
    if not INDEX.exists():
        sys.exit("找不到 site/index.html —— 請先跑 tools/build_site.py")

    manifest = SITE / "files.json"
    if not manifest.exists():
        sys.exit("找不到 site/files.json —— 請先跑 tools/build_site.py")
    files = json.loads(manifest.read_text(encoding="utf-8"))

    missing = [p for p in files if not (SITE / p).is_file()]
    # 空的 __init__.py 是正常的套件標記檔，不算錯。
    empty = [p for p in files
             if (SITE / p).is_file() and (SITE / p).stat().st_size == 0
             and not p.endswith("__init__.py")]

    for p in missing:
        print(f"缺檔: {p}")
    for p in empty:
        print(f"空檔: {p}")
    if missing or empty:
        sys.exit(f"檢查失敗：{len(missing)} 缺檔、{len(empty)} 空檔")

    if not (SITE / "solenoid_model/app.py").is_file():
        sys.exit("進入點 solenoid_model/app.py 不在 site/ 底下")

    # 反向檢查：generate_layout 會產生的圖，清單裡都必須有。上面的檢查只
    # 驗「清單列的檔案存在」，漏列的圖它看不到——layout_actuation.png 新增
    # 後就是這樣在頁面上 404，而 CI 全綠。
    missing_figs = [p for p in _declared_figures() if p not in files]
    if missing_figs:
        for p in missing_figs:
            print(f"圖檔未列入 build_site.PNGS: {p}")
        sys.exit(f"檢查失敗：{len(missing_figs)} 張圖漏列，頁面會 404")

    print(f"檢查通過：{len(files)} 個檔案齊全"
          f"（含 {len(_declared_figures())} 張 app 會顯示的圖）")


_FIG_RE = re.compile(
    r'^[A-Z_]+_OUT\s*=\s*_DOCS\s*/\s*["\']([^"\']+\.png)["\']', re.M)


def _declared_figures():
    """generate_layout 宣告為輸出的圖，轉成 repo 相對路徑。

    以正規式讀原始碼的 *_OUT 賦值，不 import 該模組：這個腳本與
    build_site.py 一樣跑在只有標準庫的 CI 步驟裡（workflow 沒有任何
    pip install），而 generate_layout 會拉進 matplotlib 與 scipy。
    先前版本直接 import，於是 CI 以 ModuleNotFoundError 中止部署。

    仍以模組本身為單一真相來源——只是改成讀而不是執行。
    """
    src = (ROOT / "solenoid_model" / "report"
           / "generate_layout.py").read_text(encoding="utf-8")
    names = _FIG_RE.findall(src)
    if not names:
        sys.exit("在 generate_layout.py 找不到任何 *_OUT 圖檔宣告——"
                 "命名或寫法改了，請同步更新 _FIG_RE")
    return sorted(f"docs/spec/{n}" for n in names)


if __name__ == "__main__":
    main()
