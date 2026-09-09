"""驗證 site/index.html 參照的每個檔案都真的存在於 site/ 底下。

在 CI 建置時擋下「頁面上線但模組缺檔」這種只有使用者才會發現的錯誤。
執行： python tools/check_site.py
"""
import json
import pathlib
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

    print(f"檢查通過：{len(files)} 個檔案齊全")


if __name__ == "__main__":
    main()
