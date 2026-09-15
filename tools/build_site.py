"""產生 site/index.html —— GitHub Pages 用的 stlite（瀏覽器內 Streamlit）包裝頁。

執行： python tools/build_site.py
程式碼有增減模組時重跑一次，檔案清單會自動更新。
"""
import json
import pathlib
import subprocess

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "site" / "index.html"

STLITE_VERSION = "0.85.1"

# app.py 只從 report 套件匯入 generate_layout；其餘 report 腳本是離線產圖用的，
# 收進來只會加大下載量。
PNGS = [
    "docs/spec/valve_schematic.png",
    "docs/spec/winding_window_schematic.png",
    "docs/spec/seal_land_schematic.png",
    "docs/spec/layout_section.png",
    "docs/spec/layout_actuation.png",
    "docs/spec/layout_architecture.png",
]


def collect_files():
    tracked = subprocess.run(
        ["git", "ls-files", "solenoid_model"],
        cwd=ROOT, capture_output=True, text=True, check=True,
    ).stdout.split()

    mods = [
        f for f in tracked
        if f.endswith(".py") and "/tests/" not in f
        and (not f.startswith("solenoid_model/report/")
             or f.endswith(("__init__.py", "generate_layout.py")))
    ]
    missing = [p for p in PNGS if not (ROOT / p).exists()]
    if missing:
        raise SystemExit(f"缺少圖檔：{missing}")
    return {p: {"url": p} for p in sorted(mods) + PNGS}


HTML = """<!doctype html>
<html lang="zh-Hant">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>電磁閥設計工具</title>
<link rel="stylesheet"
      href="https://cdn.jsdelivr.net/npm/@stlite/browser@{ver}/build/stlite.css" />
<style>
  body {{ margin: 0; font-family: system-ui, "PingFang TC", sans-serif; }}
  #boot {{
    position: fixed; inset: 0; display: flex; align-items: center;
    justify-content: center; background: #fff; z-index: 9999; text-align: center;
  }}
  #boot .inner {{ max-width: 30rem; padding: 2rem; }}
  #boot h1 {{ font-size: 1.25rem; margin: 0 0 .75rem; }}
  #boot p {{ color: #555; line-height: 1.7; margin: .4rem 0; }}
  .bar {{
    height: 4px; background: #eee; border-radius: 2px; overflow: hidden;
    margin-top: 1.25rem;
  }}
  .bar span {{
    display: block; height: 100%; width: 40%; background: #ff4b4b;
    animation: slide 1.4s ease-in-out infinite;
  }}
  @keyframes slide {{
    0% {{ transform: translateX(-100%); }}
    100% {{ transform: translateX(350%); }}
  }}
</style>
</head>
<body>
<div id="boot">
  <div class="inner">
    <h1>電磁閥設計工具</h1>
    <p>正在載入計算環境，第一次開啟約需 30 秒。</p>
    <p>程式在你的瀏覽器內執行，資料不會上傳。</p>
    <div class="bar"><span></span></div>
  </div>
</div>

<script type="module">
import {{ mount }} from
  "https://cdn.jsdelivr.net/npm/@stlite/browser@{ver}/build/stlite.js";

mount({{
  entrypoint: "solenoid_model/app.py",
  requirements: ["numpy", "scipy", "matplotlib"],
  files: {files},
  streamlitConfig: {{ "client.toolbarMode": "viewer" }},
}}, document.body);

// stlite 掛載後才移除載入畫面；否則畫面會先閃一下空白。
const boot = document.getElementById("boot");
const done = () => boot && boot.remove();
new MutationObserver((_, obs) => {{
  if (document.querySelector('[data-testid="stAppViewContainer"]')) {{
    done(); obs.disconnect();
  }}
}}).observe(document.body, {{ childList: true, subtree: true }});
setTimeout(done, 120000);  // 保險：真的卡住也不要一直蓋著畫面
</script>
</body>
</html>
"""


def main():
    files = collect_files()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    # 同時寫出清單，讓 check_site.py 不必反解 HTML。
    (OUT.parent / "files.json").write_text(
        json.dumps(sorted(files), indent=2), encoding="utf-8")
    OUT.write_text(
        HTML.format(ver=STLITE_VERSION,
                    files=json.dumps(files, indent=4, ensure_ascii=False)),
        encoding="utf-8",
    )
    n_py = sum(1 for k in files if k.endswith(".py"))
    print(f"寫入 {OUT}（{n_py} 個模組 + {len(files) - n_py} 張圖）")


if __name__ == "__main__":
    main()
