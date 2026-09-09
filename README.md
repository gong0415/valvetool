# 電磁閥設計工具

網頁版電磁閥模擬工具。開啟後在左側輸入參數，右側五個頁籤即時顯示模擬結果。

---

## 給使用者：怎麼打開

點開網址就能用，**不需要安裝任何東西**，手機和電腦的瀏覽器都可以。

> 網址：（部署完成後填在這裡）

**第一次開啟需要等約 1 分鐘**，畫面會顯示載入進度，這是在下載計算環境。
之後再開就會快很多（瀏覽器會留快取）。

操作說明看 [USAGE.md](USAGE.md)。

---

## 給管理者：怎麼把它放上網路

有兩種方式，**預設用方式一**。

### 方式一：GitHub Pages（只需要 GitHub 帳號）

程式直接在對方的瀏覽器裡執行（stlite + WebAssembly），不需要伺服器。

1. 把程式碼推到 GitHub（見下方步驟）
2. 到 repo 的 **Settings → Pages**，**Source** 選 **GitHub Actions**
3. 推上去後 Actions 會自動建置並部署，網址是
   `https://<你的帳號>.github.io/<repo名稱>/`

之後改程式只要 `git push`，網站自動重建。

實測：第一次開啟約 **70 秒**（要下載 Python 執行環境），之後改參數重算約 10 秒。
五個頁籤、所有圖表都正常。

### 方式二：Streamlit Community Cloud（速度較快）

需要多註冊一個免費帳號，但程式跑在伺服器上，開啟快、計算也快。

### 推到 GitHub

在 GitHub 開一個新的 repository，然後在專案資料夾執行：

```bash
git remote add origin https://github.com/<你的帳號>/<repo名稱>.git
git branch -M main
git push -u origin main
```

### 方式二的部署步驟

1. 到 https://share.streamlit.io 用 GitHub 帳號登入
2. 按 **New app**，選剛才那個 repository
3. **Main file path** 填 `solenoid_model/app.py`
4. 按 **Deploy**，等幾分鐘裝好套件

完成後會拿到一個 `.streamlit.app` 網址，把它貼到上面「網址」那行，再傳給對方。

改完 `git push`，網站會自動更新，網址不變。

---

## 存取權限

**兩種方式的免費方案都是公開的**，任何人拿到網址都能開啟。

- GitHub Pages：public repo 才能免費用 Pages，程式碼與網頁都公開
- Streamlit Cloud：可在 app 設定改成 Private，再逐一加入允許的帳號

若程式內容不宜外流，選 Streamlit Cloud 並設為 Private。

---

## 本機執行（開發用）

```bash
pip install -r requirements.txt
streamlit run solenoid_model/app.py
```

瀏覽器開 http://localhost:8501 。

---

## 文件

| 檔案 | 內容 |
|---|---|
| [USAGE.md](USAGE.md) | 使用說明：每個參數、每個頁籤的意義 |
| [ARCHITECTURE.md](ARCHITECTURE.md) | 程式架構與模組分工 |
| [ROADMAP.md](ROADMAP.md) | 開發進度 |

## 網頁版建置

`tools/build_site.py` 產生 `site/index.html`（stlite 包裝頁），
`tools/check_site.py` 驗證頁面參照的檔案齊全。
兩者由 `.github/workflows/pages.yml` 在每次 push 時自動執行。
新增或刪除模組後不必手動改清單，建置時會自動掃描。
