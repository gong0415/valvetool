# 電磁閥設計工具

網頁版電磁閥模擬工具。開啟後在左側輸入參數，右側五個頁籤即時顯示模擬結果。

---

## 給使用者：怎麼打開

點開網址就能用，**不需要安裝任何東西**，手機和電腦的瀏覽器都可以。

> 網址：（部署完成後填在這裡）

第一次開啟如果看到「app is waking up」，等約 30 秒讓它啟動即可。

操作說明看 [USAGE.md](USAGE.md)。

---

## 給管理者：怎麼把它放上網路

用 Streamlit Community Cloud，免費，不用自己準備伺服器。

### 步驟一：把程式碼推到 GitHub

在 GitHub 開一個新的 repository，然後在專案資料夾執行：

```bash
git remote add origin https://github.com/<你的帳號>/<repo名稱>.git
git branch -M main
git push -u origin main
```

### 步驟二：部署

1. 到 https://share.streamlit.io 用 GitHub 帳號登入
2. 按 **New app**，選剛才那個 repository
3. **Main file path** 填 `solenoid_model/app.py`
4. 按 **Deploy**，等幾分鐘裝好套件

完成後會拿到一個 `.streamlit.app` 網址，把它貼到上面「網址」那行，再傳給對方。

### 之後要改程式

改完 `git push`，網站會自動更新，網址不變。

---

## 存取權限

免費方案的 app 預設是**公開**的，任何人拿到網址都能開啟。若內容不宜公開，在 Streamlit Cloud 的 app 設定裡把它改成 Private，再逐一加入允許的 Google 或 GitHub 帳號。

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
