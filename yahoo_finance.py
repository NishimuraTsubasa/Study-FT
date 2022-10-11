
# 必要なライブラリーをインストール
from time import strftime
import matplotlib.pyplot as plt #描画ライブラリ
import pandas_datareader.data as web #データのダウンロードライブラリ
import pandas as pd
import numpy as np
import os
import datetime as date

# %%
data = ["^KLSE"]                        # Yahoo_financeの証券コード "2510.T":NOMURA BPI, "^GSPC";:S&P500
company = "yahoo"                       # Yahoo Finance からデータ取得
start = date.datetime(1983,1,4)         # データのスタート時点
end = date.datetime(2022,1,14)          # データの終了時点
today = date.date.today().strftime('%d/%m/%Y')

# %%
def get_val(data, company, start, end):
    try: 
        return ((web.DataReader(data, company, start, end)["Adj Close"].dropna()).pct_change(axis = 0)).resample(rule="W").mean() # resample で止めちゃうとエラーが出る
    except:
        pass
    return 0

df = get_val(data, company, start, today)
df.to_excel('Data_weekly_' + str(today) + '.xlsx', sheet_name = 'Data')
