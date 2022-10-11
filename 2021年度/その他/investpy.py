
"""
https://oeconomicus.jp/2021/06/python-infestpy/
"""
# %%
import investpy
from datetime import date

# 本日日付を準備
today = date.today().strftime('%d/%m/%Y')

# ETF SPYのデータ取得
df = investpy.get_etf_historical_data(index='FBMKLCI', country='Malaysia',
                                      from_date='01/01/2010',to_date=today)
# %%
# AAPLの株価データの取得
df = investpy.get_stock_historical_data(stock='AAPL', country='united states',
                                        from_date='01/01/2010',to_date='10/06/2021')
