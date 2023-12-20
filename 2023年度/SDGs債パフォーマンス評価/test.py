import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from datetime import datetime
import os

# 仮のデータフレームを作成（実際のデータに基づいて調整が必要）
# この例ではランダムなデータを使用
dates = pd.date_range(start='2022-03-31', end='2023-12-31', freq='M')
data = {
    '残存年数': pd.np.random.rand(len(dates)) * 10,  # 0〜10年
    'date': dates,
    '信用格付け': pd.np.random.choice(['AAA', 'AA', 'A', 'BBB'], len(dates)),
    '業種': pd.np.random.choice(['Sector1', 'Sector2', 'Sector3'], len(dates)),
    'T-Spread': pd.np.random.rand(len(dates)) * 200,  # 0〜200基点
    'G-Spread': pd.np.random.rand(len(dates)) * 200,  # 0〜200基点
    'OAS': pd.np.random.rand(len(dates)) * 200  # 0〜200基点
}
df = pd.DataFrame(data)

# 指定された指標に基づいて散布図を作成しPDFに出力する関数
def create_scatter_plots_by_sector_and_rating(df, metric):
    # dateでグループ化して各月末のデータを処理
    for date, group in df.groupby('date'):
        # PDFファイル名（YYYY-MM.pdf形式）
        pdf_filename = date.strftime('%Y-%m') + '.pdf'
        
        with PdfPages(pdf_filename) as pdf:
            # 各セクターごとに処理
            for sector in group['業種'].unique():
                # セクターごとのデータフレームを取得
                sector_df = group[group['業種'] == sector]

                # 信用格付けごとに散布図を作成
                plt.figure(figsize=(10, 6))
                sns.scatterplot(x='残存年数', y=metric, hue='信用格付け', data=sector_df)
                plt.title(f'Sector: {sector} - {date.strftime("%Y-%m")} - {metric}')
                plt.legend(title='信用格付け')
                plt.xlabel('残存年数')
                plt.ylabel(metric)
                
                # PDFにページを追加
                pdf.savefig()
                plt.close()

    # 完成したPDFファイルのパスを返す
    return os.listdir('.')

# 仮のデータフレームを使って関数をテスト（T-Spreadを指標とする）
generated_pdfs = create_scatter_plots_by_sector_and_rating(df, 'T-Spread')
generated_pdfs
