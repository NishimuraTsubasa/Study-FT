import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

# 整数型の日付を実際の日付に変換する関数
def convert_int_to_date(int_date):
    str_date = str(int_date)
    year = int(str_date[:4])
    month = int(str_date[4:6])
    date_obj = pd.Timestamp(year=year, month=month, day=1) - pd.Timedelta(days=1)
    return date_obj

# 散布図を作成しPDFに出力する関数
def create_scatter_plots_by_sector_and_rating(df, metric):
    pdf_files = []

    # dateでグループ化して各月末のデータを処理
    for date, group in df.groupby('date'):
        pdf_filename = f'{date.strftime("%Y-%m")}.pdf'
        pdf_files.append(pdf_filename)

        with PdfPages(pdf_filename) as pdf:
            # 各セクターごとに処理
            for sector in group['業種'].unique():
                sector_df = group[group['業種'] == sector]

                # 各信用格付けごとに散布図を作成
                for rating in sector_df['信用格付け'].unique():
                    rating_df = sector_df[sector_df['信用格付け'] == rating]

                    # データが存在する場合のみ散布図を描画
                    if not rating_df.empty:
                        plt.figure(figsize=(10, 6))
                        sns.scatterplot(x='残存年数', y=metric, data=rating_df)
                        plt.title(f'{sector} - {rating} - {date.strftime("%Y-%m")} - {metric}')
                        plt.xlabel('残存年数')
                        plt.ylabel(metric)

                        pdf.savefig()
                        plt.close()

    return pdf_files

# テスト用データフレームの作成
dates = [202203, 202204, 202205, 202206, 202207, 202208, 202209, 202210, 202211, 202212, 202301, 202302, 202303, 202304, 202305, 202306, 202307, 202308, 202309, 202310, 202311, 202312]
data = {
    '残存年数': np.random.rand(len(dates)) * 10,
    'date': dates,
    '信用格付け': np.random.choice(['AAA', 'AA', 'A', 'BBB'], len(dates)),
    '業種': np.random.choice(['Sector1', 'Sector2', 'Sector3'], len(dates)),
    'T-Spread': np.random.rand(len(dates)) * 200,
    'G-Spread': np.random.rand(len(dates)) * 200,
    'OAS': np.random.rand(len(dates)) * 200
}
df = pd.DataFrame(data)
df['date'] = df['date'].apply(convert_int_to_date)

# PDF生成の実行
generated_pdfs = create_scatter_plots_by_sector_and_rating(df, 'T-Spread')
generated_pdfs
