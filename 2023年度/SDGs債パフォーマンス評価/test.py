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


# 散布図を作成しPDFに出力する関数（カラム名を引数で指定可能）
def create_scatter_plots_by_sector_and_rating(df, metric, rating_column, sector_column, maturity_column):
    # dateでグループ化して各月末のデータを処理
    for date, group in df.groupby('date'):
        pdf_filename = f'{date.strftime("%Y-%m")}.pdf'

        with PdfPages(pdf_filename) as pdf:
            # 各セクターごとに処理
            for sector in group[sector_column].unique():
                sector_df = group[group[sector_column] == sector]

                # 各信用格付けごとに散布図を作成
                for rating in sector_df[rating_column].unique():
                    rating_df = sector_df[sector_df[rating_column] == rating]

                    # データが存在する場合のみ散布図を描画
                    if not rating_df.empty:
                        plt.figure(figsize=(10, 6))
                        sns.scatterplot(x=maturity_column, y=metric, data=rating_df)
                        plt.title(f'{sector} - {rating} - {date.strftime("%Y-%m")} - {metric}')
                        plt.xlabel('残存年数')
                        plt.ylabel(metric)

                        pdf.savefig()
                        plt.close()

# 関数のテスト実行（カラム名を引数として指定）
# ここでは、先ほどと同じテストデータを使用し、新しいカラム名指定引数を適用します。
create_scatter_plots_by_sector_and_rating(df, 'T-Spread', '信用格付け', '業種', '残存年数')



# 散布図を作成しPDFに出力する関数（修正版）
def create_scatter_plots_by_sector_and_rating(df, metric, rating_column, sector_column, maturity_column):
    # dateでグループ化して各月末のデータを処理
    for date, group in df.groupby('date'):
        pdf_filename = f'{date.strftime("%Y-%m")}.pdf'

        with PdfPages(pdf_filename) as pdf:
            # 各セクターごとに処理
            for sector in group[sector_column].unique():
                sector_df = group[group[sector_column] == sector]

                # 信用格付けのユニークな値を取得
                ratings = sector_df[rating_column].unique()

                # サブプロットの数を決定（最大6つ：3x2グリッド）
                num_plots = min(len(ratings), 6)
                cols = 3
                rows = (num_plots + 2) // cols

                # サブプロットの作成
                fig, axes = plt.subplots(rows, cols, figsize=(15, 10), squeeze=False)
                fig.suptitle(f'{sector} - {date.strftime("%Y-%m")}')

                for i, rating in enumerate(ratings):
                    if i < num_plots:
                        ax = axes[i // cols, i % cols]
                        rating_df = sector_df[sector_df[rating_column] == rating]

                        if not rating_df.empty:
                            sns.scatterplot(x=maturity_column, y=metric, data=rating_df, ax=ax)
                            ax.set_title(f'Rating: {rating}')
                            ax.set_xlabel('残存年数')
                            ax.set_ylabel(metric)

                # 空のサブプロットを非表示にする
                for j in range(i + 1, rows * cols):
                    axes[j // cols, j % cols].axis('off')

                plt.tight_layout()
                plt.subplots_adjust(top=0.9)  # タイトルのためのスペースを確保
                pdf.savefig()
                plt.close()

# 関数のテスト実行
create_scatter_plots_by_sector_and_rating(df, 'T-Spread', '信用格付け', '業種', '残存年数')






# 散布図を作成しPDFに出力する関数（修正版）
def create_scatter_plots_by_sector_rating_flag(df, metric, rating_column, sector_column, maturity_column, flag_column):
    # dateでグループ化して各月末のデータを処理
    for date, group in df.groupby('date'):
        pdf_filename = f'{date.strftime("%Y-%m")}.pdf'

        with PdfPages(pdf_filename) as pdf:
            # 各セクターごとに処理
            for sector in group[group[flag_column] == True][sector_column].unique():
                sector_df = group[group[sector_column] == sector]

                # 信用格付けを4つのグループに分割（四分位数を使用、重複削除）
                try:
                    quartiles = pd.qcut(sector_df[rating_column], 4, labels=False, duplicates='drop')
                    unique_quartiles = quartiles.unique()
                except ValueError:
                    # 四分位数の計算が不可能な場合は次のセクターへ
                    continue

                # サブプロットの数を実際に生成されたビンの数に基づいて設定
                num_plots = len(unique_quartiles)
                cols = 2
                rows = (num_plots + 1) // cols

                # サブプロットの作成
                fig, axes = plt.subplots(rows, cols, figsize=(15, 10), squeeze=False)
                fig.suptitle(f'{sector} - {date.strftime("%Y-%m")}')

                for i, q in enumerate(unique_quartiles):
                    ax = axes[i // cols, i % cols]
                    # 対応する四分位のデータをフィルタリング
                    quartile_df = sector_df[quartiles == q]

                    # フラグに基づいて色分けして散布図を描画
                    sns.scatterplot(x=maturity_column, y=metric, 
                                    hue=flag_column, palette=['blue', 'red'], 
                                    data=quartile_df, ax=ax)
                    ax.set_title(f'Quartile: {q+1}')
                    ax.set_xlabel('残存年数')
                    ax.set_ylabel(metric)
                    ax.legend(title=flag_column)

                # 空のサブプロットを非表示にする
                for j in range(i + 1, rows * cols):
                    axes[j // cols, j % cols].axis('off')

                plt.tight_layout()
                plt.subplots_adjust(top=0.9)  # タイトルのためのスペースを確保
                pdf.savefig()
                plt.close()

# 関数のテスト実行
create_scatter_plots_by_sector_rating_flag(df, 'T-Spread', '信用格付け', '業種', '残存年数', 'merged_green_flag')



# 散布図を作成しPDFに出力する関数（NaN値の処理と型の確認を追加）
def create_scatter_plots_by_sector_rating_flag(df, metric, rating_column, sector_column, maturity_column, flag_column):
    # dateでグループ化して各月末のデータを処理
    for date, group in df.groupby('date'):
        pdf_filename = f'{date.strftime("%Y-%m")}.pdf'

        with PdfPages(pdf_filename) as pdf:
            # 各セクターごとに処理
            for sector in group[sector_column].unique():
                sector_df = group[group[sector_column] == sector]

                # NaN値の行を削除し、データ型を確認
                sector_df = sector_df.dropna(subset=[rating_column])
                if not np.issubdtype(sector_df[rating_column].dtype, np.number):
                    # データ型が数値でなければスキップ
                    continue

                # 信用格付けを4つのグループに分割（四分位数を使用）
                try:
                    quartiles = pd.qcut(sector_df[rating_column], 4, labels=False, duplicates='drop')
                    unique_quartiles = quartiles.unique()
                except ValueError:
                    # 四分位数の計算が不可能な場合は次のセクターへ
                    continue

                # サブプロットの作成
                num_plots = len(unique_quartiles)
                cols = 2
                rows = (num_plots + 1) // cols
                fig, axes = plt.subplots(rows, cols, figsize=(15, 10), squeeze=False)
                fig.suptitle(f'{sector} - {date.strftime("%Y-%m")}')

                for i, q in enumerate(unique_quartiles):
                    ax = axes[i // 2, i % 2]
                    quartile_df = sector_df[quartiles == q]
                    sns.scatterplot(x=maturity_column, y=metric, 
                                    hue=flag_column, palette=['blue', 'red'], 
                                    data=quartile_df, ax=ax)
                    ax.set_title(f'Quartile: {q+1}')
                    ax.set_xlabel('残存年数')
                    ax.set_ylabel(metric)
                    ax.legend(title=flag_column)

                for j in range(i + 1, rows * cols):
                    axes[j // cols, j % cols].axis('off')

                plt.tight_layout()
                plt.subplots_adjust(top=0.9)
                pdf.savefig()
                plt.close()

# 関数のテスト実行
create_scatter_plots_by_sector_rating_flag(df, 'T-Spread', '信用格付け', '業種', '残存年数', 'merged_green_flag')
