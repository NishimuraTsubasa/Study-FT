# %%
# ライブラリーのインポート
from functools import update_wrapper
import pandas as pd
import numpy as np
import numpy.linalg as lin
import cvxopt as cvx
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import sys
import datetime
import scipy.optimize as sco
import re

# %%
# ==================================================================================
# インポート
# ==================================================================================
# # strat   : バックテスト開始日
# # end     : バックテスト終了日
# # sim_data        : シミュレーションデータ(欠損データは削除)
    # 週次データを入れる様にする（柴原さんに依頼, 現在は日時データ）
# # Int_solution    : 初期解（target_data_fundデータの初期値(Weights)） 
# # method  : 収益率の計算方法
    # log_return    : 対数収益率
    # diff_return   : 算術収益率
sim_data = pd.read_csv("C:/Users/bldyr/OneDrive/デスクトップ/自己研鑽用/02_FT勉強会/2022年度/基礎研究/for_analysis.csv",
                 encoding="shift-jis", index_col = 0) .dropna()
TDF_weight = pd.read_csv("C:/Users/bldyr/OneDrive/デスクトップ/自己研鑽用/02_FT勉強会/2022年度/基礎研究/Weight_of_TDF.csv",
                        encoding = "shift-jis", index_col = 0).dropna()
Int_solution = np.array(TDF_weight.iloc[0, :]) 
method = "log_return"   # [log_return, diff_return] 

# ==================================================================================
# 最適化における設定
# ==================================================================================
# # target_return   : 目標リターン（投資家が求める水準）
# # upper_restriction : 投資銘柄のウェイト上限
# # lower_restriction : 投資銘柄のウェイト下限
    # [国内債券，国内株式，外国債券，外国株式，その他，その他]

target_return = 0.0015   # 最適化におけるリターン
upper_restriction = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
lower_restriction = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

# %%
def calc_return(method, sim_data, i):
    """ 収益率を計算する関数
    Parameters
    ----------
    method : str
        収益率の計算方法[log_return, diff_return]
    sim_data : pd.DataFrame
        シミュレーションデータ
    Returns
    -------
    return_data : pd.DataFrame
        収益率データ
    """
    if method == "log_return":
        return_data = np.log(sim_data.iloc[0:i,:]).diff().dropna() 
    elif method == "diff_return":
        return_data = (sim_data.iloc[0:i,:].diff() / sim_data.iloc[0:i,:]).dropna() 
    else:
        print("正しい計算方法を入力してください．")
    return return_data

def calc_port_ind(return_data, i):
    """ ポートフォリオの期待収益率と分散共分散行列，相関係数行列を計算する関数
    Parameters
    ----------
    return_data : pd.DataFrame
        収益率データ
    i : int
        計算する時点の数(2以上である必要がある)
    Returns
    -------
    expected_return : np.array
        各資産の期待収益率
    cov : np.array
        ポートフォリオの分散共分散行列
    corr : np.array
        ポートフォリオの相関係数行列
    """
    return_data = return_data.iloc[0:i,:]
    expected_return = return_data.mean()
    cov = return_data.cov()
    corr = return_data.corr()
    return expected_return, cov, corr

def min_func_var(weights):
    """ 目的関数：最小化する関数
    Parameters
    ----------
    weights : np.array
        ポートフォリオウェイト
    cov : np.array
        ポートフォリオの分散共分散行列
    Returns
    -------
    vol : constant
        ポートフォリオのボラティリティ
    """
    vol = np.dot(weights.T, np.dot(cov, weights))
    return vol

def set_bnds(TDF_weight, i, upper_restriction, lower_restriction):
    """ 上下限制約を作成する関数
    Parameters
    ----------
    TDF_weights : pd.DataFrame
        ターゲットデータファンドのグライドパス
    i : int
        計算する時点の数(2以上である必要がある)
    upper_restriction : list
        各銘柄の投資ウェイト上限
    lower_restriction : list
        各銘柄の投資ウェイト下限
    Returns
    -------
    bnds : list
        ウェイトの上下限
    """

    for j in range(len(TDF_weight.columns)):
        if j == 0: 
            bnds = [(TDF_weight.iloc[i, j] - lower_restriction[j], TDF_weight.iloc[i, j] + upper_restriction[j])]
        else:
            bnds.append((TDF_weight.iloc[i, j] - lower_restriction[j], TDF_weight.iloc[i, j] + upper_restriction[j]))

    return bnds

# %%
# ==================================================================================
# 最適化の設定と実行
# ==================================================================================
# # ポートフォリオの指標計算
return_data = calc_return(method, sim_data, 1000)
expected_return, cov, corr = calc_port_ind(return_data, 1000)   # 2以上じゃないとエラーを吐くので，とりあえず2で設定（後ほどfor文で回す箇所）

# # 制約条件
# # # 'type'に'eq'を指定すると等式制約で，'fun'に左辺=0となる関数を定義
# # # 'type'に'ineq'を指定すると不等式制約で，'fun'は左辺>=0となる関数を定義
cons = [{'type': 'eq', 'fun': lambda x: np.sum(x) - 1},                                 # ウェイト合計を1に設定（x : ウェイト）
    {'type': 'ineq', 'fun': lambda x: np.sum(expected_return * x) - target_return}]     # リターンが目標リターンを超過するように設定
bnds = set_bnds(TDF_weight, 1, upper_restriction, lower_restriction)                    # 第2引数は時点情報を追加

# 最適化の実施
opts = sco.minimize(fun = min_func_var, x0 = Int_solution, method='SLSQP', bounds=bnds, constraints=cons)

# %%
# ==================================================================================
# 効率的フロンティアの作図(お遊び)
# ==================================================================================
# # 各資産のリターン最大値・最小値を計算
mean_high = expected_return.max()
mean_low = expected_return.quantile(0.05)
trets = np.linspace(mean_low, mean_high, 100)

# # 制約条件（ウェイト）
bnds = [(0, None)] * len(expected_return)

# # 最適化を複数回実施
tvols = []
weight = []
for tret in trets:
    cons = [{'type': 'eq', 'fun': lambda x: np.sum(x) - 1}, 
            {'type': 'ineq', 'fun': lambda x: np.sum(expected_return * x) - tret}]
    res = sco.minimize(fun=min_func_var, x0=Int_solution, method='SLSQP', bounds = bnds, constraints = cons)
    tvols.append(np.sqrt(res['fun']))
    weight.append(opts["jac"])
tvols = np.array(tvols)

# # 効率的フロンティアの作成
plt.figure(figsize=(10, 6))
plt.scatter(tvols, trets, c = trets / tvols, marker = 'x')
plt.grid(True)
plt.xlabel('expected volatility')
plt.ylabel('expected return')
