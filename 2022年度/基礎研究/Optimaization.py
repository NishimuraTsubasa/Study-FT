



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