# -*- coding: utf-8 -*-

hyperparam = {
    # データ前処理用
    "data_processer" : {
        "foreign_folder_path" : "", # 外債部データフォルダパス
        "domestic_folder_path" : "", #  債券部データフォルダパス
        "password" : "daiichilife", # ファイルのパスワード
        "principle_key" : ["Date, ISIN"], # 主キー
        "extension" : ".xlsx", # Excelファイルの拡張子 (.xlsx 形式)
    }
}

hyperparam = {
    # データ前処理用
    "data_processer" : {
        # 外債部データに関する処理
        "foreign_bond" : {
            "folder_path" : "", # フォルダパス
            "password" : "daiichilife", # ファイルのパスワード
            "principle_key" : ["Date, ISIN"], # 主キー
            "extension" : ".xlsx", # Excelファイルの拡張子 (.xlsx 形式)
        },
        # 債券部データに関するデータ処理
        "domestic_bond" : {
            "folder_path" : "", #  債券部データフォルダパス
            "password" : "daiichilife", # ファイルのパスワード
            "principle_key" : ["Date, ISIN"], # 主キー
            "extension" : ".xlsx", # Excelファイルの拡張子 (.xlsx 形式)          
        },
    }
}
