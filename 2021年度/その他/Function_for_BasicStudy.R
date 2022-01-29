# 【CaseStudy】
# Reirement Income Recipes in R のサンプルコードを使用して，日本データを分析
# LRPを算出・分析

library(dplyr)　# データ抽出のためのライブラリ
library(tidyverse) # データ処理や作図に使うライブラリ
library(lubridate) # 時間データ処理に使うライブラリ
# ディレクトリの指定
# setwd()

# インプットの設定
N <- 10000  # シミュレーション回数
sex <- "男性"
birthyear <- 1953 # 誕生日
thisyear <- 2018  # 現在時点
x <- thisyear - birthyear　# 年齢
maxAge <- 120

# ファイルの読み込み
setwd("C:/Users/bldyr/OneDrive/デスクトップ/自己研鑽用/データ取得")
HistData <- read.csv("data_weekly.csv", fileEncoding = "SJIS")
HistData <- select(StockData, 1, 2, 3)

# データ確認
head(HistData, 10)

bonds <- select(HistData, 1, 3)
Stock <- select(HistData, 1, 2)

####時点ごとのウェイトを決める関数#####
CalcWeight <- function(w, Tx, x){
  scock_weight <- c()
  bonds_weight <- c()
  t <- Tx-x
  for (i in 1:t){
    stock_weight[i] <- w-(i+x)
    bonds_weight[i] <- 100-stock[i]
    if (stock[i]<=0) {break}
  }
}



