
# 3.1.2�?
2*3+4
2^3
pi
sin(pi/2)
sin(0); cos(pi/2); tan(pi/2)
sin(3*pi
)
exp(1)
log(100)
log10(100); log(100,10)


# 3.1.3�?
a <- 2
a
A
abc <- 3; 身長 <- 172
abc; 身長
objects()
rm(身長)
身長
objects()
?rm
help(objects)
a*3+2
b <- a^2+2
b
3^2+2 -> c
c


# 3.1.4�?
(v <- c(2,3,6,1,4))
v[2]; v[6]  # v[6]はな�?のでnot availableである�?
length(v)
v*2+3
w <- c(1,2,3)
v+w; v
sin(v*pi/2)
1:4
5.5:1
v[1:3]


# 3.1.5�?
v <- c(1,3,5,7,9,11)
A <- matrix(v,2,3); B <- matrix(v,2,3,byrow=TRUE); C<- matrix(v,3,2)
A; B; C
A+B; A-B; A+C
A^2
dim(A)
v
dim(v) <- c(2,3)
v
c(v)
t(A)
A[2,]; A[,3]; A[2,3]
name <- c("gure", "mike"); age <- c(8,7); sex <- c("m", "f")
d <- data.frame(name, age, sex)
d
d[1,]
d[,1]


# 3.2.1�?
set.seed(60)
(test <- floor(rnorm(40, 50, 10)))
max(test); min(test)
range(test)
summary(test)
hist(test)
stem(test)
qqnorm(test)
qqline(test)
sum(test)
length(test)
sum(test)/length(test)
mean(test)
median(test)
var(test)
sum((test-mean(test))^2)/length(test); sum((test-mean(test))^2)/(length(test)-1)


# 3.2.2�?
?cars
cars
plot(cars) # �?ータフレー�?から散�?図を描�?
plot(cars[,1],cars[,2]) # 同じ長さ�?�2つのベクトルでも散�?図が描ける
cov(cars[,1],cars[,2]); cor(cars[,1],cars[,2]) # 共�?散と相関係数
kaiki <- lm(cars[,2]~cars[,1]) # 単回帰直線を求め�?
summary(kaiki) # 回帰式�?��?報を得る
abline(kaiki) # 散�?図中に回帰直線を引く


# 3.2.3�?
?iris
iris
boxplot(iris[1:50,1], iris[51:200,1]) # 箱ひげ図を描�?
t.test(iris[1:50,1], iris[51:200,1])


# 3.3節
p <- 1/6; q <- 1-p
choose(10,3)*p^3*q^7; dbinom(3,10,p) # 10回中3�?1が�?�る確�?
sum(choose(10,0:3)*p^(0:3)*q^(10-0:3)); pbinom(3,10,p) # 10回中1が�?�る�?��?3回以下�?�確�?
qbinom(0.5, 10, p) # F(x)�?0.5以上になる最小�?�x
rbinom(5,10,p) # 10回振ること�?5セ�?ト繰り返したとき�?�1の出る回数のセ�?トを表す乱数
dnorm(175, 171, 5) # x=175での�?度関数f(x)の値
pnorm(175, 171, 5) # x=175での�?�?関数F(x)の値
qnorm(0.7, 171, 5) # �?位数すなわちF(x)=0.7となるxの値
rnorm(6, 171, 5) # 乱数すなわち確�?変数の実現値�?6�?
dnorm(mean=171, 175, sd=5) # 引数の名前を指定することで�?番を変え�?
dnorm(2) # 母数を指定しな�?とき�?��?フォルトが適用され、N(0,1)の�?度関数f(x)の値
?distributions
rnorm(5)
rnorm(5)
set.seed(10); rnorm(5)
set.seed(10); rnorm(5)


# 3.4.1�?
curve(sin(x))
curve(sin(x), xlim=c(-3*pi, 3*pi))
curve(cos(x), add=T)
x <- 1:10; y <- sin(x) # x=(1,2,3,...,10), y=(sin(1),sin(2),...,sin(10))を代入
plot(x,y)
barplot(x)
hist(y)
pie(x)
plot(x)
plot(y)
plot(x,y,type="b", xlab="yoko", main="sin x", xlim=c(2,7), log="x", col="blue", pch="猫", lty=3)
plot(x,y)
oldpar <- par()	# 既存�?�グラフィ�?クスパラメータの保�?
par(pch=5)		# グラフィ�?クスパラメータを指�?
plot(x,y)
par(oldpar)		# �?のグラフィ�?クスパラメータに戻�?
plot(x,y)


# 3.4.2�?
plot(x,y)
lines(x,y)
points(x, cos(x))
text(locator(1),"y=sin x")
identify(x,y)
identify(x,y,plot="FALSE")
i <- identify(x,y)
i


# 3.5.1�?
library() # インスト�?�ルされて�?るパ�?ケージの一覧
search() # ロードされて�?るパ�?ケージの一覧
# install.packages("evd") # evdパッケージのインスト�?�ル
library(evd) # evdパッケージのロー�?
# install.packages("ismev") # evdパッケージのインスト�?�ル
library(ismev)
# install.packages("extRemes") # extRemesパッケージのインスト�?�ル
library(extRemes)
library()
search()


# 3.6.1�?
test <- floor(rnorm(40,50,10))
s <- sqrt(sum((test-mean(test))^2)/length(test))
(47-mean(test))/s*10+50
z <- (test-mean(test))/s*10+50; z; mean(z)
m <- mean(test)
hensachi <- function(x){(x-m)/s*10+50}
hensachi(47)
hensachi
hensachi2 <- function(x, group){
 # xは�?点 groupは�?団全員の得点ベクトル
 m1 <- mean(group) # �?団の平�?
 s1 <- sqrt(sum((group-m1)^2)/length(group)) # �?団の標準偏差
 (x-m1)/s1*10+50}
hensachi2(47,test)
m1; s1


# 3.6.2�?
hensachi <- function(x) {(x-m)/s*10+50}
hensachi
var
?mean
methods(mean)
mean.default
sin


# 3.6.3�?
factsum <- function(n){
  p <- 1; s <- 0 
  for (i in 1:n){ p <- p*i; s <- s+i}
  c(p,s)}
factsum(5)


# 3.6.4�?
(d <- as.data.frame(matrix(1:9,3,3)))
a <- NULL # aを空の変数とする
for (j in 1:3){a <- c(a, max(d[,j]))} # 既存�?�aにj列目の最大値を付け�?える
a
apply(d,2,max) # 列ごとの最大値�?2は列を表�?
apply(d,1,max) # 行ごとの最大値�?1は行を表�?


# 3.7.1�?
getwd() # 現在の作業�?ィレクトリを知�?
setwd("C:/Users/saigo/Documents/R") # 作業�?ィレクトリを変更する
getwd()
setwd("C:/Users/saigo/Documents/R2") # 存在しな�?フォルダに変更する


# 3.7.2�?
rm(list=ls()) # 変数がすべて消去される�?�で注�?
(x <- 3:10)
# q()
x
objects()
tenki <- read.csv("tenki.csv")
tenki
read.csv("tenki2.csv") # 表頭を�?�名とみなしてしま�?
read.csv("tenki2.csv", header=F) # 表頭を�?�名とせず、新たな列名をつける
read.table("tenki.csv") # スペ�?�ス�?タブで区�?られて�?な�?ので3�?1列とみな�?
read.table("tenki.csv", sep=",") # 表頭を�?�名とみなさず�?ータとみ�?
read.table("tenki.csv", sep=",", header=T) # 正しい読み込みとなっ�?
?BOD
BOD
write.table(BOD, "BOD.txt") # txtファイル�?が、�??容としてはssvとなって�?る�?
write.table(BOD, "BOD.csv") # csvとしたがスペ�?�ス区�?りで正しくな�?ファイルとな�?
write.table(BOD, "BOD2.csv", sep=",") # 正しいcsvファイルができる
write.csv(BOD, "BOD3.csv") # 正しいcsvファイルができる
?datasets
saikoro
source("sai.R") # プログラ�?の読み込み
saikoro
saikoro(4)
x <- 1:10; y <- 1:10
save(x,y, file="suuretsu.Rdata") # xをファイルsuuretsu.Rdataに保�?
rm(x,y)
x
y
load("suuretsu.Rdata") # ファイルsuuretsu.Rdataを読み込み
x; y
x <- 1:10
dput(x, "suuretsu.txt")
dget("suuretsu.txt")



#-----------------------------------------------------------------------------------------------------------------------------------------------

# 4�?



# 4.2節
library(evd)
curve(drweibull(x, loc=3, scale=2, shape=1.5), xlim=c(-3,3)) # ワイブルの�?度関数のグラ�?
pfrechet(5, loc=3, scale=2, shape=1.5) # フレシェの�?�?関数のx=5における値
qgumbel(0.2,loc=3,scale=2)
rgev(5, shape=-1) # 一般極値�?�?に従う乱数�?5�?



#------------------------------------------------------------------------------------------------------------------------------------------------

# 6��



# 6.2.1�?
library(evd)
?fgev # ブラウザのヘルプが立ち上がる。�?��?�ジ下からindexに移�?
library(ismev)
?ismev # ブラウザのヘルプが立ち上がる。�?��?�ジ下からindexに移�?
demo(exchange.rate) # �?モの一つを実�?



# 6.2.2�?
library(evd)
pgev(0.5)	# GEV(0,1,0)の�?�?関数のx=0.5での値
exp(-exp(-0.5)) # 上�?��?�?関数を定義通り計算す�?
curve(dnweibull(x, loc=1, scale=2, shape=3), xlim=c(-4,2)) 
# Weibull(1,2,3)�?�?の�?度関数のグラ�?
curve(pnweibull(x, loc=1, scale=2, shape=3), xlim=c(-4,2)) 
# Weibull(1,2,3)�?�?の�?�?関数のグラ�?
qgumbel(0.7, scale=2) # Gumbel(0,2)の�?�?関数F(x)=0.7となるx 
pgumbel(2.061861, scale=2) # 上�?��?位数を確�?
rfrechet(5) # Frechet(0,1,1)�?�?に従う乱数�?5�?


# 6.2.3�?
library(evd)
set.seed(5)		# 乱数の種を指�?
x <- rgev(100)	# GEV(0,1,0)の乱数�?100個発�?
plot(sort(x)) # xを�??�?に並べてグラフで様子を見る
(x.gev.evd <- fgev(x))	# GEVモ�?ルで推定す�?
fitted(x.gev.evd)
plot(x.gev.evd)	# 推定結果のEDA
# 以下�?�1つの画面に4つのプロ�?トを出す方�?
oldpar <- par() # 現在の設定を退避する 
par(mfrow=c(2,2)) # 画面�?2*2型で横を�?�に入れる
plot(x.gev.evd) # プロ�?トを描画
par(oldpar) # �?の設定に戻�?

x.prof <- profile(x.gev.evd)
oldpar <- par()
par(mfrow=c(1,3)) # 3つのグラフを並べ�?
plot(x.prof, ci=c(0.95, 0.99)) # プロファイル尤度のグラフを描く
par(oldpar)
confint(x.gev.evd, level=0.95) # ワルド�?�区間推定�?0.95のとき�?�levelを省略可
confint(x.prof) # プロファイル区間推�?

fgev(x, prob=0.01) # 100年に1回起こる
fgev(x, shape=0) # グンベルモ�?ルを使�?


# 6.2.4�?
library(ismev)
x.gev.ismev <- gev.fit(x)# GEVモ�?ルで推定す�?
gev.diag(x.gev.ismev)


# 6.3節
library(evd) # evdパッケージをロードす�?
library(ismev) # ismevパッケージをロードす�?
set.seed(100)
x <- rexp(10000)
dim(x) <- c(50,200) # 50年�?で?��年あた�?200個�?��?ータとする
x[1:5,1:5] # 行�?��?�1部がどのようなも�?�か眺める
(xmax <- apply(x, 1, max)) # �?行で最大値を取�?
(x.gev.evd <- fgev(xmax)) # evdによるGEVモ�?ル
oldpar <- par()
par(mfrow=c(2,2)) # 4つのグラフを並べ�?
plot(x.gev.evd) # そ�?�EDA
x.gev.ismev <- gev.fit(xmax) # ismevによるGEVモ�?ル
gev.diag(x.gev.ismev) # そ�?�EDA
xsort <- apply(x, 1, sort, decreasing=T) # �?行を降�??にソートす�?
xsort[1:5,1:5] # xで行だったものがxsortで列になってしまって�?�?
dim(xsort)
xsort2 <- t(xsort) # 転置を取�?
xsort2[1:5, 1:5]
x.rgev <- rlarg.fit(xsort2) # rGEVモ�?ルを当てはめる
rlarg.diag(x.rgev) # xsortの当てはまりを見る　途中で抜け�?
x.5gev <- rlarg.fit(xsort2,5) # 上�?5個とする
rlarg.diag(x.5gev)
error <- rlarg.fit(xsort)
rlarg.diag(error)
par(oldpar)

# 6.4節
library(evd)
venice2
?venice2
dim(venice2)
plot(1887:2011,venice2[,1],type="l") # 折れ線グラフで�?ータの様子を眺める
(venice.gev <- fgev(venice2[,1]))
oldpar <- par()
par(mfrow=c(2,2)) # 4つのグラフを並べ�?
plot(venice.gev)
fgev(venice2[,1], prob=0.01)
library(ismev)
venice.5gev <- rlarg.fit(venice2,5)
rlarg.diag(venice.5gev)
venice.4gev <- rlarg.fit(venice2,4)
rlarg.diag(venice.4gev)
venice.3gev <- rlarg.fit(venice2,3)
rlarg.diag(venice.3gev)
venice.2gev <- rlarg.fit(venice2,2)
rlarg.diag(venice.2gev)


fox
?fox
dim(fox)
(fox.gev <- fgev(fox[,1]))
plot(fox.gev)


data(glass)
glass
(glass.gev <- fgev(glass))
par(mfrow=c(2,2))
plot(glass.gev)

par(oldpar)

# 6.5.1�?
getwd() # 作業�?ィレクトリの確�?
Tokyo <- read.csv("TokyoTemp.csv")
head(Tokyo)
summer.n <- grep("/[789]/", Tokyo[,1]) # 7,8,9月�?�行番号を取り�?��?
length(summer.n); summer.n[1:20]
Tokyo2 <- Tokyo[summer.n,2] # 7,8,9月�?�気温を取り�?��?
head(Tokyo2)
2668/29
dim(Tokyo2) <- c(92,29) # 行�?�に整形する
colnames(Tokyo2) <- 1990:2018 # �?列に名前を付け�?
head(Tokyo2)
Tokyo3 <- apply(Tokyo2, 2, sort, decreasing=T) # �?列でソートす�?
head(Tokyo3)
tTokyo3 <- t(Tokyo3) # 関数の仕様に合わせて転置する
plot(1990:2018, tTokyo3[,1], type="l") # �?年の最高気温の折れ線グラ�?
library(ismev) # ismevパッケージをロードす�?
Tokyo_model_r1 <- rlarg.fit(tTokyo3,1) # r=1のモ�?ルを当て�?
rlarg.diag(Tokyo_model_r1) # モ�?ルの適合性を図示する
Tokyo_model_r2 <- rlarg.fit(tTokyo3,2)
rlarg.diag(Tokyo_model_r2)
library(evd)
fgev(tTokyo3[,1], prob=0.01)


# 6.5.2�?
rain <- read.csv("TateyamaAme.csv")
dim(rain)
rain2 <- na.omit(rain) # 初めの方は�?ータがな�?ようなので削除する
head(rain2)
raints <- ts(rain2[,2],start=c(1968,5),frequency=12) # 時系列�?�形式に直�?
raints
plot(raints)
wind <- read.csv("TateyamaKaze.csv")
dim(wind)
wind2 <- na.omit(wind) # 初めの方は�?ータがな�?ようなので削除する
head(wind2)
windts <- ts(wind2[,2],start=c(1968,5),frequency=12) # 時系列�?�形式に直�?
windts
plot(windts)
head(wind2)	# �?ータの先�?�
tail(wind2)	# �?ータの末尾
dim(wind2)	# �?ータの次�?
618/12
windmat <- matrix(c(rep(NA,4), wind2[,2], rep(NA,2)), 52, 12, byrow=T) # �?ータを行�?�形式にする
rownames(windmat) <- 1968:2019 # 行�?�名前を年の名前にする
(windmax <- apply(windmat, 1, max)) # �?年の最大値を取る　NAができてしま�?
(windmax <- apply(windmat, 1, max, na.rm=T)) # �?年の最大値を取り直す　NAな�?
plot(1968:2019,windmax, type="l")
length(windmax)
windmax18 <- windmax[1:51] # 2018年までの�?ータ
library(evd)
(wind.gev18 <- fgev(windmax18)) # evdによる推�?
oldpar <- par()
par(mfrow=c(2,2))
plot(wind.gev18)
str(wind.gev18)
(temp <- wind.gev18$estimate)
1-pgev(28.4, loc=temp[1], scale=temp[2], shape=temp[3]) # 上�?�確�?
wind3 <- wind2[wind2[,6]==1,]	# 6列目の�?質番号�?1の行だけ取り�?��?
head(wind3); tail(wind3); dim(wind3)
481/12
windmat2 <- matrix(c(rep(NA,4), wind3[,2], rep(NA,7)), 41, 12, byrow=T) # �?ータを行�?�形式にする
rownames(windmat2) <- 1968:2008 # 行�?�名前を年の名前にする
(windmax2 <- apply(windmat2, 1, max, na.rm=T)) # �?年の最大値を取り直す　NAな�?
(wind.gev08 <- fgev(windmax2))
plot(wind.gev08)
(wind.gev18gum <- fgev(windmax18, shape=0))
plot(wind.gev18gum)
temp <- wind.gev18gum$estimate
1-pgev(28.4, loc=temp[1], scale=temp[2], shape=0)
1/(1-pgev(28.4, loc=temp[1], scale=temp[2], shape=0))
par(oldpar)


# 6.5.3�?
olm <- read.csv("results.csv")
dim(olm)
head(olm)# 初めの6行を見る
levels(olm[,2])# 2列目が競技なのでそ�?�要�?の種類を並べ�?
men100 <- olm[olm[,2]=="100M Men",]# 男�?100m走を取り�?�す�?
dim(men100)
str(men100)# �?ータの構�?を見る
min(men100[,8])# 最速�?��?ータを知�?
(c8 <- as.numeric(as.character(men100[,8]))) # �?字�?�型を経て数値型に変換する
min(c8, na.rm=T)# 最速�?��?ータを知�?

men100[,8] <- -c8# 極値統計�?�最大値を求める方法なので最小値は�?転させ�?
men1002 <- men100[order(men100$Year),]# 年代�?に並べ�?
men100G <- men1002[men1002[,5]=="G",]# 金メダル�?け取り�?��?
men100S <- men1002[men1002[,5]=="S",]# 銀メダル�?け取り�?��?
men100B <- men1002[men1002[,5]=="B",]# �?メダル�?け取り�?��?
dim(men100G); dim(men100S); dim(men100B)# �?ータサイズの確認をする
men100B # �?�?け多いのでどこか調べ�?
# 1896年ア�?ネで2人�?と�?かる
(men100mat <- cbind(men100G[,c(4,6,8)], men100S[,c(6,8)], men100B[2:28,c(6,8)]))
# �?行が�?年に当たり�?金銀�?につ�?て氏名と記録を並べた行�?�をつくる �?を一つ外す
men100mat2 <- men100mat[,c(1,3,5,7)] # 年代と金銀�?記録のみの行�?�をつくる
plot(men100mat[,1], men100mat[,3])# 年代と金メダル記録の散�?図をつくる
# 明らかにトレンドがある
library(evd)
fgev(men100mat2[,2])# 定常モ�?ル
(men100gev <- fgev(men100mat2[,2], method="Nelder-Mead"))# ネルダー・ミ�?�ド�?
oldpar <- par()
par(mfrow=c(2,2))
plot(men100gev)
fgev(men100mat2[,2], method="Nelder-Mead", prob=0)
(men100gev.ns <- fgev(men100mat2[,2], nsloc=men100mat2[,1]))# 非定常モ�?ル
plot(men100gev.ns)# 当てはまりが悪�?
(men100gev2.ns <- fgev(men100mat2[,2],method="Nelder-Mead", nsloc=men100mat2[,1]))
plot(men100gev2.ns)
(men100gev3.ns <- fgev(men100mat2[,2], nsloc=men100mat2[,1]/100)) # 別の非定常
plot(men100gev3.ns)# 当てはまりがよい
par(oldpar)

library(ismev) 
men100gev.is <- gev.fit(men100mat2[,2])# 定常モ�?ル
gev.diag(men100gev.is)# 当てはまりがよい
# ismevで非定常モ�?ルを作る
time1 <- (men100mat2[,1]-1896)/120# 1896年から2016年�?0-1とする
time=cbind(time1, time1^2)
men100t1 <- gev.fit(men100mat2[,2], ydat=time, mul=1) # 位置母数につ�?て1次�?
gev.diag(men100t1)
men100t2 <- gev.fit(men100mat2[,2], ydat=time, mul=c(1,2))
gev.diag(men100t2)
men100gev.is$nllh*2+6	# 定常モ�?ルのAIC
men100t1$nllh*2+8		# 1次モ�?ルのAIC
men100t2$nllh*2+10		# 2次モ�?ルのAIC
men100t2$mle			# 2次モ�?ルの推定値
t1 <- (2020-1896)/120
temp <- men100t2$mle	# 何度も書く手間を省く
m100qgev <- function(year, p){
  t1 <- (year-1896)/120
  qgev(p, loc=temp[1]+temp[2]*t1+temp[3]*t1^2, scale=temp[4],
  shape=temp[5])}
m100qgev(2020,c(0.5,0.9,0.99))
plot(2000+4*(1:10), m100qgev(2000+4*(1:10), 0.9))
men100r <- rlarg.fit(men100mat2[,2:4])
rlarg.diag(men100r)
men100rt1 <- rlarg.fit(men100mat2[,2:4], ydat=time, mul=1)
rlarg.diag(men100rt1)
men100rt2 <- rlarg.fit(men100mat2[,2:4], ydat=time, mul=2)
rlarg.diag(men100rt2)
men100mat3 <- abs(100/men100mat2)	# 100/行�?��?��?成�??の行�?�をつくる
men100mat3[,1] <- men100mat2[,1]		# 1列目の年号�?け書き換�?
men100mat3
sp.gev <- gev.fit(men100mat3[,2])
gev.diag(sp.gev)
sp.gev.t1 <- gev.fit(men100mat3[,2], ydat=time, mul=1)
gev.diag(sp.gev.t1)
sp.gev.t2 <- gev.fit(men100mat3[,2], ydat=time, mul=2)
gev.diag(sp.gev.t2)



#----------------------------------------------------------------------------------------------------------------------------------------------------

# 7��

# 7.2.1��
library(evd)
dgpd(3, loc=1, scale=2, shape=0.5)	# GPの�?度関数の値
1/2*(1+0.5*(3-1)/2)^(-1/0.5-1)	# 上�?�値をg(y)に数値を代入して直接計算す�?
pgpd(5,loc=2, scale=1.5, shape=0.2)	# GPの�?�?関数の値
curve(pgpd(x,loc=2, scale=1.5, shape=0.2),xlim=c(1,10))	# �?�?関数のグラ�?
qgpd(0.8, loc=1, scale=2, shape=0.7)	# G(x)=0.8となる点?���??位数?�?
pgpd(6.957627, loc=1, scale=2, shape=0.7)	# 上�?��?位数を確�?
rgpd(5, loc=1, scale=2, shape=0.3)	# GPに従う乱数�?5�?


# 7.2.2�?
library(evd)
set.seed(100)
tmp <- rgpd(10000)	# GP(0,1,0)に従う乱数�?10000�?
hist(tmp)	# 乱数の�?�?の確�?
mrlplot(tmp)		# 標本平�?�?過�?�ロ�?�?
oldpar <- par()
par(mfrow=c(1,2))
tcplot(tmp, tlim=c(0,7))	# GPモ�?ルによる母数推定�?�ロ�?�?
par(mfcol=c(1,3))
tcplot(tmp, model="pp", tlim=c(0,7)) # PPモ�?ルによる母数推定�?�ロ�?�?
sum(tmp>7)	# 7を�?える乱数の数
(tmp.g1 <- fpot(tmp, 1)) # GPモ�?ル u=1
(tmp.g7 <- fpot(tmp, 7)) # GPモ�?ル u=7
(tmp.p1 <- fpot(tmp, model="pp", 1)) # PPモ�?ル u=1
(tmp.p7 <- fpot(tmp, model="pp", 7)) # PPモ�?ル u=7
(tmp.p7 <- fpot(tmp, model="pp", 7, method="Nelder-Mead")) # PPモ�?ル u=7
tmp.g1.pro <- profile(tmp.g1)	# GP u=1のプロファイル
plot(tmp.g1.pro)				# 上�?�プロファイルのプロ�?�?
tmp.g7.pro <- profile(tmp.g7)	# GP u=7のプロファイル　計算不可 
confint(tmp.g1)				# GP u=1のワルド区間推�?
confint(tmp.g7)				# GP u=7のワルド区間推�?
confint(tmp.p1)				# PP u=1のワルド区間推�?
confint(tmp.g1.pro)			# GP u=1のプロファイル区間推�?
par(mfrow=c(2,2))
plot(tmp.g1) # GPモ�?ル u=1 のEDA
plot(tmp.g7) # GPモ�?ル u=7 のEDA
plot(tmp.p1) # PPモ�?ル u=1 のEDA
par(oldpar)

# 7.2.4�?
library(evd)
library(ismev)
set.seed(100)
tmp <- rgpd(10000)	# GP(0,1,0)に従う乱数�?10000�?
hist(tmp)	# 乱数の�?�?の確�?
mrl.plot(tmp) # 標本平�?�?過�?�ロ�?�?
gpd.fitrange(tmp, umin=1, umax=8)	# GPモ�?ルによる母数推定�?�ロ�?�?
pp.fitrange(tmp, umin=1, umax=8)	# PPモ�?ルによる母数推定�?�ロ�?�?
tmp.g1.is <- gpd.fit(tmp,1) # GPモ�?ル u=1
tmp.g7.is <- gpd.fit(tmp,7) # GPモ�?ル u=7
tmp.p1.is <- pp.fit(tmp,1) # PPモ�?ル u=1
tmp.p7.is <- pp.fit(tmp,7) # PPモ�?ル u=7
gpd.diag(tmp.g1.is)	# GPモ�?ル u=1 のEDA
gpd.diag(tmp.g7.is)	# GPモ�?ル u=7 のEDA
pp.diag(tmp.p1.is)		# PPモ�?ル u=1 のEDA
pp.diag(tmp.p7.is)		# PPモ�?ル u=7 のEDA


# 7.3節
library(evd)
set.seed(100)
tmpn <- rnorm(2000)	# 正規乱数2000�?
hist(tmpn)	# 乱数の�?�?の確�?
mrlplot(tmpn)		# 標本平�?�?過�?�ロ�?�?
max(tmpn)	# 乱数の最大値
mrlplot(tmpn, tlim=c(-1,3))		# �?囲を指定す�?
oldpar <- par()
par(mfrow=c(2,1))
tcplot(tmpn, tlim=c(-1,3))	# GPモ�?ルによる区間推�?
par(mfrow=c(3,1))
tcplot(tmpn, model="pp", tlim=c(-1,3))	# PPモ�?ルによる区間推�?
tcplot(tmpn, model="pp", tlim=c(-1,2))	# PPモ�?ルによる区間推�? �?囲を取り直�?  
sum(tmpn>1.5)	# 1.5を�?える�?ータサイズ
(tmpn.pot <- fpot(tmp, 1.5))	# GPモ�?ル u=1.5での推�?
tmpn.prof <- profile(tmpn.pot)	# そ�?�プロファイル
confint(tmpn.pot)			# ワルド信頼区�?
par(mfrow=c(2,2))
plot(tmpn.pot)				# GPモ�?ルのEDA
(tmpn.pot.p <- fpot(tmpn, model="pp", 1.5))	# PPモ�?ル u=1.5での推�?
confint(tmpn.pot.p)						# ワルド区間推�?
plot(tmpn.pot.p)				# PPモ�?ルのEDA
set.seed(100)
tmpe <- rexp(10000)	# �?数乱数10000�?
par(oldpar)
mrlplot(tmpe)		# 標本平�?�?過�?�ロ�?�?
mrlplot(tmpe, tlim=c(3,8))	# 区間を取り直�?
sum(tmpe>4)			# 4を�?える�?ータの数
(tmpe.pot <- fpot(tmpe,4))		# GPモ�?ル u=4で推�?
confint(tmpe.pot)			# そ�?�ワルド信頼区�?
par(mfrow=c(2,2))
plot(tmpe.pot)				# GPモ�?ルのEDA
(tmpe.pot.p <- fpot(tmpe,model="pp",4))		# PPモ�?ル u=4で推�?
(tmpe.pot.p <- fpot(tmpe,model="pp",3))		# PPモ�?ルでu=3に取り直�?
confint(tmpe.pot.p)						# そ�?�ワルド信頼区�?
plot(tmpe.pot.p)					# PPモ�?ルのEDA
sum(tmpe>3)					# 3を�?える�?ータの数
gpd.fit(tmpe,4)					# ismevパッケージのGPモ�?ル
pp.fit(tmpe,3)					# ismevパッケージのPPモ�?ル
par(oldpar)


# 7.4節
library(evd)
library(ismev)
venice2
str(venice2)
head(venice2); tail(venice2)	# はじめと終わり�?��?ータを見る
ve <- c(as.matrix(venice2)) # 行�?�にした後にベクトルにする
plot(ve)
vets <- ts(ve, start=c(1887,1), end=c(2011,12), frequency=12) # 時系列にする
plot(vets) # 時系列�?�プロ�?�?
mrlplot(ve)				# 標本平�?�?過�?�ロ�?�?
tcplot(ve, tlim=c(90,140))	# 母数の推定�?�ロ�?�?
vern <- na.omit(ve)		# NAを取り除�?
sum(vern>110); sum(vern>120)	# u=110, 120を�?える�?ータサイズ
(ve.gp <- fpot(vern, 110))		# evdでGPモ�?ルによる推�?
oldpar <- par()
par(mfrow=c(2,2))
plot(ve.gp)					# EDA
(ve.pp <- fpot(vern, model="pp", 110))	# evdでPPモ�?ルによる推�?
ve.gp2 <- gpd.fit(vern, 110)	# ismevでGPモ�?ルによる推�?
gpd.diag(ve.gp2)				# EDA
ve.pp2 <- pp.fit(vern, 110)		# ismevでPPモ�?ルによる推�?
pp.diag(ve.pp2)				# EDA
plot(fpot(vern, 110, npp=12))		# 1年12月として推定�?�EDA
fpot(vern, 110, npp=12, mper=100)	# 100年再現レベル
par(oldpar)


# 7.5節
library(evd)
wind <- read.csv("TateyamaKaze.csv")	# �?ータ読み込み
wind
wind2 <- na.omit(wind) # NAの削除
head(wind2)
windts <- ts(wind2[,2], start=c(1968,5), frequency=12) # 時系列に変換
plot(windts)		# 時系列�?�グラ�?
max(wind2[,2])	# 最大の風速�?��?くら�?
tail(wind2); length(wind2[,2])	# �?ータの末尾と�?ータの長さを見る
wind3 <- wind2[1:608,2]		# 風速�?�列をとり、最後�?�10か月�?を除�?
mrlplot(wind3)				# 標本平�?�?過�?�ロ�?�?				
sum(wind3 >14); sum(wind3 >12)	# u=14, 12以上�?��?ータサイズ
tcplot(wind3, tlim=c(10,20))	# 推定�?�ロ�?�?
sum(wind3 >13)				# 13を�?える�?ータサイズ
(tate.gp <- fpot(wind3, 13, npp=12))	# 1年12月としてGPモ�?ルで推�?
oldpar <- par()
par(mfrow=c(2,2))
plot(tate.gp)						# EDA
(tate.pp <- fpot(wind3, 13, model="pp"))		# PPモ�?ルで推�?
plot(tate.pp)						# EDA
(tate.gp2 <- fpot(wind3, 13, npp=12, mper=100))	# 100年再現レベル
plot(tate.gp2)								# EDA
(temp <- tate.gp$estimate)	# 繰り返しの面倒を避ける
qgpd(28.4, scale=temp[1], shape=temp[2])	# �?位数
pgpd(28.4, scale=temp[1], shape=temp[2])	# �?�?関数
par(oldpar)


#-------------------------------------------------------------------------------------------------------------------------------------------
  
  # 8�?


# 8.4.1�?
library(evd)
portpirie
plot(portpirie, type="n")	# プロ�?ト�?��?�?�?
length(portpirie)
text(1:65, portpirie, paste(1:65))		# 番号?��時刻?��でプロ�?トす�?
abline(h=4.2)
clusters(portpirie, u=4.2, r=3) # uとrのみ設�?
clusters(portpirie, u=4.2, r=3, ulow=4, rlow=2) # ulowとrlowも設�?
clusters(portpirie, u=4.2, r=3, cmax = TRUE) # 群れ�?�最大値を取�?
clusters(portpirie, u=4.2, r=3, ulow=3.8, plot = TRUE) # 群れ�?�プロ�?トをする
tvu <- c(rep(4.2, 20), rep(4.1, 25), rep(4.2, 20))
clusters(portpirie, tvu, 3, plot = TRUE) # 閾値を時刻によって変え�?
exi(portpirie, u=4.2, r = 3, ulow = 3.8)
tvu <- c(rep(4.2, 20), rep(4.1, 25), rep(4.2, 20))
exi(portpirie, u=tvu, r = 1)
exi(portpirie, u=tvu, r = 0)
exiplot(portpirie, tlim=c(2,4.3))
exiplot(portpirie, tlim=c(2,4.3), r=3, add=T)

# 8.4.2�?
fpot(portpirie, 4.2)	# 前�?と同�?
fpot(portpirie, 4.2, r=3)	# cmax=Fなので上と同じ
fpot(portpirie, 4.2, r=3, cmax=T)		# 計算できな�?
fpot(portpirie, 4.2, r=3, cmax=T, method="Nelder-Mead")
# 群れ�?�最大値のみで推�?


# 8.4.3�?
set.seed(100)
marma(5, p=1, q=2, psi=0.2, theta=c(0.3,0.2))
mar(5, psi=0.2, init=1)
mma(5, theta=0.5, rand.gen=rexp)


# 8.5節
library(evd)
set.seed(100)
temp <- marma(500, p=1, q=2, psi=0.2, theta=c(0.3,0.2))	# MARMA時系列を発�?
plot(temp)	# 時系列�?�実現値を�?�ロ�?�?
plot(temp, type="l")		# 表示方法を変え�?
plot(temp, ylim=c(0,100))	# �?囲を変え�?
plot(temp, ylim=c(0,50))	# �?囲を変え�?
mrlplot(temp)			# 標本平�?�?過�?�ロ�?�?
tcplot(temp, c(10,30))		# 推定�?�ロ�?�?
sum(temp>13)
clusters(temp, u=13, r=3) 	# 群れを調べ�?
exi(temp, u=13, r=3)		# 極値収縮�?の推�?
fpot(temp, 13)				# �?過データをすべて使�?
fpot(temp, 13, r=3, cmax=T)	# 群れ�?�最大のみ使�?
fpot(clusters(temp, u=13, r=3, cmax=T),13) # 群れ�?�最大ベクトルをすべて使�?
fpot(temp, 13, r=3, model="pp",cmax=T)　# 群れ�?�最大のみでppモ�?ル
dim(temp) <- c(25,20) # ブロ�?ク最大で調べる準備　20ブロ�?クにする
temp2 <- apply(temp, 2, max) # �?ブロ�?クの最大を取�?
fgev(temp2) # ブロ�?ク最大のGEVモ�?ル


# 8.6節
library(evd)
library(ismev)
data(dowjones)	# �?ータをロードす�?
dowjones
plot(dowjones, type="l")	# 時系列として�?ータをグラフ化する
dj <- dowjones[,2]
dldj <- diff(log(dj))	# 対数差�?を取�?
plot(dowjones[2:1304,1], dldj, type="l")	# 対数差�?のグラ�?
max(dldj)
-min(dldj)
length(dldj)
dldjm <- matrix(dldj[1:1300], 20, 65) # 後ろ3つを捨て20*65型行�?�にする
(dldjp <- apply(dldjm, 2, max)) # �?列�?�最大値を取�?
(dldjp.gev.e <- fgev(dldjp)) # evdでブロ�?ク最大GEVモ�?ル
oldpar <- par()
par(mfrow=c(2,2))
plot(dldjp.gev.e)	# そ�?�EDA
dldjp.gev.i <-gev.fit(dldjp) # ismevでブロ�?ク最大GEVモ�?ル
gev.diag(dldjp.gev.i) # そ�?�EDA
(dldjn <- -apply(dldjm, 2, min)) # �?列�?�最小値を取り正�?を�?�れ替える
(dldjn.gev.e <- fgev(dldjn)) 	# evdでブロ�?ク最大GEVモ�?ル
plot(dldjn.gev.e)				# そ�?�EDA
dldjn.gev.i <- gev.fit(dldjn)		# ismevでブロ�?ク最大GEVモ�?ル
gev.diag(dldjn.gev.i)			# そ�?�EDA
par(oldpar)
mrlplot(dldj)					# 標本平�?�?過�?�ロ�?�?
mrlplot(dldj, tlim=c(0,0.04))
mrlplot(dldj, tlim=c(0.005,0.03))
sum(dldj>0.008)
tcplot(dldj, tlim=c(0.005,0.03))	# 推定�?�ロ�?�?
tcplot(dldj, tlim=c(0.005,0.03), method="Nelder-Mead")
sum(dldj>0.016)
(dldjp.pot <- fpot(dldj, 0.016, r=3, cmax=T))	# 母数の推�?
par(mfrow=c(2,2))
plot(dldjp.pot)							# EDA
(dldjp.pot <- fpot(dldj, 0.016, r=3, cmax=T, method="Nelder-Mead"))
plot(dldjp.pot)
par(oldpar)
mrlplot(-dldj)					# 標本平�?�?過�?�ロ�?�?
mrlplot(-dldj, tlim=c(0,0.04))
mrlplot(-dldj, tlim=c(0.015,0.03))
tcplot(-dldj, tlim=c(0.015,0.03))		# 母数の推定�?�ロ�?�?
sum(-dldj>0.018)
(dldjn.pot <- fpot(-dldj, 0.018, r=3, cmax=T, method="Nelder-Mead"))
length(grep("1999", dowjones[,1])) # 1年は何日か�?
fpot(-dldj, 0.016, r=3, cmax=T, npp=261, mper=10, method="Nelder-Mead")	# 母数の推�?
max(dldjn)
fgev(dldjn, prob=1/3650, method="Nelder-Mead")		# �?位点の推�?
# 計算できなかっ�?
(temp <- fgev(dldjn, method="Nelder-Mead")$estimate)	# 入力�?�省略のため
qgev(1-1/3650, loc=temp[1], scale=temp[2], shape=temp[3])		# �?位点の推�?


# 8.7節
nh <- read.csv("NikkeiHeikin.csv")
dim(nh)
head(nh)
plot(nh[,1], nh[,5], type="l")		# 事前に�?ータを眺める
dlnh <- diff(log(nh[,5]))			# 対数差�?
plot(nh[2:848,1], dlnh, type="l")		# 対数差�?を眺める
max(dlnh)
max(-dlnh)
nh2 <- cbind(nh, c(0,dlnh))
head(nh2)
nh3 <- nh2[9:848,]	# 初めの年の5月か�?12月を捨て�?
head(nh3)
tail(nh3)
nh.mat <- matrix(nh3[,7], 12,70)	# 1�?12ヶ�?70年�?とする
(nh.pm <- apply(nh.mat, 2, max))	# 毎年の正の最大
(nh.nm <- -apply(nh.mat, 2, min))	# 毎年の�?の最大
max(nh.pm)
max(nh.nm)
(nhp.gev <- fgev(nh.pm))	# 正につ�?てGEVモ�?ル
oldpar <- par()
par(mfrow=c(2,2))
plot(nhp.gev)			# EDA
(nhn.gev <- fgev(nh.nm))	# �?につ�?てGEVモ�?ル
plot(nhn.gev)			# EDA
mrlplot(dlnh, tlim=c(0,0.2))		# 標本平�?�?過�?�ロ�?�?
mrlplot(dlnh, tlim=c(0.03,0.15))	# �?囲を限定して詳しく見る
tcplot(dlnh, tlim=c(0.03,0.15))	# 推定�?�ロ�?�?
sum(dlnh>0.05)				# 標本サイズ
(nhp.pot1 <- fpot(dlnh, 0.05))	# GPモ�?ル
(nhp.pot1 <- fpot(dlnh, 0.05, method="Nelder-Mead"))	# 数値計算を変え�?
plot(nhp.pot1)	
(nhp.pot2 <- fpot(dlnh, 0.05, r=2, cmax=T,method="Nelder-Mead")) # 群れを入れる
plot(nhp.pot2)
(nhp.pot3 <- fpot(dlnh, 0.05, r=3, cmax=T,method="Nelder-Mead")) # 別の群�?	
plot(nhp.pot3)
mrlplot(-dlnh, tlim=c(0,0.2)) 	# �?の�?ータの標本平�?�?過�?�ロ�?�?
mrlplot(-dlnh, tlim=c(0.05,0.13)) # �?囲を限定して詳しく見る
tcplot(-dlnh, tlim=c(0.05,0.13))	# 推定�?�ロ�?�?
sum(-dlnh>0.08)				# 標本サイズ
(nhn.pot1 <- fpot(-dlnh, 0.08))				# GPモ�?ル
plot(nhn.pot1)							# EDA
(nhn.pot2 <- fpot(-dlnh, 0.08, r=2, cmax=T))	# 群れを入れる
plot(nhn.pot2)							# EDA
(nhn.pot3 <- fpot(-dlnh, 0.08, r=3, cmax=T))	# 別の群�?
plot(nhn.pot3)							# EDA
fpot(-dlnh, 0.08, r=3, cmax=T, npp=12, mper=10)	# 10年で最大下落
par(oldpar)


#----------------------------------------------------------------------------------------------------------------------------------------------
  
  # 9�?


# 9.4.1�?
t <- 2:4
dorder(2:4, distn = "norm", mean = 0.5, sd = 1.2, mlen = 5, j = 2)
4* choose(5,4)*dnorm(t,0.5,1.2)*pnorm(t,0.5,1.2)^3*(1-pnorm(t,0.5,1.2))
porder(3, distn = "exp", rate = 1.2, mlen = 3, j = 2, largest=F)
sum(choose(3, 2:3) * pexp(3, rate=1.2)^(2:3) *(1-pexp(3, rate=1.2))^(3-2:3) )
rorder(5, distn="gamma", shape = 1, mlen = 10, j = 2)
dextreme(2:4, distn = "norm", mean = 0.5, sd = 1.2, mlen = 5)
dorder(2:4, distn = "norm", mean = 0.5, sd = 1.2, mlen = 5, j=1)
pextreme(2:4, distn = "exp", rate = 1.2, mlen = 2)
porder(2:4, distn = "exp", rate = 1.2, mlen = 2, j=1)
qextreme(seq(0.9, 0.6, -0.1), distn = "exp", rate = 1.2, mlen = 2)



# 9.4.2�?
library(evd)
n <- 10; m <- 5
os <- rorder(1000, distn="exp", mlen=n, j=m, largest=F)	# �?序統計量の乱数
oldpar <- par()
par(mfrow=c(2,1))
hist(os)
xmat <- rexp(m*1000)	# m*1000個�?��?数乱数
dim(xmat) <- c(m,1000)	# m*1000型行�?�に成型
x1 <- xmat/((n-m+1):10)	# �?数乱数の割り�?
x2 <- apply(x1, 2, sum)	# �?列�?��?
hist(x2)
(x <- matrix(12, 4,4))
x/(1:4)
x/(1:5)		# 警告が出�?
par(oldpar)


# 9.4.3�?
eu <- function(n){sum(1/(1:n))-log(n)}
eu(10^8)
pi6 <- function(n){sqrt(6*sum(1/(1:n)^2))}
pi6(10^8)



# 9.5.2�?
library(evd)
# �?数�?布Fにつ�?てF^n(an x + bn)を与える関数
pnexp <- function(x,n,lambda=1){
  bn <- -1/lambda*(log(1-exp(-1/n)))
  an <- 1/(n*lambda*(exp(1/n)-1))
  pextreme(an*x+bn, distn="exp", rate=lambda, mlen=n)
}
n<- 100; x<-1:4; lambda<-2
pnexp(x,n,lambda)	# x=1:4, n=100, 平�?2として計�?
pgumbel(x)		# 収束先�?��?�?関数
# 両�?の差と上から�?�評価の比�?
abs(pnexp(x,n,lambda)-pgumbel(x)); 1/(exp(1)*(2*n-1))
# �?数�?�?の最大値のグンベル�?�?への収束
curve(pnexp(x,10,lambda), xlim=c(-3,3))
curve(pnexp(x,100,lambda), add=T)
curve(pnexp(x,1000,lambda), add=T)
curve(pgumbel, add=T)
# 両�?の差のグラ�?
curve(pnexp(x,10,lambda)-pgumbel(x), xlim=c(-3,3))
curve(pnexp(x,100,lambda)-pgumbel(x), lty=2, add=T)
curve(pnexp(x,1000,lambda)-pgumbel(x), lty=3, add=T)


#-----------------------------------------------------------------------------------------------------------------------------------------------------------
  
  # 10�?


# 10.3.2�?
# install.packages("tidyverse")
library(tidyverse)


# 10.3.3�?
?iris
head(iris)
library(ggplot2)
ggplot(data=iris) + 
  geom_point(mapping=aes(x=Sepal.Length, y=Sepal.Width)) # 散�?図
ggplot(data=iris) + geom_point(mapping=aes(x=Sepal.Length,
   y=Sepal.Width, color=Species)) # 種で色�?�?
ggplot(data=iris) + 
  geom_point(mapping=aes(x=Sepal.Length, y=Sepal.Width))+
  facet_wrap(~Species, nrow=2) # 種で散�?図を�??ける
ggplot(data=iris) + 
  geom_point(mapping=aes(x=Sepal.Length, y=Sepal.Width))+
  geom_smooth(mapping=aes(x=Sepal.Length, y=Sepal.Width)) # 平滑化する


# 10.3.4�?
library(dplyr)
head(trees)
filter(trees, Girth>10, Volume<20)
arrange(trees, Volume, Height) # Volumeを第1基準にHeightを第2基準で�?�?に並べ�?
select(trees, Girth, Volume)
mutate(trees, Volume/Height) 
summarize(group_by(iris, Species), mean(Sepal.Length)) # 種別に平�?をと�?
trees %>% filter(Girth>17) %>% select(Height) # パイプを使�?


#-----------------------------------------------------------------------------------------
  
  
  
  # 表

# anova
library(evd)
x <- rgev(100)
m1 <- fgev(x)
m2 <- fgev(x, shape=0)
anova(m1,m2)


# mtransform
library(evd)
x <- 1:10
p <- c(1,2,3)
mtransform(x,p)-log(pgev(x,loc=p[1], scale=p[2], shape=p[3]))

