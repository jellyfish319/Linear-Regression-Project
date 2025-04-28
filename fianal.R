library(tidyverse)
library(regbook)
library(psych)
library(gridExtra)

rm(list = ls())

select <- dplyr::select


earth <- read.csv("project/earthquake.csv", header = T)

tibble(earth)

earth %>% select(sig,magnitude,cdi,mmi,alert,nst,dmin,gap,depth,tsunami) -> eq01

summary(earth)
dim(earth)


eq01 %>% select(-c(alert)) -> eq02

pairs(eq02) # 산점도행렬

eq02 %>% mutate(dmin = log(dmin), gap = log(gap), depth = log(depth)) -> eq03

pairs(eq03) # dmin, gap, depth 로그변환 후 산점도 행렬

eq01 %>% select(sig,alert) %>% mutate(alert = ifelse(alert == "", "non", alert)) -> eq04

eq04 %>% ggplot() +
  aes(x = seq_along(sig), y = sig, group = alert) +
  geom_point() +
  facet_grid(~alert) +
  xlab("Index")# alert에 따른 중요도, 결측치 제거

dim(earth)

eq01 %>% select(sig,tsunami) -> eq05


eq05 %>% ggplot() +
  aes(x = tsunami, y = sig) +
  geom_point()

# EDA : 변수변환

eq01 %>% filter(!(alert == "")) %>% mutate(dmin = log(dmin), gap = log(gap), depth = log(depth)) -> eq1

dim(eq1)

# 가변수화 : red vs ~red, red vs orange vs yellow, red vs orange vs green + yellow

eq1 %>% mutate(red = ifelse(alert == "red",1,0),
               orange = ifelse(alert == "orange",1,0),
               yellow = ifelse(alert == "yellow", 1,0)
               ) %>% select(-alert) -> eq11

# 초기모형 : 수치형 변수 + 쓰나미 + red + orange + yellow

# Inf 처리

lapply(eq11, function(x){sum(is.infinite(x))})

nrow(eq11)

# Inf 제거

eq11 %>% filter(!is.infinite(gap) & !is.infinite(dmin)) -> eq12

dim(eq12)

names(eq12) # 초기모형의 설명변수

lm.eq <- lm(sig ~ . , eq12)

vif.lm <- vif(mod = lm.eq)

summary(vif.lm)


par(mfrow = c(2,2))
plot(lm.eq)
par(mfrow = c(1,1))

boxcox(lm.eq, lambda = seq(-4,-1,1/10))
hist(resid(lm.eq), main = "Histogram of Residual", xlab = "Residual")

# lambda = -2

eq12 %>% mutate(sig = (sig^(-2)-1)/(-2)) -> eq13

lm.eq <- lm(sig ~ . , eq13)

par(mfrow = c(2,2))
plot(lm.eq)
par(mfrow = c(1,1))

hist(resid(lm.eq),main = "Histogram of Residual (box-cox)", xlab = "Residual (box-cox)")

hist(eq12$sig, main = "Histogram of sig", xlab = "sig")
hist(eq13$sig, main = "Histogram of sig (box-cox)", xlab = "sig (tranformed)")



# 정규성 만족 X

# 로그변환

eq12 %>% mutate(sig = log(sig)) -> eq14

lm.eq <- lm(sig ~ . , eq14)

par(mfrow = c(2,2))
plot(lm.eq)
par(mfrow = c(1,1))

hist(resid(lm.eq), main = "Histogram of Residual (log)", xlab = "Residual (log)")

hist(eq12$sig, main = "Histogram of sig", xlab = "sig")
hist(eq13$sig, main = "Histogram of sig (box-cox)", xlab = "sig (tranformed)")

# 모형선택

eq14

step(lm.eq, direction = "back")


reg.eq <- summaryf(regsubsets(sig ~ . , eq14,nvmax = 11, nbest = 2))

# adj r^2

which(reg.eq$adjr2 == max(reg.eq$adjr2))

reg.eq$adjr2[which(reg.eq$adjr2 == max(reg.eq$adjr2))]

# mellows cp

abs(reg.eq$cp - (24-11))

abs(reg.eq$cp - (24-11))

# BIC

which(reg.eq$bic == min(reg.eq$bic))

reg.eq$bic[which(reg.eq$adjr2 == max(reg.eq$adjr2))]

# mellows cp & R^2 AIC : mag, cdi, dmic, gap, tsunami, alert

# bic : 위에서 gap tsunammi 빠짐

# 후진제거법
step(lm.eq, test = "F")

# 전진선택법
step(lm(sig ~ 1, eq14), scope = ~ magnitude + cdi + mmi + nst + dmin + gap + depth + tsunami + red + yellow + orange,direction = "forward", test = "F")

# 단계적회귀
step(lm(sig ~ 1, eq14), scope = ~ magnitude + cdi + mmi + nst + dmin + gap + depth + tsunami + red + yellow + orange,direction = "both", test = "F")

# press

pressp <- press(regsubsets(sig ~ . , eq14, nvmax = 11))

which(pressp$PRESS == min(pressp$PRESS))

lm.eq <- update(lm.eq, . ~ . -mmi -depth -nst )

summary(lm.eq)



lm.eq1 <- lm(sig ~ 1, eq14)

# p-value 이용 시 tsunami만 빠진 모형

step(lm.eq1, scope = ~ magnitude + cdi + mmi + nst + dmin + gap + depth + tsunami + red + yellow + orange,direction = "forward", test = "F")

# p-value 이용 시 tsunami만 빠진 모형

step(lm.eq1, scope = ~ magnitude + cdi + mmi + nst + dmin + gap + depth + tsunami + red + yellow + orange,direction = "both", test = "F")

# R^2, mellow cp, AIC : 15번째
# BIC : gap, tsunami 삭제
# p-value : tsunami 만 삭제

# tsunami 제거

eq14 %>% select(-tsunami) -> eq15

lm.eq.t <- lm(sig ~ ., eq15)

reg.eq <- summaryf(regsubsets(sig ~ .,eq15, nbest = 2))

# adj r^2

which(reg.eq$adjr2 == max(reg.eq$adjr2))

reg.eq$adjr2[15]

# mellows cp

which(reg.eq$cp == min(reg.eq$cp))

reg.eq$cp[which(reg.eq$cp == max(reg.eq$cp))]

# BIC

which(reg.eq$bic == min(reg.eq$bic))

reg.eq$bic[which(reg.eq$adjr2 == max(reg.eq$adjr2))]


# tsunami는 제거
# tsunami만 없는 15모형 나머지
# tsunmai, gap이 없는 15모형 BIC



lm.eq.t <- update(lm.eq.t, . ~ . -mmi -nst -depth)
lm.eq.tg <- update(lm.eq.t, . ~ . -mmi -nst -depth -gap)

anova(lm.eq.tg, lm.eq.t)
drop1(lm.eq.t)

# 최종모형 : tsunami만 뺀 15모형

# 편회귀그림

eq16 <- eq15 %>% select(-c(mmi, nst,depth))

par(mfrow = c(3,3))

# mag

y.mag <- resid(lm(sig ~ . - magnitude, eq16))
x.mag <- resid(lm(magnitude ~ . - sig, eq16))

plot(y.mag ~ x.mag, main = "Partial regression plot (mag)")
abline(lm(y.mag ~ x.mag), col = "red")
grid()

# cdi

y.cdi <- resid(lm(sig ~ . -cdi, eq16))
x.cdi <- resid(lm(cdi ~ . - sig, eq16))

plot(y.cdi ~ x.cdi, main = "Partial regression plot (cdi)")
abline(lm(y.cdi ~ x.cdi), col = "red")
grid()

# dmin

y.dmin <- resid(lm(sig ~ . -dmin, eq16))
x.dmin <- resid(lm(dmin ~ . -sig, eq16))

plot(y.dmin ~ x.dmin, main = "Partial regression plot (dmin)")
abline(lm(y.dmin ~ x.dmin), col = "red")
grid()

# gap

y.gap <- resid(lm(sig ~ . -gap, eq16))
x.gap <- resid(lm(gap ~ . - sig, eq16))

plot(y.gap ~ x.gap, main = "Partial regression plot (gap)")
abline(lm(y.gap ~ x.gap), col = "red")
grid()

# red

y.red <- resid(lm(sig ~ . -red, eq16))
x.red <- resid(lm(red ~ . - sig, eq16))

plot(y.red ~ x.red, main = "Partial regression plot (red)")
abline(lm(y.red ~ x.red), col = "red")
grid()

# orange

y.orange <- resid(lm(sig ~ . -orange, eq16))
x.orange <- resid(lm(orange ~ . -sig, eq16))

plot(y.orange ~ x.orange, main = "Partial regression plot (orange)")
abline(lm(y.orange ~ x.orange), col = "red")
grid()

# yellow

y.yellow <- resid(lm(sig ~ . -yellow, eq16))
x.yellow <- resid(lm(yellow ~ . -sig, eq16))

plot(y.yellow ~ x.yellow, main = "Partial regression plot (yellow)")
abline(lm(y.yellow ~ x.yellow), col = "red")
grid()

par(mfrow = c(1,1))

drop1(lm.eq.t)

stdcoef(lm.eq.t)

summaryf(regsubsets(sig ~ . , eq14,nvmax = 11))

# 최종모형

# gap, tsunami를 제외

lm.eq.t <- lm.eq.tg

# 이상치

inf.eq <- influence.measures(lm.eq.t)

n <- nrow(lm.eq.t$model)
p <- ncol(lm.eq.t$model)

# 쿡의 거리

plot(cooks.distance(lm.eq.t), main= "Cook's Distance", ylab = "cook's distance")
abline(h = 3.67/(n-p), col = "red", lty = 3)

which(as.vector(cooks.distance(lm.eq.t)) >= 3.67/(n-p))
sum(as.vector(cooks.distance(lm.eq.t)) >= 3.67/(n-p))

# DFFITS

plot(dffits(lm.eq.t), main = "DFFITS", ylab = 'dffits')
abline(h = 2*sqrt(p/n), col = "red", lty = 3)

which(as.vector(dffits(lm.eq.t)) > 2*sqrt(p/n))
sum(as.vector(dffits(lm.eq.t)) > 2*sqrt(p/n))

# COVRATIO

plot(covratio(lm.eq.t), main = "COVRATIO", ylab = "covraio")
abline(h = 3*p/n + 1, col ="red", lty = 3)
abline(h = -3*p/n + 1, col ="red", lty = 3)

which(as.vector(abs(covratio(lm.eq.t)-1)) > 3*p/n)
sum(as.vector(abs(covratio(lm.eq.t)-1)) > 3*p/n)

# DFBETAS

par(mfrow = c(3,2))
for (i in 2:7) {
  plot(dfbetas(lm.eq.t)[,i], ylab= names(lm.eq.t$model)[i], xlab = "Index", main = paste0("DFBETAS : ",names(lm.eq.t$model)[i]))
  abline(h = 2/sqrt(n), col = "red", lty = 3)
}
par(mfrow = c(1,1))

apply(dfbetas(lm.eq.t),2,function(x){which(as.vector(x) > 2/sqrt(n))})
apply(dfbetas(lm.eq.t),2,function(x){sum(as.vector(x) > 2/sqrt(n))})

# 지렛값

plot(hatvalues(lm.eq.t), main = "Leverage Value", ylab = "leverage value")
abline(h = 2*p/n, col ="red", lty = 3)

which(hatvalues(lm.eq.t) > 2*p/n)
sum(hatvalues(lm.eq.t) > 2*p/n)


# 내표준화잔차

plot(rstandard(lm.eq.t), main = "Internally Studentized Residual", ylab = "residual")
abline(h = 2, lty = 3, col = "red")
abline(h = -2, lty = 3, col = "red")

which(abs(rstandard(lm.eq.t)) > 2)
sum(abs(rstandard(lm.eq.t)) > 2)

# 외표준화잔차

plot(rstudent(lm.eq.t), main = "Externally Studentized Residual", ylab = "residual")
abline(h = 2, lty = 3, col = "red")
abline(h = -2, lty = 3, col = "red")

which(abs(rstudent(lm.eq.t)) > 2)
sum(abs(rstudent(lm.eq.t)) > 2)

# 이상치행렬

inf <- data.frame(cook = as.numeric(as.vector(cooks.distance(lm.eq.t)) >= 3.67/(n-p)),
                  dffits = as.numeric(as.vector(dffits(lm.eq.t)) > 2*sqrt(p/n)))

dfbetas <- apply(dfbetas(lm.eq.t),2,function(x){as.numeric(as.vector(x) > 2/sqrt(n))})
dfbetas <- as.data.frame(dfbetas)[,-1]

inf <- cbind(inf, dfbetas)

inf %>% mutate(covratio = as.numeric(as.vector(abs(covratio(lm.eq.t)-1)) > 3*p/n),
               leverage = as.numeric(hatvalues(lm.eq.t) > 2*p/n),
               rstandard = as.numeric(abs(rstandard(lm.eq.t)) > 2),
               rstudent = as.numeric(abs(rstudent(lm.eq.t)) > 2)) -> inf

inf$index <- row.names(inf)

inf %>% filter((cook != 0)|(dffits != 0)|( magnitude != 0 )|( cdi != 0 )|( dmin != 0 )|( red != 0 )| orange != 0 | yellow != 0 | covratio != 0 | leverage != 0 | rstandard != 0 | rstudent != 0) -> inf

inf %>% pivot_longer(cols = cook:rstudent, names_to = "influence", values_to = "value") -> inf


# 잔차분석

par(mfrow = c(2,2))
plot(lm.eq.t)
par(mfrow = c(1,1))

# 표준화잔차 vs 설명변수 : 설명변수는 상수 & 현 모형은 옳다, plot(rstandard(lm) ~ 설명변수)로 그림
# 표준화잔차 vs 관측순서 : 오차항에 대한 가정 epsilone ~iid N(0,sigma^2), plot(rstandard(lm))으로 그림
# 관측값 vs 적합값 : plot(관측값 ~ fitted(lm), data)
# R-F 그림 : rfs(lm) : 


# 잔차 vs 설명변수

res.plot <- list()

rstandard <- rstandard(lm.eq.t)

model <- lm.eq.t$model
model$rstandard <- rstandard
model

for (i in 3) {
  res.plot[[i-1]] <- ggplot(model, aes(x = model[,i], y = rstandard)) +
    xlab(names(lm.eq.t$coefficients)[i]) +
    ylab("Residual Standard") + 
    labs(title = paste0("Residual Standard vs ",names(lm.eq.t$coefficients)[i])) +
    geom_point()
}

res.plot[[1]] <- ggplot(model, aes(x = model[,2], y = rstandard)) +
  xlab(names(lm.eq.t$coefficients)[2]) +
  ylab("Residual Standard") + 
  labs(title = paste0("Residual Standard vs ",names(lm.eq.t$coefficients)[2])) +
  geom_point()

res.plot[[2]] <- ggplot(model, aes(x = model[,3], y = rstandard)) +
  xlab(names(lm.eq.t$coefficients)[3]) +
  ylab("Residual Standard") + 
  labs(title = paste0("Residual Standard vs ",names(lm.eq.t$coefficients)[3])) +
  geom_jitter()

res.plot[[3]] <- ggplot(model, aes(x = model[,4], y = rstandard)) +
  xlab(names(lm.eq.t$coefficients)[4]) +
  ylab("Residual Standard") + 
  labs(title = paste0("Residual Standard vs ",names(lm.eq.t$coefficients)[4])) +
  geom_point()

res.plot[[4]] <- ggplot(model, aes(x = model[,5], y = rstandard)) +
  xlab(names(lm.eq.t$coefficients)[5]) +
  ylab("Residual Standard") + 
  labs(title = paste0("Residual Standard vs ",names(lm.eq.t$coefficients)[5])) +
  geom_jitter()

res.plot[[5]] <- ggplot(model, aes(x = model[,6], y = rstandard)) +
  xlab(names(lm.eq.t$coefficients)[6]) +
  ylab("Residual Standard") + 
  labs(title = paste0("Residual Standard vs ",names(lm.eq.t$coefficients)[6])) +
  geom_jitter()

res.plot[[6]] <- ggplot(model, aes(x = model[,7], y = rstandard)) +
  xlab(names(lm.eq.t$coefficients)[7]) +
  ylab("Residual Standard") + 
  labs(title = paste0("Residual Standard vs ",names(lm.eq.t$coefficients)[7])) +
  geom_jitter()

do.call(grid.arrange,res.plot)

plot(rstandard(lm.eq.t), ylab = "Residual Standard", main = "Residual Standard vs Index")




