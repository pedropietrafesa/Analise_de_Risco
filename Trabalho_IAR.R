---
title: "Trabalho de Introdução à Análise de Risco"
author: "Jaqueline Rodrigues de Souza Gentil e Pedro Araujo Pietrafesa"
date: "18/09/2022"
---



library(reshape2) 
library(magrittr)
library(tidyverse)



## Baixando o banco de dados




# Os dados selecionados foram o dados15 
library(readr)
dados <- read_delim("~/Downloads/dados15.csv", 
                    delim = ";", escape_double = FALSE, locale = locale(decimal_mark = ",",grouping_mark = "."), 
                    trim_ws = TRUE)




## Organizando o banco de dados 




# nomeando a primeira coluna que estava sem identificação inicial

dados %<>% set_colnames(c('ID','VALOR_PARCELA',paste0('P',1:15)))
head(dados)

# Ordenando por linha o banco de dados para criar a Variável X e N
d <- dados %>% reshape2::melt(id = c('ID','VALOR_PARCELA'))

# Estruturando o banco de dados pelas observações individuais
d_ind <- d %>% mutate(valor_da_parcela = ifelse(value == 1, VALOR_PARCELA,0), 
                      sinistro = ifelse(value == 1,1,0)) %>% na.omit(d_ind) %>%
  group_by(ID) %>% summarise(areceber = sum(valor_da_parcela),
                             n_de_sinistro = sum(sinistro)) %>% as.data.frame() %>% 
  set_colnames(c('id','valor_a_receber','n_de_sinistro'))








# Selecionando as observações que possuem sinistro (parcelas não pagas)

bd <- d_ind %>% filter(n_de_sinistro != 0)

bd %<>% set_colnames(c('ID','X','N'))




# Prêmio puro, Probabilidade de Ruina para um tempo finito - testes iniciais das distribuições e ajuste do modelo para Garantir a solvência da empresa para 1 ano.



# E[S_col] para Binomial Negativa Composta 
# E[S_col] = r*(q/p) E[X]
# V[S_col] para Binomial Negativa Composta 
# V[S_col] = r*(q/p)*E[X^2] + r*(q^2/p^2)*E[X]^2

bd1 <- bd  %>% filter(N < 13)


# Distribuição de N

dispersion_test <- function(x) 
{
  res <- 1-2 * abs((1 - pchisq((sum((x - mean(x))^2)/mean(x)), length(x) - 1))-0.5)
  
  cat("Teste de dispersão de dados de contagem:\n",
      length(x), " observações.\n",
      "Média: ",mean(x),"\n",
      "Variância: ",var(x),"\n",
      "Probabilidade de seguir uma distribuição de Poisson: ", 
      round(res, 3),"\n", sep = "")
  
  invisible(res)
}


ggplot(bd1, aes(x = bd1$N)) + 
  geom_histogram(aes(y = ..density..)) +
  geom_density()

dispersion_test(bd1$N)

# Com a análise dos dados sobre o número de sinistros (Não pagamento de parcelas) foi observado que Var[N] > E[N]. Desta forma, segundo Ferreira (2010), a distribuição de Poisson para N não é adequada. A alternativa seria utilizar a Binoial Negativa


# Os parâmetros da distribuição de N
f1 <- fitdist(bd1$N, "nbinom")
f1

# p será a probabilidade do cliente não pagar a prestação
# p = size/(size+mu)
p1 <- (3.256527)/(3.256527+4.570958)
q1 <- 1-p1

# O parâmetro r será o output size do software R 
r <- 3.256527

# Distribuição de X

w1 <- ggplot(bd1, aes(x = bd1$X)) + 
  geom_histogram(aes(y = ..density..)) +
  geom_density()


w2 <- ggplot(bd1, aes(x = ID , y = X, fill = X)) + 
  geom_violin()

ggpubr::ggarrange(w1,w2, 
                  ncol = 2, nrow = 1,
                  labels = c("a) Histograma", "b) Boxplot no formato de Violino"))

descdist(bd1$X, discrete=FALSE, boot=5000)
fitdist(bd1$X, "gamma", method = "mme")


summary(bd1$X)


#Como o plano de resseguros é de excesso de danos, então,
# Vamos supor Bd1$X ~ gama(0.695631324, 0.001158683)

# a seguradora pega até LT. Se passar de LT a 

# Assim, 

#E[X^k] = int^{LT}_{0} x^k f_{X}(x)dx + LT^k (1-F_{X}(LT))


x <- grouped.data(bd1$X)

hist(x)

Fnt <- ogive(x)

knots(Fnt)

Fnt(knots(Fnt))
plot(Fnt)


# Teste inicial do cálculo do Limite Técnico e Probabilidade de Ruína.

lt <- 40
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret <- mean(bd1$N)*ex_ret
dp <- sqrt((r*(q1/p1)*ex2_ret) + (r*((q1^2)/(p1^2))*ex_ret^2))

# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - e_col_ret) / dp

# Probabilidade de ruina 
1 - pnorm(z)




# Simulação da Probabilidade de Ruína com a Reserva de Risco igual a 0



# Mu = 0 e LT = 200 

lt <- 200
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret1 <- i$value + mais


# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret1 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] e DP[S_ret]

e_col_ret1 <- mean(bd1$N)*ex_ret1
dp1 <- sqrt((r*(q1/p1)*ex2_ret1) + (r*((q1^2)/(p1^2))*ex_ret1^2))

# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret1 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret1) / sqrt((r*(q1/p1)*ex2_ret1) + (r*((q1^2)/(p1^2))*ex_ret1^2))

# Probabilidade de ruina 
ruina1 <- 1 - pnorm(z)


# Mu = 0 e LT = 400 

lt <- 400
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret2 <- i$value + mais


# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret2 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] e DP[S_ret]

e_col_ret2 <- mean(bd1$N)*ex_ret2
dp2 <- sqrt((r*(q1/p1)*ex2_ret2) + (r*((q1^2)/(p1^2))*ex_ret2^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret2 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret2) / sqrt((r*(q1/p1)*ex2_ret2) + (r*((q1^2)/(p1^2))*ex_ret2^2))

# Probabilidade de ruina 
ruina2 <- 1 - pnorm(z)


# Mu = 0 e LT = 600 

lt <- 600
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret3 <- i$value + mais



# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret3 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret3 <- mean(bd1$N)*ex_ret3
dp3 <- sqrt((r*(q1/p1)*ex_ret3) + (r*((q1^2)/(p1^2))*ex_ret3^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret3 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret3) / sqrt((r*(q1/p1)*ex2_ret3) + (r*((q1^2)/(p1^2))*ex_ret3^2))

# Probabilidade de ruina 
ruina3 <- 1 - pnorm(z)


# Mu = 0 e LT = 800 

lt <- 800
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret4 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret4 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret4 <- mean(bd1$N)*ex_ret4
dp4 <- sqrt((r*(q1/p1)*ex2_ret4) + (r*((q1^2)/(p1^2))*ex_ret4^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret4 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret4) / sqrt((r*(q1/p1)*ex2_ret4) + (r*((q1^2)/(p1^2))*ex_ret4^2))

# Probabilidade de ruina 
ruina4 <- 1 - pnorm(z)



# Mu = 0 e LT = 1000 

lt <- 1000
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret5 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret5 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret5 <- mean(bd1$N)*ex_ret5
dp5 <- sqrt((r*(q1/p1)*ex2_ret5) + (r*((q1^2)/(p1^2))*ex_ret5^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret5 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret5) / sqrt((r*(q1/p1)*ex2_ret5) + (r*((q1^2)/(p1^2))*ex_ret5^2))

# Probabilidade de ruina 
ruina5 <- 1 - pnorm(z)


# Mu = 0 e LT = 1200 

lt <- 1200
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret6 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret6 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret6 <- mean(bd1$N)*ex_ret6
dp6 <- sqrt((r*(q1/p1)*ex2_ret6) + (r*((q1^2)/(p1^2))*ex_ret6^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret6 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret6) / sqrt((r*(q1/p1)*ex2_ret6) + (r*((q1^2)/(p1^2))*ex_ret6^2))

# Probabilidade de ruina 
ruina6 <- 1 - pnorm(z)

# Mu = 0 e LT = 1400 

lt <- 1400
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret7 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret7 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret7 <- mean(bd1$N)*ex_ret7
dp7 <- sqrt((r*(q1/p1)*ex2_ret7) + (r*((q1^2)/(p1^2))*ex_ret7^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret7 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret7) / sqrt((r*(q1/p1)*ex2_ret7) + (r*((q1^2)/(p1^2))*ex_ret7^2))

# Probabilidade de ruina 
ruina7 <- 1 - pnorm(z)


# Mu = 0 e LT = 1600 

lt <- 1600
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret8 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret8 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret8 <- mean(bd1$N)*ex_ret8
dp8 <- sqrt((r*(q1/p1)*ex2_ret8) + (r*((q1^2)/(p1^2))*ex_ret8^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret8 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret8) / sqrt((r*(q1/p1)*ex2_ret8) + (r*((q1^2)/(p1^2))*ex_ret8^2))

# Probabilidade de ruina 
ruina8 <- 1 - pnorm(z)


# Mu = 0 e LT = 1800 

lt <- 1800
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret9 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret9 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret9 <- mean(bd1$N)*ex_ret9
dp9 <- sqrt((r*(q1/p1)*ex2_ret9) + (r*((q1^2)/(p1^2))*ex_ret9^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret9 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret9) / sqrt((r*(q1/p1)*ex2_ret9) + (r*((q1^2)/(p1^2))*ex_ret9^2))

# Probabilidade de ruina 
ruina9 <- 1 - pnorm(z)

# Mu = 0 e LT = 2000 

lt <- 2000
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret10 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret10 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret10 <- mean(bd1$N)*ex_ret10
dp10 <- sqrt((r*(q1/p1)*ex2_ret10) + (r*((q1^2)/(p1^2))*ex_ret10^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret10 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret10) / sqrt((r*(q1/p1)*ex2_ret10) + (r*((q1^2)/(p1^2))*ex_ret10^2))

# Probabilidade de ruina 
ruina10 <- 1 - pnorm(z)


########################## Data frame com os resultados  ##################################
#


r1 <- data.frame(lt = c(200,400,600,800,1000,1200,1400,1600,1800,
                        2000),
                 es = c(e_col_ret1, e_col_ret2, e_col_ret3, e_col_ret4, 
                        e_col_ret5,e_col_ret6,e_col_ret7,e_col_ret8,
                        e_col_ret9, e_col_ret10),
                 ruina = c(ruina1, ruina2, ruina3,ruina4,
                           ruina5, ruina6, ruina7, ruina8, 
                           ruina9, ruina10),
                 dp = c(dp1, dp2, dp3, dp4, dp5, dp6, dp7, dp8, dp9, dp10))







# Simulação da Probabilidade de Ruína com a Reserva de Risco Igual a 664





# Mu = 0 e LT = 200 

lt <- 200
mu <- 664

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret1 <- i$value + mais


# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret1 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] e DP[S_ret]

e_col_ret1 <- mean(bd1$N)*ex_ret1
dp1 <- sqrt((r*(q1/p1)*ex2_ret1) + (r*((q1^2)/(p1^2))*ex_ret1^2))

# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret1 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret1) / sqrt((r*(q1/p1)*ex2_ret1) + (r*((q1^2)/(p1^2))*ex_ret1^2))

# Probabilidade de ruina 
ruina1 <- 1 - pnorm(z)


# Mu = 0 e LT = 400 

lt <- 400
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret2 <- i$value + mais


# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret2 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] e DP[S_ret]

e_col_ret2 <- mean(bd1$N)*ex_ret2
dp2 <- sqrt((r*(q1/p1)*ex2_ret2) + (r*((q1^2)/(p1^2))*ex_ret2^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret2 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret2) / sqrt((r*(q1/p1)*ex2_ret2) + (r*((q1^2)/(p1^2))*ex_ret2^2))

# Probabilidade de ruina 
ruina2 <- 1 - pnorm(z)


# Mu = 0 e LT = 600 

lt <- 600
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret3 <- i$value + mais



# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret3 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret3 <- mean(bd1$N)*ex_ret3
dp3 <- sqrt((r*(q1/p1)*ex_ret3) + (r*((q1^2)/(p1^2))*ex_ret3^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret3 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret3) / sqrt((r*(q1/p1)*ex2_ret3) + (r*((q1^2)/(p1^2))*ex_ret3^2))

# Probabilidade de ruina 
ruina3 <- 1 - pnorm(z)


# Mu = 0 e LT = 800 

lt <- 800
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret4 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret4 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret4 <- mean(bd1$N)*ex_ret4
dp4 <- sqrt((r*(q1/p1)*ex2_ret4) + (r*((q1^2)/(p1^2))*ex_ret4^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret4 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret4) / sqrt((r*(q1/p1)*ex2_ret4) + (r*((q1^2)/(p1^2))*ex_ret4^2))

# Probabilidade de ruina 
ruina4 <- 1 - pnorm(z)



# Mu = 0 e LT = 1000 

lt <- 1000
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret5 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret5 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret5 <- mean(bd1$N)*ex_ret5
dp5 <- sqrt((r*(q1/p1)*ex2_ret5) + (r*((q1^2)/(p1^2))*ex_ret5^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret5 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret5) / sqrt((r*(q1/p1)*ex2_ret5) + (r*((q1^2)/(p1^2))*ex_ret5^2))

# Probabilidade de ruina 
ruina5 <- 1 - pnorm(z)


# Mu = 0 e LT = 1200 

lt <- 1200
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret6 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret6 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret6 <- mean(bd1$N)*ex_ret6
dp6 <- sqrt((r*(q1/p1)*ex2_ret6) + (r*((q1^2)/(p1^2))*ex_ret6^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret6 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret6) / sqrt((r*(q1/p1)*ex2_ret6) + (r*((q1^2)/(p1^2))*ex_ret6^2))

# Probabilidade de ruina 
ruina6 <- 1 - pnorm(z)

# Mu = 0 e LT = 1400 

lt <- 1400
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret7 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret7 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret7 <- mean(bd1$N)*ex_ret7
dp7 <- sqrt((r*(q1/p1)*ex2_ret7) + (r*((q1^2)/(p1^2))*ex_ret7^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret7 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret7) / sqrt((r*(q1/p1)*ex2_ret7) + (r*((q1^2)/(p1^2))*ex_ret7^2))

# Probabilidade de ruina 
ruina7 <- 1 - pnorm(z)


# Mu = 0 e LT = 1600 

lt <- 1600
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret8 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret8 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret8 <- mean(bd1$N)*ex_ret8
dp8 <- sqrt((r*(q1/p1)*ex2_ret8) + (r*((q1^2)/(p1^2))*ex_ret8^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret8 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret8) / sqrt((r*(q1/p1)*ex2_ret8) + (r*((q1^2)/(p1^2))*ex_ret8^2))

# Probabilidade de ruina 
ruina8 <- 1 - pnorm(z)


# Mu = 0 e LT = 1800 

lt <- 1800
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret9 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret9 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret9 <- mean(bd1$N)*ex_ret9
dp9 <- sqrt((r*(q1/p1)*ex2_ret9) + (r*((q1^2)/(p1^2))*ex_ret9^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret9 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret9) / sqrt((r*(q1/p1)*ex2_ret9) + (r*((q1^2)/(p1^2))*ex_ret9^2))

# Probabilidade de ruina 
ruina9 <- 1 - pnorm(z)

# Mu = 0 e LT = 2000 

lt <- 2000
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret10 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret10 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret10 <- mean(bd1$N)*ex_ret10
dp10 <- sqrt((r*(q1/p1)*ex2_ret10) + (r*((q1^2)/(p1^2))*ex_ret10^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret10 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret10) / sqrt((r*(q1/p1)*ex2_ret10) + (r*((q1^2)/(p1^2))*ex_ret10^2))

# Probabilidade de ruina 
ruina10 <- 1 - pnorm(z)


########################## Data frame com os resultados  ##################################
#


r2 <- data.frame(lt = c(200,400,600,800,1000,1200,1400,1600,1800,
                        2000),
                 es = c(e_col_ret1, e_col_ret2, e_col_ret3, e_col_ret4, 
                        e_col_ret5,e_col_ret6,e_col_ret7,e_col_ret8,
                        e_col_ret9, e_col_ret10),
                 ruina = c(ruina1, ruina2, ruina3,ruina4,
                           ruina5, ruina6, ruina7, ruina8, 
                           ruina9, ruina10),
                 dp = c(dp1, dp2, dp3, dp4, dp5, dp6, dp7, dp8, dp9, dp10))






# Simulação da Probabilidade de Ruína com a Reserva de Risco Igual a 1328



# Mu = 0 e LT = 200 

lt <- 200
mu <- 1328

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret1 <- i$value + mais


# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret1 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] e DP[S_ret]

e_col_ret1 <- mean(bd1$N)*ex_ret1
dp1 <- sqrt((r*(q1/p1)*ex2_ret1) + (r*((q1^2)/(p1^2))*ex_ret1^2))

# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret1 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret1) / sqrt((r*(q1/p1)*ex2_ret1) + (r*((q1^2)/(p1^2))*ex_ret1^2))

# Probabilidade de ruina 
ruina1 <- 1 - pnorm(z)


# Mu = 0 e LT = 400 

lt <- 400
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret2 <- i$value + mais


# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret2 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] e DP[S_ret]

e_col_ret2 <- mean(bd1$N)*ex_ret2
dp2 <- sqrt((r*(q1/p1)*ex2_ret2) + (r*((q1^2)/(p1^2))*ex_ret2^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret2 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret2) / sqrt((r*(q1/p1)*ex2_ret2) + (r*((q1^2)/(p1^2))*ex_ret2^2))

# Probabilidade de ruina 
ruina2 <- 1 - pnorm(z)


# Mu = 0 e LT = 600 

lt <- 600
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret3 <- i$value + mais



# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret3 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret3 <- mean(bd1$N)*ex_ret3
dp3 <- sqrt((r*(q1/p1)*ex_ret3) + (r*((q1^2)/(p1^2))*ex_ret3^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret3 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret3) / sqrt((r*(q1/p1)*ex2_ret3) + (r*((q1^2)/(p1^2))*ex_ret3^2))

# Probabilidade de ruina 
ruina3 <- 1 - pnorm(z)


# Mu = 0 e LT = 800 

lt <- 800
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret4 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret4 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret4 <- mean(bd1$N)*ex_ret4
dp4 <- sqrt((r*(q1/p1)*ex2_ret4) + (r*((q1^2)/(p1^2))*ex_ret4^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret4 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret4) / sqrt((r*(q1/p1)*ex2_ret4) + (r*((q1^2)/(p1^2))*ex_ret4^2))

# Probabilidade de ruina 
ruina4 <- 1 - pnorm(z)



# Mu = 0 e LT = 1000 

lt <- 1000
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret5 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret5 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret5 <- mean(bd1$N)*ex_ret5
dp5 <- sqrt((r*(q1/p1)*ex2_ret5) + (r*((q1^2)/(p1^2))*ex_ret5^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret5 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret5) / sqrt((r*(q1/p1)*ex2_ret5) + (r*((q1^2)/(p1^2))*ex_ret5^2))

# Probabilidade de ruina 
ruina5 <- 1 - pnorm(z)


# Mu = 0 e LT = 1200 

lt <- 1200
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret6 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret6 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret6 <- mean(bd1$N)*ex_ret6
dp6 <- sqrt((r*(q1/p1)*ex2_ret6) + (r*((q1^2)/(p1^2))*ex_ret6^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret6 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret6) / sqrt((r*(q1/p1)*ex2_ret6) + (r*((q1^2)/(p1^2))*ex_ret6^2))

# Probabilidade de ruina 
ruina6 <- 1 - pnorm(z)

# Mu = 0 e LT = 1400 

lt <- 1400
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret7 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret7 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret7 <- mean(bd1$N)*ex_ret7
dp7 <- sqrt((r*(q1/p1)*ex2_ret7) + (r*((q1^2)/(p1^2))*ex_ret7^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret7 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret7) / sqrt((r*(q1/p1)*ex2_ret7) + (r*((q1^2)/(p1^2))*ex_ret7^2))

# Probabilidade de ruina 
ruina7 <- 1 - pnorm(z)


# Mu = 0 e LT = 1600 

lt <- 1600
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret8 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret8 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret8 <- mean(bd1$N)*ex_ret8
dp8 <- sqrt((r*(q1/p1)*ex2_ret8) + (r*((q1^2)/(p1^2))*ex_ret8^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret8 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret8) / sqrt((r*(q1/p1)*ex2_ret8) + (r*((q1^2)/(p1^2))*ex_ret8^2))

# Probabilidade de ruina 
ruina8 <- 1 - pnorm(z)


# Mu = 0 e LT = 1800 

lt <- 1800
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret9 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret9 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret9 <- mean(bd1$N)*ex_ret9
dp9 <- sqrt((r*(q1/p1)*ex2_ret9) + (r*((q1^2)/(p1^2))*ex_ret9^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret9 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret9) / sqrt((r*(q1/p1)*ex2_ret9) + (r*((q1^2)/(p1^2))*ex_ret9^2))

# Probabilidade de ruina 
ruina9 <- 1 - pnorm(z)

# Mu = 0 e LT = 2000 

lt <- 2000
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret10 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret10 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret10 <- mean(bd1$N)*ex_ret10
dp10 <- sqrt((r*(q1/p1)*ex2_ret10) + (r*((q1^2)/(p1^2))*ex_ret10^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret10 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret10) / sqrt((r*(q1/p1)*ex2_ret10) + (r*((q1^2)/(p1^2))*ex_ret10^2))

# Probabilidade de ruina 
ruina10 <- 1 - pnorm(z)


########################## Data frame com os resultados  ##################################
#


r3 <- data.frame(lt = c(200,400,600,800,1000,1200,1400,1600,1800,
                        2000),
                 es = c(e_col_ret1, e_col_ret2, e_col_ret3, e_col_ret4, 
                        e_col_ret5,e_col_ret6,e_col_ret7,e_col_ret8,
                        e_col_ret9, e_col_ret10),
                 ruina = c(ruina1, ruina2, ruina3,ruina4,
                           ruina5, ruina6, ruina7, ruina8, 
                           ruina9, ruina10),
                 dp = c(dp1, dp2, dp3, dp4, dp5, dp6, dp7, dp8, dp9, dp10))






# Simulação da Probabilidade de Ruína com a Reserva de Risco Igual a 1992


# Mu = 0 e LT = 200 

lt <- 200
mu <- 1992

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret1 <- i$value + mais


# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret1 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] e DP[S_ret]

e_col_ret1 <- mean(bd1$N)*ex_ret1
dp1 <- sqrt((r*(q1/p1)*ex2_ret1) + (r*((q1^2)/(p1^2))*ex_ret1^2))

# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret1 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret1) / sqrt((r*(q1/p1)*ex2_ret1) + (r*((q1^2)/(p1^2))*ex_ret1^2))

# Probabilidade de ruina 
ruina1 <- 1 - pnorm(z)


# Mu = 0 e LT = 400 

lt <- 400
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret2 <- i$value + mais


# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret2 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] e DP[S_ret]

e_col_ret2 <- mean(bd1$N)*ex_ret2
dp2 <- sqrt((r*(q1/p1)*ex2_ret2) + (r*((q1^2)/(p1^2))*ex_ret2^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret2 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret2) / sqrt((r*(q1/p1)*ex2_ret2) + (r*((q1^2)/(p1^2))*ex_ret2^2))

# Probabilidade de ruina 
ruina2 <- 1 - pnorm(z)


# Mu = 0 e LT = 600 

lt <- 600
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret3 <- i$value + mais



# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret3 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret3 <- mean(bd1$N)*ex_ret3
dp3 <- sqrt((r*(q1/p1)*ex_ret3) + (r*((q1^2)/(p1^2))*ex_ret3^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret3 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret3) / sqrt((r*(q1/p1)*ex2_ret3) + (r*((q1^2)/(p1^2))*ex_ret3^2))

# Probabilidade de ruina 
ruina3 <- 1 - pnorm(z)


# Mu = 0 e LT = 800 

lt <- 800
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret4 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret4 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret4 <- mean(bd1$N)*ex_ret4
dp4 <- sqrt((r*(q1/p1)*ex2_ret4) + (r*((q1^2)/(p1^2))*ex_ret4^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret4 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret4) / sqrt((r*(q1/p1)*ex2_ret4) + (r*((q1^2)/(p1^2))*ex_ret4^2))

# Probabilidade de ruina 
ruina4 <- 1 - pnorm(z)



# Mu = 0 e LT = 1000 

lt <- 1000
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret5 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret5 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret5 <- mean(bd1$N)*ex_ret5
dp5 <- sqrt((r*(q1/p1)*ex2_ret5) + (r*((q1^2)/(p1^2))*ex_ret5^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret5 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret5) / sqrt((r*(q1/p1)*ex2_ret5) + (r*((q1^2)/(p1^2))*ex_ret5^2))

# Probabilidade de ruina 
ruina5 <- 1 - pnorm(z)


# Mu = 0 e LT = 1200 

lt <- 1200
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret6 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret6 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret6 <- mean(bd1$N)*ex_ret6
dp6 <- sqrt((r*(q1/p1)*ex2_ret6) + (r*((q1^2)/(p1^2))*ex_ret6^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret6 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret6) / sqrt((r*(q1/p1)*ex2_ret6) + (r*((q1^2)/(p1^2))*ex_ret6^2))

# Probabilidade de ruina 
ruina6 <- 1 - pnorm(z)

# Mu = 0 e LT = 1400 

lt <- 1400
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret7 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret7 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret7 <- mean(bd1$N)*ex_ret7
dp7 <- sqrt((r*(q1/p1)*ex2_ret7) + (r*((q1^2)/(p1^2))*ex_ret7^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret7 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret7) / sqrt((r*(q1/p1)*ex2_ret7) + (r*((q1^2)/(p1^2))*ex_ret7^2))

# Probabilidade de ruina 
ruina7 <- 1 - pnorm(z)


# Mu = 0 e LT = 1600 

lt <- 1600
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret8 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret8 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret8 <- mean(bd1$N)*ex_ret8
dp8 <- sqrt((r*(q1/p1)*ex2_ret8) + (r*((q1^2)/(p1^2))*ex_ret8^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret8 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret8) / sqrt((r*(q1/p1)*ex2_ret8) + (r*((q1^2)/(p1^2))*ex_ret8^2))

# Probabilidade de ruina 
ruina8 <- 1 - pnorm(z)


# Mu = 0 e LT = 1800 

lt <- 1800
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret9 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret9 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret9 <- mean(bd1$N)*ex_ret9
dp9 <- sqrt((r*(q1/p1)*ex2_ret9) + (r*((q1^2)/(p1^2))*ex_ret9^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret9 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret9) / sqrt((r*(q1/p1)*ex2_ret9) + (r*((q1^2)/(p1^2))*ex_ret9^2))

# Probabilidade de ruina 
ruina9 <- 1 - pnorm(z)

# Mu = 0 e LT = 2000 

lt <- 2000
mu <- 0

# Encontrando E[X_ret], E[X] = int^{LT}_{0} x f_{X}(x)dx + LT (1-F_{X}(LT))
i <- integrate(function(x){x*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais <- lt*(1-pgamma(lt,0.695631324, 0.001158683))

ex_ret10 <- i$value + mais

i$value

# Encontrando E[X^2_ret], E[X^2] = int{LT}_{0} x^2 f_{X}(x)dx + LT^2 (1-F_{X}(LT))
i2 <- integrate(function(x){x^2*dgamma(x,0.695631324, 0.001158683)}, 0, lt) 
mais2 <- lt^2*(1-pgamma(lt,0.695631324, 0.001158683))

ex2_ret10 <- i2$value + mais2


# O prêmio puro retido é calculado em função do LT

# Carregamento de Segurança
#E[S_col] = E[N]E[X]
e_col <- mean(bd1$N)*mean(bd1$X)

# Variância de S_coletivo
#V[S_col] = E[N]V[X]+E[X^2]V[N]
x2 <- (bd1$X)^2  
v_col <- (mean(bd1$N)*var(bd1$X))+(mean(x2)*var(bd1$N))

# carregamento de segurança com Z a 95% de confiança
teta <- (1.645*sqrt(v_col))/e_col


# E[S_ret] 

e_col_ret10 <- mean(bd1$N)*ex_ret10
dp10 <- sqrt((r*(q1/p1)*ex2_ret10) + (r*((q1^2)/(p1^2))*ex_ret10^2))


# Prêmio Puro em função do LT - Prêmio_ret = E[S_ret](1 + theta)

premio_ret = e_col_ret10 * (1 +teta)


#Calculo otimizado de LT

z <- (mu + premio_ret - (r*(q1/p1))*ex_ret10) / sqrt((r*(q1/p1)*ex2_ret10) + (r*((q1^2)/(p1^2))*ex_ret10^2))

# Probabilidade de ruina 
ruina10 <- 1 - pnorm(z)


########################## Data frame com os resultados  ##################################
#


r4 <- data.frame(lt = c(200,400,600,800,1000,1200,1400,1600,1800,
                        2000),
                 es = c(e_col_ret1, e_col_ret2, e_col_ret3, e_col_ret4, 
                        e_col_ret5,e_col_ret6,e_col_ret7,e_col_ret8,
                        e_col_ret9, e_col_ret10),
                 ruina = c(ruina1, ruina2, ruina3,ruina4,
                           ruina5, ruina6, ruina7, ruina8, 
                           ruina9, ruina10),
                 dp = c(dp1, dp2, dp3, dp4, dp5, dp6, dp7, dp8, dp9, dp10))



















