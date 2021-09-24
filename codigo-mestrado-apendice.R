#################################################
##      Carregando os pacotes do R: library      ##
###################################################
rm(list=ls(all=TRUE))

library(evd)
library(POT)
library(evdbayes)


#analise dados Petropolis-Teresopolis-JF###########################
dados=read.table("max_anual_Tere_Petr_JF.txt", dec=".",h=T)                                        
dados

                  ## Armazenando por mês de cada ano do estudo em x
y1=dados$ Max_Petr_ano[1:78] ;y1
y2=dados$max_Tere_ano[1:78] ;y2
y3=dados$Max_jf_ano[1:78] ;y3
x=dados$Ano;x 
#################Análise descritiva dos dados################
summary(y1)
summary(y2)
summary(y3)
###############Cálculo do desvio Padrão###################
sd(y1)
sd(y2)
sd(y3)
##################Cálculo do coeficiente de variação###########################
CV1 = 100*sd(y1)/mean(y1)
CV2 = 100*sd(y2)/mean(y2)
CV3 = 100*sd(y3)/mean(y3)
CV1
CV2
CV3
##############################################
##    Verificando a independencia
###############################################
Box.test(y1,type=("Ljung-Box"))   #teste de independencia
Box.test(y2,type=("Ljung-Box"))   #teste de independencia
Box.test(y3,type=("Ljung-Box"))   #teste de independencia
#################################################
# Verificando Tendência Teste Mannkendall #######
#################################################
require(Kendall)
library(Kendall)

MannKendall(y1)
MannKendall(y2)
MannKendall(y3)

##########################################################
##               Análise gráfica preliminar             ##
##########################################################

par(mfrow = c(3,1))

plot(x,y1,main='Petrópolis (RJ) ',ylim=c(32,186),xlim=c(1943,2020),col="black",pch =19, xlab='Anos', ylab='PMA (em mm)')
plot(x,y2,main='Teresópolis (RJ) ',ylim=c(32,186),xlim=c(1943,2020),col="black",pch =19, xlab='Anos', ylab='PMA (em mm)')
plot(x,y3,main=' Juiz de Fora (MG) ',ylim=c(32,186),xlim=c(1943,2020),col="black",pch =19, xlab='Anos', ylab='PMA (em mm)')

#########AS obsevações foram divididas para a análise do EMP######################
######################################################################################

x=dados$Ano;x                      ## Armazenando por mês de cada ano do estudo em x


yc1=dados$Max_Petr_ano[1:38] ;yc1
yc2=dados$max_Tere_ano[1:38] ;yc2
yc3=dados$Max_jf_ano[1:38] ;yc3

xc=dados$Ano[1:38] ;xc          
summary(yc1)
summary(yc2)
summary(yc3)

##############################################
##    Verificando a independencia
###############################################


Box.test(yc1,type=("Ljung-Box"))   #teste de independencia
Box.test(yc2,type=("Ljung-Box"))   #teste de independencia
Box.test(yc3,type=("Ljung-Box"))   #teste de independencia

###############################################
##    Verificando a aleatoriedade dos dados  ##
##                  Runs test                ##
###############################################
# Mediana
require(snpar)
library(snpar)

runs.test(yc1)
length(yc1)
length(yc2)
length(yc3)

runs.test(yc1)
runs.test(yc2)

runs.test(yc3)

#################################################
# Verificando Tendência Teste Mannkendall #######
#################################################
require(Kendall)
library(Kendall)

MannKendall(yc1)
MannKendall(yc2)
MannKendall(yc3)

##################################
############Cálculo das EMV da GEV#########################################

fgev(yc1)
fgev(yc2)
fgev(yc3)
confint(fgev(yc1),level=0.95) # Intervalo de confiança aprox pela normal
confint(fgev(yc2),level=0.95) # Intervalo de confiança aprox pela normal
confint(fgev(yc3),level=0.95) # Intervalo de confiança aprox pela normal


##############################################################################
## Avaliando ajuste com priori Normal Trivariada NÃO informativa ##
###################################################################


vmed <-c(0,0,0)  ## Hiperparâmetros para a média
mvar <- diag(c(10000000, 100000, 1000)) ## Hiperparâmetros para matriz de variância
p_nao_inf <- prior.norm(mean =vmed , cov = mvar)  # priori Normal Trivariada 

#
#sd=sqrt(100000);sd
#hist(rnorm(1000,0,sd))

t_cadeia=320000
t_queima=20000
salto=60

T_fin=(t_cadeia-t_queima)/salto;T_fin #deixar >4000

posteNI<-posterior(t_cadeia, init = c(50, 30,0.1), prior = p_nao_inf, lh = "gev",
data =yc1 , psd = c(1,1,.1),burn=t_queima,thin=salto)


#####################################################################
#    UTILIZANDO O PACOTE CODA para analisar congergência das cadeias
###################################################################
 
library("coda")# Pacote para analisar as cadeias a posteriori

param.mcmcNI <- as.mcmc(posteNI)

summary(param.mcmcNI)

raftery.diag(param.mcmcNI)
geweke.diag(param.mcmcNI)
heidel.diag(param.mcmcNI)
HPDinterval(param.mcmcNI)


# Análise gráfica de convergência
#plot(param.mcmcNI)

#======================================================================
## Salvando as estimatinas médias a posteriori
#======================================================================

mu2<-mean(param.mcmcNI[,1]); mu2    ## Estimativa do parâmetro posição (Gev)
sig2<-mean(param.mcmcNI[,2]);sig2   ## Estimativa do parâmetro escala (GEV)
xi2<-mean(param.mcmcNI[,3]);xi2    ## Estimativa do parâmetro Forma(GEV)

#======================================================================
#   Teste de Kolmogorov Smirnov e A-D
#======================================================================

library(goftest)
ks.test(yc1,"pgev",mu2,sig2,xi2,alternative="greater") # Dist GEV

ad.test(yc1, "pgev",loc=mu2,scale=sig2,shape=xi2, estimated=TRUE)

#======================================================================
# Calcula DIC gev Não Informativa
#=========================================

logLikelihood = function(yc1, theta) {
  #Get the individual parameters out of theta.
  mu = theta[1]
  sigma = theta[2]
qsi = theta[3]
  
  #sum of log likelihoods = log of product of likelihoods
  sum( dgev(yc1, mu, sigma, qsi, log=TRUE) )
}

calculateDIC_ni = function(yc1, theta_post, llFun) {
  #Calculate L
  theta_hat = apply(theta_post, 2, mean)
  L = llFun(yc1, theta_hat)

  #Calculate P
  S = nrow(theta_post) #S = number of iterations
  #Add up the log likelihoods of each iteration
  llSum = 0
  for (s in 1:S) {
    theta_s = theta_post[s,]
    llSum = llSum + llFun(yc1, theta_s)
  }
  P = 2 * (L - (1 / S * llSum))
  
  #Calculate DIC
  DIC = -2 * (L - P)
  
  #Return the results
  list(DIC=DIC, P=P, L=L)
}

calculateDIC_ni(yc1,param.mcmcNI, llFun=logLikelihood)


#======================================================================
#           Rotina R para PREDITIVA para a NÃO INFORMATIVA
#======================================================================

k<-length(posteNI[,1])
 
loc <- numeric(k) 
escala <- numeric(k) 
forma <- numeric(k) 

TR=c(10,20,30,40) #INDICANDO tempos de retorno, 
prob=  1/TR

np<-length(prob) 

q1 <- numeric(k) 
q2 <- numeric(k)
q3 <- numeric(k)
q4 <- numeric(k) 
for (i in 1:k)
{ 
          loc[i] <- posteNI[i,1] 
          escala[i] <- posteNI[i,2] 
          forma[i] <- posteNI[i,3] 
q1[i]<-(loc[i]+((escala[i]/forma[i])*(((-log(1-prob[1]))^ 
(-forma[i]))-1))) 
q2[i]<-(loc[i]+((escala[i]/forma[i])*(((-log(1-prob[2]))^ 
(-forma[i]))-1))) 
q3[i]<-(loc[i]+((escala[i]/forma[i])*(((-log(1-prob[3]))^ 
(-forma[i]))-1))) 
q4[i]<-(loc[i]+((escala[i]/forma[i])*(((-log(1-prob[4]))^ 
(-forma[i]))-1))) 
} 
 

nit<-length(q1)  
qoNI<-matrix(0,nit,4)  
qoNI[,1]<-q1 
qoNI[,2]<-q2 
qoNI[,3]<-q3
qoNI[,4]<-q4

#======================================================================
#     Análise dos níveis de retorno
#======================================================================
nivel_ret_mcmcNI <- as.mcmc(qoNI)

summary(nivel_ret_mcmcNI)

HPDinterval(nivel_ret_mcmcNI)

raftery.diag(nivel_ret_mcmcNI)
geweke.diag(nivel_ret_mcmcNI)
heidel.diag(nivel_ret_mcmcNI)

# plot(nivel_ret_mcmcNI)
#120.1 185.8 185.8 185.8

#======================================================================
#     CALCULANDO ERRO MÉDIO DE  PREDIÇÃO 
#======================================================================
#dados
yobs=dados$Max_Petr_ano[39:78];yobs

o1<- max(yobs[1:10]);o1
o2<- max(yobs[1:20]);o2
o3<- max(yobs[1:30]);o3
o4<- max(yobs[1:40]);o4
Vo = c(o1,o2,o3,o4);Vo   # Valores observados


q10ni<-mean(nivel_ret_mcmcNI[,1])
q20ni<-mean(nivel_ret_mcmcNI[,2])
q30ni<-mean(nivel_ret_mcmcNI[,3])
q40ni<-mean(nivel_ret_mcmcNI[,4])

pred_bayes_NI=c(q10ni,q20ni,q30ni,q40ni);pred_bayes_NI

EpB_NI= abs((Vo-pred_bayes_NI)/Vo);EpB_NI

round(mean(EpB_NI)*100,2)

#=================================================
#  fAZENDO AMP MÉDIA
#===================================================
HP=HPDinterval(nivel_ret_mcmcNI);HP
A10=HP[5]-HP[1];A10
A20=HP[6]-HP[2];A20
A30=HP[7]-HP[3];A30
A40=HP[8]-HP[4];A40

amps=c(A10,A20,A30 ,A40);amps
Amp_media=mean(amps);Amp_media



#####################################################################
#####################################################################
#
#          INFORMATIVA
#
#####################################################################
#####################################################################

##############################################################################
## Estimação Bayesiana com priori Normal Trivariada INFORMATIVA ##
###################################################################
## Informação Teresópolis0 ##
       

yc2


Box.test(yc2,type=("Ljung-Box")) 
#-------------------------------
require(snpar)
library(snpar)
runs.test(yc2)

#--------------------------
###############Cálculo dos parâmetros da Matriz de Variância e Covariâncias#################################
modj<-fgev(yc2);modj

mvarI<-modj$var.cov 

   mvarI
  e_mu<-(modj)[[1]][1:1] ;e_mu         ## Estimativa do parâmetro posição (Gumbel)
  e_sig<-(modj)[[1]][2:2];e_sig
e_xi<-(modj)[[1]][3:3];e_xi
  var_mu<-modj$var.cov[1] ; var_mu       ## Estimativa do  variancia do parâmetro posição (Gumbel)
  var_sig<-modj$var.cov[5]; var_sig
 var_xi<-modj$var.cov[9]; var_xi        ## Estimativa da variancia do parametro escala (Gumbel)

 cov_m_s=modj$var.cov[2]; cov_m_s
 cov_m_x=modj$var.cov[3]; cov_m_x
cov_s_x=modj$var.cov[6]; cov_s_x

#Informação a priori

vmed1<-c(e_mu,e_sig,e_xi);vmed1## Hiperparâmetros para a média## vetor de média para primeiro mês (janeiro)
mvar1<-matrix(0,3,3)

mvar1<-mvarI; mvar1

p_inf <- prior.norm(mean =vmed1 , cov = mvar1)  # priori Normal Trivariada 

#mvarD <- diag(c(var_mu, var_sig, var_xi))
#mvarD
#pinf_D<-prior.norm(mean =vmed1 , cov = mvarD)

####################################################################
## Estimação Bayesiana com priori Normal Trivariada INFORMATIVA ##
###################################################################

posteI<-posterior(120000, init = c(80,30,0.2), prior = p_inf, lh = "gev",
data =yc1 , psd = c(1,1,.1),burn=20000,thin=20)

#####################################################################
#    UTILIZANDO O PACOTE CODA para analisar congergência das cadeias
###################################################################
 
library("coda")# Pacote para analisar as cadeias a posteriori

param.mcmcI <- as.mcmc(posteI)


 
raftery.diag(param.mcmcI)
geweke.diag(param.mcmcI)
heidel.diag(param.mcmcI)
HPDinterval(param.mcmcI)
summary(param.mcmcI)

# Análise gráfica de convergência
 #plot(param.mcmcI)

HPDinterval(param.mcmcI4)

summary(param.mcmcI4)

# Análise gráfica de convergência
 #plot(param.mcmcI4)
