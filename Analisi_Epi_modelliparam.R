library(survival)
library(MASS)
e <-read.table("epi.dat", header=TRUE)

# fattori
e$ad_chemio <- as.factor(e$ad_chemio)
e$ad_hormon <- as.factor(e$ad_hormon)
e$cur_radiation <- as.factor(e$cur_radiation)
e$pre_hormon <- as.factor(e$pre_hormon)
e$pre_chemio <- as.factor(e$pre_chemio)
e$fam_heart <- as.factor(e$fam_heart)
e$n_sites <- as.factor(e$n_sites)
e$performance <- as.factor(e$performance)
e$cardiotox <- as.factor(e$cardiotox)
e$ad_radiation <- as.factor(e$ad_radiation)
e$rig_radiation <- as.factor(e$rig_radiation)
e$medi_radiation <- as.factor(e$medi_radiation)
e$lef_radiation <- as.factor(e$lef_radiation)
attach(e)

# verifico quale modello parametrico puo essere plausibile. 
# verifica sui rischi  proporzionali e sull'ipotesi distributiva (esponenziale o Weibull)
# in realta una verifica grafica della proporzionalita dei rischi, � gia stata fatta tramite stimatori non parametrici
# curve log(-log(S(t)))
S_hat = survfit(Surv(Time, status) ~ 1, data = e)
plot(S_hat)

#  MODELLO WEIBULL ed esponenziale
# SELEZIONE DEL MODELLO 
e1 <-na.omit(e) 
e1$perform <- e1$performance
e1$cumdose <-cut(e1$epi_m2,breaks=c(0,470,870,970,1600))
levels(e1$perform) <- c("1","2","2")

fp <- survreg(Surv(lifetime,event)~ ad_chemio +cur_radiation +pre_hormon +pre_chemio + n_sites+
              surface + perform +ad_radiation + medi_radiation + lef_radiation + eta + cumdose, data=e1,dist="weibull")
summary(fp)
stepAIC(fp,direction="both")
?stepAIC
# The final model is:
# survreg(formula = Surv(lifetime, event) ~ ad_chemio + pre_hormon + 
#           n_sites + perform + medi_radiation + eta + cumdose, data = e1, 
#         dist = "weibull")
# 
# Coefficients:
#   (Intercept)           ad_chemio2          pre_hormon2             n_sites2             n_sites3 
# 3.352309385          0.218317865          0.193849656         -0.666491336         -0.754478340 
# n_sites4             perform2      medi_radiation2                  eta     cumdose(470,870] 
# -0.749932304         -0.353804435         -0.251727316         -0.009560036          0.831216187 
# cumdose(870,970] cumdose(970,1.6e+03] 
# 1.176645594          1.181830036 
# 
# Scale= 0.8898227 
# 
# Loglik(model)= -3599.2   Loglik(intercept only)= -3823.2
# Chisq= 447.92 on 11 degrees of freedom, p= 0 



# Se non consideriamo la variabile cumdose in partenza, otteniamo:
fp1 <- survreg(Surv(lifetime,event)~ ad_chemio +cur_radiation +pre_hormon +pre_chemio + n_sites+
                surface + perform +ad_radiation + medi_radiation + lef_radiation + eta , data=e1,dist="weibull")
summary(fp1)
stepAIC(fp1,direction="both")

# survreg(formula = Surv(lifetime, event) ~ ad_chemio + pre_hormon + 
#           pre_chemio + n_sites + perform + eta, data = e1, dist = "weibull")
# 
# Coefficients:
#   (Intercept)  ad_chemio2 pre_hormon2 pre_chemio2    n_sites2    n_sites3    n_sites4    perform2         eta 
# 3.73529231  0.39213600  0.21299272  0.36897793 -0.72462566 -0.81871502 -0.84841159 -0.57413418 -0.01144286 
# 
# Scale= 0.9549786 
# 
# Loglik(model)= -3686.6   Loglik(intercept only)= -3823.2
# Chisq= 273.24 on 8 degrees of freedom, p= 0

# Osserviamo che il modello contiene variabili diverse rispetto a quello di Cox.
# Provando la procedura backward conil criterio del p-value, si ottiene lo stesso modello :
fp1 <- survreg(Surv(lifetime,event)~ ad_chemio +pre_hormon +pre_chemio + n_sites+
                  perform + eta , data=e1,dist="weibull")
summary(fp1)


# l'AIC è calcolato come segue:
AICw1 <- -2*fp1$loglik[2] +2*(8+2)      #  8 coeff di regressione, due parametri (mu, sigma) del modello Weibull
AICw1


# Confrontiamo il modello Weibull con quello esponenziale, sia tramite AIC sia tramite il test del rapporto delle log-verosimiglianze
fp2 <- survreg(Surv(lifetime,event)~ ad_chemio +pre_hormon +pre_chemio + n_sites+
                 perform + eta , data=e1,scale=1)
summary(fp2)
AICw2 <- -2*fp2$loglik[2] +2*(8+1)      #  8 coeff di regtessione, due parametri (mu, sigma) del modello Weibull
AICw2

# È più basso l'AIC del modello Weibull, quindi scegliamo quest'ultimo
# Secondo il test del rapporto delle log-verosimiglianze, otteniamo un risultato simile
# (per poco non significativo, ma comunque vicino al borderline)
LLR <- 2*( fp1$loglik[2] - fp2$loglik[2])
1-pchisq(LLR, df=1)
#0.06940241
# Concluderemmo a favore del modello Weibull





############   ADEGUATEZZA DEL MODELLO    ###############


# dovremmo verificare l'adeguatezza del modello rispetto all'ipotesi di proporzionalit� dei rischi e 
# all'ipotesi distributiva Weibull.
# La prima assunzione è in realtà già stata testata quando trattavamo il modello di Cox.
# Per la seconda assunzione (distr. Weibull), possiamo testare globalmente il modello, oppure confrontare le curve di sopravvivenza sotto l'ipotesi Weibull
# con le curve stimate non parametricamente:


# ad_chemio
s1 <- survfit(Surv(lifetime,event)~ ad_chemio)
plot(s1,conf.int=F,col=c(1,2), lty=2, ylab="Survival probability", xlab="Tempo (mesi)")
library(survminer)
ggsurvplot(s1, data = e, pval = TRUE, pval.method = TRUE, log.rank.weights = "1", conf.int = FALSE)
w1 <- survreg(Surv(lifetime,event)~ ad_chemio, data=e1,dist="weibull")
mu <-  w1$coefficients[1]
sigma <- w1$scale
beta <-  w1$coefficients[-1]
alpha <- 1/sigma            # stime (da qua in giu) col "ripar"
lambda <- exp(-mu/sigma)
gamma <- - beta /sigma
HR <- exp(gamma)

# osservo il min e max di sort(s1$time)
tt <- seq(0,235,by=0.01)
S0 <- exp(-lambda* tt^alpha)
S1 <- S0^HR
lines(tt,S0,lty=2, col=1, type="l")
lines(tt,S1, lty=2, col=2, type="l")

# confronto log(-log S(t))
plot(s1,conf.int=F,fun="cloglog", col=c(1,2), lty=2, ylab="Survival probability", xlab="Tempo (mesi)")
lines(tt,log(-log(S0)),lty=2, col=1, type="l")
lines(tt,log(-log(S1)), lty=2, col=2, type="l")



## pre_hormon
s4 <- survfit(Surv(lifetime,event)~ pre_hormon)
plot(s4,conf.int=F,col=c(1,2), lty=2, ylab="Survival probability", xlab="Tempo (mesi)")

w4 <- survreg(Surv(lifetime,event)~ pre_hormon, data=e1,dist="weibull")
mu <-  w4$coefficients[1]
sigma <- w4$scale
beta <-  w4$coefficients[-1]
alpha <- 1/sigma
lambda <- exp(-mu/sigma)
gamma <- - beta /sigma
HR <- exp(gamma)
# osservo il min e max di sort(s4$time)
tt <- seq(0,235,by=0.01)
S0 <- exp(-lambda* tt^alpha)
S1 <- S0^HR
lines(tt,S0,lty=2, col=1, type="l")
lines(tt,S1, lty=2, col=2, type="l")

# confronto log(-log S(t))
plot(s4,conf.int=F,fun="cloglog", col=c(1,2), lty=2, ylab="Survival probability", xlab="Tempo (mesi)")
lines(tt,log(-log(S0)),lty=2, col=1, type="l")
lines(tt,log(-log(S1)), lty=2, col=2, type="l")




#############   VERIFICA GLOBALE DEL MODELLO WEIBULL - RESIDUI COX-SNELL
#########   riparametrizziamo, passando alla rappresentazione a rischi moltiplicativi del modello fp1
mu <-  fp1$coefficients[1]
sigma <- fp1$scale
beta <- fp1$coefficients[-1]
alpha <- 1/sigma
lambda <- exp(-mu/sigma)
gamma <- - fp1$coefficients[-1] /sigma

HR <- exp(gamma) 

#######  svolgiamo un'analisi globale del modello Weibull usando i residui di Cox-Snell
stand_resid <- (log(e1$lifetime)- fp1$linear)/sigma   # r_i^S
# per passare ai residui di Cox-Snell, uso la funzione di rischio cumulato della W ~ valori estremi standard
#  S(t) = exp(-exp(w))
#  H(t) = - log(exp(-exp(w))) = exp(w)
cox_snell <- exp(stand_resid)   # r_i
# proseguo come fatto per il modello di Cox, costruendo il grafico della f. di rischio cumulato degli r_i,  verso  t_i
fitr <-survfit(Surv(cox_snell,e1$event)~1)
plot(fitr$time, -log(fitr$surv),type="s",xlab="Residui di Cox-Snell",ylab="Funzione di rischio cumulato stimata")
abline(0,1,col="red",lty=2)
# chiaramente il modello Weibull si adatta molto male ai dati,
# e sembra di gran lunga peggiore rispetto al modello di Cox.


 
