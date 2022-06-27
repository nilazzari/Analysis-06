library(survival)
library(MASS)
library(ggsurvplot)
e <-read.table("epi.dat", header=TRUE)

# factors
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
# in realta una verifica grafica della proporzionalita dei rischi, è gia stata fatta tramite stimatori non parametrici
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

table(event)
#  0    1 
#  56 1041 

table(cardiotox)   # 1= si
#   1   2 
#  125 962 

table(ad_chemio)   # 1=si
summary(lifetime)
summary(cardiotime)
# etc. per ogni variabile esplicativa
# censura tipo I generico e censura casuale.
# Per i modelli di regressione, assumiamo che la censura sia di tipo 'indipendente'.

# NONPARAMETRIC ANALYSIS
s1 <- survfit(Surv(lifetime,event)~ ad_chemio)
s1[1]$surv  # stime di sopravvivenza per gruppo ad_chemio=1
s1[2]$surv  # stime di sopravvivenza per gruppo ad_chemio=2
plot(s1,conf.int=F,col=c(1,2), ylab="Survival probability", xlab="Tempo (mesi)", lwd = 2)
ggsurvplot(s1, data = e, pval = TRUE, conf.int = FALSE, pval.method = TRUE, log.rank.weights = "1")
title("ad_chemio")
# nero: si, rosso: no
legend(200,1, c("si","no"), lty=c(1,1), col=c("black","red"))
# test dei ranghi logaritmici
survdiff(Surv(lifetime,event)~ ad_chemio)
# Chisq= 6.7  on 1 degrees of freedom, p= 0.00968 
# la differenza tra le curve stimate di sopravvivenza nei due gruppi
# (con ad_chemio, senza ad_chemio) è tale da spiegare che esiste una 
# relazione importante tra la variabile ad_chemio ed la sopravvivenza (sopravvivenza signif. più bassa per i pazienti sottoposti a chemioterapia adiuvante prima dell'inizio dello studio)
# notiamo curve proporzionali 

s2 <- survfit(Surv(lifetime,event) ~ ad_hormon)
ggsurvplot(s2, data = e, pval = TRUE, conf.int = FALSE, pval.method = TRUE, log.rank.weights = "1")
survdiff(Surv(lifetime,event)~ ad_hormon)
# la variabile non sembra rilevante per la sopravvivenza
# notiamo curve non proporzionali
survdiff(Surv(lifetime,event)~ ad_hormon, rho = 1)

s3 <-survfit(Surv(lifetime,event)~ cur_radiation)
ggsurvplot(s3, data = e, pval = TRUE, conf.int = FALSE, pval.method = TRUE, log.rank.weights = "1")
title("cur_radiation")
survdiff(Surv(lifetime,event)~ cur_radiation)
# la variabile non sembra molto rilevante per la sopravvivenza

s4 <- survfit(Surv(lifetime,event)~ pre_hormon)
ggsurvplot(s4, data = e, pval = TRUE, conf.int = FALSE, pval.method = TRUE, log.rank.weights = "1")
survdiff(Surv(lifetime,event)~ pre_hormon)
# la variabile sembra importante 

s5 <- survfit(Surv(lifetime,event)~ pre_chemio)
ggsurvplot(s5, data = e, pval = TRUE, conf.int = FALSE, pval.method = TRUE, log.rank.weights = "1")
survdiff(Surv(lifetime,event)~ pre_chemio)
# la variabile sembra importante 
# difficoltà nel valutare la proporzionalità tra le curve

s6 <- survfit(Surv(lifetime,event) ~ fam_heart)
ggsurvplot(s6, data = e, pval = TRUE, conf.int = FALSE, pval.method = TRUE, log.rank.weights = "1")
survdiff(Surv(lifetime,event)~ fam_heart)
# la variabile non sembra rilevante per la sopravvivenza

s7 <- survfit(Surv(lifetime,event)~ n_sites)
ggsurvplot(s7, data = e, pval = TRUE, conf.int = FALSE, pval.method = TRUE, log.rank.weights = "1")
survdiff(Surv(lifetime,event)~ n_sites)
# la variabile sembra molto importante: la sopravvivenza peggiora all'aumentare del numero di siti del tumore

summary(surface)
surf <-cut(surface,breaks=c(0,1.6,1.7,1.8,2.3))
table(surf)
s8 <- survfit(Surv(lifetime,event)~surf)
ggsurvplot(s8, data = e, pval = TRUE, conf.int = FALSE, pval.method = TRUE, log.rank.weights = "1")
survdiff(Surv(lifetime,event) ~ surf)
# le curve si incrociano e non sembrano proporzionali:
survdiff(Surv(lifetime,event)~surf, rho=1)
# al diminuire della superficie corporea, diminuisce la sopravvivenza

s9 <- survfit(Surv(lifetime,event)~ performance)
ggsurvplot(s9, data = e, pval = TRUE, conf.int = TRUE, pval.method = TRUE, log.rank.weights = "1")
survdiff(Surv(lifetime,event)~ performance)
#  all'aumentare dello stato, peggiora la sopravvivenza
# il terzo livello del fattore contiene poche unità (11), lo accorpiamo al secondo
perform <- performance
levels(perform) <- c("1","2","2")
s9b <-survfit(Surv(lifetime,event)~ perform)
plot(s9b,conf.int=F,col=c(1,2,3), ylab="Survival probability", xlab="Tempo (mesi)")
survdiff(Surv(lifetime,event)~ perform) # il pvalue non cambia di molto

s10 <- survfit(Surv(lifetime,event)~ ad_radiation)
ggsurvplot(s10, data = e, pval = TRUE, conf.int = FALSE, pval.method = TRUE, log.rank.weights = "1")
survdiff(Surv(lifetime,event)~ ad_radiation)
# le curve si incrociano e non sembrano proporzionali:
survdiff(Surv(lifetime,event)~ ad_radiation,rho=1)
# la variabile sembra rilevante per la sopravvivenza

s10 <- survfit(Surv(lifetime,event)~ rig_radiation)
plot(s10,conf.int=F,col=c(1,2), ylab="Survival probability", xlab="Tempo (mesi)")
survdiff(Surv(lifetime,event)~ rig_radiation)
# le curve si incrociano e non sembrano proporzionali:
survdiff(Surv(lifetime,event)~ rig_radiation,rho=2)
# la variabile non sembra rilevante per la sopravvivenza

s11 <- survfit(Surv(lifetime,event)~ medi_radiation)
plot(s11,conf.int=F,col=c(1,2), ylab="Survival probability", xlab="Tempo (mesi)")
survdiff(Surv(lifetime,event)~ medi_radiation)
# le curve si incrociano e non sembrano proporzionali:
survdiff(Surv(lifetime,event)~ medi_radiation,rho=1)

s12 <- survfit(Surv(lifetime,event)~ lef_radiation)
plot(s12,conf.int=F,col=c(1,2), ylab="Survival probability", xlab="Tempo (mesi)")
survdiff(Surv(lifetime,event)~ lef_radiation)
# le curve si incrociano e non sembrano proporzionali:
survdiff(Surv(lifetime,event)~ lef_radiation, rho=1)

summary(eta)
etac <-cut(eta,breaks=c(0,49,56,63,79))
table(etac)
s13 <- survfit(Surv(lifetime,event)~ factor(etac))
plot(s13,conf.int=F,col=c(1,2,3,4), ylab="Survival probability", xlab="Tempo (mesi)")
survdiff(Surv(lifetime,event)~etac)
# i pazienti più anziani hanno una sopravvivenza signif. minore rispetto ai pzienti più giovani.

# cardiotox è una variabile tempo-dipendente (si verifica durante lo studio).
# Consapevoli di ciò, studiamo ugualmente le curve di sopravvivenza nei due gruppi con i metodi che conosciamo  
s14 <- survfit(Surv(lifetime,event)~ cardiotox)
plot(s14,conf.int=F,col=c(1,2), ylab="Survival probability", xlab="Tempo (mesi)")
title("cardiotox")
# nero: si, rosso: no
legend(200,1, c("si","no"), lty=c(1,1), col=c("black","red"))
survdiff(Surv(lifetime,event)~ cardiotox)
# Chisq= 9.2  on 1 degrees of freedom, p= 0.00239
# C'è una differenza significativa nella sopravvivenza dei pazienti senza e 
# con cardiotossicità. Attenzione: cardiotox e mortalità sono 'rischi competitivi' 
# in prima battuta, sembrerebbe che chi ha avuto l'evento di cardiotox sopravvive più a lungo (ci aspettiamo in realtà il risultato opposto!! ). 
# Osserviamo ciò a causa della presenza dei rischi competitivi. Il decesso 
# (precoce, si guardi infatti il primo tratto della curva rossa) del paziente 
# preclude per lui l'osservazione di un eventuale evento di cardiotossicità.
# Proviamo a studiare il tempo di attesa, arrestato all'evento di cardiotossicità (quando si verifica)
s14c <- survfit(Surv(cardiotime,event)~ cardiotox)
plot(s14c,conf.int=F,col=c(1,2), ylab="Survival probability", xlab="Tempo (mesi)")
title("cardiotox")
survdiff(Surv(cardiotime,event)~ cardiotox)

# la dose totale cumulata è una variabile tempo-dipendente
summary(epi_m2)
cumdose <-cut(epi_m2,breaks=c(0,500,1000,1600))
table(cumdose)
s15 <- survfit(Surv(lifetime,event)~cumdose)
plot(s15,conf.int=F,col=c(1,2,3), ylab="Survival probability", xlab="Tempo (mesi)")
title("cumdose")
# nero: si, rosso: no
legend(150,1, c("(0,500]","(500,1000]",">1000"), lty=c(1,1,1), col=c("black","red","green"))
survdiff(Surv(lifetime,event)~cumdose)
#  Chisq= 337  on 2 degrees of freedom, p= 0 

hist(epi_m2)
cumdose <-cut(epi_m2,breaks=c(0,470,870,970,1600))
table(cumdose)
s15 <- survfit(Surv(lifetime,event)~cumdose)
plot(s15,conf.int=F,col=c("black","red","green","blue"), ylab="Survival probability", xlab="Tempo (mesi)")
title("cumdose")
# nero: si, rosso: no
legend(150,1, c("(0,470]","(470,870]","(870,970]",">970"), lty=c(1,1,1,1), col=c("black","red","green","blue"))
survdiff(Surv(lifetime,event)~cumdose)
# Le variabili che mostrano differenze significative nella funzione di sopravvivenza sono:
# ad_chemio (s1), cur_radiation (?) (s3), pre_hormon (s4), pre_chemio (s5), n_sites (s7), surface (s8), 
# perform (s9b), ad_radiation (s10), medi_radiation_(?) (s11), left_radiation (?) (s12), eta (s13), cumdose (s15)
# cardiotox è anche fortemente significativa, ma non verrà inclusa nel modello perchè violerebbe l'assunzione di censura indipendente.

#  verifico quali modelli possono essere plausibili, tramite una verifica grafica delle assunzioni 
# 1) verifica sui rischi proporzionali 
plot(s1,col=1:2,fun="cloglog",xlab="log(t)",ylab="log(-log(S(t)))")  # dubbio
plot(s3,col=1:2,fun="cloglog",xlab="log(t)",ylab="log(-log(S(t)))")   # OK
plot(s4,col=1:2,fun="cloglog",xlab="log(t)",ylab="log(-log(S(t)))")  # OK
plot(s5,col=1:2,fun="cloglog",xlab="log(t)",ylab="log(-log(S(t)))")    # dubbio
plot(s7,col=1:4,fun="cloglog",xlab="log(t)",ylab="log(-log(S(t)))")    # no
plot(s8,col=1:4,fun="cloglog",xlab="log(t)",ylab="log(-log(S(t)))")    # no
plot(s9b,col=1:2,fun="cloglog",xlab="log(t)",ylab="log(-log(S(t)))")     # no
plot(s10,col=1:2,fun="cloglog",xlab="log(t)",ylab="log(-log(S(t)))")    # no
plot(s11,col=1:2,fun="cloglog",xlab="log(t)",ylab="log(-log(S(t)))")    # no
plot(s12,col=1:2,fun="cloglog",xlab="log(t)",ylab="log(-log(S(t)))")   # no
plot(s13,col=1:4,fun="cloglog",xlab="log(t)",ylab="log(-log(S(t)))")   # dubbio
plot(s15,col=1:4,fun="cloglog",xlab="log(t)",ylab="log(-log(S(t)))")   # NO!

#  MODELLO SEMI-PARAMETRICO DI COX
f1 <- coxph(Surv(lifetime,event)~ ad_chemio, ties="breslow",data=e)
summary(f1)  # si
# Akaike information criterion: AIC = -2 log L + k* p
#  k=2 sempre 
# p = numero di coeff regressione nel modello
AIC1 <- -2*f1$loglik +2*c(0,1)
AIC1
lrt1 <- 2*(f1$loglik[2] - f1$loglik[1])
lrt1
f2 <- coxph(Surv(lifetime,event)~ ad_chemio +pre_hormon, ties="breslow",data=e)
summary(f2)  # si
AIC2 <- -2*f2$loglik +2*c(0,2)
AIC2
lrt2 <- 2*(f2$loglik[2] - f2$loglik[1])
lrt2
## local test per beta_2 =0
lrt2 - lrt1   #  20.3568
1 -pchisq( 20.3568,df=1)

f3 <- coxph(Surv(lifetime,event)~ ad_chemio +pre_hormon +n_sites, ties="breslow",data=e)
summary(f3)  # si
AIC3 <- -2*f3$loglik +2*c(0,5)
AIC3

# metodo veloce che ci calcola l'AIC per tutti i possibili sottogruppi di variabili esplicative
# ed effettua una procedura stepwise per la selezione delle variabili nel modello 
e1 <-na.omit(e) 
e1$perform <- e1$performance
e1$cumdose <-cut(e1$epi_m2,breaks=c(0,470,870,970,1600))
levels(e1$perform) <- c("1","2","2")
ff <- coxph(Surv(lifetime,event)~ ad_chemio +cur_radiation +pre_hormon +pre_chemio + n_sites+
       surface + perform +ad_radiation + medi_radiation + lef_radiation + eta + cumdose, ties="breslow",data=e1)
stepAIC(ff)

# coef exp(coef) se(coef)      z       p
# ad_chemio2           -0.15577   0.85575  0.07299  -2.13  0.0328
# n_sites2              0.59091   1.80564  0.08527   6.93 4.2e-12
# n_sites3              0.71622   2.04669  0.09686   7.39 1.4e-13
# n_sites4              0.74519   2.10684  0.11719   6.36 2.0e-10
# surface              -0.43619   0.64649  0.20290  -2.15  0.0316
# perform2              0.44881   1.56645  0.08243   5.44 5.2e-08
# medi_radiation2       0.31684   1.37278  0.11721   2.70  0.0069
# eta                   0.01008   1.01013  0.00371   2.72  0.0066
# cumdose(470,870]     -0.95574   0.38453  0.09742  -9.81 < 2e-16
# cumdose(870,970]     -1.41261   0.24351  0.09998 -14.13 < 2e-16
# cumdose(970,1.6e+03] -1.38001   0.25158  0.10395 -13.28 < 2e-16

# Criterio backward del p-value:
f1 <- coxph(Surv(lifetime,event)~ ad_chemio +cur_radiation +pre_hormon +pre_chemio + n_sites+
              surface + perform +ad_radiation + medi_radiation + lef_radiation + eta + cumdose, ties="breslow",data=e1)
summary(f1)
# via ad_radiation
# via cur_radiation
# via pre_chemio
# via lef_radiation
# via pre_hormon
# Modello finale
f2 <- coxph(Surv(lifetime,event)~ ad_chemio + n_sites + surface + perform + 
              medi_radiation + eta + cumdose, ties="breslow",data=e1)
summary(f2)

#                           coef exp(coef)  se(coef)       z Pr(>|z|)    
# ad_chemio2             -0.155771  0.855755  0.072990  -2.134  0.03283 *  
#   n_sites2              0.590913  1.805636  0.085272   6.930 4.22e-12 ***
#   n_sites3              0.716222  2.046686  0.096864   7.394 1.42e-13 ***
#   n_sites4              0.745190  2.106842  0.117186   6.359 2.03e-10 ***
#   surface              -0.436194  0.646493  0.202905  -2.150  0.03158 *  
#   perform2              0.448810  1.566448  0.082433   5.445 5.19e-08 ***
#   medi_radiation2       0.316840  1.372783  0.117213   2.703  0.00687 ** 
#   eta                   0.010080  1.010131  0.003711   2.716  0.00661 ** 
#   cumdose(470,870]     -0.955738  0.384528  0.097417  -9.811  < 2e-16 ***
#   cumdose(870,970]     -1.412611  0.243507  0.099976 -14.129  < 2e-16 ***
#   cumdose(970,1.6e+03] -1.380013  0.251575  0.103952 -13.275  < 2e-16 ***

# stesso risultato della procedura stepwise basata sull'AIC

dropterm(f3,test = "Chisq")
f3 <- coxph(Surv(lifetime,event)~ ad_chemio  +n_sites+
              surface + perform + medi_radiation +  eta + cumdose, ties="breslow",data=e1)
summary(f3)
# si arriva allo stesso modello
anova(f1,f2)  # OK
# f2 è il modello finale, per il quale occore studiare il buon adattamento ai dati, e la verifica delle asunzioni sottostanti
#  NOTA: la selezione delle varibili è avvenuta eliminando tutti i valori mancanti per le var esplicative --> risultati diversi
# NOTA:  cumdose spiega buona parte del modello perchè è chiaro che chi ha avuto dosi maggiori di chemioterapia, 
# è vissuto più a lungo!
# NOTA: dalle analisi univariate (stime nonparam. delle curve di sopravvivenza, oppure modello Cox con una variabile)
# erano risultate significative anche le variabili: pre_hormon, pre_chemio, ad_radiation
# Perciò vediamo che succede se togliamo "cumdose" dal modello di partenza:
ff1 <- coxph(Surv(lifetime,event)~ ad_chemio +cur_radiation +pre_hormon +pre_chemio + n_sites+
              surface + perform +ad_radiation + medi_radiation + lef_radiation + eta , ties="breslow",data=e1)
stepAIC(ff1)
#                   coef exp(coef) se(coef)     z       p
# ad_chemio2    -0.31920   0.72673  0.07591 -4.21 2.6e-05
# pre_chemio2   -0.41059   0.66326  0.10885 -3.77 0.00016
# n_sites2       0.60806   1.83687  0.08463  7.19 6.7e-13
# n_sites3       0.74410   2.10454  0.09604  7.75 9.3e-15
# n_sites4       0.80105   2.22788  0.11716  6.84 8.1e-12
# surface       -0.42649   0.65279  0.20539 -2.08 0.03785
# perform2       0.63875   1.89410  0.08100  7.89 3.1e-15
# ad_radiation2 -0.12443   0.88300  0.07550 -1.65 0.09935
# eta            0.01197   1.01204  0.00359  3.34 0.00085
# 
# Likelihood ratio test=227  on 9 df, p=0
# Osserviamo che ora  pre_chemio, ad_radiation sono entrate nel modello!
# La sopravvivenza dipende dalla dose assunta, e viceversa (non è chiara la relazione di causa-effetto!),
# tuttavia, sia la dose totale assunta sia la sopravvivenza dipendno anche dalla cardiotossicità (eventi cardiotox)  

# VALIDAZIONE
# verifichiamo globalmente il modello tramite i residui di Cox-snell:
# consideriamo il dataset originario senza dati mancanti
f2 <- coxph(Surv(lifetime,event)~ ad_chemio + n_sites + surface + perform + 
              medi_radiation + eta + cumdose, ties="breslow",data=e1)
summary(f2)
?residuals.coxph
# molti tipi di residui: # martingale, deviance, score, schoenfeld, dfbeta, dfbetas, scaledsch, partial
# residui Cox-snell:
mres <- resid(f2, type="martingale")
csres <- e1$event-mres
r.surv <- survfit(Surv(csres,event) ~1, type="fleming-harrington",data=e1)
plot(r.surv$time, -log(r.surv$surv), type="s", xlab="Cox-snell Residual",
     ylab="Estimated Cum Hazards")
lines(c(0,10),c(0,10),col=2)
# chiaramente il modello di Cox non si adatta bene
# Qualche assunzione sottostante il modello non rispetta l'andamento spiegato dai dati
# il più importante aspetto da verificare è l'assunzione di proporzionalità dei rischi.
# abbiamo già verificato con il grafico non parametrico log(-log) ciascuna variabile singolarmente, capendo 
# che diverse variabili violano l'assunzione di proporzionalità dei rischi.
# usiamo anche la strategia della stratificazione: consideriamo il modello finale dove la variabile da testare è
# usata come variabile di stratificazione:
a1 <- coxph(Surv(lifetime,event)~ strata(ad_chemio) + n_sites + surface + perform + 
              medi_radiation + eta + cumdose, ties="breslow",data=e1)
summary(a1)
ch <- basehaz(a1, centered=F)
ch1 <- ch[ch$strata=="1",]
ch2 <- ch[ch$strata=="2",]
plot(log(ch1$time), log(ch1$hazard), type="s")
lines(log(ch2$time), log(ch2$hazard), type="s",col=2)
# rette appross parallele: no violazione causata da ad_chemio

a2 <- coxph(Surv(lifetime,event)~ ad_chemio + strata(n_sites) + surface + perform + 
              medi_radiation + eta + cumdose, ties="breslow",data=e1)
summary(a2)
ch <- basehaz(a2, centered=F)
ch1 <- ch[ch$strata=="1",]
ch2 <- ch[ch$strata=="2",]
ch3 <- ch[ch$strata=="3",]
ch4 <- ch[ch$strata=="4",]
plot(log(ch1$time), log(ch1$hazard), type="s")
lines(log(ch2$time), log(ch2$hazard), type="s",col=2)
lines(log(ch3$time), log(ch3$hazard), type="s",col=3)
lines(log(ch4$time), log(ch4$hazard), type="s",col=4)
# curve appross parallele, ad eccezione della rossa: dubbio per n_sites
......
#  così via per ciascuna variabile (le variabili quantitative devo essere categorizzate!)
......
e1$surf <- cut(e1$surface,breaks=c(0,1.6,1.7,1.8,2.3))
levels(e1$surf) <- c(1:4)
a2 <- coxph(Surv(lifetime,event)~ ad_chemio + n_sites + strata(surf) + perform + 
              medi_radiation + eta + cumdose, ties="breslow",data=e1)
summary(a2)
ch <- basehaz(a2, centered=F)
ch1 <- ch[ch$strata=="1",]
ch2 <- ch[ch$strata=="2",]
ch3 <- ch[ch$strata=="3",]
ch4 <- ch[ch$strata=="4",]
plot(log(ch1$time), log(ch1$hazard), type="s")
lines(log(ch2$time), log(ch2$hazard), type="s",col=2)
lines(log(ch3$time), log(ch3$hazard), type="s",col=3)
lines(log(ch4$time), log(ch4$hazard), type="s",col=4)
# rette appross parallele: no violazione causata da surf
# rette non sono parallele per perform, medi_radiation, cumdose.
# Queste variabili sono responsabili della violazione dell'assunzione di prop rischi 

# SCHOENFELD RESIDUALS ANALYSIS
# Per valutare la proporzionalità, effettuiamo anche un'analisi dei residui di Schoenfeld:
cox.zph(f2, global=TRUE)
# cambiamo la trasformazione g(t)= theta * t (forma lineare):
cox.zph(f2, transform="identity", global=TRUE)
# cambiamo la trasformazione g(t)= theta * log(t) (forma logaritmica):
cox.zph(f2, transform="log", global=TRUE)
# ora ad_chemio sembra anche violare lássunto di proporzionalità
# Il test fornisce risultato significativo al 5%, per le variabili ad_chemio, perform, medi_radiation, cumdose.
# Il test globale indica che il modello non si adatta bene ai dati

# residui di Schoenfeld:
s <- cox.zph(f2, transform="log", global=TRUE)
View(s$y)  # matrice D*p dei residui
s*x  # tempi trasformati = g(t)
cor(s$y,s$x)   # prima colonna (rho) dell'output

# Unendo tutte le informazioni fin qui raccolte con i vari metodi di verifica sul modello, decidiamo di 
# (1) eliminare cumdose dal modello, essendo la principale causa di violazione. Essa è una variabile tempo-dipendente
# Il modello finale ottenuto senza cumdose, dalla procedura stepwise con AIC era:
f3 <- coxph(Surv(lifetime,event)~ pre_chemio + ad_radiation + n_sites + surface + perform + ad_radiation + eta, ties="breslow",data=e1)
summary(f3)
cox.zph(f3, transform="log", global=TRUE)

# (2) la variabile maggiormente responsabile della violazione è perform. 
# Questa variabile indica lo stato di gravità dei soggetti. Decidiamo di usarla come fattore di stratificazione:
# eliminiamo medi_radiation dal modello
f4 <- coxph(Surv(lifetime,event)~ ad_chemio + pre_chemio + n_sites + surface + strata(perform) + ad_radiation + eta, ties="breslow",data=e1)
summary(f4)
cox.zph(f4, transform="log", global=TRUE)
# la variabile surface è ora l'unica responsabile della violazione del modello.
# Tuttavia notiamo che il test globale non è significativo. Ciò indica che possiamo globalmente accettare il modello
# come adeguato:
#                   rho  chisq      p
# ad_chemio2    -0.0391  1.381 0.2399
# pre_chemio2   -0.0156  0.219 0.6395
# n_sites2       0.0148  0.195 0.6585
# n_sites3      -0.0302  0.811 0.3680
# n_sites4      -0.0193  0.339 0.5605
# surface        0.0868  6.566 0.0104
# ad_radiation2  0.0363  1.202 0.2728
# eta            0.0153  0.199 0.6554
# GLOBAL             NA 10.562 0.2278

#  Potremmo poi stimare le curve di sopravvivenza per varie tipologie di pazienti:
mydata <- data.frame(ad_chemio="2",pre_chemio="1",n_sites="1", surface=1.7,  ad_radiation="2", eta=55)
mydata1 <- data.frame(ad_chemio="1",pre_chemio="1",n_sites="1", surface=1.7, ad_radiation="2", eta=55)
fit1 <- survfit(f4, newdata= mydata)
fit2 <- survfit(f4, newdata= mydata1)
## perform==1
plot(fit1[1]$time, fit1[1]$surv, type="s", ylab="Survival curve",xlab="Time (months)")  # senza chemio adiuvante
lines(fit2[1]$time, fit2[1]$surv, col="red")   # con chemio adiuvante

## perform==2
lines(fit1[2]$time, fit1[2]$surv, lty=2, type="s", ylab="Survival curve",xlab="Time (months)")  # senza chemio adiuvante
lines(fit2[2]$time, fit2[2]$surv,  lty=2, col="red")   # con chemio adiuvante

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

# è più basso l'AIC del modello Weibull, quindi scegliamo quest'ultimo
# Secondo il test del rapporto delle log-verosimiglianze, otteniamo un risultato simile
# (per poco non significativo, ma comunque vicino al borderline)
LLR <- 2*( fp1$loglik[2] - fp2$loglik[2])
1-pchisq(LLR, df=1)
#0.06940241
# Concluderemmo a favore del modello Weibull

# VALIDAZIONE
# dovremmo verificare l'adeguatezza del modello rispetto all'ipotesi di proporzionalità dei rischi e 
# all'ipotesi distributiva Weibull.
# La prima assunzione Ã¨ in realtÃ  giÃ  stata testata quando trattavamo il modello di Cox.
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

# pre_hormon
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

# VERIFICA GLOBALE DEL MODELLO WEIBULL - RESIDUI COX-SNELL
# riparametrizziamo, passando alla rappresentazione a rischi moltiplicativi del modello fp1
mu <-  fp1$coefficients[1]
sigma <- fp1$scale
beta <- fp1$coefficients[-1]
alpha <- 1/sigma
lambda <- exp(-mu/sigma)
gamma <- - fp1$coefficients[-1] /sigma
HR <- exp(gamma) 

# svolgiamo un'analisi globale del modello Weibull usando i residui di Cox-Snell
stand_resid <- (log(e1$lifetime)- fp1$linear)/sigma   # r_i^S
# per passare ai residui di Cox-Snell, uso la funzione di rischio cumulato della W ~ valori estremi standard
#  S(t) = exp(-exp(w))
#  H(t) = - log(exp(-exp(w))) = exp(w)
cox_snell <- exp(stand_resid)   # r_i
# proseguo come fatto per il modello di Cox, costruendo il grafico della f. di rischio cumulato degli r_i,  verso  t_i
fitr <-survfit(Surv(cox_snell,e1$event)~1)
plot(fitr$time, -log(fitr$surv),type="s",xlab="Residui di Cox-Snell",ylab="Funzione di rischio cumulato stimata")
abline(0,1,col="red",lty=2)
# Il modello di Cox esibisce un adattamento ampiamente migliore
