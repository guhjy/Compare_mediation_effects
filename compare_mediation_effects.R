#################################################################################
#  compare_mediation_effects.R                                                  # 
#                                                                               #
#  Simulation code for                                                          #
# "Understanding interventional effects: a more natural approach to mediation   #
#  analysis?"                                                                   #
#                                                                               #
#  Margarita Moreno-Betancur, 18 December 2017                                  #
#                                                                               #
#################################################################################

library(boot)


### Set parameters

set.seed(14187)
n<-1000000
pp0<-0.2  #(Prevalence of M if A=0 and L=0)=(Prevalence of L if A=0)=(Prevalence of Y if A=L=M=0)
pp1<-0.3  #Prevalence of M if A=1 and L=0


tot<-data.frame()
ind<-data.frame()
dir<-data.frame()

pL<-NULL
pM<-NULL
pout<-NULL
rem<-NULL

 for(beta in seq(-3,3,1))
 {

te<-vector()
ie<-vector()
de<-vector()

### Simulate trial arm (G), exposure (A), first mediator (L), second mediator (M) and outcome (Y)
### according to characteristics of each arm in Table 1 of main text.

G<-apply(rmultinom(n,1,prob=rep(1/7,7))==1,2,which)
A<-ifelse(G==1|G==7,0,1)
L<-rbinom(n,1,inv.logit(logit(pp0)+beta*(!G%in%c(1,3,4,7)))) 
M<-rbinom(n,1,ifelse(G%in%c(1,3),inv.logit(logit(pp0)+beta*L),
              ifelse(G%in%c(4,5,7), inv.logit(logit(pp0)+beta)*pp0+pp0*(1-pp0),
              ifelse(G%in%c(2),inv.logit(logit(pp1)+beta*L),
                                 inv.logit(logit(pp1)+beta)*inv.logit(logit(pp0)+beta)+
                                 inv.logit(logit(pp1))*(1-inv.logit(logit(pp0)+beta))))))
                     
Y<-rbinom(n,1,inv.logit(logit(pp0)-2*A+3*M+beta*L))


### Calculate probabilities of mediators and outcome for each arm

pL<-cbind(pL,(table(G,L)/cbind(table(G),table(G)))[,2])
pM<-cbind(pM,(table(G,M)/cbind(table(G),table(G)))[,2])
pr<-(table(G,Y)/cbind(table(G),table(G)))[,2]
pout<-cbind(pout,pr)

### Calculate/estimate effects

# Calculate total causal effect
TCE<-pr[2]-pr[1]

# Calculate M-M interventional effects 
IDE<-pr[3]-pr[1]
IIEL<-pr[5]-pr[4]
IIEM<-pr[6]-pr[5]
REM<-(pr[2]-pr[6])-(pr[3]-pr[4]) 

te<-c(te,sum(IDE,IIEL,IIEM,REM))
de<-c(de,IDE+IIEL)
ie<-c(ie,IIEM+REM)
rem<-c(rem,REM)

# Calculate S-B interventional effects
IDE_S<-pr[5]-pr[7]
IIE_S<-pr[6]-pr[5]
te<-c(te,sum(IDE_S,IIE_S))
de<-c(de,IDE_S)
ie<-c(ie,IIE_S)

# Estimate natural effects with g-computation (ignoring L)

  # fit outcome and mediator models based on data from 
fit<-glm(Y~A+M,family=binomial, data=data.frame(Y,A,M)[G%in%c(1,2),])
fit1<-glm(M~A,family=binomial, data=data.frame(A,M)[G%in%c(1,2),])

NDE<-    (predict(fit,newdata=data.frame(A=1,M=1),type="response")-
          predict(fit,newdata=data.frame(A=0,M=1),type="response"))*
          predict(fit1,newdata=data.frame(A=0),type="response")+
          (predict(fit,newdata=data.frame(A=1,M=0),type="response")-
          predict(fit,newdata=data.frame(A=0,M=0),type="response"))*
          (1-predict(fit1,newdata=data.frame(A=0),type="response"))

NIE<-    (predict(fit,newdata=data.frame(A=1,M=1),type="response")*
         (predict(fit1,newdata=data.frame(A=1),type="response")-
          predict(fit1,newdata=data.frame(A=0),type="response")))+
       (predict(fit,newdata=data.frame(A=1,M=0),type="response")*
       ((1-predict(fit1,newdata=data.frame(A=1),type="response"))-
        (1-predict(fit1,newdata=data.frame(A=0),type="response"))))

  # # Equivalent:
  #  library(mediation)
  #  med<-mediate(model.m=fit1, model.y=fit, treat = "A", mediator = "M")

te<-c(te,sum(NDE,NIE))
de<-c(de,NDE)
ie<-c(ie,NIE)

tot<-rbind(tot,te)
ind<-rbind(ind,ie)
dir<-rbind(dir,de)
 }

### Plot Figure Supp1

png("FigureSupp1.png",width=1000,height=400)
yran<-c(-0.25,0.35)
par(mfrow=c(1,3),cex=1.2)
plot(seq(-3,3),ind[,1],type="b",lty=1,main="Indirect effect (via M)",ylim=yran,ylab="Risk difference",
     xlab="",pch=20)#A-L, L-M and L-Y association (log OR)")
lines(seq(-3,3),ind[,2],type="b",lty=5,pch=4,cex=.7)
lines(seq(-3,3),ind[,3],type="b",lty=3,pch=2,cex=.7)

legend("bottomright",legend=c("Interventional (estimand, M-M)",
                              "Interventional (estimand, S-B)",
                              "Natural (estimate, ignoring L)"),
       lty=c(1,5,3),bty="n",cex=0.8,pch=c(20,4,2),pt.cex=c(1,.7,.7))

plot(seq(-3,3),dir[,1],type="b",lty=1,main="Direct effect (not via M)",ylim=yran,ylab="Risk difference",
     xlab="",pch=20)
lines(seq(-3,3),dir[,2],type="b",lty=5,pch=4,cex=.7)
lines(seq(-3,3),dir[,3],type="b",lty=3,pch=2,cex=.7)



plot(seq(-3,3),tot[,1],type="b",lty=1,main="Total effect (sum)",ylim=yran,ylab="Risk difference",
     xlab="",pch=20)
lines(seq(-3,3),tot[,2],type="b",lty=5,pch=4,cex=.7)
lines(seq(-3,3),tot[,3],type="b",lty=3,pch=2,cex=.7)

mtext("A-L, L-M and L-Y association (log OR)",side=1,outer=T,line=-1.5,cex=1.3)

dev.off()

### Create Table Supp1

res<-apply(pL,1,function(x)return(ifelse(round(range(x)[1],2)==round(range(x)[2],2),
                                    round(range(x)[1],2),paste(format(round(range(x),2),trim=T),collapse="-"))))

res<-cbind(res,apply(pM,1,function(x)return(ifelse(round(range(x)[1],2)==round(range(x)[2],2),
                                            round(range(x)[1],2),paste(format(round(range(x),2),trim=T),collapse="-")))))

res<-cbind(res,apply(pout,1,function(x)return(ifelse(round(range(x)[1],2)==round(range(x)[2],2),
                                            round(range(x)[1],2),paste(format(round(range(x),2),trim=T),collapse="-")))))
                      

write.csv(res,"TableSupp1.csv")

### Check values of remainder

range(rem)