#  fluswall - code for example 6.16 using all parameters 

# read data
y=scan("/mydata/flu.dat")
num=length(y);


# initial estimates 
nstate=4;
alpha1=1.4;  alpha2=-.5;  # eqn 6.155
beta0=.3; beta1=.1;         # eqn 6.156
sQ1=.1; sQ2=.1;           # stdev for w1, w2
sR=.1;                      # stdev for obs eqn
initpar=c(alpha1,alpha2,beta0,beta1,sQ1,sQ2,sR)


# -- Function to Calculate Likelihood --
Linn=function(para){
  alpha1=para[1]; alpha2=para[2]
  beta0=para[3]; beta1=para[4]      
  sQ1=para[5];  sQ2=para[6]; 
  sR=para[7] 
  y=as.matrix(y)
  M1=as.matrix(cbind(1,0,0,1))        # obs matrix normal
  M2=as.matrix(cbind(1,0,1,1))        # obs matrix flu epi
  xf=matrix(0, nstate, 1)  # x filter
  xp=matrix(0, nstate, 1)  # x pred
  Pf=diag(.1, nstate)      # filter cov
  Pp=diag(.1, nstate)      # prec cov
  pi11=.75   # initial switch probs
  pi12=.25
  pi22=.75
  pi21=.25
  pif1=.5; pif2=.5   # pi_j(t|t) - filter
  #  
  phi=matrix(0,nstate,nstate)
  phi[1,1]=alpha1; phi[1,2]=alpha2; phi[3,3]=beta1 
  phi[2,1]=1; phi[4,4]=1 
  Gam=as.matrix(rbind(0,0,beta0,0))
  Q=matrix(0,nstate,nstate)
  Q[1,1]=sQ1^2; Q[3,3]=sQ2^2; 
  R=sR^2
  #-----------begin filter------------
  like=0
  for(i in 1:num){
    xp = phi%*%xf + Gam 
    Pp = phi%*%Pf%*%t(phi) + Q
    sig1=as.real(M1%*%Pp%*%t(M1) + R)          # innov var
    sig2=as.real(M2%*%Pp%*%t(M2) + R)
    k1=Pp%*%t(M1)/sig1
    k2=Pp%*%t(M2)/sig2 
    e1=y[i]-M1%*%xp                   # innov
    e2=y[i]-M2%*%xp
    #
    pip1=pif1*pi11 + pif2*pi21;       # pi_j(t|t-1) - predict
    pip2=pif1*pi12 + pif2*pi22;
    den1= (1/sqrt(sig1))*exp(-.5*e1^2/sig1);
    den2= (1/sqrt(sig2))*exp(-.5*e2^2/sig2);
    denom=pip1*den1+pip2*den2;
    pif1=pip1*den1/denom;
    pif2=pip2*den2/denom;
    #
    pif1=as.real(pif1); pif2=as.real(pif2)
    e1=as.real(e1); e2=as.real(e2)
    xf = xp + pif1*k1*e1 + pif2*k2*e2
    eye=diag(1,nstate)
    Pf = pif1*(eye-k1%*%M1)%*%Pp + pif2*(eye-k2%*%M2)%*%Pp 
    like=like-log(pip1*den1 + pip2*den2)
  }
  return(like)
}
#---- end Like function

# -- Estimation --
est=optim(initpar,Linn,NULL,method="BFGS",hessian=TRUE,control=list(trace=1,REPORT=1))
stderr=sqrt(diag(solve(est$hessian)))
est      
estim=cbind(est$par,stderr)
rownames(estim)=c("alpha1","alpha2","beta0","beta1","sQ1","sQ2","sR")
estim


# -- Final run and plots --
prob=matrix(0,num,1)                # pi_1(t|t)
y=as.matrix(y)
yp=y                                # y(t|t-1)
xfil=matrix(0,num,nstate)           # x(t|t) for plotting
alpha1=est$par[1]; alpha2=est$par[2]
beta0=est$par[3]; beta1=est$par[4]      
sQ1=est$par[5];  sQ2=est$par[6]
sR=est$par[7]
M1=as.matrix(cbind(1,0,0,1))        # obs matrix normal
M2=as.matrix(cbind(1,0,1,1))        # obs matrix flu epi
xf=matrix(0, nstate, 1)  # x filter
xp=matrix(0, nstate, 1)  # x pred
Pf=diag(.1, nstate)  # filter cov
Pp=diag(.1, nstate)  # prec cov
pi11=.75   # initial switch probs
pi12=.25
pi22=.75
pi21=.25
pif1=.5; pif2=.5   # pi_j(t|t) - filter
#  
phi=matrix(0,nstate,nstate)
phi[1,1]=alpha1; phi[1,2]=alpha2;  
phi[2,1]=1; phi[3,3]=beta1; phi[4,4]=1 
Gam=as.matrix(rbind(0,0,beta0,0))
Q=matrix(0,nstate,nstate)
Q[1,1]=sQ1^2; Q[3,3]=sQ2^2; 
R=sR^2
#-----------begin filter------------
for(i in 1:num){
  xp = phi%*%xf + Gam 
  Pp = phi%*%Pf%*%t(phi) + Q
  sig1=as.real(M1%*%Pp%*%t(M1) + R)          # innov var
  sig2=as.real(M2%*%Pp%*%t(M2) + R)
  k1=Pp%*%t(M1)/sig1
  k2=Pp%*%t(M2)/sig2 
  e1=y[i]-M1%*%xp                   # innov
  e2=y[i]-M2%*%xp
  #
  pip1=pif1*pi11 + pif2*pi21;       # pi_j(t|t-1) - predict
  pip2=pif1*pi12 + pif2*pi22;
  den1= (1/sqrt(sig1))*exp(-.5*e1^2/sig1);
  den2= (1/sqrt(sig2))*exp(-.5*e2^2/sig2);
  denom=pip1*den1+pip2*den2;
  pif1=pip1*den1/denom; 
  pif2=pip2*den2/denom;
  pif1=as.real(pif1); pif2=as.real(pif2)
  e1=as.real(e1); e2=as.real(e2)
  xf = xp + pif1*k1*e1 + pif2*k2*e2
  eye=diag(1,nstate)
  Pf = pif1*(eye-k1%*%M1)%*%Pp + pif2*(eye-k2%*%M2)%*%Pp
  #
  xfil[i,]=xf
  prob[i]=pif1
  if (pip1 > pip2) yp[i]=M1%*%xp
  else yp[i]=M2%*%xp  
}
#
time=seq(1969,1978.99,by=1/12)
par(mfrow=c(3,1))
plot(time, y[13:132], type="o", ylim=c(0,1))
lines(time, prob[13:132], col=4, lty="dashed")
#  
plot(time, xfil[13:132,1], type="l", ylim=c(-.5,.8))
lines(time, xfil[13:132,3], type="l", col=2)
lines(time, xfil[13:132,4], type="l", col=4)
#
plot(time, y[13:132], type="p", ylim=c(0,1))
lines(time, yp[13:132], type="l")