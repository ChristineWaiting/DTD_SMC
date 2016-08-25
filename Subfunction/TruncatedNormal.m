function X=TruncatedNormal(mu,sigma,LB,UB,Nsim)

u=rand(Nsim,1);

X=norminv(normcdf(LB,mu,sigma)+u*(normcdf(UB,mu,sigma)-normcdf(LB,mu,sigma)),mu,sigma);




