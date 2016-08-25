function ESS=TempFun(gamma_inc,lnw,lnew,ESS_bound)
%returns ESS(gamma)-(ESS_bound-1)
incrementalw=bsxfun(@times,exp(gamma_inc),sum(lnew,2));

j=find(isnan(incrementalw) | imag(incrementalw)~=0);

incrementalw(j)=-inf;

lnw=bsxfun(@plus,lnw,incrementalw); %%  

W=exp(bsxfun(@minus,lnw,max(lnw)));

ESS=sum(W).^2./sum(W.^2);


