function sample = my_mvnrnd(muProposal, varProposal, Nparam1, Nparam2)

K = length(varProposal);
muM = repmat(muProposal, Nparam1, 1);
[U, Lambda] = eig(varProposal);
ld = diag(Lambda);
j = Lambda < max(ld) * 1e-12;
Lambda(j) = 0;
Sig = U * Lambda .^ .5;
sample = muM + (Sig * randn(K, Nparam2))';

end
