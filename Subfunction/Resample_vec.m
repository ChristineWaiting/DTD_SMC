function [Y,bin] = Resample_vec(X, weights, N)

[PH, bin] = histc(rand(N,1), cumsum([0; weights]));
Y = reRankData( X, bin);
% Y.m = X.m(bin,:);
% Y.beta = X.beta(bin,:);
% Y.sigma = X.sigma(bin,:);