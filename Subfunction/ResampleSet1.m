function [theta,lnw,l,lnpdf_int]=ResampleSet1(theta,lnw,l,lnpdf_int)
%this function resamples the fixed parameter set

Nparam=size(l,1);

cumw = [0; cumsum(exp(lnw - max(lnw)))];
cumw = cumw/cumw(end);

%use stratified resampling
[PH, bin] = histc(((0:Nparam-1)+rand(1,Nparam))/Nparam,cumw);

theta = reRankData( theta, bin);

l = l(bin,:);

lnpdf_int = lnpdf_int(bin,:);

lnw = zeros(Nparam, 1);
