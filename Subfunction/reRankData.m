function [ rankData ] = reRankData( Data, order)
%ExtractData : only output partial Data from beg_i to end_i
% Input structure Data
% row numbers are the same for each field
Fname = fieldnames(Data);
N = length(Fname);
rankData = Data;
for i=1:N
    tmp=getfield(rankData, Fname{i});
    rankData = setfield(rankData, Fname{i},tmp(order,:));
end
end

