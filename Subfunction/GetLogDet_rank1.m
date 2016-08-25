function lnDet=GetLogDet_rank1(A,U)
%compute log determinant of Sig_i=diag(A(i,:))+U(i,:)'*U(i,:)
%consequence of matrix determinant lemma
[Nparticles,N]=size(A);

if all(size(U)~=size(A))
    error('A and U need to have the same size!');
end

lnDet=sum(log(A),2)+log(1+sum(U.^2./A,2));


