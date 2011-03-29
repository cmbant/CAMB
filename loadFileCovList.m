%load list of covariance matrices for each ell
function [n,ls,covs]=loadFileCovList(varargin)

d=load(varargin{1});

if (nargin>1)
lmax=varargin{2};
else
lmax=d(size(d,1),1);
end;
imax=lmax-1;
ls=d(:,1);
n=sqrt(size(d,2)-1);

covs=cell(imax,1);
for i=1:imax
 L=ls(i);
 cov=zeros(n,n);
 for j=1:n
  cov(:,j)=d(i,1+(j-1)*n+1:1+j*n);
 end;

 covs{L}=cov;
 
end;

ls=ls(1:imax);
