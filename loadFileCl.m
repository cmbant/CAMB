%Load an cross or auto spectrum from a CAMB sources scalCovCls file
%e.g. loadFileCl('test_scalCovCls.dat',1,2,1000) reads the cross-spectrum of the 
%first window function with the second, up to lmax=1000
function [ls,cl]=loadFileCl(varargin)

if (nargin>3)
[~,ls,covs]=loadFileCovList(varargin{1},varargin{4});
else
[~,ls,covs]=loadFileCovList(varargin{1});    
end

imax=size(ls,1);
cl=zeros(imax,1);
for i=1:imax
    L=ls(i);
    cov=covs{L};
    cl(i)=cov(3+varargin{2},3+varargin{3});
end;
