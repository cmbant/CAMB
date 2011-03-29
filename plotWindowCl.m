%e.g. plotWindowCl('test_scalCovCls.dat',1,1,100);
%or to maxlimum ell in the file
%e.g. plotWindowCl('test_scalCovCls.dat',1,1);
function plotWindowCl(fname, i,j,varargin)

clf;

if (size(varargin,1)>0)
 lmax=varargin{1};
 [ls,cl]=loadFileCl(fname,i,j,lmax);
else
 [ls,cl]=loadFileCl(fname,i,j);
 lmax=ls(size(ls,1));
end

plot(ls,cl);
xlim([2,lmax]);

xlabel('$l$','interpreter','latex');
ylabel('$l(l+1)C_l/2\pi$','interpreter','latex');