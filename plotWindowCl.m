%e.g. plotWindowCl('test_scalCovCls.dat',1,1);
%or to maxlimum ell in the file in red
%e.g. plotWindowCl('test_scalCovCls.dat',1,1,100,'-r');
%or plot multiple files
%plotWindowCl({'test_scalCovCls.dat','test2_scalCovCls.dat'},1,1);

function plotWindowCl(fname, i,j,varargin)

lmax=0;
if (iscell(fname))
    clf;
    ncl=size(fname,2);
else
    ncl=1;
end

linespecs={'-b','--r','-.b','-g','--c','-.m'};

for ix=1:ncl
     
if (iscell(fname))
    %multiple files
filename=fname{ix};
else
filename=fname;
end
   
if (size(varargin,1)>0)
 lmax=varargin{1};
 [ls,cl]=loadFileCl(filename,i,j,lmax);
else
 [ls,cl]=loadFileCl(filename,i,j);
 lmax=max(lmax,ls(size(ls,1)));
end
if (ncl==1)
 plot(ls,cl,varargin{2:end});
else
    h1=subplot(2,1,1);
    plot(ls,cl,linespecs{ix});
    hold on;
    if(ix==1)
     comparecl=cl;
     else
         h2=subplot(2,1,2);
         plot(ls,(cl-comparecl)./comparecl,linespecs{ix});
         hold on;
    end
end

end; % loop over files

if (ncl>1)  
%    jointfig(gcf,2,1);
    xlim([2,lmax]);  
    xlabel('$l$','interpreter','latex');
    ylabel('$\Delta C_l /C_l$');
    set(gcf,'CurrentAxes',h1)
    legend(fname,'interpreter','none','location','southeast','fontsize',10);
    hold off;
end
xlim([2,lmax]);
xlabel('$l$','interpreter','latex');
ylabel('$l(l+1)C_l/2\pi$','interpreter','latex');