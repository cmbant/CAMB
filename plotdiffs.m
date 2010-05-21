%use plotdiffs('scalCls.dat','scalCls2.dat','scalCls3.dat') to plot fractional differences

function plotdiffs(varargin)
global l_label

clf;
ticks = [2 10 40 100 200 400 700 1000 1500];

colours = 'krbcmgy';
l_label = 'l';

for i=1:nargin

dats{i}=load(varargin{i});
x=size(dats{i});
colnum(i)=x(2);
end

lmin=2;
lmax=1;

plotx=1;ploty=4;

for i=1:nargin
x=dats{i};

col=colours(i);
ls=x(:,1);

TT=x(:,2);
EE=x(:,3);
noB = colnum(i)==4 || colnum(i)==6;
if (noB)
 ploty=3;
 TE=x(:,4);
 BB=0;
else
 BB=x(:,4);
 TE=x(:,5);
end 

if i==1    
 compTT=TT;
 compEE=EE;
 compTE=TE;
 if ~noB
  compBB=BB;
 end
end

lmin=min(lmin,min(ls));
if i==1
 lmax = max(ls);
end

if i>1   

  lmax=min(lmax,max(ls));
  subplot(ploty,plotx,1);

  dm = min(size(TT),size(compTT));
  cmp= (TT(1:dm) - compTT(1:dm)) ./ compTT(1:dm);

  plot(sqrt(ls(1:dm)),cmp,col);
  hold on;
  ylabel('{\Delta C^{TT} / C^{TT}}');
  xlabel(l_label);

  if i==2
   if (lmax >= 2000) ticks = [ticks 2000];end
   if (lmax >= 3000) ticks = [ticks 3000];end
  end

  set(gca,'XTick',sqrt(ticks));
  set(gca,'XTickLabel',ticks);

  set(gca,'XLim',[sqrt(lmin) sqrt(lmax)]);

  subplot(ploty,plotx,2);

  cmp= (TE(1:dm) - compTE(1:dm)) ./ sqrt(compEE(1:dm).*compTT(1:dm));

  plot(sqrt(ls(1:dm)),cmp,col);
  hold on;

  ylabel('\Delta C^{TE}/(C^{EE} C^{TT})^{1/2}');
  xlabel(l_label);

  set(gca,'XTick',sqrt(ticks));
  set(gca,'XTickLabel',ticks);

  set(gca,'XLim',[sqrt(lmin) sqrt(lmax)]);

  subplot(ploty,plotx,3);

  cmp= (EE(1:dm) - compEE(1:dm)) ./ compEE(1:dm);
  maxcmp = max(cmp);

  plot(sqrt(ls(1:dm)),cmp,col);
  hold on;
  ylabel('{\Delta C^{EE} / C^{EE}}');
  xlabel(l_label);

  set(gca,'XTick',sqrt(ticks));
  set(gca,'XTickLabel',ticks);
  set(gca,'XLim',[sqrt(lmin) sqrt(lmax)]);

  if ~noB
  subplot(ploty,plotx,4);

  cmp= (BB(1:dm) - compBB(1:dm)) ./ compBB(1:dm);
  maxcmp = max(cmp);

  plot(sqrt(ls(1:dm)),cmp,col);
  hold on;
  ylabel('{\Delta C^{BB} / C^{BB}}');
  xlabel(l_label);

  set(gca,'XTick',sqrt(ticks));
  set(gca,'XTickLabel',ticks);
  set(gca,'XLim',[sqrt(lmin) sqrt(lmax)]);
  end

end


end;

%print -dpsc2 plotdiffcls.ps;

