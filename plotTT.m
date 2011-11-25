%Usage: e.g. plotTT('scalCls.dat'), or plotTT('scalCls.dat','scalCls2.dat')
function plotTT(varargin)
global l_label

hin=ishold;
h=false;
colours = 'krbcmgy';
l_label = 'l';

for i=1:nargin

dats{i}=load(varargin{i});
x=size(dats{i});
colnum(i)=x(2);
end

lmin=2;
lmax=1;
minTE=1e30;minEE=1e30;minTT=1e30;
maxTE=-1e30;maxEE=-1e30;maxTT=-1e30;
plotx=1;ploty=1;
if nargin>1
 ploty=2;
end

for i=1:nargin
x=dats{i};

col=colours(i);
ls=x(:,1);


TT=x(:,2);

if i==1    
 compTT=TT;
end

minTT = min(minTT,min(TT(:)));
maxTT = max(maxTT,max(TT(:)));

lmin=min(lmin,min(ls));
lmax=min(3000,max(lmax,max(ls)));
subplot(ploty,plotx,1);
seth(h);
plot(sqrt(ls),TT,':k','LineWidth',2);
setaxes(minTT*0.9,maxTT*1.1,1,sqrt(lmax));
if (lmax>2500)
ticks = [2 10 40 100 200 400 700 1000 1500 2000 2500];
else
ticks = [2 10 40 100 200 400 700 1000 1500 2000];
end
set(gca,'XTick',sqrt(ticks));
set(gca,'XTickLabel',ticks);

if i>1   
  print i
  subplot(ploty,plotx,2);
  if i==2
   seth(hin);
  else
   seth(h)
  end
  dm = min(size(TT),size(compTT));
  dm = min(3000,dm);
  cmp= (TT(1:dm) - compTT(1:dm)) ./ compTT(1:dm);

  plot(sqrt(ls(1:dm)),cmp,col);
  hold on;
  plot([1 sqrt(lmax)],[0 0],':k');
  ylabel('\Delta C_l / C_l');
  xlabel(l_label);
  set(gca,'XLim',[2 sqrt(lmax)]);
set(gca,'XTick',sqrt(ticks));
set(gca,'XTickLabel',ticks);

  hold on;
end

h=true;
end;

seth(hin);
% print -dpsc2 plotcls.ps;

function seth(h)

if h
 hold on
else
 hold off
end

function setaxes(a,b,c,d)
global l_label

set(gca,'YLim',[a,b],'XLim',[c,d]);
xlabel(l_label);ylabel('l(l+1) C_l / 2\pi');

