
function plotcls(varargin)
global l_label


clf;
hin=ishold;
h=false;
colours = 'krbcmgy';
l_label = 'l';


ticks = [2 40 200 600 1200];


for i=1:nargin

dats{i}=load(varargin{i});
x=size(dats{i});
colnum(i)=x(2);
end

lmin=2;
lmax=1;
minTE=1e30;minEE=1e30;minTT=1e30;
maxTE=-1e30;maxEE=-1e30;maxTT=-1e30;
plotx=2;ploty=2;
if nargin>1
 ploty=3;
end

for i=1:nargin
x=dats{i};

col=colours(i);
ls=x(:,1);

TT=x(:,2);
EE=x(:,3);
noB = colnum(i)==4 || colnum(i)==6;
if (noB)
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
end

minTT = min(minTT,min(TT(:)));
maxTT = max(maxTT,max(TT(:)));
maxTE = max(maxTE,max(TE(:)));
maxEE = max(maxEE,max(EE(:)));
minTE = min(minTE,min(TE(:)));
minEE = min(minEE,min(EE(:)));

lmin=min(lmin,min(ls));
lmax=max(lmax,max(ls));

if i==1
 if (lmax >= 2000) ticks = [ticks 2000];end
 if (lmax >= 3000) ticks = [ticks 3000];end
end

subplot(ploty,plotx,1);
seth(h);
loglog(ls,TT,col);
setaxes(minTT*0.9,maxTT*1.1,lmin,lmax);

subplot(ploty,plotx,3);
seth(h);
plot(sqrt(ls),TT,col);
setaxes(0,maxTT*1.1,sqrt(lmin),sqrt(lmax));
set(gca,'XTick',sqrt(ticks));
set(gca,'XTickLabel',ticks);


subplot(ploty,plotx,2);

seth(h);
e=minEE*0.5;
t=max(maxEE,maxTE)*2;
loglog(ls,EE,['-' col]);
setaxes(e,t, lmin,lmax);

hold on;
loglog(ls,TE,['-' col]);
loglog(ls,-TE,[':' col]);

if ~noB
 loglog(ls,BB,['--' col]);
end

subplot(ploty,plotx,4);
seth(h);
plot(sqrt(ls),TE,['-' col]);
hold on;
plot(sqrt(ls),EE,['-' col]);
setaxes(minTE*1.1,maxTE*1.1,sqrt(lmin),sqrt(lmax));

set(gca,'XTick',sqrt(ticks));
set(gca,'XTickLabel',ticks);


if i>1   
  subplot(ploty,plotx,5);
  if i==2
   seth(hin);
  else
   seth(h)
  end

  dm = min(size(TT),size(compTT));
  cmp= (TT(1:dm) - compTT(1:dm)) ./ compTT(1:dm);
  plot(sqrt(ls(1:dm)),cmp,col);
  set(gca,'XTick',sqrt(ticks));
  set(gca,'XTickLabel',ticks);
  set(gca,'XLim',[sqrt(lmin) sqrt(lmax)]);

  ylabel('{\Delta C_l / C_l} (TT)');
  xlabel(l_label);
  hold on;
  subplot(ploty,plotx,6);
  if i==2
   seth(hin);
  else
   seth(h)
  end
  cmp= (EE(1:dm) - compEE(1:dm)) ./ compEE(1:dm);
  maxcmp = max(cmp);
  plot(sqrt(ls(1:dm)),cmp,col);
  hold on;  
  for j=1:dm
    cmp= (TE(1:dm) - compTE(1:dm)) ./ sqrt(compEE(1:dm).*compTT(1:dm));
  end
  plot(sqrt(ls(1:dm)),cmp,[':' col]);
  set(gca,'XTick',sqrt(ticks));
  set(gca,'XTickLabel',ticks);
  set(gca,'XLim',[sqrt(lmin) sqrt(lmax)]);
   
  ylabel('{{\Delta C_l / C_l}}(EE/TE)');
  xlabel(l_label);
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