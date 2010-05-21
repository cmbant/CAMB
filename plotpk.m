
function plotpk(varargin);

colours = 'krbcmgy';
mink = 1e-4;
clf;
maxk=0;

for i=1:nargin

x=load(varargin{i});
col=colours(i);
kh=x(:,1);
sz = size(kh);
maxk=max(kh(sz(1)),maxk);
loglog(kh,x(:,2),['-' col]);
hold on;

end

xlim([mink maxk]);
xlabel('k/[h^{-1} Mpc]');
ylabel('P_k/(h Mpc^{-1})^3');
hold off;

%print -depsc2 pk.eps
end
