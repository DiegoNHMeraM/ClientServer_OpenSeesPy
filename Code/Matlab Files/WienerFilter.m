function w_stat=WienerFilter(x,par)
% Computes regularized Wiener Filter coefficients w_stat.
% x: data to estimate Wiener
% par: filter parameters [L,N,l]
% L: filter length, N : number of points, l: regularization factor
% Example: w_stat = WienerFilter(x,[L,l,N])
if nargin==1
    L=8; l=1e-10; N=2000; % default parameters
else
    L=par(1); l=par(2); N=par(3);
end
Rxd=zeros(L,1);Rx=zeros(L); gd=0; 
for i=L:N
    Rxd=Rxd+x(i:-1:i-L+1)*x(i+1);
    Rx=Rx+x(i:-1:i-L+1)*x(i:-1:i-L+1)';
    gd=gd+x(i+1)^2;
end
Rxd=Rxd/(N-L+1);
Rx=Rx/(N-L+1);
w_stat=(Rx+l*eye(L))\Rxd;
gd=gd/(N-L+1);
cost=gd-Rxd'*((Rx+l*eye(L))\Rxd);
c1=cond(Rx);
fprintf('Regularized Wiener Filter with L= %i, lambda= %8.3e, cost %8.3e, cond %8.3e\n',L,l,cost,c1);