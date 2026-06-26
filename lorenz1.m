function doall=lorenz1
clear all; close all;
% Eli Tziperman: example for course in nonlinear dynamics and
% chaos.

% Solve Lorenz system using the "modified Euler scheme" (a second
% order, constant step scheme), see Strogatz page 33.

% the equation:
%  dx/dt = sigma*(y-x)
%  dy/dt = r*x-y-x*z
%  dz/dt = x*y-b*z

% define parameters:
% ------------------
% parameters:
sigma=10.0;
b=8.0/3.0;
r=28;

delta_t=0.001;

%print (and save) every nprint steps:
nprint=1;
% integration time:
T=200.0;

% number of steps:
N=floor(T/delta_t);
NP=floor(N/nprint);

% initial conditions:
n=1;
x=zeros(3,NP);
n=1;
x(:,n)=[0.1,0,0]';

% integrate:
for n=1:1:NP-1
  xp=x(:,n);
  for j=1:1:nprint
    x_tilde=xp+delta_t*rhs(xp,sigma,b,r);
    xp=xp+0.5*(rhs(xp,sigma,b,r)  ...
	       + rhs(x_tilde,sigma,b,r))*delta_t;
  end
  x(:,n+1)=xp;
end

% plot simple time series:
beg=floor(NP*9/10);
end1=NP;
time=delta_t*[1:1:NP];
plot(time(beg:end1),x(1,beg:end1))
xlabel('time')
ylabel('x')

% plot phase space:
figure
plot(x(1,beg:end1),x(3,beg:end1))
xlabel('x')
ylabel('z')

% plot phase space:
figure
b=floor(NP/2);
e=NP;
plot(x(1,beg:end1),x(2,beg:end1))
xlabel('x')
ylabel('y')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ff=rhs(x,sigma,b,r)
% return the rhs of the dynamical system:

ff=[sigma*(x(2)-x(1)),r*x(1)-x(2)-x(1)*x(3),x(1)*x(2)-b*x(3)]';
