% Ari Cooper
% 1D GUI SIM
% 12/28/2016

clc, close all, clear all;

%% Run GUI
% OneDimGUI

%% WPorbit Test Code

%%Implicit Differncing Method
% I = eye(N);
% K = full(gallery('tridiag',N,-1,2,-1))*dt;
% b = zeros(N,1);
% b(1) = -1; 
% b(end) = -1;
% b = b*dt;

% Initial Conditions
N = 500; % Num of points
r = 100;
x = linspace(r,0,N); % Coordinate vector
dx(1:N) = x(1)- x(2);
dx(1:N) = dx.*(1-2./x);
nm = 1:N;
dt = 1e-3;
% dt(1:N) = 1e-3*(1-2./x); % time step
tf = 100;  % final time
v = (dx/dt).*nm;
a(1:N) = .1 + f(x);
% a(1:N) = 1;
u = zeros(100,N);
u(1,:) = ic(x,r);
% ut = u(1,:)+dt*g(x);

% for i = 2:N-1
% u(2,i+1) = .5*(a(i)*u(1,i+1)+2*(1-a(i))*u(1,i)+a(i)*u(1,i-1))+dt*g(i,r);
% end

u(2,1:end-1) = u(1,1:end-1);

ax = gca;
tn = 2;
while (a(tn)< inf)
% M = ((4./a).*I + K)^(-1)*(2*((4./a).*I - K).*u - ((4./a).*I + K).*u0 - 4*b); %implicit method
% u0 = u;
% u = M;
tn = tn +1;
for i = 2:N-1
u(tn,i) = a(i)*u(tn-1,i+1)+2*(1-a(i))*u(tn-1,i)+a(i)*u(tn-1,i-1)-u(tn-2,i);
end

WP = conj(u).*u;

plot(x,WP(tn,:),'b',x,WP(1,:),'k');  % Plot Wave Packet Solutions
title(sprintf('t =%0.2f',tn*dt))
drawnow limitrate
% pause(.1)
end

function y = ic(x,r)
% initial condition function
y = exp(-.25*(x-r).^2);
% y = max(exp(-(x/8e-6).^2),r);
% y = max(exp(1-16*(x+r).^2).*sin(w.*x),r);
end


function y = f(x)
% metric component
 k = (1-(2./x));
y = 1-(.0001+(2*(x(1)-x(2))))*k;


end
