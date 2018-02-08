function [ WP, x ] = ComputeWPorbits( hObject )
%COMPUTEWPORBITS Summary of this function goes here
%   Detailed explanation goes here

%% 1 Dim wavepacket propagation (in t)
handles = guidata(hObject);

% Initial Conditions

% r = get(handles.sliderHandles(1), 'value');
r= 40; %% only r that works well right now
guidata(hObject, handles);

%% Central differencing method
% Initial Conditions
N = 500; % Num of points
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
WP = u;
u(1,:) = ic(x,r);
% ut = u(1,:)+dt*g(x);

% for i = 2:N-1
% u(2,i+1) = .5*(a(i)*u(1,i+1)+2*(1-a(i))*u(1,i)+a(i)*u(1,i-1))+dt*g(i,r);
% end

u(2,1:end-1) = u(1,1:end-1);


ax = gca;
cla(ax, 'reset');
tn = 2;
while (a(tn) < 100)
% M = ((4./a).*I + K)^(-1)*(2*((4./a).*I - K).*u - ((4./a).*I + K).*u0 - 4*b); %implicit method
% u0 = u;
% u = M;
tn = tn +1;
for i = 2:N-1
u(tn,i) = a(i)*u(tn-1,i+1)+2*(1-a(i))*u(tn-1,i)+a(i)*u(tn-1,i-1)-u(tn-2,i);
end

WP = conj(u).*u;
guidata(hObject, handles);
plot(ax,x,WP(tn,:),'b',x,WP(1,:),'k');  % Plot Wave Packet Solutions
drawnow limitrate
title(sprintf('Wavepacket Propagation : t =%0.2f',tn*dt))
xlim(ax,[0, r]);
ylim(ax,[0, 1]);
end

plot(ax,x,WP(end, :), 'b');

function y = ic(x,r)
% initial condition function
y = exp(-.5*(x-r).^2);
% y = max(exp(-(x/8e-6).^2),r);
% y = max(exp(1-16*(x+r).^2).*sin(w.*x),r);
end

function y = f(x)
% metric compenent in a
 k = (1-(2./x));
y = 1-(.0001+(2*(x(1)-x(2))))*k;
end



%% failed attempt
% %Implicit Differncing Method
% % W0 = 8e-6;
% W0 = 1;
% u0 = 0;
% x = linspace(handles.r, 0, s);
% u = exp (-(x/W0).^2); %initial Guassian wp
% I = eye(s);
% K = full(gallery('tridiag',s,-1,2,-1));
% b = zeros(s,1);
% b(1) = -1; 
% b(end) = -1;
% r=.5;
% % r = -(1-(2/handles.r));
% 
% 
% while(handles.simSize < s)
%     
% handles.simSize = handles.simSize + 1;    
% nu = ((4/r^2).*I + K)^(-1)*(2*((4/r^2).*I - K).*u - ((4/r^2).*I + K).*u0 - 4.*b);
% WP = conj(nu).*nu;
% 
% end

% Plot Wave Packet Solutions
% figure
% plot(x,WP);



end