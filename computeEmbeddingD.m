function [fig] =  computeEmbeddingD( handles )
%COMPUTEEMBEDDINGD Generates the 3D Embedding Diagram and Infalling
%particle Geodesics
%   3D Embedding diagram of equitorial planar spacial curvature of space
%   time with the color representing the curvature of the time curvature
%   over that space. Geodesics particles falling at rest only fall radially
%   inward due to the time component of their four-velocity displaying the
%   effect of the time curvature of spacetime.

 %% Schwarzschild Embedding Diagram
c = 50;    % number of points
% The effect of spatial curvature
S = zeros(4);
theta = sym('theta');
r = sym('r');
theta = (pi/2);


S = subs(handles.M, [handles.theta, handles.r], [theta, r]);
 % constant time
 % constant theta @ pi/2

% Generate z(r) using r component of metric
syms x z(r) Dz;
Dz = simplify(solve((x)^2 + 1 == S(2,2)));
z(r) = simplify(int(Dz(1),r));    % integrate to get z(r)

theta = linspace(-2*pi,2*pi,c);
rho = linspace(2,100,c);
[x,y] = meshgrid(theta,rho);
[X,Y] = pol2cart(x,y);

% Z = simplify(subs(z(r), {r}, {sqrt(X.^2 + Y.^2)}));
Z = sqrt(8*(sqrt(X.^2 + Y.^2) - 2));

% the effect of time curvature

syms y z(t) dz
dz = simplify(solve((y)^2 - 1 == S(1,1)));
z(t) = simplify(int(dz(1),r));   % integrate to get z(r) for t component
% T = simplify(subs(z(t), {r}, {sqrt(X.^2 + Y.^2)}));
T = sqrt(8*sqrt(X.^2+Y.^2));


az = 100;
el = 45;
% G = subs(handles.G, {handles.r}, {sqrt(X.^2+Y.^2)});
P = (-1 ./ (r^2));
lG = 100;
for rho = 2:100
G(2,rho) = lG - (subs(P, {r}, {rho}));
lG = G(2,rho);
G(1,rho) = 0;
end

colorbar;
rotate3d on;



colormap(jet);
axes(handles.axes1);   
surf(X,Y,real(Z),real(T));
title('Flamms Paraboloid')
xlabel('r in (GM)')
ylabel('r in (GM)')
zlabel('r in (GM)')
view(az, el);

fig = figure;
colormap(jet);
contourf(X,Y,real(T));
title('Infalling Particle at Rest Geodesic')
hold on
[x,y] = pol2cart(G(1,:),G(2,:));
plot(x,y, 'LineWidth', 2);

% figure;
% contour(X,Y,real(T));


% axes(handles.axes2);    surf(X,Y,real(T),real(Z));
% title('Embedding Diagram of t')
% xlabel('r in (GM)')
% ylabel('r in (GM)')
% zlabel('r in (GM)')
% view(az, el);



end

