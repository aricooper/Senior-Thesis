function [ U,V,W ] = generatePoints( handles )
%GENERATEPOINTS Creates arrays U,V,W for streamline function
%   Detailed explanation goes here

% Generate points (X = U, Y = V, Z = W)

for r = 2:5
    for theta = 1:4
    X(r,theta) = subs(handles.G(2), [handles.r, handles.theta], [r,sym(2*pi/theta)]);   %r
    Y(r,theta) = subs(handles.G(3), [handles.r, handles.theta], [r,sym(2*pi/theta)]);   %theta
    Z(r,theta) = r^2*sin(theta)^2;      %phi (constant)
    end
end

% Convert to cartesian (x,y,z)

[U,V,W] = sph2cart(X,Y,Z);


end

