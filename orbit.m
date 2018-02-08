function [ Output ] = orbit(handles)
%ORBIT generates path outputting matrix of r and phi
%   Detailed explanation goes here
R = zeros(1, handles.SimSize);
Phi = zeros(1,handles.SimSize);

for n = 1:handles.SimSize
    
    R(n) = 2*handles.r - handles.rp + (handles.t^2)*( -(handles.GM/handles.r^2) + ...
    (handles.l^2/handles.r^3) - ((3*handles.GM*handles.l^2)/handles.r^4) );

    Phi(n) = handles.phi + handles.t*( handles.l/ ((1/2)*(R(n)+handles.r))^2 );
    
      handles.phi = Phi(n);
      handles.rp = handles.r;
      handles.r = R(n);
      handles.l = handles.l*diff(handles.phi);
end
% size(R)
% size(Phi)
Output(1,:) = R;
Output(2,:) = Phi;

end

