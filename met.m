function  [X,coords] = met(handles)
%FUNC Summary of this function goes here
%   Detailed explanation goes here
coords = [handles.t,handles.r,handles.theta,handles.phi];
rs = 2*handles.GM;
X = [-(1 - (rs/handles.r)); 1/(1 - (rs/handles.r)); handles.r^2; (handles.r^2)*sin(handles.theta)^2];


