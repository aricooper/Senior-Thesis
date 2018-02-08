function [GP] =  changeConditions( hObject, handles, n )
%CHANGECONDITIONS Summary of this function goes here
%   Detailed explanation goes here
handles.T = 1/(n*n);
handles.r = 10*n/2;
handles.u = 1*n;
size = 100 + n*10;
handles.simSize = 0;
[GP] = computeTrips(size, hObject, handles);
end

