function [GP, hObject] = ComputePath( hObject,handles, c )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%% Compute Path
if (c == 1)
   handles = guidata(hObject);
end

   handles.nur = 2*handles.r1 - handles.r + (handles.T^2 * ...
       (-(handles.GM/handles.r1^2) + (handles.l^2/handles.r1^3) - (3*handles.GM*handles.l^2/handles.r1^4)));
   guidata(hObject, handles);

   handles.nuphi = handles.nuphi + handles.T * (handles.l/(.5*(handles.nur+handles.r1)^2));
   guidata(hObject, handles);
   
   % create array for r and phi
   GP = [handles.nur, handles.nuphi];

   % reset new variables
   handles.phi = handles.nuphi;
   handles.r = handles.r1;
   handles.r1 = handles.nur;
   guidata(hObject, handles);




end

