function [ GP ] = computeTrips(c, hObject, handles )
%COMPUTETRIPS Summary of this function goes here
%   Detailed explanation goes here
size = c;
handles.l = handles.r^2 * handles.u;
handles.nuphi = .5*handles.u*handles.T;
handles.r1 = handles.r;
handles.nur = handles.r;
guidata(hObject, handles);

%% Compute Path

while(handles.simSize < size)
   handles.simSize = handles.simSize + 1;
   handles.nur = 2*handles.r1 -  handles.r + handles.T^2 * ...
       (-(handles.GM/handles.r1^2) + (handles.l^2/handles.r1^3) - (3*handles.GM*handles.l^2/handles.r1^4));
   guidata(hObject, handles);

   handles.nuphi = handles.nuphi + handles.T * (handles.l/(.5*(handles.nur+handles.r1)^2));
   guidata(hObject, handles);
   
   % create array for r and phi
   GP(handles.simSize,:) = [handles.nur; handles.nuphi];
   
   % reset new variables
   handles.phi = handles.nuphi;
   handles.r = handles.r1;
   handles.r1 = handles.nur;
   guidata(hObject, handles);

end

handles.GP = GP;
guidata(hObject, handles);




end

