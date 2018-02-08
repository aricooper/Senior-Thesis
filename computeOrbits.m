function [GP, hObject] = computeOrbits( hObject )
%COMPUTEORBITS Generates schwarzschild planar orbits around mass
%   Allows for input of initial conditions and updated calculations for
%   different trajectories and animation toggle.


%% Initial Conditions
handles = guidata(hObject);
handles.GM = .9 + (get(handles.sliderHandles(4), 'value')/10);
handles.r = get(handles.sliderHandles(1), 'value');
u = get(handles.sliderHandles(2), 'value'); % initial angular speed
s = get(handles.sliderHandles(3), 'value')*1000;
handles.simSize = 0;

guidata(hObject, handles);

% handles.l = handles.r^2 * (handles.u);
handles.l = handles.r/sqrt(handles.r - 3)-(u/50);
Tc = 2*pi*(handles.r^2/handles.l); % proper time orbit period
handles.T = Tc/500;
handles.u = handles.l / handles.r^2;
handles.size = 1000+u^2+s;


handles.nuphi = .5*handles.u*handles.T;
handles.r1 = handles.r;
handles.nur = handles.r1-5;
guidata(hObject, handles);
h = hObject;
%% Compute Path
while(handles.simSize < handles.size)
    
handles.simSize = handles.simSize + 1;    
[G, hObject] = ComputePath(h,handles, 1);
h = hObject;
GP(handles.simSize,:) = [G(1), G(2)];

end
guidata(hObject, handles);




end

