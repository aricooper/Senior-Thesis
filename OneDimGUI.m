function varargout = OneDimGUI(varargin)
% ONEDIMGUI MATLAB code for OneDimGUI.fig
%      ONEDIMGUI, by itself, creates a new ONEDIMGUI or raises the existing
%      singleton*.
%
%      H = ONEDIMGUI returns the handle to a new ONEDIMGUI or the handle to
%      the existing singleton*.
%
%      ONEDIMGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ONEDIMGUI.M with the given input arguments.
%
%      ONEDIMGUI('Property','Value',...) creates a new ONEDIMGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before OneDimGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to OneDimGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help OneDimGUI

% Last Modified by GUIDE v2.5 02-May-2017 17:18:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OneDimGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @OneDimGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
end
% --- Executes just before OneDimGUI is made visible.
function OneDimGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to OneDimGUI (see VARARGIN)
% Choose default command line output for OneDimGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes OneDimGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
%% Initialization
setGlobalx(1);
handles = guidata(hObject);
handles.t = sym('t');
handles.r = sym('r');
handles.theta = sym('theta');
% handles.theta = sym(pi/2);
handles.phi = sym('phi');
handles.propt = sym('T');
handles.GM = 1;
% handles.GM = 1;
[handles.M, handles.IM, handles.coords, handles.G, handles.cs] = Initialization(handles);
guidata(hObject, handles);
movegui('north');
[handles.fig1] = computeEmbeddingD(handles);
[handles.sliderHandles] = showControls(hObject, handles);
% handles.fig2 = figure();
guidata(hObject, handles);
% implay('trippy.avi',5);
% AcidComputer(hObject, handles);
end
% --- Outputs from this function are returned to the command line.
function varargout = OneDimGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

% --- Executes on button press in updateButton.
function updateButton_Callback(hObject, eventdata, handles)
% hObject    handle to updateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Generating path of Geodesics and Wavepacket

handles = guidata(hObject);
if ishandle(handles.fig1)
close(handles.fig1);
end
ax = handles.axes1;
children = get(ax, 'children');
delete(children(:));
% cla(ax,'reset');
 flat = get(handles.flatspaceButton,'Value');
 curved = get(handles.curvedspaceButton,'Value');
 orbits = get(handles.orbitsphereButton,'Value');
switch get(handles.popupmenu,'Value')   
    case 1
      geodesics = true;
    case 2
      geodesics = false;
    otherwise
        geodesics = true;
 end  


if (geodesics)
[GP, hObject] = computeOrbits(hObject); %% Geodesics
handles = guidata(hObject);
end

if (~geodesics && flat)
% Wave Packet Solutions
ComputeWPorbits(hObject); %% Wave Solutions

grid on;
end
%%
if (geodesics && curved)
% Embedding Diagram
theta = linspace(-2*pi,2*pi,100);
rho = linspace(15,1000,100);
[x,y] = meshgrid(theta,rho);
[X,Y] = pol2cart(x,y);

% the effect of space curvature

Z = sqrt(8*(sqrt(X.^2 + Y.^2) - 2));

% the effect of time curvature

T = sqrt(8*sqrt(X.^2+Y.^2));
colormap(jet);
surf(ax,X,Y,real(Z),real(T), 'EdgeColor', 'none','FaceLighting','flat');
hold on
[x,y,z] = sphere;
surf(ax,2*x,2*y,2*z,'FaceColor','k');
hold on

% animated line plot
[X,Y] = pol2cart(GP(:,2),GP(:,1));
T = sqrt(8*sqrt(GP(:,1)));
Z = sqrt(8*GP(:,1) - 2) + T;
maxR = max(GP(:,1));
handles.maxZ = max(Z);
guidata(hObject, handles);

zlim(ax,[-10, handles.maxZ+handles.maxZ*.5]);
xlim(ax,[-maxR, maxR]);
ylim(ax,[-maxR, maxR]);
xlabel('r in (GM)')
ylabel('r in (GM)')
zlabel('r in (GM)')
title('3D Planar Orbit')
view(100, 45);
grid on
rotate3d on

comet3(ax,X,Y,Z);
% h = animatedline(ax, 'Color','r', 'LineWidth',.001);
% for k = 1:length(GP(:,1))
%     addpoints(h,X(k),Y(k),Z(k));
%     drawnow limitrate
% end

end

%%
if (geodesics && flat)
ax = gca;
cla(ax,'reset');
maxR = max(GP(:,1))+10;
xlim(ax,[-maxR, maxR]);
ylim(ax,[-maxR, maxR]);
hold(ax)
% black circle
r = 2;
theta = linspace(0,2*pi);
x = r*cos(theta);
y = r*sin(theta);
fill(x,y,'k')
hold on
% 2D polar plot
[X,Y] = pol2cart(GP(:,2),GP(:,1));
comet(ax,X,Y);
view(180, 90); 
title('Flatspace Diagram')
end

%%
if (geodesics && orbits)
% 3D trajectories
[X,Y] = pol2cart(GP(:,2),GP(:,1));

T = sqrt(8*sqrt(GP(:,1)));
Z = sqrt(8*GP(:,1) - 2) + 2*T;
v = [1 1 1];
d = diag(v);
d(4,:) = [1 1 0];
d(5,:) = [1 0 1];
d(6,:) = [0 1 1];
rotate3d on
for i = 1:6
    for k = linspace(0,360,72)
   p = line(X,Y,Z);
   rotate(p,d(i,:),k*5,[0 0 0]);
   plot3(ax,p.XData,p.YData,p.ZData, 'LineWidth', .5);
   maxR = max(GP(:,1))+10;
xlim(ax,[-maxR, maxR]);
ylim(ax,[-maxR, maxR]);
zlim(ax,[-maxR, maxR]);
view(100, 45);
   title('3D Rotated Orbits')
   hold on
    end
end
end


guidata(hObject, handles);

end


% --- Executes on button press in Animatetoggle.
function Animatetoggle_Callback(hObject, eventdata, handles)
% hObject    handle to Animatetoggle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of Animatetoggle
%% Compute Path

flat = get(handles.flatspaceButton,'Value');
curved = get(handles.curvedspaceButton,'Value');
orbits = get(handles.orbitsphereButton,'Value');
v = [1 1 1];
d = diag(v);
d(4,:) = [1 1 0];
d(5,:) = [1 0 1];
d(6,:) = [0 1 1];


ax = handles.axes1;
maxR = handles.r*2;
xlim(ax,[-maxR, maxR]);
ylim(ax,[-maxR, maxR]);
if (curved)
zlim(ax,[-10, handles.maxZ+handles.maxZ*.5]);
end
if(get(hObject,'Value'))
children = get(ax, 'children');
delete(children(1));
delete(children(2));
delete(children(3));

end
l = animatedline(ax, 'Color','r', 'LineWidth',.001);
h = animatedline(ax, 'Color','r', 'LineWidth',.001);
output = hObject;
while(get(hObject,'Value'))

    [G, hObject] = ComputePath(output,handles, 2);
    handles = guidata(hObject);
    output = hObject;
    [X,Y] = pol2cart(G(2),G(1));
    T = sqrt(8*sqrt(G(1)));
    Z = sqrt(8*G(1) - 2) + T;
    
    if (curved)
    children = get(ax, 'children');
    delete(children(1));
    [x,y,z] = sphere;
    surf(ax,X+x*handles.GM,Y+y*handles.GM,Z+z*handles.GM,'FaceColor','k');
    addpoints(l,X,Y,Z);
    drawnow limitrate

    end
    
    if (flat)
    children = get(ax, 'children');
    delete(children(1));
    r = 5;
    theta = linspace(0,2*pi);
    x = r*cos(theta);
    y = r*sin(theta);
    fill(ax,x+X*handles.GM,y+Y*handles.GM,'k');
    addpoints(h,X,Y);
    drawnow limitrate

    end
%     % 3D trajectories
%     if (orbits)
%         cla(ax)
% for i = 1:6
%     for k = linspace(0,360,72)
%     l = animatedline(ax,'LineWidth',.001);
%     rotate(l,d(i,:),k,[0 0 0]);
%     addpoints(l,X,Y,Z);
%     hold on
%     end
% end
%     end
    xlabel('r in (GM)')
    ylabel('r in (GM)')
    zlabel('r in (GM)')
    title('Animation')

end
end



% --- Executes on selection change in popupmenu.
function popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu
end

% --- Executes during object creation, after setting all properties.
function popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function setGlobalx(val)
global x
x = val;
end

function r = getGlobalx
global x
r = x;
end



% --- Executes on button press in saveButton.
function saveButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

x = getGlobalx;
setGlobalx(x+1);
outputString = sprintf('save_%d.png', x);
saveas(gca, outputString);

end
