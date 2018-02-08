function [sliderHandle] =  showControls(hObject, handles)
%SHOWCONTROLS Summary of this function goes here
%   Detailed explanation goes here

%% 3) Create multiple sliderPanels as children of a UIPANEL.

figure;
movegui('northwest');
h = uipanel(gcf,'title','MAIN','units','normalized','pos',[0.2 0.1 0.6 0.8]);
PnlOpt.title = 'Controls';
PnlOpt.bordertype = 'none';
PnlOpt.titleposition = 'centertop';
PnlOpt.fontweight = 'bold';
EditOpts = {'fontsize',10};
LabelOpts = {'fontsize',9,'fontweight','b'};
numFormat = '%0.0f';
titleStrings = {'Initial r', 'Angular Velocity ', 'Simulation Duration','Mass'};
sliderTags = {'rSlider', 'uSlider ', 'simSlider','mSlider'};

startPos = {
	        [0.1 0.28 0.8 0.21];
	        [0.1 0.51 0.8 0.21];
	        [0.1 0.74 0.8 0.21]
            [0.1 0.05 0.8 0.21]};
        
sldrCallbacks = {'1';'2';'3';'4'};
        
setMin = [5.0; 0; 0.0; 1];
setMax = [500; 100; 100; 10];
setInitial = [100; 0; 0; 1];

for ii = 1:4
    
	PnlOpt.position = startPos{ii};
	PnlOpt.title = titleStrings{ii};
	SldrOpt.callback = sldrCallbacks{ii};
    SldrOpt.min = setMin(ii);
    SldrOpt.max = setMax(ii);
    SldrOpt.value = setInitial(ii);
    SldrOpt.string = sliderTags(ii);
	sliderHandle(ii) = sliderPanel(hObject, handles, h,PnlOpt,SldrOpt,EditOpts,LabelOpts,numFormat);

end

end

