function varargout = ROI_selection_gui_PCAICA(varargin)
% ROI_selection_gui_PCAICA MATLAB code for ROI_selection_gui_PCAICA.fig
%      ROI_selection_gui_PCAICA, by itself, creates a new ROI_selection_gui_PCAICA or raises the existing
%      singleton*.
%
%      H = ROI_selection_gui_PCAICA returns the handle to a new ROI_selection_gui_PCAICA or the handle to
%      the existing singleton*.
%
%      ROI_selection_gui_PCAICA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ROI_selection_gui_PCAICA.M with the given input arguments.
%
%      ROI_selection_gui_PCAICA('Property','Value',...) creates a new ROI_selection_gui_PCAICA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ROI_selection_gui_PCAICA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ROI_selection_gui_PCAICA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ROI_selection_gui_PCAICA

% Last Modified by GUIDE v2.5 18-May-2018 15:19:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ROI_selection_gui_PCAICA_OpeningFcn, ...
                   'gui_OutputFcn',  @ROI_selection_gui_PCAICA_OutputFcn, ...
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

% --- Executes just before ROI_selection_gui_PCAICA is made visible.
function ROI_selection_gui_PCAICA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ROI_selection_gui_PCAICA (see VARARGIN)

% Choose default command line output for ROI_selection_gui_PCAICA
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Set keyboard shortcuts
set(handles.figure1, 'UserData', handles);
set(handles.figure1, 'KeyPressFcn', @keyboardShortcuts);


% set(gcf, 'WindowButtonDownFcn', @selectRoi);
% set(src, 'Pointer', 'crosshair'); % Optional

global premasks
global viewed_ROIs
if sum(premasks.ROI_list) > 0
    viewed_ROIs = double(premasks.ROI_list);
    curr_nROI = sum(premasks.ROI_list>0 & premasks.ROI_list~=5);
    set(handles.textCurrnROI, 'String', num2str(curr_nROI));
else
    viewed_ROIs = zeros(length(premasks.ica),1);
end
% determine if jump to next cell automatically
premasks.JumpNext = 1;
premasks.image_to_show = premasks.movm;
premasks.group_to_use = 1;
% set default erosion threshold for masks
premasks.erode_step = 5; 
% check to see if a dff_images exists in that directory
try tmp_path = premasks.fName;
    ind_slash = find(tmp_path == '\');
    path_dff_images = [tmp_path(1:ind_slash(end)) 'dFF_images'];
    if exist(path_dff_images)
        tmpDirR = dir(path_dff_images)
        for i = 1:length(tmpDirR)
            if strfind(tmpDirR(i).name,'dff_mean')
                set(handles.pushbuttonmaxOri, 'Enable', 'on');
                RR = PPPack.hf.readtiff([path_dff_images '\' tmpDirR(i).name]);
                RR = nanmax(RR,[],3);
                % HAVE TO REMOVE 2* what normally would remove from buffer
                % bc already spatially downsampled
                if size(RR,1) < 500 % if using old dFF images
                    cols_remove = 32;
                    rows_remove = 10;
                    RR(:,[1:cols_remove end-cols_remove+1:end],:) = [];
                    RR([1:rows_remove end-rows_remove+1:end],:,:) = [];
                else
                    RR = PPPack.hf.bin_mov_xyt(RR,2,1,0);
                    cols_remove = 32;
                    rows_remove = 10;
                    RR(:,[1:cols_remove end-cols_remove+1:end],:) = [];
                    RR([1:rows_remove end-rows_remove+1:end],:,:) = [];
                end
                premasks.max_mov_Ori = RR;
            end
        end
    else
        set(handles.pushbuttonmaxOri, 'Enable', 'off');
    end
catch err
    set(handles.pushbuttonmaxOri, 'Enable', 'off');
end    

% This sets up the initial plot - do when we are invisible
% so window can get raised using ROI_selection_gui_PCAICA.
updateICAplots(handles, 1, 0);
set(handles.pushbutton_commit, 'Enable', 'off');
set(handles.pushbutton_remove_roi, 'Enable', 'off');

% --- Outputs from this function are returned to the command line.
function varargout = ROI_selection_gui_PCAICA_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton_next.
function pushbutton_next_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global premasks
% first update last ROI viewed selection
last_ROI = get(handles.popupmenu_icanumber, 'Value');
set(handles.textLastROIn, 'String', num2str(last_ROI));

popup_sel_index = get(handles.popupmenu_icanumber, 'Value')+1;
if popup_sel_index <= length(premasks.ica)
	set(handles.popupmenu_icanumber, 'Value', popup_sel_index);
	popupmenu_icanumber_Callback(hObject, eventdata, handles);
end

% --- Executes on button press in pushbutton_back.
function pushbutton_back_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
popup_sel_index = get(handles.popupmenu_icanumber, 'Value')-1;
if popup_sel_index >= 1
	set(handles.popupmenu_icanumber, 'Value', popup_sel_index);
	popupmenu_icanumber_Callback(hObject, eventdata, handles);
end
updateICAplots(handles, popup_sel_index, 0);

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu_icanumber.
function popupmenu_icanumber_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_icanumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu_icanumber contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_icanumber
global premasks
currentICA = get(handles.popupmenu_icanumber, 'Value');
updateICAplots(handles, currentICA);

% Make pushbutton_remove_roi available iff a number has been assigned
% already
% 
% if  ismember(currentICA, find(arrayfun(@(s) ~isempty(s.roiNum), premasks.ica)))
% 	set(handles.pushbutton_remove_roi, 'Enable', 'on');
% 	set(handles.pushbutton_keep_roi, 'Enable', 'off');
% else
% 	set(handles.pushbutton_remove_roi, 'Enable', 'off');
% 	set(handles.pushbutton_keep_roi, 'Enable', 'on');
% end


% --- Executes during object creation, after setting all properties.
function popupmenu_icanumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_icanumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

global premasks
set(hObject, 'String', [{}, cellstr(num2str([1:length(premasks.ica)]'))]);

% --- Executes on button press in pushbutton_keep_roi.
function pushbutton_keep_roi_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_keep_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global premasks
% global JumpNext
currentICAnumber = get(handles.popupmenu_icanumber, 'Value');
% old way
lowest = min(setdiff(1:length(premasks.ica), [premasks.ica.roiNum]));
premasks.ica(currentICAnumber).roiNum = lowest;
% new way
premasks.ROI_list(currentICAnumber) = premasks.group_to_use;

% update masks variable
premasks.masks = premasks.masks + ...
    currentICAnumber*single(premasks.ica(currentICAnumber).filter>0);
% masks with group number
premasks.masks_with_group = premasks.masks_with_group + ...
    premasks.group_to_use*single(premasks.ica(currentICAnumber).filter>0);
% can't have over 5 groups
premasks.masks_with_group(premasks.masks_with_group > 5) = 5;

% updateAvailableRoiNumbers(handles);
curr_nROI = sum(premasks.ROI_list>0 & premasks.ROI_list~=5);
set(handles.textCurrnROI, 'String', num2str(curr_nROI));
updateICAplots(handles, currentICAnumber, currentICAnumber);
set(handles.pushbutton_remove_roi, 'Enable', 'on');


%Automatically advance to next segment:
if premasks.JumpNext == 1;
    pushbutton_next_Callback(hObject, eventdata, handles)
end


% --- Executes on button press in pushbutton_load_corr.
function pushbutton_load_corr_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load_corr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global premasks

premasks.image_to_show = premasks.movcorr;
curr_ROI = get(handles.popupmenu_icanumber, 'Value');
updateICAplots(handles, curr_ROI, curr_ROI);

% --- Executes on button press in Load_mean_image.
function Load_mean_image_Callback(hObject, eventdata, handles)
% hObject    handle to Load_mean_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global premasks
% if isfield(premasks,'backup_movm')
%     premasks.movm = premasks.backup_movm;
%     updateICAplots(handles, 1, 0);
%     premasks.use_other_image_tag = 0;
% end
premasks.image_to_show = premasks.movm;
curr_ROI = get(handles.popupmenu_icanumber, 'Value');
updateICAplots(handles, curr_ROI, curr_ROI);

% --- Executes on button press in pushbutton_remove_roi.
function pushbutton_remove_roi_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_remove_roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global premasks
currentICAnumber = get(handles.popupmenu_icanumber, 'Value');
% first remove from list
premasks.ROI_list(currentICAnumber) = 0;

% now remake the masks and masks with group variable
% masks with group number
premasks.masks_with_group  = zeros(size(premasks.masks_with_group));
premasks.masks  = zeros(size(premasks.masks_with_group));
for curr_ROI = 1:length(premasks.ROI_list)
    if premasks.ROI_list(curr_ROI) ~= 0
        % for masks with group variable
        premasks.masks_with_group = premasks.masks_with_group + ...
            premasks.ROI_list(curr_ROI)*single(premasks.ica(curr_ROI).filter>0);
        % for matrix with ROI number
        premasks.masks = premasks.masks + ...
            curr_ROI*single(premasks.ica(curr_ROI).filter>0);
    end
end
premasks.masks_with_group(premasks.masks_with_group > 5) = 5;



% % now remove from both mask images
% ind_to_remove = find(premasks.ica(currentICAnumber).filter > 0);
% premasks.masks(ind_to_remove) = 0;
% premasks.masks_with_group(ind_to_remove) = 0;

% set(handles.pushbutton_remove_roi, 'Enable', 'off');
set(handles.pushbutton_keep_roi,'Enable','on');

curr_nROI = sum(premasks.ROI_list>0);
set(handles.textCurrnROI, 'String', num2str(curr_nROI));

updateICAplots(handles, currentICAnumber, 0);
popupmenu_icanumber_Callback(hObject, eventdata, handles);



% --- Executes on button press in checkbox_confirmcommit.
function checkbox_confirmcommit_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_confirmcommit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_confirmcommit

if get(hObject,'Value')
	set(handles.pushbutton_commit, 'Enable', 'on');
else
	set(handles.pushbutton_commit, 'Enable', 'off');
end
	
% --- Executes on button press in pushbutton_commit.
% Parses data into a format usable by the movie object.
function pushbutton_commit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_commit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global premasks

premasks = PPPack.hf.build_masks_for_extraction(premasks);

% Return to the function that called the GUI
close(handles.figure1);
return


% --- Executes on selection change in popupmenu_mousefunction.
function popupmenu_mousefunction_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_mousefunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_mousefunction contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_mousefunction

% The handle of meanimg is not available globally, so...
meanimg_handle = findobj(handles.axes_meanimg, 'Type', 'image');
contents = cellstr(get(hObject,'String'));	
selection = contents(get(hObject,'Value'));
switch selection{1}
    case 'Choose option'
        
	case 'Zoom' % Zoom
		axes(handles.axes_meanimg);
		zoom on
		set(meanimg_handle, 'ButtonDownFcn', []);
		
	case 'Select ROI' % Select Roi
		axes(handles.axes_meanimg);
		zoom off
% 		set(meanimg_handle, 'ButtonDownFcn', @selectRoi);
        set(gcf, 'WindowButtonDownFcn', @selectRoi);
		set(gca, 'UserData', handles);
    case 'ROI in area' % Get ROI in region
		axes(handles.axes_meanimg);
		zoom off
        set(gcf, 'WindowButtonDownFcn', '');
		set(meanimg_handle, 'ButtonDownFcn', @chooseROIRegion);
		set(gca, 'UserData', handles);
    case 'Erode mask'
		axes(handles.axes_meanimg);
		zoom off
        set(gcf, 'WindowButtonDownFcn', @selectRoi);
		set(gcf, 'WindowScrollWheelFcn', @ScrollWheelThresh);
% 		set(gca, 'UserData', handles);   
    case 'Draw ROI'
		axes(handles.axes_meanimg);
		zoom off
        set(gcf, 'WindowButtonDownFcn', '');
        addinROIs(handles)
	otherwise
		error('Dont recognize menu item -- weird.')
end
		

% --- Executes during object creation, after setting all properties.
function popupmenu_mousefunction_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_mousefunction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject, 'String', [{'Choose option'},{'Select ROI'},{'Zoom'},{'ROI in area'},{'Erode mask'},{'Draw ROI'}]);

function addinROIs(handles)
persistent h
global premasks
global viewed_ROIs

meanimg_handle = findobj(handles.axes_meanimg, 'Type', 'image');
mm = imfreehand(handles.axes_meanimg);
newmask = createMask(mm,meanimg_handle(2));
nROIs_ica = length(premasks.ica)
curr_fields = fieldnames(premasks.ica(end));
% now preallocate new mask
for i = 1:length(curr_fields)
    premasks.ica(nROIs_ica+1).(curr_fields{i}) = [];
end
% now stick stuff in to fake this mask
premasks.ica(nROIs_ica+1).filter = newmask;
premasks.ica(nROIs_ica+1).trace = nan(size(premasks.ica(nROIs_ica).trace));
tmpR = regionprops(newmask);
premasks.ica(nROIs_ica+1).centroid = sqrt(tmpR.Centroid(1).^2+tmpR.Centroid(2).^2);
[tmpidx(:,1),tmpidx(:,2)] = find(newmask);
premasks.ica(nROIs_ica+1).idx = tmpidx;

% now lets go to that ROI
set(handles.popupmenu_icanumber, 'String', [{}, cellstr(num2str([1:length(premasks.ica)]'))]);
% popupmenu_icanumber_CreateFcn([],[],handles)
set(handles.popupmenu_icanumber, 'Value', nROIs_ica+1);
popupmenu_icanumber_Callback([], [], handles);
% protect agains trying to draw ROI weird place
set(handles.popupmenu_mousefunction,'Value',1)

function updateICAplots(handles, icanumber, ics_to_draw)
persistent h
global premasks
global viewed_ROIs

if nargin < 3
	full_redraw = 0;
	ics_to_draw = find(arrayfun(@(s) ~isempty(s.roiNum), premasks.ica));
	ics_to_draw = max(ics_to_draw);
end

if ics_to_draw == 0
	full_redraw = 1;
else
	full_redraw = 0;
end

% try clearing java
if mod(icanumber,10) == 0
    java.lang.System.gc()
    display('Clearing java')
end

% Assign colors depending on whether an ICA has been assigned an ROI
% number:
contourcolor = [0 250 154]./255;
textcolor = [255 0 0]./255;
numstring = num2str(icanumber);

axes(handles.axes_meanimg);
hold on;

% draw image
cla;
colormap('Gray');
% new version
img = premasks.image_to_show;
img(img<0) = 0;
img(isnan(img)) = 0;
img = sqrt(img);
img = img/max(img(:));
img = adapthisteq(img);
tmpRRR = median(double(nonzeros(img)));
h.movm_image = imagesc(img,[0 tmpRRR*2.5]);
xlim([0 size(premasks.movm, 2)])
ylim([0 size(premasks.movm, 1)])
set(gca, 'XTick', [], ...
        'YTick', [], ...
        'YDir','reverse',...
        'DataAspect', [1 1 1]);


% Delete existing plot elements
try delete(h.current_contour), catch, end
try delete(h.current_num), catch, end

% mask group colors
mask_group_colors = [0 255 255;... % light blue
                        255 140 0;... % orange
                        50 205 50;... % lime green
                        128 0 0;... % dark salmon
                        0 0 128;... % add many lines of navy, so that if there are too many overlapping regions the GUI doesnt crash
                        0 0 128;... % again, navy
                        0 0 128;... % navy
                        0 0 128;... % naby
                        0 0 128]./255;  % navy
                    
effective_group_number = unique(premasks.masks_with_group(:));
effective_group_number(effective_group_number == 0) = [];

% plot masks that have been accepted
if ~isempty(effective_group_number)
    for curr_group = effective_group_number'
    %     mask_img = premasks.masks;
    %     mask_img = mask_img > 0;
        mask_img = premasks.masks_with_group == curr_group;
    %     mask_color = cat(3, 0*ones(size(mask_img)), 255*ones(size(mask_img)), 255*ones(size(mask_img)))./255;
        mask_color = cat(3, mask_group_colors(curr_group,1)*ones(size(mask_img)), ...
            mask_group_colors(curr_group,2)*ones(size(mask_img)), ...
            mask_group_colors(curr_group,3)*ones(size(mask_img)));
        hold on
        m = imagesc(mask_color);
        hold off
        transparency_factor = 0.3;
        set(m,'AlphaData',mask_img.*transparency_factor)
    end
end

% plot current mask
curr_mask_img = single(premasks.ica(icanumber).filter>0);
curr_mask_color = cat(3, 255*ones(size(curr_mask_img)), 0*ones(size(curr_mask_img)), 255*ones(size(curr_mask_img)))./255;
hold on
m = imagesc(curr_mask_color);
hold off
transparency_factor = 0.15;
set(m,'AlphaData',curr_mask_img.*transparency_factor)

% write text
ind_to_use_text = [premasks.ica(icanumber).idx(end,2),premasks.ica(icanumber).idx(end,1)]; 
try h.current_num = text(ind_to_use_text(1)+3,ind_to_use_text(2)+3, numstring,...
        'Color', textcolor,...
        'FontSize', 12,...
        'FontWeight', 'bold',...
        'HorizontalAlignment', 'center');
catch err
    h.current_num = text(ind_to_use_text(1),ind_to_use_text(2), numstring,...
        'Color', textcolor,...
        'FontSize', 12,...
        'FontWeight', 'bold',...
        'HorizontalAlignment', 'center');
end

% update viewed ROIs variable
viewed_ROIs(icanumber) = 1;
updateViewedROIs(handles);

hold off

% Plot for bottom 2 axes
axes(handles.axes_aux2);
imagesc(img,[0 tmpRRR*2.5]);
hold on
m = imagesc(curr_mask_color);
hold off
transparency_factor2 = 0.25;
set(m,'AlphaData',curr_mask_img.*transparency_factor2)
curr_centroid = [nanmean(premasks.ica(icanumber).idx(:,1)) nanmean(premasks.ica(icanumber).idx(:,2))];
zoom_factor = 35;
axis([curr_centroid(2)-zoom_factor curr_centroid(2)+zoom_factor curr_centroid(1)-zoom_factor curr_centroid(1)+zoom_factor])
set(gca, 'XTick', [], 'YTick', []);
% now for just filter
axes(handles.axes_aux1);
rr = colormap('gray');
rr(1,:) = ([211 211 211])./255;
tmp_mask = premasks.ica(icanumber).filter;
tmp_mask(tmp_mask<0) = 0;
imagesc(tmp_mask);colormap(rr)
axis([curr_centroid(2)-zoom_factor curr_centroid(2)+zoom_factor curr_centroid(1)-zoom_factor curr_centroid(1)+zoom_factor])
set(gca, 'XTick', [], 'YTick', []);
colorbar;


axes(handles.axes_trace);
plot((1:numel(premasks.ica(1).trace))./6.2,premasks.ica(icanumber).trace,'k','Linewidth',1);
set(gca, 'YTick', []);
xlabel('Time (s, assuming 31 fps acquisition and 5 times temporal downsampling for ICA)');

% Make sure some things are up to date:
if full_redraw
	popupmenu_mousefunction_Callback(handles.popupmenu_mousefunction, [], handles)
end

function updateViewedROIs(handles)
global premasks
global viewed_ROIs
% availableRoiNums = setdiff(1:length(premasks.ica), [premasks.ica.roiNum]);

%Update ICA number list to also show ROI numbers
ICAlist = (1:length(premasks.ica))';
ICAlabels = cell(size(ICAlist));
for i = 1:length(ICAlabels)
	if viewed_ROIs(i) == 0
		ICAlabels(i) = {num2str(i)};
    else        
		ICAlabels(i) = {[num2str(i) ' --']};
	end
end
set(handles.popupmenu_icanumber, 'String', ICAlabels);

% Callback for click on ROI in mean img:
function selectRoi(hObject, eventdata, ~)
global premasks
% handles = get(gca, 'Userdata');
handles = guidata(hObject);

coords = get(gca, 'Currentpoint');
coords = round(coords([1 3]));
selected_roi = 0;

% if outside range then ignore
xLimits = get(handles.axes_meanimg, 'xlim');
yLimits = get(handles.axes_meanimg, 'ylim');

if coords(1) > xLimits(2) | coords(2) > yLimits(2)
    return
end


% Find ROI which matches the clicked point:
for i = 1:numel(premasks.ica)
	if premasks.ica(i).filter(coords(2), coords(1))
		selected_roi = i;
		break
	end
end

if selected_roi
    % first update last ROI viewed selection
    last_ROI = get(handles.popupmenu_icanumber, 'Value');
    set(handles.textLastROIn, 'String', num2str(last_ROI));
    
	set(handles.popupmenu_icanumber, 'Value', selected_roi);
	popupmenu_icanumber_Callback(hObject, eventdata, handles);
end

function chooseROIRegion(hObject, eventdata, ~)
global premasks

% handles = guidata(hObject);

try handles = get(gca, 'Userdata');
    coords_box = imrect(handles.axes_meanimg);
    position = round(getPosition(coords_box));
    delete(coords_box);


    % iterate through each ROI and find the ROIs in this box
    ROI_in_selected_region = [];
    size_each_ROI = [];
    for curr_ROI = 1:length(premasks.ica)
        curr_mask = premasks.ica(curr_ROI).filter > 0;
        curr_square = curr_mask(position(2):position(2)+position(4),position(1):position(1)+position(3));
        if sum(curr_square(:)) > 1
            ROI_in_selected_region = [ROI_in_selected_region curr_ROI];
            size_each_ROI = [size_each_ROI sum(curr_square(:))];
        end
    end
    % now sort by size
    [~,ind_sort] = sort(size_each_ROI,'descend');
    ROI_in_selected_region = ROI_in_selected_region(ind_sort);

    % set select menu
    if ~isempty(ROI_in_selected_region)
        listboxROIArea_CreateFcn([],eventdata,handles,num2cell(ROI_in_selected_region))
    else
        listboxROIArea_CreateFcn([],eventdata,handles,{'No ROI'})
    end
catch err
end
    
function ScrollWheelThresh(hObject, eventdata, handles)
global premasks
% handles = get(gca, 'Userdata');
handles = guidata(hObject);
curr_ROI = get(handles.popupmenu_icanumber, 'Value');

% pre-allocate the threshold for that filter and save backup filter
if ~isfield(premasks.ica(curr_ROI),'backup_filter')
    curr_mask = premasks.ica(curr_ROI).filter;
    premasks.ica(curr_ROI).backup_filter = curr_mask;
    premasks.ica(curr_ROI).moving_threshold = 0;
elseif isempty(premasks.ica(curr_ROI).backup_filter)
    curr_mask = premasks.ica(curr_ROI).filter;
    premasks.ica(curr_ROI).backup_filter = curr_mask;
    premasks.ica(curr_ROI).moving_threshold = 0;
else
    curr_mask = premasks.ica(curr_ROI).backup_filter;
end
%normalize maskl
curr_mask_norm = curr_mask./max(curr_mask(:));
baseline_threshold = 0;
nonzero_values = (curr_mask_norm(curr_mask_norm > 0));
[nonzero_values_sort sort_ind] = sort(nonzero_values,'ascend');

% get the amount that scrolled
step_size = premasks.erode_step;
if eventdata.VerticalScrollCount > 0
    increment_threshold = -1*step_size;
else
    increment_threshold = 1*step_size;
end
curr_threshold = premasks.ica(curr_ROI).moving_threshold + increment_threshold;
if curr_threshold < 0
    curr_threshold = 0;
end
% reset threshold
premasks.ica(curr_ROI).moving_threshold  = curr_threshold;


value_of_nonzero_values_to_toss = nonzero_values_sort(1:min([curr_threshold length(nonzero_values_sort)]));
curr_mask(ismember(curr_mask_norm,value_of_nonzero_values_to_toss)) = 0;


    % additional step where ask if not connected to highest point then chuck
    % image
    % this is the largest point
try [idx,idy] = find(curr_mask == max(curr_mask(:)));
    idx = idx(1);idy = idy(1);
    % now see what are connected things
    bwc = bwconncomp(curr_mask>0,8);
    % make image of non connected things
    labmat = labelmatrix(bwc);
    % this is the value in the labelmatrix of the highest weight
    curr_labmat_highestpx = double(labmat(idx,idy));
    new_idx_to_toss = find(~ismember(labmat,[0 curr_labmat_highestpx]));
    curr_mask(new_idx_to_toss) = 0;
catch err
    warning('Struggle with highest point erosion')
end

premasks.ica(curr_ROI).filter = curr_mask;
updateICAplots(handles, curr_ROI, curr_ROI)

function keyboardShortcuts(hObject, eventdata)
global premasks
global viewed_ROIs
handles = get(hObject, 'UserData');
switch eventdata.Key
	case 'k' % Keep
		pushbutton_keep_roi_Callback(hObject, eventdata, handles);
		
	case 'n' % Next
		pushbutton_next_Callback(hObject, eventdata, handles);
		
	case 'b' % Back
		pushbutton_back_Callback(hObject, eventdata, handles)
        
    case 't' % Trash cell
        % set popupmenu
        curr_group = get(handles.popupmenuROIGroup, 'Value');
        set(handles.popupmenuROIGroup, 'Value', 5);
        premasks.group_to_use = 5;
        pushbutton_keep_roi_Callback(hObject, eventdata, handles);
        set(handles.popupmenuROIGroup, 'Value', curr_group);
        premasks.group_to_use = curr_group;
    case 'h'
        % find next ROI that hasn't been viewed
        % this is curr ROI
        curr_ROI = get(handles.popupmenu_icanumber, 'Value');
        % set last ROI string handle
        set(handles.textLastROIn, 'String', num2str(curr_ROI));
        % next ROI that is larger than curr ROI
        unviewed_ROIs = find(viewed_ROIs == 0);
        next_ROI = unviewed_ROIs(find(unviewed_ROIs > curr_ROI,1,'first'));
        % set next ROI
        try set(handles.popupmenu_icanumber, 'Value', next_ROI);
            popupmenu_icanumber_Callback(hObject, eventdata, handles);
        catch err 
        end
        
        
	otherwise
		fprintf('No function is assigned to %s.\n', eventdata.Key);
end

% --- Executes on button press in pushbutton_auto.
function pushbutton_auto_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_auto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global premasks

if ~isfield(premasks, 'tuning')
	return
end

% Go through all ROIs
for i = 1:numel(premasks.ica)
	set(handles.popupmenu_icanumber, 'Value', i);
	popupmenu_icanumber_Callback(hObject, eventdata, handles);
	
	if premasks.ica(i).p_val < 0.1
		pushbutton_keep_roi_Callback(hObject, eventdata, handles);
	end
	pause(0.5)
end



function Enter_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to Enter_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Enter_ROI as text
%        str2double(get(hObject,'String')) returns contents of Enter_ROI as a double

global premasks
desired_ROI = str2double(get(hObject,'String')); % returns contents of Enter_ROI as a double
updateICAplots(handles, desired_ROI);
set(handles.popupmenu_icanumber, 'Value', desired_ROI);
popupmenu_icanumber_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function Enter_ROI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Enter_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listboxROIArea.
function listboxROIArea_Callback(hObject, eventdata, handles)
% hObject    handle to listboxROIArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxROIArea contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxROIArea

% get click
contents = cellstr(get(hObject,'String'));
ROI_selection = contents{get(hObject,'Value')};

set(handles.popupmenu_icanumber, 'Value', str2num(ROI_selection));
popupmenu_icanumber_Callback(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function listboxROIArea_CreateFcn(hObject, eventdata, handles, ROI_to_input)
% hObject    handle to listboxROIArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
try if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
    end

    if nargin == 4
        set(handles.listboxROIArea,'string',{ROI_to_input{:}});
    end
catch err
end


% --- Executes on selection change in popupmenuROIGroup.
function popupmenuROIGroup_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenuROIGroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenuROIGroup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenuROIGroup
global premasks
% get selection
contents = cellstr(get(hObject,'String'));
selection = contents(get(hObject,'Value'));
switch selection{1}
	case 'Group 1'
        premasks.group_to_use = 1;
	case 'Group 2'
        premasks.group_to_use = 2;
    case 'Group 3'
        premasks.group_to_use = 3;
    case 'Group 4'
        premasks.group_to_use = 4;    
    case 'Non ROI'
        premasks.group_to_use = 5;
	otherwise
		error('Dont recognize menu item -- weird.')
end

% --- Executes during object creation, after setting all properties.
function popupmenuROIGroup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenuROIGroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% possible group options
set(hObject, 'String', [{'Group 1'},{'Group 2'},{'Group 3'},{'Group 4'},{'Non ROI'}]);


% --- Executes on button press in pushbuttonmaxOri.
function pushbuttonmaxOri_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonmaxOri (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global premasks

premasks.image_to_show = premasks.max_mov_Ori;
curr_ROI = get(handles.popupmenu_icanumber, 'Value');
updateICAplots(handles, curr_ROI, curr_ROI);




% --- Executes on button press in checkboxNext.
function checkboxNext_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxNext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxNext

global premasks

if get(hObject,'Value')
    premasks.JumpNext = 1;
else
    premasks.JumpNext = 0;
end


% --- Executes on slider movement.
function sliderThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to sliderThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% range that I will allow is from 1-60
global premasks
curr_erode_step_value = round(get(hObject,'Value'))
premasks.erode_step =  (curr_erode_step_value);

% --- Executes during object creation, after setting all properties.
function sliderThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
