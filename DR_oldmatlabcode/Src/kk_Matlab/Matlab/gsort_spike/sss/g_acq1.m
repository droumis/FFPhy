function varargout = g_acq1(varargin)
% G_ACQ1 M-file for g_acq1.fig
%      G_ACQ1, by itself, creates a new G_ACQ1 or raises the existing
%      singleton*.
%
%      H = G_ACQ1 returns the handle to a new G_ACQ1 or the handle to
%      the existing singleton*.
%
%      G_ACQ1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in G_ACQ1.M with the given input arguments.
%
%      G_ACQ1('Property','Value',...) creates a new G_ACQ1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before g_acq1_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to g_acq1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help g_acq1

% Last Modified by GUIDE v2.5 30-Jan-2006 20:04:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @g_acq1_OpeningFcn, ...
                   'gui_OutputFcn',  @g_acq1_OutputFcn, ...
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


% --- Executes just before g_acq1 is made visible.
function g_acq1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to g_acq1 (see VARARGIN)

% Choose default command line output for g_acq1
handles.output = hObject;

warning off;

%%%%%%%%%%%%%%%%%%FIGURE PROPERTIES%%%%%%%%%%%%%%%%
handles.guifighandle=get(0,'CurrentFigure');
set(handles.guifighandle, 'HandleVisibility', 'Off');
handles.oriunits = get(handles.guifighandle,'Units');
handles.oriposition = get(handles.guifighandle,'Position');
set(0,'Units','Characters');
handles.screensize = get (0, 'Screensize');
set(0,'Units','pixels');
set(handles.guifighandle,'Position',[1.4 5 handles.screensize(3)/5 handles.screensize(4)-10 ]);
handles.position = get(handles.guifighandle,'Position');
% a= get (0, 'Screensize');
% set(gcf, 'Position', [a(3)*0.3 a(4)*(0.51) a(3)*0.68 a(4)*0.42])
handles.parameterfile = [];
handles.clr = ['y' 'm' 'b' 'c' 'r' 'y' 'm' 'b' 'c' 'r' 'y' 'm' 'b' 'c' 'r'];

%%%%%%%%%%%%%%Initialize AI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
handles.N_chan=0:3;
handles.gain=1000; 
handles.maxsignal=500;
handles.srate=30000;
handles.time=5;
handles.Savefile='test';
handles.bufferblocks=14;
handles.buffersize=131072;

%%%%%%%%%%%%%%VARIABLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

handles.load_nblocks=1; handles.load_nsecs=[0,5]; handles.load_ch=handles.N_chan;
handles.loaddata=[]; handles.loadtime=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Update handles structure
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% UIWAIT makes g_acq1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = g_acq1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  CALLBACKS  %%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in make_ai.
function make_ai_Callback(hObject, eventdata, handles)
% hObject    handle to make_ai (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gain=handles.gain; maxsignal=handles.maxsignal;  %in uV
time=handles.time;		% FOR TIMER FUNCTION, WHICH WILL STOP AI
ai=analoginput('nidaq',1);
setverify(ai,'InputType','NonReferencedSingleEnded');
addchannel(ai,handles.N_chan);

setverify(ai,'LoggingMode','Disk');
setverify(ai,'LogToDiskMode','Index');
setverify(ai,'LogFileName',[handles.Savefile]);

%%%% UNITS AND GAIN
%scaled value = (A/D value)(units range)/(sensor range)
setverify(ai.Channel,'Units','uV');	%get(ai.Channel,'UnitsRange')
%setverify(ai.Channel,'SensorRange',[-5 5]*1000000/gain)		%in uV

maxsignal_V = maxsignal/(10^6); ip = maxsignal_V*gain;
setverify(ai.Channel,'InputRange',[-ip ip]);   % in V: Any signal above this will be clipped
setverify(ai.Channel,'UnitsRange',[-maxsignal maxsignal]);    % Set units range (/sensor range) to be same, in uV, as of input range (eg [-1000 1000]uV units range for [-1000 1000]uV input range)
setverify(ai.Channel,'SensorRange',[-ip ip]);     %% corresponding to units range in V
%%%  (In volts: setverify(ai.Channel,'SensorRange',[-5 5]/gain)   )
setverify(ai,'SampleRate',handles.srate);
%setverify(ai,'SamplesPerTrigger',32000)   % 1 sec of data
setverify(ai,'SamplesPerTrigger',Inf); 	% Continuous function till timer stop is

%setverify(ai,'TimerPeriod',time);  ai.TimerFcn = {'stop'};	%setverify(ai,'TimerFcn','timerstop.m')
actualsamprate=get(ai,'SampleRate');
setverify(ai,'SamplesAcquiredFcnCount',round(time*actualsamprate));
ai.SamplesAcquiredFcn = {'stop'};


ai
disp('SampleRate'); handles.srate
disp('Gain'); handles.gain
disp('maxsignal'); handles.maxsignal 
disp('Filename'); get(ai,'LogFileName')
disp('Time'); handles.time

handles.ai=ai; save ai ai;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function set_channels_Callback(hObject, eventdata, handles)
% hObject    handle to set_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_channels as text
%        str2double(get(hObject,'String')) returns contents of set_channels as a double

% Get the new value for the Kf Gain
NewStrVal = get(hObject, 'String');
NewVal = str2num(NewStrVal);
% Check that the entered value falls within the allowable range
if  isempty(NewVal) | (length(NewVal)< 1) | (length(NewVal)>8),
    disp('Error: Check N_channels');
else, % Use new Chan value
    %handles.N_chan=0:NewVal-1;
    handles.N_chan=NewVal;
    handles.load_ch=handles.N_chan;
end
%addchannel(handles.ai,handles.N_chan);
disp('Channels'); handles.N_chan
disp('!!!!!Re-make AI!!!!!')
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function set_channels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to set_channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Set_gain_Callback(hObject, eventdata, handles)
% hObject    handle to Set_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Set_gain as text
%        str2double(get(hObject,'String')) returns contents of Set_gain as a double

NewStrVal = get(hObject, 'String');
NewVal = str2double(NewStrVal);
% Check that the entered value falls within the allowable range
if  isempty(NewVal) | (NewVal< 1) | (NewVal>10000),
    disp('Error: Check Gain');    
else, % Use new Chan value
    handles.gain=NewVal;
    %set(handles.gain_slider,'Value',NewVal); 
end

disp('Gain'); handles.gain
maxsignal_V = handles.maxsignal/(10^6); ip = maxsignal_V*handles.gain;
setverify(handles.ai.Channel,'InputRange',[-ip ip]);   % in V: Any signal above this will be clipped
setverify(handles.ai.Channel,'UnitsRange',[-handles.maxsignal handles.maxsignal]);    % Set units range (/sensor range) to be same, in uV, as of input range (eg [-1000 1000]uV units range for [-1000 1000]uV input range)
setverify(handles.ai.Channel,'SensorRange',[-ip ip]);     %% corresponding to un

disp('Gain (and all ranges) Updated');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Set_gain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Set_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Set_maxsignal_Callback(hObject, eventdata, handles)
% hObject    handle to Set_maxsignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Set_maxsignal as text
%        str2double(get(hObject,'String')) returns contents of Set_maxsignal as a double

NewStrVal = get(hObject, 'String');
NewVal = str2double(NewStrVal);
% Check that the entered value falls within the allowable range
if  isempty(NewVal) | (NewVal< 1) | (NewVal>5*(10^6)),
    disp('Error: Check MaxSignal');    
else, % Use new Chan value
    handles.maxsignal=NewVal;
end

disp('Max Signal'); handles.maxsignal
maxsignal_V = handles.maxsignal/(10^6); ip = maxsignal_V*handles.gain;
setverify(handles.ai.Channel,'InputRange',[-ip ip]);   % in V: Any signal above this will be clipped
setverify(handles.ai.Channel,'UnitsRange',[-handles.maxsignal handles.maxsignal]);    % Set units range (/sensor range) to be same, in uV, as of input range (eg [-1000 1000]uV units range for [-1000 1000]uV input range)
setverify(handles.ai.Channel,'SensorRange',[-ip ip]);     %% corresponding to un

disp('MaxSignal (and all ranges) Updated');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Set_maxsignal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Set_maxsignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Set_Rate_Callback(hObject, eventdata, handles)
% hObject    handle to Set_Rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Set_Rate as text
%        str2double(get(hObject,'String')) returns contents of Set_Rate as a double

NewStrVal = get(hObject, 'String');
NewVal = str2double(NewStrVal);
% Check that the entered value falls within the allowable range
if  isempty(NewVal) | (NewVal< 20000) | (NewVal>31000),
    disp('Error: Check Rate');   
else, % Use new Chan value
    handles.srate=NewVal;
end

disp('SampleRate'); handles.srate
setverify(handles.ai,'SampleRate',handles.srate);

actualsamprate=get(handles.ai,'SampleRate');
setverify(handles.ai,'SamplesAcquiredFcnCount',round(handles.time*actualsamprate));
handles.ai.SamplesAcquiredFcn = {'stop'};

disp('Sample Rate, and samples to get Updated')
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Set_Rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Set_Rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Set_Time_Callback(hObject, eventdata, handles)
% hObject    handle to Set_Time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Set_Time as text
%        str2double(get(hObject,'String')) returns contents of Set_Time as a double

NewStrVal = get(hObject, 'String');
NewVal = str2double(NewStrVal);
% Check that the entered value falls within the allowable range
if  isempty(NewVal) | (NewVal< 0.1) | (NewVal>Inf),
    disp('Error: Check Time');   
else, % Use new Chan value
    handles.time=NewVal;
end

disp('Time'); handles.time
%setverify(handles.ai,'TimerPeriod',handles.time),

actualsamprate=get(handles.ai,'SampleRate');
setverify(handles.ai,'SamplesAcquiredFcnCount',round(handles.time*actualsamprate));
handles.ai.SamplesAcquiredFcn = {'stop'};


disp('Time, and samples to get Updated')
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Set_Time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Set_Time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function buffer_size_Callback(hObject, eventdata, handles)
% hObject    handle to buffer_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of buffer_size as text
%        str2double(get(hObject,'String')) returns contents of buffer_size as a double

NewStrVal = get(hObject, 'String');
NewVal = str2double(NewStrVal);
% Check that the entered value falls within the allowable range

if  isempty(NewVal) | (NewVal< 1000) | (NewVal>200000),
    disp('Error: Check Buffer Size');   
else, % Use new value
    handles.buffersize=NewVal;
end
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function buffer_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to buffer_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function buffer_blocks_Callback(hObject, eventdata, handles)
% hObject    handle to buffer_blocks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of buffer_blocks as text
%        str2double(get(hObject,'String')) returns contents of buffer_blocks as a double

NewStrVal = get(hObject, 'String');
NewVal = str2double(NewStrVal);
% Check that the entered value falls within the allowable range

if  isempty(NewVal) | (NewVal< 1) | (NewVal>50),
    disp('Error: Check Buffer Blocks');   
else, % Use new value
    handles.bufferblocks=NewVal;
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function buffer_blocks_CreateFcn(hObject, eventdata, handles)
% hObject    handle to buffer_blocks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in buffer_do.
function buffer_do_Callback(hObject, eventdata, handles)
% hObject    handle to buffer_do (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

setverify(handles.ai,'BufferingMode','Manual');
setverify(handles.ai,'BufferingConfig',[handles.buffersize handles.bufferblocks]);
disp('Buffer Updated: Size and Blocks'); 
handles.buffersize, handles.bufferblocks,
guidata(hObject, handles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function set_file_Callback(hObject, eventdata, handles)
% hObject    handle to set_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of set_file as text
%        str2double(get(hObject,'String')) returns contents of set_file as a double

[Savefile,Pathsavefile] = uiputfile('*.*','Save file name','test');
if isequal(Savefile,0)
    disp('User selected Cancel')
else
    handles.Savefile = Savefile;
end
disp(['User selected', fullfile(Pathsavefile, Savefile)])
setverify(handles.ai,'LogFileName',[handles.Savefile]);

disp('Savefile Updated')
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in fileindex.
function fileindex_Callback(hObject, eventdata, handles)
% hObject    handle to fileindex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns fileindex contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fileindex

list1 = get(handles.fileindex,'String');
index1 = get(handles.fileindex,'Value');
handles.index = handles.fileindex;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function fileindex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileindex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% --- Executes on button press in start_ai.
function start_ai_Callback(hObject, eventdata, handles)
% hObject    handle to start_ai (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

start(handles.ai);
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in stop_ai.
function stop_ai_Callback(hObject, eventdata, handles)
% hObject    handle to stop_ai (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

stop(handles.ai);
guidata(hObject, handles);

% --- Executes on button press in load_ai.
function load_ai_Callback(hObject, eventdata, handles)
% hObject    handle to load_ai (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ai=handles.ai;
save ai ai
%load ai
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in load_daqfile.
function load_daqfile_Callback(hObject, eventdata, handles)
% hObject    handle to load_daqfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[Filename,Pathname] = uigetfile('*.daq','Select the daq-file');
if isequal(Filename,0)
    disp('User selected Cancel')
else
    disp(['User selected', fullfile(Pathname, Filename)])
    handles.daqfile = Filename;
    handles.pathname = Pathname;
    if length(handles.load_nsecs)==0
        [handles.loaddata handles.loadtime]=daqread([handles.daqfile], 'Channels', handles.load_ch+1);
    else
        [handles.loaddata handles.loadtime]=daqread([handles.daqfile],'Time',handles.load_nsecs, 'Channels', handles.load_ch+1);
    end
end
guidata(hObject, handles);

% --- Executes on button press in Which_Channels.
function Which_Channels_Callback(hObject, eventdata, handles)
% hObject    handle to Which_Channels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tlines={'Channels - Enter channels to load with ; separator'};
defent={''}; % default entries
infonames= {'ch'};
info = inputdlg(Tlines, 'For all ch, all or no change', 1, defent); 

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    ch = str2num(info.ch);   %convert string to number
    handles.load_ch = ch;
    disp('User chose following channels')
    handles.load_ch 
else
    disp('User selected Cancel')
end

guidata(hObject, handles);

% --- Executes on button press in load_nsecs.
function load_nsecs_Callback(hObject, eventdata, handles)
% hObject    handle to load_nsecs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tlines={'Time - Enter time RANGE to load'};
defent={'0, 5'}; % default entries
infonames= {'nsecs'};
info = inputdlg(Tlines, 'For whole file: all', 1, defent); 

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    nsecs = str2num(info.nsecs);   %convert string to number
    handles.load_nsecs = nsecs;
    disp('User chose following nsecs')
    handles.load_nsecs 
else
    disp('User selected Cancel')
end

guidata(hObject, handles);

% --- Executes on button press in load_nblocks.
function load_nblocks_Callback(hObject, eventdata, handles)
% hObject    handle to load_nblocks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tlines={'Blocks - Enter n_blocks to load'};
defent={''}; % default entries
infonames= {'nblocks'};
info = inputdlg(Tlines, 'nblocks', 1, defent); 

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    nblocks = str2num(info.nblocks);   %convert string to number
    handles.load_nblocks = nblocks;
    disp('User chose following nblocks')
    handles.load_nblocks 
else
    disp('User selected Cancel')
end

guidata(hObject, handles);

% --- Executes on button press in Plot_Ch.
function Plot_Ch_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_Ch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tlines={'Enter Loaded Channels to plot'};
defent={''}; % default entries
infonames= {'ch'};
info = inputdlg(Tlines, 'ch', 1, defent); 

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    ch = str2num(info.ch);   %convert string to number
    disp('User chose following channels')
    ch
    disp('Converting to ...')
    ch=ch+1
else
    disp('User selected Cancel')
end

user_response = questdlg('Really Plot?','Plot?','Yes','No','Yes')
switch user_response
case {'No'}
    % take no action
    disp('User selected Cancel')
case 'Yes'
    % 
    figure(100); 
    for i=1:length(ch), subplot(length(ch),1,i); plot(handles.loadtime, handles.loaddata(:,i)); axis tight; end
    figure(200); hold on;
    for i=1:length(ch), plot(handles.loadtime, handles.loaddata(:,i),[handles.clr(i)]); end
    
end

guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%MENU%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%DISPLAY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------------------------------------------------------------------

function display_Callback(hObject, eventdata, handles)
% hObject    handle to display (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --------------------------------------------------------------------

function display_ai_Callback(hObject, eventdata, handles)
% hObject    handle to display_ai (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp('Current AI object')
handles.ai
guidata(hObject, handles);

% --------------------------------------------------------------------

% --------------------------------------------------------------------
function BufferDisp_Callback(hObject, eventdata, handles)
% hObject    handle to BufferDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp('Buffer Size and Blocks')
get(handles.ai,'BufferingMode'),  get(handles.ai,'BufferingConfig')

guidata(hObject, handles);

% --------------------------------------------------------------------


% --------------------------------------------------------------------
function Refresh_Callback(hObject, eventdata, handles)
% hObject    handle to Refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%%%%%%%%%%%%%%%%%%%%%SOFTSCOPE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in softscope.
function softscope_Callback(hObject, eventdata, handles)
% hObject    handle to softscope (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ai=handles.ai;
save ai ai
%load ai

[Filename,Pathname] = uigetfile('*.si','Select the Osc-setting-file');
if isequal(Filename,0)
    disp('User selected Cancel')
else
    disp(['User selected', fullfile(Pathname, Filename)])
    handles.setfile = Filename;
    handles.pathname = Pathname;
end
guidata(hObject, handles);


softscope(handles.ai)
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



