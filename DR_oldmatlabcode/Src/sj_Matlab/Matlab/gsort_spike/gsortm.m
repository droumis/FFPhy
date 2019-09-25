function varargout = gsortm(varargin)
% GUI_SPKSORT M-file for gui_spksort.fig
%      GUI_SPKSORT, by itself, creates a new GUI_SPKSORT or raises the existing
%      singleton*.
%
%      H = GUI_SPKSORT returns the handle to a new GUI_SPKSORT or the
%      handle to
%      the existing singleton*.
%
%      GUI_SPKSORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SPKSORT.M with the given input arguments.
%
%      GUI_SPKSORT('Property','Value',...) creates a new GUI_SPKSORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_spksort_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property
%      application
%      stop.  All inputs are passed to gui_spksort_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_spksort

% Last Modified by GUIDE v2.5 15-Apr-2011 11:54:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gsortm_OpeningFcn, ...
    'gui_OutputFcn',  @gsortm_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% ------------------------------------------------------------------
% OPENING FUNCTION

% --- Executes just before gui_spksort is made visible.
function gsortm_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_spksort (see VARARGIN)

warning off;

% Choose default command line output for gui_spksort
handles.output = hObject;

% Update handles structure
handles.spkfile = [];
handles.spikes = [];
handles.Spiketype = 1;
handles.elec='1';
handles.Subset=[];
handles.draw=[];
handles.parameters.outliers_a = [];
handles.parameters.kmeans_options.divisions = [];
handles.parameters.reint_out = []; handles.parameters.tmin = []; handles.parameters.tref=[]; handles.parameters.cutoff=[];;
handles.rep=0; handles.repNF=0;


handles.Savefile =[]; 
handles.spikes1=[]; handles.spikes2=[]; handles.spikes3=[]; handles.spikes4=[]; handles.spikes5=[];
handles.parameters1=[]; handles.parameters2=[]; handles.parameters3=[]; handles.parameters4=[]; handles.parameters5=[];
handles.plot_sd=0;


handles.tmin=0.001;         % TMIN WILL BE IN SPIKES.PARAMETERS IN FUTURE VERSIONS
handles.arfp = handles.tmin*1000;

%%% FOR AUTOSPIKECUT PARAMETERS
handles.spkcut = [];

% SOME PLOTTING VARAIBLES
handles.plot_ampl = 400; handles.plot_type = 1; handles.nch=4;

% DECLARING VARIABLES MAKES IT EASIER TO REFERESH

handles.guifighandle=get(0,'CurrentFigure');
set(handles.guifighandle, 'HandleVisibility', 'Off');

handles.oriunits = get(handles.guifighandle,'Units');
handles.oriposition = get(handles.guifighandle,'Position');
set(0,'Units','Characters');
handles.screensize = get (0, 'Screensize');
set(0,'Units','pixels');

set(handles.guifighandle,'Position',[1.4 5.9 handles.screensize(3)/5 handles.screensize(4)-10 ]);
handles.position = get(handles.guifighandle,'Position');




% a= get (0, 'Screensize');
% set(gcf, 'Position', [a(3)*0.3 a(4)*(0.51) a(3)*0.68 a(4)*0.42])


handles.parameterfile = [];
handles.clr = ['y' 'm' 'b' 'c' 'r' 'y' 'm' 'b' 'c' 'r' 'y' 'm' 'b' 'c' 'r'];


% if exist('default_parameters.mat','file')==2
%     load default_parameters;
%  
%     supplied = lower(fieldnames(parameters));   % which options did the user specify?
%     for op = 1:length(supplied)              % copy those over the parameters
%         if (version('-release') < 13)        % annoyingly, pre-R13 Matlab doesn't do dynamic field names, ...
%             opts = setfield(opts, supplied(op), getfield(options, supplied(op)));  % so we use an older syntax
%         else
%             handles.parameters.(supplied{op}) = parameters.(supplied{op});  % this is the preferred syntax as of R13 --
%         end                                                %   we include it b/splice_waves 'setfield' is deprecated
%     end
% end


guidata(hObject, handles);
% set(hObject, 'HandleVisibility', 'Off');



% UIWAIT makes gui_spksort wait for user response (see UIRESUME)
% uiwait(handles.Main_GUI);

% ------------------------------------------------------------------
% OUTPUT FUNCTION

% --- Outputs from this function are returned to the command line.
function varargout = gsortm_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.spikes;

% ------------------------------------------------------------------

function Main_GUI_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.guifighandle=get(0,'CurrentFigure');
handles.position = get(handles.guifighandle,'Position');
guidata(hObject, handles);

%--------------------------------------------------------------------


%--------------------START CUSTOM FUNCTIONS--------------------------

%************************************************************************
% MENU ITEMS
%************************************************************************




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. SPIKECUT FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Shantanu's AUTOSPIKECUT Functions
function Spikecut_Callback(hObject, eventdata, handles)
% hObject    handle to Spikecut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);

% --------------------------------------------------------------------
function spkcut_params_Callback(hObject, eventdata, handles)
% hObject    handle to spkcut_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%handles.spkcut = [];  %%  REFRESH INFO BEFORE START: dont: use refresh button   
disp('Choose params')

%Tlines={['Std Dev above Thrs (1 value only)'],['StimCode-Please enter' num2str(size(ascii_cell,1))...
Tlines={['Choose channels to cut'],['Sampling Freq'],['Std Dev above Thrs (1 value only)'],...
    ['Threshold: +ve spikes,thr=0; -ve spikes, thr=1'],['StimCode-Please enter no_of_files' ...
    'stimcodes as an array'],['Push Sweeps?(default = 1:yes)'], ['Sweepall: sweep range'], ['Sweepstim:'...
    'Cell array of sw ranges: 1 for each file'],...
    ['p.stimonset(default=200ms)'],['p.bckwindow(default=100ms)'],['p.sweepd(default=1000ms)'] };

infonames = {'Ch','Fs','sdcuts', 'threshold', 'stimcode','consweep','sweepall', 'sweepstim','stimonset','bckwindow','sweepd'};

%defent={'1:4','30000','5', '0','[0 1 2]','1','1:1000','{[1:100 201:300]}; [101:200]','200','100','1000'}; % default entries

defent={'1:4','30000','5', '0','[0]','1','','','','',''}; % default entries


info = inputdlg(Tlines, 'Defaults: Empty or no change', 1, defent); 
if ~isempty(info)              %see if user hit cancel
    str = cell2struct(info,infonames);
    handles.spkcut.Ch = str2num(str.Ch); 
    handles.spkcut.Fs = str2num(str.Fs); 
    handles.spkcut.sdcuts = str2num(str.sdcuts);   %convert string to number
    handles.spkcut.threshold = str2num(str.threshold); 
    handles.spkcut.stimcode = str2num(str.stimcode);   
    handles.spkcut.consweep = str2num(str.consweep);
    handles.spkcut.sweepall = str2num(str.sweepall);
    handles.spkcut.sweepstim = str2num(str.sweepstim);
    handles.spkcut.params.stimonset = str2num(str.stimonset);
    handles.spkcut.params.bckwindow = str2num(str.bckwindow);
    handles.spkcut.params.sweepd = str2num(str.sweepd);
else 
    handles.spkcut = [];
end

disp('Params chosen');
handles.spkcut
guidata(hObject, handles);

% --------------------------------------------------------------------

% --------------------------------------------------------------------
function matlab_cont3_Callback(hObject, eventdata, handles)
% hObject    handle to matlab_cont3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ascii_cell]=uigetfile('*.daq*','Select the .daq Files','MultiSelect','on');
if iscell(ascii_cell), handles.spkcut.files=ascii_cell'; else, handles.spkcut.files={ascii_cell}; end
disp('Files Chosen')
handles.spkcut.files'
disp('SpikeCut in Progress.....')
if size(handles.spkcut.files)~=0
    [spikes,aspkfile]=sss_spkcut_matlabcont3(handles.spkcut.files,...
        handles.spkcut.Ch, handles.spkcut.Fs, handles.spkcut.sdcuts,...
        handles.spkcut.threshold, handles.spkcut.stimcode, handles.spkcut.params,...
        handles.spkcut.consweep, handles.spkcut.sweepall, handles.spkcut.sweepstim);
    %aspkfile = uiputfile('*.mat','Save SPKFILE name');
    handles.spikes=spikes; handles.spkfile = aspkfile;
else
    disp('Choose files and params first!!!')
end

guidata(hObject, handles);

% --------------------------------------------------------------------
function labview_cont3_Callback(hObject, eventdata, handles)
% hObject    handle to labview_cont3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ascii_cell]=uigetfile('*.*','Select the labview Files','MultiSelect','on');
if iscell(ascii_cell), handles.spkcut.files=ascii_cell'; else, handles.spkcut.files={ascii_cell}; end
disp('Files Chosen')
handles.spkcut.files'
disp('SpikeCut in Progress.....')
if size(handles.spkcut.files)~=0
    [spikes,aspkfile]=sss_spkcut_labviewcont3(handles.spkcut.files,...
        handles.spkcut.Ch, handles.spkcut.Fs, handles.spkcut.sdcuts,...
        handles.spkcut.threshold, handles.spkcut.stimcode, handles.spkcut.params,...
        handles.spkcut.consweep, handles.spkcut.sweepall, handles.spkcut.sweepstim);
    %aspkfile = uiputfile('*.mat','Save SPKFILE name');
    handles.spikes=spikes; handles.spkfile = aspkfile;
else
    disp('Choose files and params first!!!')
end

guidata(hObject, handles);


% --------------------------------------------------------------------

function basic_1ch_1ms_Callback(hObject, eventdata, handles)
% hObject    handle to basic_1ch_1ms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ascii_cell]=uigetfile('*.swp*','Select the ascii .swp Files','MultiSelect','on');
if iscell(ascii_cell), handles.spkcut.files=ascii_cell'; else, handles.spkcut.files={ascii_cell}; end
disp('Files Chosen')
handles.spkcut.files'
disp('SpikeCut in Progress.....')
if size(handles.spkcut.files)~=0
    [spikes,aspkfile]=sss_spkcut_1ms_basic(handles.spkcut.files,...
        handles.spkcut.sdcuts, handles.spkcut.stimcode, handles.spkcut.params,...
        handles.spkcut.consweep, handles.spkcut.sweepall, handles.spkcut.sweepstim);
    %aspkfile = uiputfile('*.mat','Save SPKFILE name');
    handles.spikes=spikes; handles.spkfile = aspkfile;
else
    disp('Choose files and params first!!!')
end

guidata(hObject, handles);


% --------------------------------------------------------------------
function basic_1ch_1halfms_Callback(hObject, eventdata, handles)
% hObject    handle to basic_1ch_1halfms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ascii_cell]=uigetfile('*.swp*','Select the ascii .swp Files','MultiSelect','on');
if iscell(ascii_cell), handles.spkcut.files=ascii_cell'; else, handles.spkcut.files={ascii_cell}; end
disp('Files Chosen')
handles.spkcut.files'
disp('SpikeCut in Progress....')
if size(handles.spkcut)~=0
    [spikes,aspkfile]=sss_spkcut_1halfms_basic(handles.spkcut.files,...
        handles.spkcut.sdcuts, handles.spkcut.stimcode, handles.spkcut.params,...
        handles.spkcut.consweep, handles.spkcut.sweepall, handles.spkcut.sweepstim);
    %aspkfile = uiputfile('*.mat','Save SPKFILE name');
    handles.spikes=spikes; handles.spkfile = aspkfile;
else
    disp('Choose files and params first!!!')
end

guidata(hObject, handles);

% --------------------------------------------------------------------
function basic_2ch_Callback(hObject, eventdata, handles)
% hObject    handle to basic_2ch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%handles.spkcut = [];  %%  REFRESH INFO BEFORE START: Parameters chosen before. Dont clear   
[files1]=uigetfile('*.swp*','First Select Ch1 the ascii .swp Files','MultiSelect','on');
if iscell(files1), handles.spkcut.files1=files1'; else, handles.spkcut.files1={files1}; end
[files2]=uigetfile('*.swp*','Now Select the Ch2 ascii .swp Files','MultiSelect','on');
if iscell(files2), handles.spkcut.files2=files2'; else, handles.spkcut.files2={files2}; end

disp('Files Chosen:Ch1'), handles.spkcut.files1', disp('Files Chosen:Ch2'), handles.spkcut.files2'
disp('SpikeCut in Progress.....')
if size(handles.spkcut.files1)~=0
    [spikes,aspkfile]=sss_spkcut_2ch_basic(handles.spkcut.files1, handles.spkcut.files2,...
        handles.spkcut.sdcuts, handles.spkcut.stimcode, handles.spkcut.params,...
        handles.spkcut.consweep, handles.spkcut.sweepall, handles.spkcut.sweepstim);
    %aspkfile = uiputfile('*.mat','Save SPKFILE name');
    handles.spikes=spikes; handles.spkfile = aspkfile;
else
    disp('Choose files and params first!!!')
end

guidata(hObject, handles);

% --------------------------------------------------------------------
function basic_4ch_Callback(hObject, eventdata, handles)
% hObject    handle to basic_4ch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%handles.spkcut = [];  %%  REFRESH INFO BEFORE START   
[files1]=uigetfile('*.swp*','First Select Ch1 the ascii .swp Files','MultiSelect','on');
if iscell(files1), handles.spkcut.files1=files1'; else, handles.spkcut.files1={files1}; end
[files2]=uigetfile('*.swp*','Now Select the Ch2 ascii .swp Files','MultiSelect','on');
if iscell(files2), handles.spkcut.files2=files2'; else, handles.spkcut.files2={files2}; end
[files3]=uigetfile('*.swp*','Now Select the Ch3 ascii .swp Files','MultiSelect','on');
if iscell(files3), handles.spkcut.files3=files3'; else, handles.spkcut.files3={files3}; end
[files4]=uigetfile('*.swp*','Now Select the Ch4 ascii .swp Files','MultiSelect','on');
if iscell(files4), handles.spkcut.files4=files4'; else, handles.spkcut.files4={files4}; end

disp('Files Chosen:Ch1'), handles.spkcut.files1', disp('Files Chosen:Ch2'), handles.spkcut.files2'
disp('Files Chosen:Ch3'), handles.spkcut.files3', disp('Files Chosen:Ch4'), handles.spkcut.files4'
disp('SpikeCut in Progress.....')

if size(handles.spkcut.files)~=0
    [spikes,aspkfile]=sss_spkcut_4ch_basic(handles.spkcut.files1,handles.spkcut.files2,...
        handles.spkcut.files3,handles.spkcut.files4,...
        handles.spkcut.sdcuts, handles.spkcut.stimcode, handles.spkcut.params,...
        handles.spkcut.consweep, handles.spkcut.sweepall, handles.spkcut.sweepstim);
    %aspkfile = uiputfile('*.mat','Save SPKFILE name');
    handles.spikes=spikes; handles.spkfile = aspkfile;
else
    disp('Choose files and params first!!!')
end

guidata(hObject, handles);

% --------------------------------------------------------------------

% --------------------------------------------------------------------
function clear_spkcut_params_Callback(hObject, eventdata, handles)
% hObject    handle to clear_spkcut_params (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.spkcut = [];

guidata(hObject, handles);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2. RAWDATAPLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------------------------------------------------------------------
function rawdataplot_Callback(hObject, eventdata, handles)
% hObject    handle to rawdataplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);


% --------------------------------------------------------------------
function rawload_matlab_Callback(hObject, eventdata, handles)
% hObject    handle to rawload_matlab (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp('Choose Channels and time to load')
Tlines={['Choose channels to cut'],['Time Range in sec']};
infonames = {'Ch','time'};
defent={'1:4','[0,300]'}; % default entries

info = inputdlg(Tlines, 'Defaults: Empty or no change', 1, defent);
if ~isempty(info)              %see if user hit cancel
    str = cell2struct(info,infonames);
    Ch = str2num(str.Ch);
    time = str2num(str.time);

    disp('Raw matlab Params chosen');
    Ch
    time

    [ascii_cell]=uigetfile('*.daq*','Select the ascii .daq File','MultiSelect','on');
    if iscell(ascii_cell), files=ascii_cell'; else, files={ascii_cell}; end
    disp('Files Chosen')
    files'
    disp('Loading in Progress.....')
    if size(files)~=0
        [handles.rawdata handles.rawtime handles.abstime handles.Events handles.daqinfo] =...
            daqread(files,'Channels',Ch,'Time',time);
    else
        disp('Choose files!!!')
    end

else
    disp('Choose Params!');
end

guidata(hObject, handles);


% --------------------------------------------------------------------
function rawload_labview_Callback(hObject, eventdata, handles)
% hObject    handle to rawload_labview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp('Choose Fs and time to load')
Tlines={['Choose Fs'],['Time Range in sec']};
infonames = {'Fs','time'};
defent={'30000','[0,300]'}; % default entries

info = inputdlg(Tlines, 'Defaults: Empty or no change', 1, defent);
if ~isempty(info)              %see if user hit cancel
    str = cell2struct(info,infonames);
    Fs = str2num(str.Fs);
    time = str2num(str.time);

    disp('Raw Labview Params chosen');
    Fs
    time

    [ascii_cell]=uigetfile('*.*','Select the labview File','MultiSelect','on');
    if iscell(ascii_cell), files=ascii_cell'; else, files={ascii_cell}; end
    disp('Files Chosen')
    files'
    disp('Loading in Progress.....')
    if size(files)~=0
        fn=files;
        handles.rawdata=labview_loadspike4(fn,8,Fs,time(1),time(2));  %% NEED TO READ ALL 8 channels
    else
        disp('Choose files!!!')
    end

else
    disp('Choose Params!');
end

guidata(hObject, handles);


% --------------------------------------------------------------------
function rawplot_Callback(hObject, eventdata, handles)
% hObject    handle to rawplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);


% --------------------------------------------------------------------
function rawplot_thrs_Callback(hObject, eventdata, handles)
% hObject    handle to rawplot_thrs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


guidata(hObject, handles);

%------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. REFRESH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in Refresh.
function Refresh_Callback(hObject, eventdata, handles)
% hObject    handle to Refresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% USE TO REFRESH STUFF BEFORE START OF EVERY FILE

guidata(hObject, handles);



% --------------------------------------------------------------------
function Refresh_Spikes_Callback(hObject, eventdata, handles)
% hObject    handle to Refresh_Spikes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.spkfile = [];
handles.spikes = [];
handles.draw=[];
handles.rep=0; handles.repNF=0;
handles.Savefile =[]; 
handles.spikes1=[]; handles.spikes2=[]; handles.spikes3=[]; handles.spikes4=[]; handles.spikes5=[];
handles.parameters1=[]; handles.parameters2=[]; handles.parameters3=[]; handles.parameters4=[]; handles.parameters5=[];
handles.show =[];
msgbox('Cleared!','Refresh','warn')

guidata(hObject, handles);

% --------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. SET PARAMETERS MENU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Parameters_Callback(hObject, eventdata, handles)
% hObject    handle to Parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Display_Parameters_Callback(hObject, eventdata, handles)
% hObject    handle to Display_Parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.spikes.parameters=handles.parameters; handles.spikes.nch=handles.nch;
handles.parameters
guidata(hObject, handles);

% --------------------------------------------------------------------
function Set_Parameter_File_Callback(hObject, eventdata, handles)
% hObject    handle to Set_Parameter_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[ParFile,ParPath] = uigetfile('*.mat','Select the Mat-file with Parameters');
if isequal(ParFile,0)
    disp('User selected Cancel')
else
    disp(['User selected', fullfile(ParPath, ParFile)])
    handles.parameterfile = ParFile;


    load(handles.parameterfile);

    supplied = lower(fieldnames(parameters));   % which options did the user specify?
    for op = 1:length(supplied)              % copy those over the defaults
        if length(parameters.(supplied{op}))~=0
            if (version('-release') < 13)        % annoyingly, pre-R13 Matlab doesn't do dynamic field names, ...
                opts = setfield(opts, supplied(op), getfield(options, supplied(op)));  % so we use an older syntax
            else
                handles.parameters.(supplied{op}) = parameters.(supplied{op});  % this is the preferred syntax as of R13 --
            end                                                %   we include it b/c 'setfield' is deprecated
        end
    end
end

handles.spikes.parameters=handles.parameters;  handles.spikes.nch=handles.nch;

handles.parameters
guidata(hObject, handles);



% --------------------------------------------------------------------
function Clear_Parameter_File_Callback(hObject, eventdata, handles)
% hObject    handle to Clear_Parameter_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.parameterfile = [];
guidata(hObject, handles);


% --------------------------------------------------------------------
function Choose_Parameters_Callback(hObject, eventdata, handles)
% hObject    handle to Choose_Parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

params = choose_params2(hObject,handles);

if length(params)==0
    disp('User selected Cancel: Defaults not chosen')
else

    supplied = lower(fieldnames(params.parameters));   % which options did the user specify?
    for op = 1:length(supplied)              % copy those over the defaults
        if length(params.parameters.(supplied{op}))~=0
            if (version('-release') < 13)        % annoyingly, pre-R13 Matlab doesn't do dynamic field names, ...
                opts = setfield(opts, supplied(op), getfield(options, supplied(op)));  % so we use an older syntax
            else
                handles.parameters.(supplied{op}) = params.parameters.(supplied{op});  % this is the preferred syntax as of R13 --
            end
        end

     % clustnum = str2num(info.Clusters);  %   we include it b/splice_waves 'setfield' is deprecated
     % disp(['User selected', supplied{op} params.parameters.(supplied{op}) ]);
    
    end
end

handles.spikes.parameters=handles.parameters; handles.spikes.nch=handles.nch;
handles.parameters
guidata(hObject, handles);


% --------------------------------------------------------------------

function Clear_Parameters_Callback(hObject, eventdata, handles)
% hObject    handle to Clear_Parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

supplied = lower(fieldnames(handles.parameters));  
for op = 1:length(supplied)  
    handles.parameters.(supplied{op})=[];
end

guidata(hObject, handles);



% --------------------------------------------------------------------
function Make_Parameter_File_Callback(hObject, eventdata, handles)
% hObject    handle to Make_Parameter_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

params = make_params(hObject,handles);
parameters = params.parameters;

if length(params)==0
    disp('User selected Cancel: Defaults File not made')
else

    defname = ['default_parameters'];
    [Savedeffile,Pathsavedeffile] = uiputfile('*.mat','Save file name',[defname]);
    handles.parameterfile = Savedeffile;
    save (handles.parameterfile, 'parameters');
end

guidata(hObject, handles);

% --------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Noise Plots MENU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --------------------------------------------------------------------
function Noise_Callback(hObject, eventdata, handles)
% hObject    handle to Noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);

% --------------------------------------------------------------------

% --------------------------------------------------------------------
function Plot_Den_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_Den (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure; colormap hot; hold on;
[n,x,y] = hist2d(handles.spikes.waveforms);
imagesc(x,y,n); axis xy; colormap hot;
axis tight;
plot (mean (handles.spikes.waveforms), 'y', 'Linewidth', 2);
ylabel ('Voltage (uV)');
xlabel ('Sample (32 pts/ms)')
text (0.8*size(handles.spikes.waveforms,2), 0.8*max(max(handles.spikes.waveforms)), ['#Sp: ' num2str(length(handles.spikes.waveforms))]);

% --------------------------------------------------------------------
function Plot_Spikes_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_Spikes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure; colormap hot; hold on;
if (size (handles.spikes.waveforms,1) > 1000)
    idxs = randperm (size (handles.spikes.waveforms,1));
    unit_sm = handles.spikes.waveforms (idxs(1:1000), :);
else
    unit_sm=handles.spikes.waveforms;
end
plot(unit_sm','b.-'); hold on;
plot (mean (handles.spikes.waveforms), 'y', 'Linewidth', 2);
text (0.8*size(handles.spikes.waveforms,2), 0.8*max(max(unit_sm)), ['#Sp: ' num2str(length(handles.spikes.spiketimes))]);
h = gca;


% --------------------------------------------------------------------
function Manual_Removal_Callback(hObject, eventdata, handles)
% hObject    handle to Manual_Removal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cmd=sprintf('handles.spikes%d = handles.spikes;',handles.repNF); eval(cmd);
handles.repNF = handles.repNF + 1;
handles.spikes = sss_windowcut(handles.spikes);

string=strtok(handles.spkfile,'.');
filename1 = [string '_NF'];
spikes=handles.spikes;
handles.spikes=handles.spikes;

[SaveNFfile,PathsaveNFfile] = uiputfile('*.mat','Save file name',[filename1]);

%handles.spkNFfile = Savedeffile;
if isequal(SaveNFfile,0)
    disp('User selected Cancel')
else
    save (SaveNFfile, 'spikes');
end

guidata(hObject, handles);

% --------------------------------------------------------------------

function Tetrode_Noise_Cut_Callback(hObject, eventdata, handles)
% hObject    handle to Tetrode_Noise_Cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cmd=sprintf('handles.spikes%d = handles.spikes;',handles.repNF); eval(cmd);
handles.repNF = handles.repNF + 1;
handles.spikes = sss_windowcut_tet(handles.spikes, handles.nch);

string=strtok(handles.spkfile,'.');
filename1 = [string '_tetNF'];
spikes=handles.spikes;
handles.spikes=handles.spikes;

[SaveNFfile,PathsaveNFfile] = uiputfile('*.mat','Save file name',[filename1]);

%handles.spkNFfile = Savedeffile;
if isequal(SaveNFfile,0)
    disp('User selected Cancel')
else
    save (SaveNFfile, 'spikes');
end

guidata(hObject, handles);

% --------------------------------------------------------------------

function hist_tetrode_noise_cut_Callback(hObject, eventdata, handles)
% hObject    handle to hist_tetrode_noise_cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cmd=sprintf('handles.spikes%d = handles.spikes;',handles.repNF); eval(cmd);
handles.repNF = handles.repNF + 1;
handles.spikes = sss_windowcut_tet_hist(handles.spikes, handles.nch);

string=strtok(handles.spkfile,'.');
filename1 = [string '_tetNF'];
spikes=handles.spikes;
handles.spikes=handles.spikes;

[SaveNFfile,PathsaveNFfile] = uiputfile('*.mat','Save file name',[filename1]);

%handles.spkNFfile = Savedeffile;
if isequal(SaveNFfile,0)
    disp('User selected Cancel')
else
    save (SaveNFfile, 'spikes');
end

guidata(hObject, handles);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     PANELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%************************************************************************
% LOAD PANEL
%************************************************************************


% ------------------------------------------------------------------
% 1. LOAD SPKFILE 
%-------------------------------------------------------------------

% --- Executes on button press in Load_file.
function Load_file_Callback(hObject, eventdata, handles)
% hObject    handle to Load_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% d = dir('*.mat');
% str = {d.name};
% [s,v] = listdlg('PromptString','Select a file:','SelectionMode','single','ListString',str);
% handles.spkfile = str{s};


%% First get channels %%
Tlines={'How many channels?'};
defent={'4'}; % default entries
infonames= {'Channels'};
info = inputdlg(Tlines, 'N_channels', 1, defent);

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    nch = str2num(info.Channels);   %convert string to number
    handles.nch=nch;
end
disp('Nchannels')
handles.nch



%% Now get file
[Filename,Pathname] = uigetfile('*.mat','Select the Mat-file with spikes');
if isequal(Filename,0)
    disp('User selected Cancel')
else
    disp(['User selected', fullfile(Pathname, Filename)])
    handles.spkfile = Filename;
    handles.pathname = Pathname;

    x = load(handles.spkfile,'spikes');
    
    % If there is no waveform field, generate it from individual channel
    % waveforms
    if (~isfield(x.spikes, 'waveforms')), 
        x.spikes.waveforms = [];
        for i=1:handles.nch,
            cmd = sprintf('x.spikes.waveforms = [x.spikes.waveforms, x.spikes.waveforms_ch%d];',i); eval(cmd);
        end
    end
        
    handles.spikes = x.spikes;
    
    % Raw spikes does not need individual channel info - just get rid of this
%     for i=1:handles.nch,
%             cmd = sprintf('x.spikes.waveforms_ch%d=[];',i); eval(cmd);
%     end
%     handles.rawspikes = x.spikes; 
    clear x;
    
end
guidata(hObject, handles);


% ------------------------------------------------------------------
% 2. Electrode Number - Default 1
% ------------------------------------------------------------------


% --- Executes on button press in Elec_No.
function Elec_No_Callback(hObject, eventdata, handles)
% hObject    handle to Elec_No (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tlines={'Enter Electrode No'};
defent={'1'}; % default entries
infonames= {'Elec'};
info = inputdlg(Tlines, 'Elec', 1, defent);

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    elec = info.Elec;   % keep as string 
    handles.elec=elec;
end
disp('Electrode No')
handles.elec


%-----------------------------------------------------------------------
% Spiketype - +ve or -ve Threshold for Spikes needed for Dejittering 
%----------------------------------------------------------------------

function Spiketype_Callback(hObject, eventdata, handles)
% hObject    handle to Spiketype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tlines={'Spiketype: positive thrs (1)  or negative thrs (0)?'};
defent={'1'}; % default entries
infonames= {'Spiketype'};
info = inputdlg(Tlines, 'Spiketype', 1, defent);

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    Spiketype = str2num(info.Spiketype);   %convert string to number
    handles.Spiketype=Spiketype;
end
disp('Spiketype')
handles.Spiketype
guidata(hObject, handles);



%************************************************************************
% Pre-Sort Panel - Do PCA of data. Plot Amplitude or PCA before clustering
%************************************************************************


% ---------------------------------------------------------------------
% 1a.) Do PCA and save as spikes.pca
% ---------------------------------------------------------------------

% --- Executes on button press in do_pca.
function do_pca_Callback(hObject, eventdata, handles)
% hObject    handle to do_pca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% SVD the data because we want to use the PCs as axes.
%[pca.scores,pca.u,pca.s,pca.v] = pcasvd(handles.spikes.waveforms);


if isfield(handles.spikes,'pca')
    disp('PCA Field already exists in spike structure')
    disp('Hit Cancel if you do not want to re-calculate PCA')
end
Tlines={'How many channels on electrode?',['Use L2 norm(1) or simple svd (0)']};
defent={'4','1'}; % default entries
infonames= {'Channels','PCAType'};
info = inputdlg(Tlines, 'N_channels and PCAType', 1, defent);

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    nch = str2num(info.Channels);   %convert string to number
    pcatype = str2num(info.PCAType);   %convert string to number
    handles.nch=nch;
end
disp('Nchannels')
handles.nch

%% Now do PCA %%
[handles.spikes.pcadata] = sj_pcasvd(handles.spikes.waveforms, handles.nch, pcatype);
disp('PCA Done');
%handles.spikes.pca = pcadata;
%data = pca.scores;
guidata(hObject, handles);

% ---------------------------------------------------------------------
% 1a.) DoAmpl - and save as spikes.ampl
% ---------------------------------------------------------------------


% --- Executes on button press in do_ampl.
function do_ampl_Callback(hObject, eventdata, handles)
% hObject    handle to do_ampl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles.spikes,'ampl')
    disp('Ampl Field already exists in spike structure')
else
    
    %% Now do ampl %%
    [handles.spikes.ampl] = sj_do_ampl(handles.spikes, handles.nch);
    disp('Ampl Done');
end
%handles.spikes.pca = pcadata;
%data = pca.scores;
guidata(hObject, handles);



% ---------------------------------------------------------------------
% 2.) pAmpPl - Presort Amplitude Plots 
% ---------------------------------------------------------------------

% --- Executes on button press in presort_amppl.
function presort_amppl_Callback(hObject, eventdata, handles)
% hObject    handle to presort_amppl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in .

Tlines={'Plot PCA? 1 for yes, 0 for no'};
defent={'1'}; % default entries
infonames= {'plotpca'};
info = inputdlg(Tlines, 'Plot PCA?', 1, defent);

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    plotpca = str2num(info.plotpca);   %convert string to number
else
    dopca=0;
end

sss_presort_ampplot(handles.spikes,handles.nch,plotpca);    
guidata(hObject, handles);

%--------------------------------------------------------------------------

% ---------------------------------------------------------------------
% 3.) pAmp2D - 2d Amplitude Plot 
% ---------------------------------------------------------------------

% --- Executes on button press in pAmp2D.
function pAmp2D_Callback(hObject, eventdata, handles)
% hObject    handle to pAmp2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 sss_presort_amp(handles.spikes,'xy', handles.nch);
 guidata(hObject, handles);

% ---------------------------------------------------------------------
% 4.) pAmp3D - 3d Amplitude Plot 
% ---------------------------------------------------------------------

 
% --- Executes on button press in pAmp3D.
function pAmp3D_Callback(hObject, eventdata, handles)
% hObject    handle to pAmp3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sss_presort_amp(handles.spikes,'xyz', handles.nch);
guidata(hObject, handles);

% %------------------------------------------------------------------------
% --

% ---------------------------------------------------------------------
% 5.) pPCA2D - 2d PCA Plot 
% ---------------------------------------------------------------------


% --- Executes on button press in presort_PCA2d.
function presort_PCA2d_Callback(hObject, eventdata, handles)
% hObject    handle to presort_PCA2d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.assigns =[]; handles.z= 0; %handles.pc =1;

%ssgtest(handles.spikes,handles.assigns, handles.z);
sss_presort_pca(handles.spikes,'xy',handles.nch);
guidata(hObject, handles);

% ---------------------------------------------------------------------
% 6.) pPCA3D - 3d PCA Plot 
% ---------------------------------------------------------------------


% --- Executes on button press in presort_PCA3d.
function presort_PCA3d_Callback(hObject, eventdata, handles)
% hObject    handle to presort_PCA3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.assigns =[]; handles.z= 0; %handles.pc =1;
guidata(hObject, handles);

%ssgtest(handles.spikes,handles.assigns, handles.z);
sss_presort_pca(handles.spikes,'xyz',handles.nch);
guidata(hObject, handles);






%************************************************************************
%. 1ch Sort PANEL
%************************************************************************


% ------------------------------------------------------------------
% 1. RUN_SORT - All steps of algorithm for single-channel data 
% ------------------------------------------------------------------

% --- Executes on button press in Run_Sort.
function Run_Sort_Callback(hObject, eventdata, handles)
% hObject    handle to Run_Sort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% MAIN SORTING FUNCTION

handles.rep = handles.rep + 1;

handles.spikes.parameters=handles.parameters; handles.spikes.nch=handles.nch;
[handles.spikes,draw] = sss_spksort(handles.spikes, handles.parameters,handles.spkfile,handles.Spiketype);

handles.draw = draw; handles.tempspikes = handles.spikes;
cmd=sprintf('handles.spikes%d = handles.spikes;',handles.rep); eval(cmd);
handles.spikes.parameters = handles.parameters;
cmd=sprintf('handles.parameters%d = handles.parameters;',handles.rep); eval(cmd);
handles.spikes.draw=draw;
handles.show = unique(handles.spikes.hierarchy.assigns);

guidata(hObject, handles);



% ---------------------------------------------------------------------
% 2. Save - Choose File Name for Saving Sort Results
% ---------------------------------------------------------------------

% --- Executes on button press in SaveSort_Yes.
function SaveSort_Yes_Callback(hObject, eventdata, handles)
% hObject    handle to SaveSort_Yes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set(hObject, 'HandleVisibility', 'Off');
set(handles.guifighandle, 'HandleVisibility', 'Off');
close all;
%cmd=sprintf('handles.spikes = handles.spikes%d;',handles.rep); eval(cmd);
string=strtok(handles.spkfile,'.');
filename1 = [string '_sort'];
spikes=handles.spikes;


[Savefile,Pathsavefile] = uiputfile('*.mat','Save file name',[filename1]);
if isequal(Savefile,0)
    disp('User selected Cancel')
else
    handles.Savefile = Savefile;
    save (handles.Savefile, 'spikes');
end

guidata(hObject, handles);


%----------------------------------------------------------------------
% 3. Repeat - Does not Save current Data. Can Run algorithm again
%----------------------------------------------------------------------

% --- Executes on button press in Repeat_Sort.
function Repeat_Sort_Callback(hObject, eventdata, handles)
% hObject    handle to Repeat_Sort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% NO NEED TO WRITE LOOP FOR SAME FILE, SINCE HANDLES.REP UPDATED EACH TIME
% RUN_SORT IS CALLED.

% set(hObject, 'HandleVisibility', 'Off');
set(handles.guifighandle, 'HandleVisibility', 'Off');
close all;

msgbox('Sort File Not Saved', 'Re-define params/etc & Re-run');

guidata(hObject, handles);

%------------------------------------------------------------------------




% *********************************************************************
% Ntrode Sorting - Sorting Tetrode or Any multiple Channel Data
% *********************************************************************


%---------------------------------------------------------------------
% 1.) Tet_Sort - Executes All Functions of Algorithm. De-Jitter has to be
% working for your data to use this
%---------------------------------------------------------------------


% --- Executes on button press in Tet_Sort.
function Tet_Sort_Callback(hObject, eventdata, handles)
% hObject    handle to Tet_Sort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.spikes.parameters=handles.parameters; handles.spikes.nch=handles.nch;
[handles.spikes,draw] = ssm_tetsort(handles.spkfile, handles.parameters, handles.nch, handles.Spiketype, handles.elec, handles.Subset);

handles.draw = draw;      % handles.tempspikes = tempspikes;
cmd=sprintf('handles.spikes%d = handles.spikes;',handles.rep); eval(cmd);
handles.spikes.parameters = handles.parameters;
cmd=sprintf('handles.parameters%d = handles.parameters;',handles.rep); eval(cmd);
handles.spikes.draw=draw;
handles.show = unique(handles.spikes.hierarchy.assigns);

guidata(hObject, handles);


%---------------------------------------------------------------------
% STEPWISE SORT FUNCTIONS FOR TETRODE
%--------------------------------------------------------------------

%---------------------------------------------------------------------
% 2.) Dejitter - Spiketype has to be set properly for De-Jitterring to 
% work. Skip This Step if it is erroring out. 
%---------------------------------------------------------------------


% --- Executes on button press in DeJitter.
function DeJitter_Callback(hObject, eventdata, handles)
% hObject    handle to DeJitter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Tlines={['Enter maxshift']};
% ch = inputdlg(Tlines, 'Maxshift', 1, {'3'}); 

maxshift=3;

% figure; colormap hot; 
% if (size (handles.spikes.waveforms,1) > 1000)
%     idxs = randperm (size (handles.spikes.waveforms,1));
%     unit_sm = handles.spikes.waveforms (idxs(1:1000), :);
% else
%     unit_sm=handles.spikes.waveforms;
% end
% subplot(2,1,1); plot(unit_sm','b.-'); axis tight; title('Original Data Subset');

handles.spikes = ssm_dejitter_fortet1(handles.spikes,maxshift,handles.nch);

if (size (handles.spikes.waveforms,1) > 1000)
    idxs = randperm (size (handles.spikes.waveforms,1));
    unit_sm = handles.spikes.waveforms (idxs(1:1000), :);
else
    unit_sm=handles.spikes.waveforms;
end

%figure; hold on; plot(unit_sm','b.-'); axis tight; title('Centered Data Subset');
%plot (mean (handles.spikes.waveforms), 'y', 'Linewidth', 2);
%text (0.8*size(handles.spikes.waveforms,2),
%0.8*max(max(handles.spikes.waveforms)), ['#Sp: '
%num2str(size(handles.spikes.waveforms,1))]);
%save([handles.spkfile '-dejitter'], 'handles.spikes');

guidata(hObject, handles);

%---------------------------------------------------------------------
% N.) Splice Waves in Out-dated
%---------------------------------------------------------------------


% % --- Executes on button press in Splice_Waves.
% function Splice_Waves_Callback(hObject, eventdata, handles)
% % hObject    handle to Splice_Waves (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% Tlines={'Channel No - Enter channel nos with ; separator'};
% defent={'1','2','3','4'}; % default entries
% infonames= {'Channels'};
% info = inputdlg(Tlines, 'Choose channels', 1, defent);
% 
% if ~isempty(info)              %see if user hit cancel
%     info = cell2struct(info,infonames);
%     channelnum = str2num(info.Channels);   %convert string to number
%     handles.show = channelnum;
%     handles.spikes.oriwaveforms=handles.spikes.waveforms;
%     spikes=handles.spikes; spikes.waveforms=[];
%     for ch=1:length(channelnum)
%         cmd=sprintf('swaves = spikes.waveforms_ch%d;',ch); eval(cmd);
%         spikes.waveforms=[spikes.waveforms swaves];
%     end
%     handles.spikes=spikes;
% end
% 
% handles.draw.dejitter = handles.spikes.waveforms;
% disp('User chose following channels')
% handles.show
% guidata(hObject, handles);


%---------------------------------------------------------------------
% 3.) Remove Outliers 
%---------------------------------------------------------------------


% --- Executes on button press in Outliers.
function Outliers_Callback(hObject, eventdata, handles)
% hObject    handle to Outliers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.spikes = ssm_outliers(handles.spikes, handles.parameters.outliers_a);
handles.draw.nooutliers = handles.spikes.waveforms;
handles.draw.outliers = handles.spikes.outliers.waveforms;
handles.spikes.parameters=handles.parameters; handles.spikes.nch=handles.nch;

%save([handles.spkfile '-out'],'handles.spikes');

guidata(hObject, handles);


%---------------------------------------------------------------------
% 4.) KMeans Step
%---------------------------------------------------------------------


% --- Executes on button press in KMeans.
function KMeans_Callback(hObject, eventdata, handles)
% hObject    handle to KMeans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.spikes = sss_kmeans(handles.spikes, handles.parameters.kmeans_options);
guidata(hObject, handles);


%---------------------------------------------------------------------
% 5.) Energy Step
%---------------------------------------------------------------------


% --- Executes on button press in Energy.
function Energy_Callback(hObject, eventdata, handles)
% hObject    handle to Energy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.spikes = ssm_energy_subset(handles.spikes, handles.nch);

%save([handles.spkfile '-autoenergy'],'handles');
guidata(hObject, handles);


%---------------------------------------------------------------------
% 6.) Aggr (Aggregation) Step
%---------------------------------------------------------------------


% --- Executes on button press in Aggr.
function Aggr_Callback(hObject, eventdata, handles)
% hObject    handle to Aggr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

agg.reint_out=handles.parameters.reint_out; agg.tmin = handles.parameters.tmin; agg.tref = handles.parameters.tref; agg.cutoff=handles.parameters.cutoff

if handles.rep==0, handles.aggspikes=handles.spikes; end

tempspikes = ssm_aggregate(handles.aggspikes, agg);
% ALWAYS AGGREGATE THE 1ST SPIKES, WITH ALL THE ENERGY INFO, AND NO
%%% AGGREGATION
handles.rep = handles.rep + 1;
handles.tempspikes = tempspikes;
cmd=sprintf('handles.spikes%d = tempspikes;',handles.rep); eval(cmd);
tempspikes.parameters = handles.parameters;
cmd=sprintf('handles.parameters%d = tempspikes.parameters;',handles.rep); eval(cmd);
handles.spikes =  tempspikes;
%%% SAVE TO HANDLES.SPIKES FOR TEMP PLOTTING, AND IF WANTED, SAVING
handles.show = unique(handles.spikes.hierarchy.assigns);

handles.draw.sorted = handles.spikes.waveforms;
guidata(hObject, handles);





%************************************************************************
%   FINALIZE SORT
%************************************************************************


% --- Executes on button press in Load_Spikes.
function Load_Spikes_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Spikes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[SortFile,SortPath] = uigetfile('*sort*.mat','Select the Sorted Mat-file');
if isequal(SortFile,0)
    disp('User selected Cancel')
else
    disp(['User selected', fullfile(SortPath, SortFile)])
    x = load(SortFile,'spikes');
    handles.spikes = x.spikes;
    handles.sortfile=SortFile;
    handles.spkfile = SortFile;
end

if exist('handles.spikes.draw')
    handles.draw=handles.spikes.draw;
end

disp('Current clusters')
unique(handles.spikes.hierarchy.assigns)
handles.show = unique(handles.spikes.hierarchy.assigns);

clear x;

guidata(hObject, handles);

% --- Executes on button press in Raster_ISI.
function Raster_ISI_Callback(hObject, eventdata, handles)
% hObject    handle to Raster_ISI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

spikes = handles.spikes;

x_isi=20; tmin=0.001;
arfp = tmin*1000;
skip = [1:4:41]; skipp= skip+1;
clr = ['y' 'm' 'b' 'c' 'r' 'y' 'm' 'b' 'c' 'r' 'y' 'm' 'b' 'c' 'r'];

sss_rawrastersisi (spikes.waveforms, spikes.fstimes, spikes.hierarchy.assigns, clr, arfp, skip, skipp, spikes.stimonset, x_isi, 7, spikes.sweepd, spikes.nsweeps, handles.show);

% --- Executes on button press in Plot_Type.
function Plot_Type_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_Type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tlines={'Plot Type'};
defent={'1'}; % default entries
infonames= {'plot_type'};
info = inputdlg(Tlines,'Plot Type: 0 for normal, 1 for density plot', 1, defent); 

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    handles.plot_type = str2num(info.plot_type);   %convert string to number
end

guidata(hObject, handles);

% --- Executes on button press in plot_ampl.
function plot_ampl_Callback(hObject, eventdata, handles)
% hObject    handle to plot_ampl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tlines={'Choose Max AMPL for spike waveform plots'};
defent={'400'}; % default entries
infonames= {'plot_ampl'};
info = inputdlg(Tlines, 'Choose PLOT_AMPL', 1, defent); 

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    handles.plot_ampl = str2num(info.plot_ampl);   %convert string to number
end

guidata(hObject, handles);

% --- Executes on button press in Plot_SD.
function Plot_SD_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_SD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tlines={'Plot SD?'};
defent={num2str(0)}; % default entries
infonames= {'plot_sd'};
info = inputdlg(Tlines,'Plot SD? 0 for no, 1 for yes', 1, defent); 

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    handles.plot_sd = str2num(info.plot_sd);   %convert string to number
end

guidata(hObject, handles);


% --- Executes on button press in Plot_All.
function Plot_All_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_All (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

spikes = handles.spikes;
clusters = handles.spikes.hierarchy.assigns;
unqcs =  unique(clusters);

sss_sdss2arcplots(spikes.waveforms, spikes.fstimes ,clusters, unqcs, spikes.nsweeps, spikes.sweepd,...
    spikes.stimonset, spikes.window, handles.arfp, handles.plot_ampl, handles.plot_sd, handles.plot_type);

%disp(['Clusters: ' num2str(unique(handles.spikes.hierarchy.assigns))]);
disp('Current clusters')
unique(handles.spikes.hierarchy.assigns)
%guidata(hObject, handles);


% --- Executes on button press in choose_clu.
function Choose_Clu_Callback(hObject, eventdata, handles)
% hObject    handle to choose_clu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

unique(handles.spikes.hierarchy.assigns)
Tlines={'Clusters - Enter cluster nos with ; separator'};
defent={'', '','',''}; % default entries
infonames= {'Clusters'};
info = inputdlg(Tlines, 'Choose clusters (Current cluster nos. in main window)', 1, defent); 

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    clustnum = str2num(info.Clusters);   %convert string to number
    handles.show = clustnum;
end

disp('User chose following clusters')
handles.show
guidata(hObject, handles);



% --- Executes on button press in Agg_Tree.
function Agg_Tree_Callback(hObject, eventdata, handles)
% hObject    handle to Agg_Tree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.hdraw(11)=figure; aggtreen(handles.spikes); title('Aggregation Tree');
guidata(hObject, handles);


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.hdraw(find(handles.hdraw==0))=[];
close(handles.hdraw);
guidata(hObject, handles);





% --- Executes on button press in igorstim_TC.
function igorstim_TC_Callback(hObject, eventdata, handles)
% hObject    handle to igorstim_TC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

spikes = handles.spikes;
show_clu =  handles.show;

sss_igorTC(spikes, show_clu, handles.plot_ampl, handles.arfp);
guidata(hObject, handles);

% --- Executes on button press in stim_TC.
function stim_TC_Callback(hObject, eventdata, handles)
% hObject    handle to stim_TC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

spikes = handles.spikes;
show_clu =  handles.show;

sss_stimTC(spikes, show_clu, handles.plot_ampl, handles.arfp);
guidata(hObject, handles);

% --- Executes on button press in Plot_Sel_Clu.
function Plot_Sel_Clu_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_Sel_Clu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

spikes = handles.spikes;
clusters = handles.spikes.hierarchy.assigns;
unqcs =  handles.show;

sss_sdss2arcplots(spikes.waveforms, spikes.fstimes ,clusters, unqcs, spikes.nsweeps, spikes.sweepd,...
    spikes.stimonset, spikes.window, handles.arfp,handles.plot_ampl, handles.plot_sd, handles.plot_type);
disp('Current clusters')
unique(handles.spikes.hierarchy.assigns)

% --- Executes on button press in Just_Waves.
function Just_Waves_Callback(hObject, eventdata, handles)
% hObject    handle to Just_Waves (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

spikes = handles.spikes;
clusters = handles.spikes.hierarchy.assigns;
unqcs =  handles.show;
gui_waveplot(spikes,clusters,unqcs,handles.arfp, handles.plot_ampl, handles.plot_sd, handles.plot_type);



% --- Executes on button press in Compare.
function Compare_Callback(hObject, eventdata, handles)
% hObject    handle to Compare (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tlines={'Clusters - Please enter only cluster numbers with ; in between', ['Cluster ' 'color code= combine(0) or no combine(1)']};
defent={'', '1'}; % default entries
infonames= {'Clusters', 'Coloropt'};
info = inputdlg(Tlines, 'Which clusters you want to plot overlaid', 1, defent); 

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    clustnum = str2num(info.Clusters);   %convert string to number
    Coloropt= str2num(info.Coloropt);
    sss_compare (handles.spikes,clustnum, handles.arfp,Coloropt, handles.parameters);
end

%disp(['Clusters: ' num2str(unique(handles.spikes.hierarchy.assigns))]);
disp('Current clusters')
unique(handles.spikes.hierarchy.assigns)
guidata(hObject, handles);


% --- Executes on button press in compare_part.
function compare_part_Callback(hObject, eventdata, handles)
% hObject    handle to compare_part (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


Tlines={'Clusters - Please enter only cluster numbers with ; in between', ['Cluster ' 'color code= combine(0) or no combine(1)']};
defent={'', '1'}; % default entries
infonames= {'Clusters', 'Coloropt'};
info = inputdlg(Tlines, 'Which clusters you want to plot overlaid', 1, defent); 

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    clustnum = str2num(info.Clusters);   %convert string to number
    Coloropt= str2num(info.Coloropt);
    sss_compare_part (handles.spikes,clustnum, handles.arfp,Coloropt, handles.parameters);
end

%disp(['Clusters: ' num2str(unique(handles.spikes.hierarchy.assigns))]);
disp('Current clusters')
unique(handles.spikes.hierarchy.assigns)
guidata(hObject, handles);



% --- Executes on button press in plot_isiStats.
function plot_isiStats_Callback(hObject, eventdata, handles)
% hObject    handle to plot_isiStats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tlines={'Clusters - Please enter TWO cluster numbers with ; in between'};
defent={'', '1'}; % default entries
infonames= {'Clusters'};
info = inputdlg(Tlines, 'Which clusters do you want to test?', 1, defent); 

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    clustnum = str2num(info.Clusters);   %convert string to number
    plot_isiStats(handles.spikes,clustnum , handles.parameters);
end

%disp(['Clusters: ' num2str(unique(handles.spikes.hierarchy.assigns))]);
disp('Current clusters')
unique(handles.spikes.hierarchy.assigns)
guidata(hObject, handles);




% --- Executes on button press in Combine.
function Combine_Callback(hObject, eventdata, handles)
% hObject    handle to Combine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


Tlines={['Please enter Cluster Numbers with ; in between']};
combinecell = inputdlg(Tlines, 'Enter Clusters to Combine', 1, {''}); 
if ~isempty(combinecell)              %see if user hit cancel
    combinestr = cell2struct(combinecell,'combine');
    combine = str2num(combinestr.combine);   %convert string to number
else 
    combine=[];
end

if size(combine,1)~=0
    for lpnc=2:length(combine)
        handles.spikes = merge_clusters(handles.spikes, combine(1), combine(lpnc));
    end
   
end

%disp(['Clusters: ' num2str(unique(handles.spikes.hierarchy.assigns))]);
handles.show = unique(handles.spikes.hierarchy.assigns);
disp('Current clusters')
unique(handles.spikes.hierarchy.assigns)
guidata(hObject, handles);


% --- Executes on button press in Split.
function Split_Callback(hObject, eventdata, handles)
% hObject    handle to Split (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


Tlines={['Please enter Only 1 Cluster Number']};
splitcell = inputdlg(Tlines, 'Enter Cluster to Split', 1, {''}); 
if ~isempty(splitcell)              %see if user hit cancel
    split = splitcell{1};
    split=str2num(split);
else 
    split=[];
end
if size(split,1)~=0
    spikes = sss_split_cluster(handles.spikes, split);
    handles.spikes = spikes;
end

%disp(['Clusters: ' num2str(unique(handles.spikes.hierarchy.assigns))]);
handles.show = unique(handles.spikes.hierarchy.assigns);
disp('Current clusters')
unique(handles.spikes.hierarchy.assigns)
guidata(hObject, handles);


% --- Executes on button press in List_Energy.
function List_Energy_Callback(hObject, eventdata, handles)
% hObject    handle to List_Energy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sss_list_intenergy (handles.spikes, handles.show, handles.parameters);
guidata(hObject, handles);


% --- Executes on button press in Update_energy.
function Update_energy_Callback(hObject, eventdata, handles)
% hObject    handle to Update_energy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

spikes=sss_update_intenergy (handles.spikes);
%handles.spikes.hierarchy.interface_energyn=new_intenergy;
handles.spikes = spikes;
disp('Update Energy: Done')
guidata(hObject, handles);


% --- Executes on button press in add_to_outliers.
function add_to_outliers_Callback(hObject, eventdata, handles)
% hObject    handle to add_to_outliers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tlines={['Please enter Only 1 Cluster Number']};
splitcell = inputdlg(Tlines, 'Enter Cluster to Add', 1, {''}); 
if ~isempty(splitcell)              %see if user hit cancel
    out = splitcell{1};
    out=str2num(out);
    handles.spikes.hierarchy.assigns(find(handles.spikes.hierarchy.assigns==out))=0;
else
    disp('User hit Cancel')
end

handles.show = unique(handles.spikes.hierarchy.assigns);
disp('Current clusters')
unique(handles.spikes.hierarchy.assigns)
guidata(hObject, handles);


% --- Executes on button press in rename_clu.
function rename_clu_Callback(hObject, eventdata, handles)
% hObject    handle to rename_clu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


Tlines={'Original Cluster Name', 'Destination CLuster Name'};
defent={'',''}; % default entries
infonames= {'sCluster', 'dCluster'};
info = inputdlg(Tlines, 'Enter Soure and Destination Clusters', 1, defent); 

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    sCluster = str2num(info.sCluster);   %convert string to number
    dCluster = str2num(info.dCluster);
    handles.spikes.hierarchy.assigns(find(handles.spikes.hierarchy.assigns==sCluster))=dCluster;
else
    disp('User Selcted Cancel')
end

%disp(['Clusters: ' num2str(unique(handles.spikes.hierarchy.assigns))]);
disp('Current clusters')
unique(handles.spikes.hierarchy.assigns)
guidata(hObject, handles);


% --- Executes on button press in Remove_Clu.
function Remove_Clu_Callback(hObject, eventdata, handles)
% hObject    handle to Remove_Clu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tlines={'Cluster Name to remove'};
defent={''}; % default entries
infonames= {'rCluster'};
info = inputdlg(Tlines, 'Enter Cluster Name to Remove', 1, defent); 

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    rCluster = str2num(info.rCluster);   %convert string to number
    spikes = sss_split_cluster(handles.spikes, rCluster,handles.nch);
    handles.spikes = spikes;
else
    disp('User Selcted Cancel')
end

%disp(['Clusters: ' num2str(unique(handles.spikes.hierarchy.assigns))]);
disp('Current clusters')
unique(handles.spikes.hierarchy.assigns)
guidata(hObject, handles);




% --- Executes on button press in save_spktime.
function save_spktime_Callback(hObject, eventdata, handles)
% hObject    handle to save_spktime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tlines={['Enter Cluster Numbers for Spiketimes']};
splitcell = inputdlg(Tlines, 'Enter Cluster to save', 1, {''}); 
if ~isempty(splitcell)              %see if user hit cancel
    clusters = splitcell{1};
    clusters=str2num(clusters);
    
    for clnum=1:length(clusters)
        spktimes=handles.spikes.fstimes(find(handles.spikes.hierarchy.assigns==clusters(clnum)));
        spiketimes{clnum}=spktimes;
        clufilename = ['Cluster' num2str(clusters(clnum)) '.txt'];
        save(clufilename, 'spktimes', '-ascii');
        clear spktimes
    end   
    
    name=strtok(handles.sortfile,'.');
    filename=[name '_spktimes'];
    save(filename, 'spiketimes');
    
    disp('Saved Clusters')
    clusters
    
else
    disp('User hit Cancel')
end

guidata(hObject, handles);


% --- Executes on button press in save_waveform.
function save_waveform_Callback(hObject, eventdata, handles)
% hObject    handle to save_waveform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

unique(handles.spikes.hierarchy.assigns)
Tlines={'Clusters - Enter cluster nos with ; separator'};
defent={'', '','',''}; % default entries
infonames= {'Clusters'};
info = inputdlg(Tlines, 'Choose clusters (Current cluster nos. in main window)', 1, defent); 

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    clustnum = str2num(info.Clusters);   %convert string to number
    gui_meanwaves(handles.spikes,clustnum,handles.sortfile)
end

disp('User chose following clusters')
clustnum
guidata(hObject, handles);


% --- Executes on button press in Save_Mclust.
function Save_Mclust_Callback(hObject, eventdata, handles)
% hObject    handle to Save_Mclust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tlines={['Enter Number of Spks']};
splitcell = inputdlg(Tlines, 'Enter Number of Spks to save', 1, {'20000'}); 
if ~isempty(splitcell)              %see if user hit cancel
    nspks = splitcell{1};
    nspks=str2num(nspks);
    %sss_savespksforMclust(handles.spikes,nspks, handles.nch)
    disp('Saved TT_tetNF for MClust')
else
    disp('User hit Cancel')
end

guidata(hObject, handles);


% --- Executes on button press in Saveclu.
function Saveclu_Callback(hObject, eventdata, handles)
% hObject    handle to Saveclu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clufile = [handles.sortfile '.clu.' handles.elec];

assignments = handles.spikes.hierarchy.assigns;
nclu = length(unique(assignments)); % Including 0
fid = fopen(clufile, 'wt');
fprintf(fid, '%3.0f\n', nclu);
fprintf(fid, '%3.0f\n', assignments);
fclose(fid);

disp('Saved .clu file');
guidata(hObject, handles);


% --- Executes on button press in save_autocorr.
function save_autocorr_Callback(hObject, eventdata, handles)
% hObject    handle to save_autocorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tlines={['Enter Cluster Numbers to Save in File']};
splitcell = inputdlg(Tlines, 'Enter Clusters to save', 1, {''}); 
if ~isempty(splitcell)              %see if user hit cancel
    clusters = splitcell{1};
    clusters=str2num(clusters);
    sss_savespksforautocorr(handles.spikes, clusters);
    disp('Saved Sort_autocorr for autocorr')
else
    disp('User hit Cancel')
end


guidata(hObject, handles);


% --- Executes on button press in save_finalize.
function save_finalize_Callback(hObject, eventdata, handles)
% hObject    handle to save_finalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

string=strtok(handles.spkfile,'.');
%filenamef = [string '_sortf'];
filenamef = [string];

[SortfFile,SortfPath] = uiputfile('*sortf*.mat','Save the final Sortf Mat-file',[filenamef]);
if isequal(SortfFile,0)
    disp('User selected Cancel')
else
    disp(['User selected', fullfile(SortfPath, SortfFile)])
    handles.sortffile = SortfFile;
    spikes = handles.spikes;
    save (handles.sortffile, 'spikes');

end

%disp(['Clusters: ' num2str(unique(handles.spikes.hierarchy.assigns))]);
disp('Saved Current clusters')
unique(handles.spikes.hierarchy.assigns)
guidata(hObject, handles);


% --- Executes on button press in Cl_Clu.
function Cl_Clu_Callback(hObject, eventdata, handles)
% hObject    handle to Cl_Clu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

string=strtok(handles.spkfile,'.');
%filenamef = [string '_sortf'];
filenamef = [string '-cl'];

[SortfFile,SortfPath] = uiputfile('MIC*.mat','Save the final Sortf Mat-file',[filenamef]);
if isequal(SortfFile,0)
    disp('User selected Cancel')
else
    disp(['User selected', fullfile(SortfPath, SortfFile)])
    handles.sortffile = SortfFile;
    spikes = sss_cleanclustersfinal(handles.spikes, handles.nch)    
    save (handles.sortffile, 'spikes');
    handles.spikes = spikes;

end

%disp(['Clusters: ' num2str(unique(handles.spikes.hierarchy.assigns))]);
disp('Saved Current clusters')
unique(handles.spikes.hierarchy.assigns)
guidata(hObject, handles);






%***************************************************************
% Algo Figs Panel - Plots Intermediate Steps of Algorithm
%***************************************************************


%------------------------------------------------------------------
% 1. De-Jitterred Waveforms
%------------------------------------------------------------------


% --- Executes on button press in De_Jitter.
function De_Jitter_Callback(hObject, eventdata, handles)
% hObject    handle to De_Jitter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.hdraw(4)=figure; colormap hot;
if (size (handles.spikes.waveforms,1) > 5000)
    idxs = randperm (size (handles.spikes.waveforms,1));
    unit_sm = handles.spikes.waveforms (idxs(1:5000), :);
else
    unit_sm=handles.spikes.waveforms;
end
plot(unit_sm','b.-'); hold on;
plot (mean (handles.spikes.waveforms), 'y', 'Linewidth', 2); axis tight; title('Centered Data - Subset 5000');
text (0.8*size(handles.spikes.waveforms,2), 0.8*max(max(unit_sm)), ['#Sp: ' num2str(length(handles.spikes.spiketimes))]);

guidata(hObject, handles);


%------------------------------------------------------------------
% 2. Waveforms with outliers removed
%------------------------------------------------------------------


% --- Executes on button press in No_Outliers.
function No_Outliers_Callback(hObject, eventdata, handles)
% hObject    handle to No_Outliers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.hdraw(5)=figure;
nooutidx = find(handles.spikes.hierarchy.assigns~=0);
plot(handles.spikes.waveforms(nooutidx,:)','b'); axis tight; title('Centered Data w/ Outliers Removed');
text (0.8*size(handles.spikes.waveforms(nooutidx,:),2), 0.8*max(max(handles.spikes.waveforms(nooutidx,:))), ['#Sp: ' num2str(length(nooutidx))]);

guidata(hObject, handles);


%------------------------------------------------------------------
% 3. Outlier Waveforms 
%------------------------------------------------------------------


% --- Executes on button press in Outliers_Only.
function Outliers_Only_Callback(hObject, eventdata, handles)
% hObject    handle to Outliers_Only (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles.spikes,'outliers')
    handles.hdraw(6)=figure; plot(handles.spikes.outliers.waveforms','b'); axis tight; title(' Outliers ');
    text (0.8*size(handles.spikes.draw.outliers,2), 0.8*max(max(handles.spikes.draw.nooutliers)), ['#Sp: ' num2str(size(handles.spikes.draw.outliers,1))]);
else
    disp('Cant find Outliers Field - Showing Cluster 0');
    outidx = find(handles.spikes.hierarchy.assigns==0);
    handles.hdraw(6)=figure; plot(handles.spikes.waveforms(outidx,:)'); axis tight; title(' Outliers-Cluster 0 ');
    text (0.8*size(handles.spikes.waveforms(outidx,:),2), 0.8*max(max(handles.spikes.waveforms(outidx,:))), ['#Sp: ' num2str(length(outidx))]);
end

guidata(hObject, handles);


%------------------------------------------------------------------
% 4. Waveforms when Over-Clustered
%------------------------------------------------------------------


% --- Executes on button press in Over_Clusters.
function Over_Clusters_Callback(hObject, eventdata, handles)
% hObject    handle to Over_Clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.hdraw(7)=figure;  set(gcf, 'Renderer', 'OpenGL');
clustshow(handles.spikes, handles.spikes.overcluster.assigns, [1:32],handles.nch); subplot(2,1,1); title('Local Clusters');

guidata(hObject, handles);



%------------------------------------------------------------------
% 5. Raw Waveforms - Not Used Anymore. Use Noise Menu for Plotting
%--------------------------------------------------------------


% --- Executes on button press in Raw_Waveforms.
function Raw_Waveforms_Callback(hObject, eventdata, handles)
% hObject    handle to Raw_Waveforms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.hdraw(1)=figure; colormap hot;
if isfield(handles.spikes,'oriwaveforms')
    subplot(2,1,1); 
    if (size (handles.spikes.oriwaveforms,1) > 1000)
        idxs = randperm (size (handles.spikes.oriwaveforms,1));
        unit_sm = handles.spikes.oriwaveforms (idxs(1:1000), :);
    else
        unit_sm=handles.spikes.oriwaveforms;
    end
    plot(unit_sm','b.-'); hold on;
    axis tight; title('Raw Data');
    text (0.8*size(handles.spikes.oriwaveforms,2), 0.8*max(max(handles.spikes.oriwaveforms)), ['#Sp: ' num2str(length(handles.spikes.spiketimes))]);
    title(['Ori Waveforms: Subset']);
    subplot(2,1,2); hist2d(handles.spikes.oriwaveforms); h = gca;
    text (0.8*size(handles.spikes.oriwaveforms,2), 0.8*max(max(handles.spikes.oriwaveforms)), ['#Sp: ' num2str(length(handles.spikes.spiketimes))]);
   
else
    subplot(2,1,1); 
    if (size (handles.spikes.waveforms,1) > 1000)
        idxs = randperm (size (handles.spikes.waveforms,1));
        unit_sm = handles.spikes.waveforms (idxs(1:1000), :);
    else
        unit_sm=handles.spikes.waveforms;
    end
    plot(unit_sm','b.-'); hold on;
    plot (mean (handles.spikes.waveforms), 'y', 'Linewidth', 2);
    text (0.8*size(handles.spikes.waveforms,2), 0.8*max(max(unit_sm)), ['#Sp: ' num2str(length(handles.spikes.spiketimes))]);
    title(['Dejittered Waveforms: Subset']);
    subplot(2,1,2); hist2d(unit_sm); 
    text (0.8*size(unit_sm,2), 0.8*max(max(unit_sm)), ['#Sp: ' num2str(length(handles.spikes.spiketimes))]);
end


guidata(hObject, handles);




%*******************************************************************
% Cluster Prop Panel
%*******************************************************************


%-----------------------------------------------------------------
%1. AmplPl - Amplitude Plots for CLustered Data
%-----------------------------------------------------------------


% --- Executes on button press in ampplot.
function ampplot_Callback(hObject, eventdata, handles)
% hObject    handle to ampplot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


Tlines={'Plot PCA? 1 for yes, 0 for no'};
defent={'1'}; % default entries
infonames= {'dopca'};
info = inputdlg(Tlines, 'Plot PCA?', 1, defent);

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    dopca = str2num(info.dopca);   %convert string to number
else
    dopca=0;
end

sss_ampplot(handles.spikes,handles.spikes.hierarchy.assigns, handles.show,handles.nch, dopca);    
guidata(hObject, handles);



%------------------------------------------------------------------
% 2. Amplitude Drift With Time for chosen clusters
%------------------------------------------------------------------


% --- Executes on button press in Drift_draw.
function Drift_draw_Callback(hObject, eventdata, handles)
% hObject    handle to Drift_draw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% If multiple channles, get channel to plot %%
if handles.nch>1
    Tlines={'Which channel?'};
    defent={'1'}; % default entries
    infonames= {'Channel'};
    info = inputdlg(Tlines, 'use_channel', 1, defent);
    
    if ~isempty(info)              %see if user hit cancel
        info = cell2struct(info,infonames);
        usech = str2num(info.Channel);   %convert string to number
        %handles.nch=nch;
    end
    disp('Using channel')
    usech
else
    usech = 1;
end

spikes = handles.spikes;
%if isempty(handles.spikes), spikes = handles.spikes; end
sss_drift(spikes,spikes.hierarchy.assigns, handles.show,handles.nch, usech)

guidata(hObject, handles);


%------------------------------------------------------------------
% 3. 2D AMP - 2-d Amplitude Plots 
%------------------------------------------------------------------


% --- Executes on button press in amp_2d.
function amp_2d_Callback(hObject, eventdata, handles)
% hObject    handle to amp_2d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sss_amptest(handles.spikes,handles.spikes.hierarchy.assigns, handles.show, 'xy',handles.nch);
guidata(hObject, handles);


%------------------------------------------------------------------
% 4. 3D AMP - 3-d Amplitude Plots 
%------------------------------------------------------------------


% --- Executes on button press in amp_3d.
function amp_3d_Callback(hObject, eventdata, handles)
% hObject    handle to amp_3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


sss_amptest(handles.spikes,handles.spikes.hierarchy.assigns, handles.show, 'xyz',handles.nch);
guidata(hObject, handles);



%------------------------------------------------------------------
% 5. 2D PCA - 2-d PCA Plots 
%------------------------------------------------------------------

% --- Executes on button press in Final_2D.
function Final_2D_Callback(hObject, eventdata, handles)
% hObject    handle to Final_2D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.assigns =[]; handles.z= 0; %handles.pc =1;
guidata(hObject, handles);

%ssgtest(handles.spikes,handles.assigns, handles.z);
ssgtest(handles.spikes,handles.spikes.hierarchy.assigns, handles.show, 'xy',handles.nch);



%------------------------------------------------------------------
% 6. 3D PCA - 3-d PCA Plots 
%------------------------------------------------------------------


% --- Executes on button press in Final_3D.
function Final_3D_Callback(hObject, eventdata, handles)
% hObject    handle to Final_3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.assigns =[]; handles.z = 1; %handles.pc =1;
guidata(hObject, handles);

%ssgtest(handles.spikes,handles.assigns, handles.z);
ssgtest(handles.spikes,handles.spikes.hierarchy.assigns, handles.show, 'xyz',handles.nch);


%------------------------------------------------------------------
% 7. Clust Show - Overlaid Waveforms and Width-Height Plots using
% thresholded peaks
%------------------------------------------------------------------



% --- Executes on button press in Clust_Show.
function Clust_Show_Callback(hObject, eventdata, handles)
% hObject    handle to Clust_Show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.hdraw(9)=figure;clustshow(handles.spikes, handles.spikes.hierarchy.assigns, handles.show, handles.nch);

guidata(hObject, handles);


%------------------------------------------------------------------
% 8. WavePlot - Waveform Plots + ISI 
%------------------------------------------------------------------

% --- Executes on button press in wave_plot.
function wave_plot_Callback(hObject, eventdata, handles)
% hObject    handle to wave_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.hdraw(10)=figure; redimscreen_halfvert(0); waveplotn(handles.spikes, handles.spikes.hierarchy.assigns, handles.show);
guidata(hObject, handles);


%------------------------------------------------------------------
% 9. DenPlot - Waveform Density Plots + ISI + Log(ISI)
%------------------------------------------------------------------



% --- Executes on button press in Den_Plot.
function Den_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to Den_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.hdraw(10)=figure; sss_densityplot(handles.spikes, handles.spikes.hierarchy.assigns, handles.show, handles.plot_ampl);
guidata(hObject, handles);


%------------------------------------------------------------------
% 9. FinPlot - Waveform Density Plots + ISI + AutoCorr + Log ISI 
%------------------------------------------------------------------

% --- Executes on button press in Fin_Plot.
function Fin_Plot_Callback(hObject, eventdata, handles)
% hObject    handle to Fin_Plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.hdraw(10)=figure; 
sss_finalplot(handles.spikes, handles.spikes.hierarchy.assigns, handles.show, handles.plot_ampl);
guidata(hObject, handles);






















% --- Executes on button press in DynParamPlot.
function DynParamPlot_Callback(hObject, eventdata, handles)
% hObject    handle to DynParamPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sss_matclusttest(handles.spikes,handles.spikes.hierarchy.assigns, handles.show, 'xy',handles.nch);
guidata(hObject, handles);


% --- Executes on button press in ISI_Viol.
function ISI_Viol_Callback(hObject, eventdata, handles)
% hObject    handle to ISI_Viol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tlines={'Please enter Only 1 Cluster Number','ISI'};
defent={'','1'};
infonames={'Cluster','ISI'};
info = inputdlg(Tlines, 'Enter Cluster No', 1, defent); 
if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    Cluster = str2num(info.Cluster);
    ISI = str2num(info.ISI); 
    handles.hdraw(12)=figure; 
    handles.spikes = sss_isi_viol(handles.spikes,handles.spikes.hierarchy.assigns, Cluster, ISI, handles.plot_ampl);
else
    disp('User hit Cancel')
end

guidata(hObject, handles);


% --- Executes on button press in TwoParamPlot.
function TwoParamPlot_Callback(hObject, eventdata, handles)
% hObject    handle to TwoParamPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


Tlines={'Enter Parameter Idxs (1=Time, 2=Amp1, 3=Amp3, 4=Amp4, 5=Amp5, 6=HtCh, 7=WidCh, 8=Xpos, 9=Ypos)',['Param2']};
defent={'2', '3'}; % default entries
infonames= {'Param1', 'Param2'};
info = inputdlg(Tlines, 'Which parameters do you want to plot', 1, defent); 

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    Param1 = str2num(info.Param1);   %convert string to number
    Param2 = str2num(info.Param2);
    sss_twoparamplot(handles.spikes, Param1, Param2, handles.spikes.hierarchy.assigns, handles.show, handles.nch)
end

guidata(hObject, handles);



% --- Executes on button press in CorrPlot.
function CorrPlot_Callback(hObject, eventdata, handles)
% hObject    handle to CorrPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.hdraw(13)=figure; 
sss_corrplot(handles.spikes, handles.spikes.hierarchy.assigns, handles.show, handles.plot_ampl, handles.nch);
guidata(hObject, handles);




% --- Executes on button press in SaveGclust.
function SaveGclust_Callback(hObject, eventdata, handles)
% hObject    handle to SaveGclust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Tlines={'Enter Day','Enter Date','Enter Tet No'};
defent={'01', '032510','5'}; % default entries
infonames= {'Day', 'Date', 'Tet'};
info = inputdlg(Tlines, 'Enter Info to save', 1, defent); 

if ~isempty(info)              %see if user hit cancel
    info = cell2struct(info,infonames);
    day = info.Day;
    date = info.Date;
    tetno = info.Tet;
    sss_save_gclust(handles.spikes, tetno, day, date);
    disp('Saved gclust file')
else
    disp('User hit Cancel')
end

guidata(hObject, handles);







%***********************************************************************
% BATCH MODE PANEL
%***********************************************************************


%----------------------------------------------------------------------
% 1. BATCH MODE
%----------------------------------------------------------------------

% --- Executes on button press in Batch_Mode.
function Batch_Mode_Callback(hObject, eventdata, handles)
% hObject    handle to Batch_Mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[sort_cell]=uigetfile('*.mat','Select the MatFiles with spikes/Read Spikes','MultiSelect','on');

if ~isequal(sort_cell,0)

    button = questdlg('Really? All Those Files?','Batch-Mode');
    if strcmp(button,'Yes')

        for files = 1:size(sort_cell,2)
            handles.spkfile = sort_cell{files};
            
            if handles.nch ==1
                [tempspikes,draw] = sss_spksort(handles.spkfile, handles.parameters);
            else
                [tempspikes,draw] = ssm_tetsort(handles.spkfile, handles.parameters, handles.nch, handles.Spiketype, handles.elec, handles.Subset);
            end

            handles.rep = handles.rep + 1;
            handles.draw = draw; handles.tempspikes = tempspikes;
            cmd=sprintf('handles.spikes%d = tempspikes;',handles.rep); eval(cmd);
            tempspikes.parameters = handles.parameters;
            cmd=sprintf('handles.parameters%d = tempspikes.parameters;',handles.rep); eval(cmd);
            handles.spikes =  tempspikes;

            set(handles.guifighandle, 'HandleVisibility', 'Off');
            close all;

            string=strtok(handles.spkfile,'.');
            filename1 = [string '_sort'];
            spikes=handles.spikes;

            handles.Savefile = filename1;
            save (handles.Savefile, 'spikes','handles');
        end

    end
end

guidata(hObject, handles);






