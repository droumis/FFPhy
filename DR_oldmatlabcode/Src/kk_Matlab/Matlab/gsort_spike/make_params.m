function varargout = make_params(varargin)
% make_PARAMS M-file for make_params.fig
%      make_PARAMS, by itself, creates a new make_PARAMS or raises the existing
%      singleton*.
%
%      H = make_PARAMS returns the handle to a new make_PARAMS or the handle to
%      the existing singleton*.
%
%      make_PARAMS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in make_PARAMS.M with the given input arguments.
%
%      make_PARAMS('Property','Value',...) creates a new make_PARAMS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before make_params_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to make_params_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  make "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help make_params

% Last Modified by GUIDE v2.5 18-May-2005 21:29:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @make_params_OpeningFcn, ...
                   'gui_OutputFcn',  @make_params_OutputFcn, ...
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


% --- Executes just before make_params is made visible.
function make_params_OpeningFcn(hmakep, makep_eventdata, handles_makep, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to make_params (see VARARGIN)

% make default command line output for make_params
handles_makep.output = hmakep;

handles_makep.input_handles = varargin{2};  % HANDLES FROM MAIN GUI
handles_makep.input_hObject = varargin{1};

handles_makep.guifighandle=get(0,'CurrentFigure');
% Update handles structure
guidata(hmakep, handles_makep);

% UIWAIT makes make_params wait for user response (see UIRESUME)
uiwait(handles_makep.guifighandle);


% --- Outputs from this function are returned to the command line.
function varargout = make_params_OutputFcn(hmakep, makep_eventdata, handles_makep)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles_makep.output;

close(handles_makep.guifighandle);


%--------------------------------------------------

% INI_CLUSTERS: SELECTING INITIAL NUMBER OF CLUSTERS

% --- Executes during object creation, after setting all properties.
function Ini_Clusters_CreateFcn(hmakep, makep_eventdata, handles_makep)
% hObject    handle to Ini_Clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hmakep,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hmakep,'BackgroundColor','white');
end

% --- Executes on selection change in Ini_Clusters.
function Ini_Clusters_Callback(hmakep, makep_eventdata, handles_makep)
% hObject    handle to Ini_Clusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Ini_Clusters contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Ini_Clusters


% --- Executes during object creation, after setting all properties.
function Outlier_a_CreateFcn(hmakep, makep_eventdata, handles_makep)
% hObject    handle to Outlier_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hmakep,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hmakep,'BackgroundColor','white');
end

% --- Executes on selection change in Outlier_a.
function Outlier_a_Callback(hmakep, makep_eventdata, handles_makep)
% hObject    handle to Outlier_a (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Outlier_a contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Outlier_a



%---------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function Reint_Out_CreateFcn(hmakep, makep_eventdata, handles_makep)
% hObject    handle to Reint_Out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hmakep,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hmakep,'BackgroundColor','white');
end


% --- Executes on selection change in Reint_Out.
function Reint_Out_Callback(hmakep, makep_eventdata, handles_makep)
% hObject    handle to Reint_Out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Reint_Out contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Reint_Out


%---------------------------------------------------------------------

% --- Executes on selection change in tmin.
function tmin_Callback(hmakep, makep_eventdata, handles_makep)
% hObject    handle to tmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns tmin contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tmin


% --- Executes during object creation, after setting all properties.
function tmin_CreateFcn(hmakep, makep_eventdata, handles_makep)
% hObject    handle to tmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hmakep,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hmakep,'BackgroundColor','white');
end

%---------------------------------------------------------------------

% --- Executes on selection change in tref.
function tref_Callback(hmakep, makep_eventdata, handles_makep)
% hObject    handle to tref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns tref contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tref


% --- Executes during object creation, after setting all properties.
function tref_CreateFcn(hmakep, makep_eventdata, handles_makep)
% hObject    handle to tref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hmakep,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hmakep,'BackgroundColor','white');
end

%---------------------------------------------------------------------


% --- Executes on selection change in Cutoff.
function Cutoff_Callback(hmakep, makep_eventdata, handles_makep)
% hObject    handle to Cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Cutoff contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Cutoff


% --- Executes during object creation, after setting all properties.
function Cutoff_CreateFcn(hmakep, makep_eventdata, handles_makep)
% hObject    handle to Cutoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hmakep,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hmakep,'BackgroundColor','white');
end



%*********************************************************************

% --- Executes on button press in Param_Done.
function Param_Done_Callback(hmakep, makep_eventdata, handles_makep)
% hObject    handle to Param_Done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

list1 = get(handles_makep.Ini_Clusters,'String');
index1 = get(handles_makep.Ini_Clusters,'Value');
handles_makep.input_handles.parameters.kmeans_options.divisions = log2(str2num(list1{index1}));
handles_makep.output.parameters.kmeans_options.divisions = log2(str2num(list1{index1}));

list2 = get(handles_makep.Outlier_a,'String');
index2 = get(handles_makep.Outlier_a,'Value');
handles_makep.input_handles.parameters.outliers_a = str2num(list2{index2});
handles_makep.output.parameters.outliers_a = str2num(list2{index2});

list3 = get(handles_makep.Reint_Out,'String');
index3 = get(handles_makep.Reint_Out,'Value');
handles_makep.input_handles.reint_out = str2num(list3{index3});
handles_makep.output.parameters.reint_out = str2num(list3{index3});

list4 = get(handles_makep.tmin,'String');
index4 = get(handles_makep.tmin,'Value');
handles_makep.input_handles.parameters.tmin = str2num(list4{index4});
handles_makep.output.parameters.tmin = str2num(list4{index4});

list5 = get(handles_makep.tref,'String');
index5 = get(handles_makep.tref,'Value');
handles_makep.input_handles.parameters.tref = str2num(list5{index5});
handles_makep.output.parameters.tref = str2num(list5{index5});

list6 = get(handles_makep.Cutoff,'String');
index6 = get(handles_makep.Cutoff,'Value');
handles_makep.input_handles.parameters.cutoff = str2num(list6{index6});
handles_makep.output.parameters.cutoff = str2num(list6{index6});



%guidata(handles_makep.input_hObject, handles_makep.input_handles);
guidata(hmakep, handles_makep);

uiresume(handles_makep.guifighandle);
%



















