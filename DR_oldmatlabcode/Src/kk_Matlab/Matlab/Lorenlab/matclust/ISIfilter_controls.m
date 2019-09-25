function ISIfilter_controls(varargin)
%These functions are the gui controls for the filebrowser function 

feval(varargin{:});
%--------------------------------------------
function CloseRequestFnc(handles)

handles.filterindex = [];
guidata(handles.fighandle, handles);
uiresume(handles.fighandle);
%--------------------------------------------
function cancelbutton_Callback(handles)

handles.filterindex = [];
guidata(handles.fighandle, handles);
uiresume(handles.fighandle);
%----------------------------------------------
function okbutton_Callback(handles)
%the user clicked the ok button.  We do some checks before we allow the
%figure to close


%all checks cleared, so we release the hold on the figure
uiresume(handles.fighandle);
%---------------------------------------------
function slider1_Callback(handles)
global clustdata;
global clustattrib;
value = get(handles.slider1,'Value');
sign = get(handles.signbutton,'UserData');
cutoff = 10000*value;
set(handles.text1,'String',[num2str(round(cutoff/10)), ' ms']);
if length(handles.index) > 2
    timediff = diff(clustdata.params(handles.index,1));
    if (sign == 1)
        handles.filterindex = handles.index(find((timediff> cutoff))+1);
    else
        handles.filterindex = handles.index(find((timediff< cutoff))+1);
    end
end


set(handles.p(1),'XData',clustdata.params(handles.filterindex,2));
set(handles.p(1),'YData',clustdata.params(handles.filterindex,3));   

set(handles.p(2),'XData',clustdata.params(handles.filterindex,2));    
set(handles.p(2),'YData',clustdata.params(handles.filterindex,4));  

set(handles.p(3),'XData',clustdata.params(handles.filterindex,2));    
set(handles.p(3),'YData',clustdata.params(handles.filterindex,5));  

set(handles.p(4),'XData',clustdata.params(handles.filterindex,3));    
set(handles.p(4),'YData',clustdata.params(handles.filterindex,4));  

set(handles.p(5),'XData',clustdata.params(handles.filterindex,3));    
set(handles.p(5),'YData',clustdata.params(handles.filterindex,5));  

set(handles.p(6),'XData',clustdata.params(handles.filterindex,4));    
set(handles.p(6),'YData',clustdata.params(handles.filterindex,5));  

guidata(handles.fighandle, handles);

%--------------------------------------------
function signbutton_Callback(handles)

sign = get(handles.signbutton,'UserData');
if (sign == 1)
    set(handles.signbutton,'UserData',2);
    set(handles.signbutton,'String','<');
    newsign = 2;
elseif (sign == 2)
    set(handles.signbutton,'UserData',1);
    set(handles.signbutton,'String','>');
    newsign = 1;
end

slider1_Callback(handles);
    
