
function gui_meanwaves(spikes,clusters,sortfile)
%%%% TO SAVE  MEANS OF WAVEFORMS
%--------------------------------------------------------------------------

assigns=spikes.hierarchy.assigns;
nclu=length(clusters);

for n=1:nclu
    currclu = clusters(n);
    sp_idx=find (assigns==currclu);

    m = mean(spikes.waveforms(sp_idx,:));
    s = std(spikes.waveforms(sp_idx,:));
    wave(n).mean=m;
    wave(n).std=s;
end

savename= strtok(sortfile,'.');
savename=[savename '-WAVE'];

save (savename, 'wave');


% % --- Executes on button press in choose_clu.
% function Choose_Clu_Callback(hObject, eventdata, handles)
% % hObject    handle to choose_clu (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% unique(handles.spikes.hierarchy.assigns)
% Tlines={'Clusters - Enter cluster nos with ; separator'};
% defent={'', '','',''}; % default entries
% infonames= {'Clusters'};
% info = inputdlg(Tlines, 'Choose clusters (Current cluster nos. in main window)', 1, defent); 
% 
% if ~isempty(info)              %see if user hit cancel
%     info = cell2struct(info,infonames);
%     clustnum = str2num(info.Clusters);   %convert string to number
%     handles.show = clustnum;
% end
% 
% disp('User chose following clusters')
% handles.show
% guidata(hObject, handles);


