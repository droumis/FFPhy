
gain=1; maxsignal=5*(10^6);  %in uV
time=5		% FOR TIMER FUNCTION, WHICH WILL STOP AI
ai=analoginput('nidaq',1)
setverify(ai,'InputType','NonReferencedSingleEnded');
addchannel(ai,0:7)

setverify(ai,'LoggingMode','Disk')
setverify(ai,'LogToDiskMode','Index')
setverify(ai,'LogFileName','test')

setverify(ai,'TimerPeriod',time),  ai.TimerFcn = {'stop'};	%setverify(ai,'TimerFcn','timerstop.m')

%%%% UNITS AND GAIN
%scaled value = (A/D value)(units range)/(sensor range)
setverify(ai.Channel,'Units','uV');	%get(ai.Channel,'UnitsRange')
%setverify(ai.Channel,'SensorRange',[-5 5]*1000000/gain)		%in uV

maxsignal_V = maxsignal/(10^6); ip = maxsignal_V*gain;
setverify(ai.Channel,'InputRange',[-ip ip]);   % in V: Any signal above this will be clipped
setverify(ai.Channel,'UnitsRange',[-maxsignal maxsignal]);    % Set units range (/sensor range) to be same, in uV, as of input range (eg [-1000 1000]uV units range for [-1000 1000]uV input range)
setverify(ai.Channel,'SensorRange',[-ip ip])     %% corresponding to units range in V
ai

%%%  (In volts: setverify(ai.Channel,'SensorRange',[-5 5]/gain)   )

setverify(ai,'SampleRate',30000)
%setverify(ai,'SamplesPerTrigger',32000)   % 1 sec of data
setverify(ai,'SamplesPerTrigger',Inf) 	% Continuous function till timer stop is called

%%%%%%%%%% TO ADD %%%%%%%%%%%%%%

%get(ai,'TriggerType')
%get(ai,'ManualTriggerHwOn')
%get(ai,'TriggerFcn')
%get(ai,'StopFcn')
%get(ai,'StartFcn')
%get(ai,'BufferingMode'),  get(ai,'BufferingConfig')          %% (blocksize*n_blocks)*n_channels_bytes/channel
%setverify(ai,'BufferingMode','Manual')
%setverify(ai,'BufferingConfig',[4000 20])


%%%%%%%%%%%%%IMP IMP %%%%%%%%%%%%%%%%%%%%%%%%%%
%daqpropedit(ai)   %GRAPHICALLY SETUP ai
% data = getdata(ai)

%You can replay data saved to disk with the daqread function. Refer to
%Logging Information to Disk for more information about LoggingMode and daqread

%%GETDATA%% GET Data from the engine(memory, not disk file)

% data = getdata(obj)
% data = getdata(obj,samples)
% data = getdata(obj,'type')
% data = getdata(obj,samples,'type')
% [data,time] = getdata(...)
% [data,time,abstime] = getdata(...)
% [data,time,abstime,events] = getdata(...)

%%%%DAQREAD
% daqread Read a Data Acquisition Toolbox (.daq) file Syntaxdata = daqread('file')
% data = daqread('file','PropertyName',PropertyValue,...)
% [data,time] = daqread('test.daq');
% [data,time,abstime] = daqread(...)
% [data,time,abstime,events] = daqread(...)
% [data,time,abstime,events,daqinfo] = daqread(...)
% daqinfo = daqread('file','info')

%%%PEEKDATA%%%
%data = peekdata(obj,samples)
%More About Using peekdataUnlike getdata, peekdata is a nonblocking
%function that immediately returns control to MATLAB. Because peekdata
%does not block execution control, data might be missed or repeated.
%peekdata takes a "snapshot" of the most recent acquired data and does not
%remove samples from the data acquisition engine. Therefore, the 
%SamplesAvailable property value is not affected when peekdata is called
%%%Rules for Using peekdata
%You can call peekdata before a trigger executes. 
%Therefore, peekdata is useful for previewing data before it is
%logged to the engine or to a disk file. In most cases, you will
%call peekdata while the device object is running. However, you
%can call peekdata once after the device object stops running.
%If samples is greater than the number of samples currently acquired,\
%all available samples are returned
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%get(ai.Channel,'InputRange')
%setverify(ai,'TransferMode','SingleDMA')
%daqmem


%out = daqhwinfo(ai);
%out.NativeDataType
%ans =
%int16

%start(ai)
%nativedata = getdata(ai,1000,'native');
%scaling = get(ai.Channel(1),'NativeScaling');
%offset = get(ai.Channel(1),'NativeOffset');
%data = double(nativedata)*scaling + offset;