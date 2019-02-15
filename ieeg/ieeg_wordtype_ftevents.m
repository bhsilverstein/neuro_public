function [trl] = ieeg_wordtype_ftevents(cfg)

%%% Import events into fieldtrip format.
% send it a cfg.eventfile and cfg.dname with the full path and filename of the event
% file.

% Event codes for Word Type event files.
% 101: Who onset    301: What onset     401: Where onset    501: When onset     601 Filler word onset     701: Filler word onset      
% 102: Who offset   302: What offset    402: Where offset   502: When offset    602 Filler word offset    702: Filler word offset
% 201: Baseline


%% Read in event data
% Open file
fid=fopen([cfg.dname cfg.eventfile]);

% Skip header lines.
% dataformat='%d %d %d %s';
% nextline=textscan(fid,'%d %d %d %s','HeaderLines',1);
tmp = fgetl(fid);
nextline = strsplit(tmp);
chk = str2double(nextline{1});
while isnan(chk) || chk==0
    tmp=fgetl(fid);
    nextline = strsplit(tmp);
    chk = str2double(nextline{1});
end

% Read event data
idx=1;
eventdata = nextline;
while ~isempty(nextline{1})
    idx=idx+1;
%     nextline=textscan(fid,'%d %d %d %s','HeaderLines',1);
    tmp = fgetl(fid);
    if tmp==-1
        break
    end
    nextline = strsplit(tmp);
    if size(nextline{1},1)>1
        for c = 1:4
            eventdata{idx,c}=nextline{c}(1);
        end
        idx=idx+1;
        for c = 1:4
            eventdata{idx,c}=nextline{c}(2);
        end
    else
        for c = 1:4
            eventdata{idx,c}=nextline{c};
        end
    end
end
fclose(fid);

%% Grab event type from third column
if ischar(eventdata{1,3})
    type = cellfun(@str2double,eventdata(:,3));    
else
    type = cell2mat(eventdata(:,3));
end
% Find bad event flags
bad = type==0;
type(bad)=[];
types=unique(type); 
nevent = numel(type);

%% Grab event latencies from first column and convert to samples
if ischar(eventdata{1,1})
    latencies = cellfun(@str2double,eventdata(:,1));
else
    latencies = double(cell2mat(eventdata(:,1)));   
end
latencies = 1E-3 * latencies;

% Drop bad events
latencies(bad)=[];

%% Event windows
eventwins(:,1) = ones(nevent,1) .* -1999;
eventwins(:,2) = ones(nevent,1) .* 1000;
% Different window for baseline events
eventwins(type == 201,1) = -599;
eventwins(type == 201,2) = -200;

%% Offset vector
offset = zeros(nevent,1); stopind = offset; duration=offset;
for t = 1:numel(types)
    if types(t) == 201
        offset(type==types(t),1)= -599;
        stopind(type==types(t),1)= -200;
        duration(type==types(t),1)=400;
    else
        offset(type==types(t),1)= -1999;
        stopind(type==types(t),1)= 1000;
        duration(type==types(t),1)=3000;
    end
    
end       

stopind = stopind + latencies;

%% Store in fieldtrip format
% Type
temp = num2cell(type);
[event(1:length(temp)).type] = temp{:};
% Sample
temp = num2cell(latencies);
[event.sample] = temp{:};
% Offset
temp = num2cell(offset);
[event.offset] = temp{:};
% Duration
temp = num2cell(duration);
[event.duration]=temp{:};

% trl
trl=[latencies + offset stopind offset type];

% Drop trials that still have type=0;
drop = trl(:,4)==0;
trl(drop,:)=[];
event(drop)=[];











































