%% Analyze Ship
%
% Martin Gassmann, 2 May 2016, SIO-MPL/UCSD
%
% * Analyzes a ship's AIS and acoustic data from CINMS_30_B array 
%
% Version:
%
% 2: Written for analyzing 2nd ship 229989000
%    Faster Spectral Averaging
% 3: Finds correct acoustic data files automatically
%    Can read in data segment over two consecutive SoundTrap data files
% 4: Uses ship's reference point according to ANSI/ASA S12.64-2009 and ISO
%    17208-1:2016(E) based on location of GPS antenna from AIS type 5 messages
%    Processes df20 data (rather than df100)
%    Added Anglegram (angle,freq, and SL)
% 5: Finds correct now also df20 data files automatically
%    Generates beampattern filenames automatically
%    Filters out times with electronic noise
%    Fixed Megan computation: now uses mean instead of sum
%
% 6: Fixed averaging error (forgot to devide by 3 in ANSI eq. (8)
%    Uses Ship Info from Lloyd (Ship_CIB30_SCI_CntShpInfoLloyd7.mat )
% 7: Fixed noise excluder for HARP 0530
%    Uses theo-corrected preamp 780 transfer function
%
%
% TO DO: - filter/smooth AIS track for better x and y speed
%        - track acoustically ship
%        - add SoundTraps and Bprobes
%        - check if file start time in file name is true for SoundTraps
%        - postcalibrate HARP 0545 more precise: e.g. w/ ambient noise data
%% Initialize 

clear all
close all

StartTime = datenum('01/19/2016 08:40:00');
EndTime = datenum('01/19/2016 09:10:00');

MMSI = '477637600';                       %MMSI number of the ship
Stern2EngRoom = 60;                       %Distance from Stern to engine room [m]
ShipStandard = 'ANSI'                     %'ISO', 'ANSI' (requires Stern2EngRoom) or 'GPS'
load Ship_CIB30_SCI_CntShpInfoLloyd7.mat; %load processed AIS data
load CINMS_B_30_ArrayDefinition.mat       %load array definition of CINMS_B_30
SelRec = [1 1 1 0 0 1 0 0 1 1 1 1];       %selected Recorders (1=YES | 0=NO)
RefRec = 'CINMS_B_30_03';                 %defines (0,0,0) center of array
Num1secSeg = 5;                           %number of segments for averaging
TimeIntSec = 3;                           %time-difference [sec] btw. beginning of segments
fstartind = 10;                           %10Hz for computing received levels
fstopind = 1001;                          %1000Hz for computing received levels
c = 1484;                                 %sound speed [m/s]
ResultFileNamePrefix = ['MMSI' MMSI];     %prefix of file name to save results
NumRec = length(find(SelRec));            %number of selected recorders
load AutoMinusTheoGainCorrdB780.mat       %Correction for LF part of 780 preamp

%sampling freq. of data files [Hz]
%Fs = 1e3*[2 2 2 2.4 2.4 2 2.048 2.048 2 2 2 2]; %df100
Fs = 1e3*[10 10 10 2.4 2.4 10 5.048 2.048 10 10 10 10]; %df20 

%% Acoustic data files

% % d100 file path for each acoustic recorder
% filepath = {...
% 'J:\CINMS_B_30_0115\CINMS_B_30_0115_decimated\CINMS_B_30_0115_disks01-08_df100\';
% 'J:\CINMS_B_30_0130\CINMS_B_30_0130_decimated\CINMS_B_30_0130_disks01-08_df100\';
% 'J:\CINMS_B_30_0145\CINMS_B_30_0145_decimated\CINMS_B_30_0145_disks01-07_df100\';
% 'I:\CINMS_B_30_0260\SoundTrap1611153475\decimated_df40\';
% 'I:\CINMS_B_30_0275\SoundTrap805629958\decimated_df40\';
% 'K:\CINMS_B_30_03_decimated\CINMS_B_30_03_disks01-07_df100\';
% 'I:\CINMS_B_30_0460\CIB460_30_A_B007\';
% 'I:\CINMS_B_30_0475\CIB475_30_A_B019\';
% 'L:\CINMS_B_30_0515\CINMS_B_30_0515_decimated\CINMS_B_30_0515_disks01-06_df100\';
% 'L:\CINMS_B_30_0530\CINMS_B_30_0530_decimated\CINMS_B_30_0530_disks01-15_df100\';
% 'L:\CINMS_B_30_0545\CINMS_B_30_0545_decimated\CINMS_B_30_0545_disks01-07_df100\';
% }

% d20 file path for each acoustic recorder
filepath = {...
'H:\CINMS_B_30_0115\CINMS_B_30_0115_decimated\';
'H:\CINMS_B_30_0130\CINMS_B_30_0130_decimated\';
'H:\CINMS_B_30_0145\CINMS_B_30_0145_decimated\';
'H:\CINMS_B_30_0260\SoundTrap1611153475\decimated_df40\';
'H:\CINMS_B_30_0275\SoundTrap805629958\decimated_df40\';
'H:\CINMS_B_30_03_decimated\';
'H:\CINMS_B_30_0460\CIB460_30_A_B007\';
'H:\CINMS_B_30_0475\CIB475_30_A_B019\';
'H:\CINMS_B_30_0515\CINMS_B_30_0515_decimated\';
'H:\CINMS_B_30_0530\CINMS_B_30_0530_decimated\';
'H:\CINMS_B_30_0545\CINMS_B_30_0545_decimated\';
'H:\CINMS_B_30_00\';
}

%use only filepaths of selected acoustic recorders
filepath = filepath(find(SelRec));

%build cell array with filenames
filename = {};
% for RecInd = 1:NumRec
%     [~,dummy]=dos(['dir ' filepath{RecInd} '*.wav /s']);
%     if isempty(dummy)
%         filename{RecInd} = '';
%     else    
%         for FileNameIndex = 1:length(dummy)
%             filename{RecInd}{FileNameIndex} = ...
%                 [filepath{RecInd} dummy(FileNameIndex).name];
%         end    
%     end    
% end    

filename = {};
for RecInd = 1:NumRec
    dummy=dos(['dir ' filepath{RecInd} '*.wav /b /s > df20filelist.txt']);
    if isempty(dummy)
        filename{RecInd} = '';
    else        
        shellcmd=['findstr /C:"df20" df20filelist.txt'];
        [status,cmdout] = system(shellcmd);
        filename{RecInd} = strsplit(cmdout,'\n');
        filename{RecInd}(end) =[];
    end    
end
%% Pre-process input parameters              

%select desired recorders
SelRec = find(SelRec);
CINMS_B_30 = CINMS_B_30(SelRec);

%index for ref. recorder
RefRecInd = find(strcmp(RefRec,{CINMS_B_30.label})); 
if isempty(RefRecInd); 
    error('RefRec: Invalid input for reference recorder.'); 
end

Fs = Fs(SelRec);
SegLengthSamples = 1*Fs; % for 1 Hz bins
DataLengthSamples = Num1secSeg*SegLengthSamples;

switch ShipStandard 
    case 'ANSI'
        ShipStandardLabel = 'ANSI S12.64-2009';
    case 'ISO'
        ShipStandardLabel = 'ISO 17208-1:2016';
end 

%% Transfer function for each recorder
% one-sided to save memory

G = cell(NumRec,1);
fRspectra = cell(NumRec,1);

for RecInd=1:NumRec
    %frequency vector for each recorder
    fRspectra{RecInd} =  Fs(RecInd)/2*linspace(0,1,SegLengthSamples(RecInd)/2+1); %frequencies [Hz]
    if strcmp(CINMS_B_30(RecInd).recType,'HARP')
        %preamp 780 was recalibrated with theoratical curve        
        if strcmp(CINMS_B_30(RecInd).preAmp,'780')         
            G{RecInd}(:,1) = 10.^(interp1(CINMS_B_30(RecInd).TF{1}(:,1),CINMS_B_30(RecInd).TF{1}(:,2),fRspectra{RecInd},'spline')/20);
            G{RecInd}(2:22,1) = 10.^((20*log10(G{RecInd}(2:22,1))+AutoMinusTheoGainCorrdB780(2,:).').'/20);
        else
        for findex=1:length(fRspectra{RecInd})
            [~,minindex] = ...
                min(abs(CINMS_B_30(RecInd).TF{1}(:,1)-abs(fRspectra{RecInd}(findex))));
            G{RecInd}(findex,1) = 10.^(CINMS_B_30(RecInd).TF{1}(minindex,2)/20); %uPa
        end
        end
    elseif strcmp(CINMS_B_30(RecInd).recType,'SoundTrap')    
        if strcmp(CINMS_B_30(RecInd).DL_ID,'805629958')
            G{RecInd} = 10^(186/20)*ones(length(fRspectra{RecInd}),1);
        
        elseif strcmp(CINMS_B_30(RecInd).DL_ID,'1611153475')
            G{RecInd} = 10^(180.2/20)*ones(length(fRspectra{RecInd}),1);
        else
            error('Unrecognized SoundTrap.')
        end        
    elseif strcmp(CINMS_B_30(RecInd).recType,'Bprobe')
        G{RecInd} = ones(length(fRspectra{RecInd}),1);
    else
        error(['Unknown Instrument Type: ' CINMS_B_30(RecInd).recType '.'])
    end
end

%% Display

DataIDstring = ['CINMS_B_30 | Ship MMSI' MMSI];
display(['======== ' DataIDstring ' =========='])
display(['Selected Recorders:'])
display([{CINMS_B_30.label}.'])

%% Get File Start Times from File Headers or File Name
for RecInd=1:NumRec    
    %Get file Start Times
    for FileIndex=1:length(filename{RecInd})         
        switch(CINMS_B_30(RecInd).recType)               
            case('HARP')
                %Get Timestamp of first sample from xwavheader
                %==============================================
                % timestamp of first sample = timestamp of first raw file    

                %open file
                fid = fopen(filename{RecInd}{FileIndex},'r');
                %go to time stamp place in xwav header
                fseek(fid,100,'bof');
                %read data
                xwavStartTimeVector(1) = fread(fid,1,'uchar') + 2000;                         % Year
                xwavStartTimeVector(2) = fread(fid,1,'uchar');                                % Month
                xwavStartTimeVector(3) = fread(fid,1,'uchar');                                % Day
                xwavStartTimeVector(4) = fread(fid,1,'uchar');                                % Hour
                xwavStartTimeVector(5) = fread(fid,1,'uchar');                                % Minute
                xwavStartTimeVector(6) = fread(fid,1,'uchar')+ fread(fid,1,'uint16')/1000;    % Decimal seconds    
                %close file
                fclose(fid);
                %convert to serial date in decimal days
                xwavReferenceTime{RecInd}(FileIndex) = datenum(xwavStartTimeVector);
            
            case('SoundTrap')                
                %get file start time from file name + 8h to go from local
                %to GMT in decimal days
                xwavReferenceTime{RecInd}(FileIndex) = ...
                    datenum(filename{RecInd}{FileIndex}(end-19:end-8),...
                    'yymmddHHMMSS') + (8/24); %decimal days in GMT
            case('Bprobe')    
                [~,~,info] = MTRead(filename{RecInd}{FileIndex},1,1,'p');
                xwavReferenceTime{RecInd}(FileIndex) = info.datenumber;
        end
    end
end

%% Time stamps (decimal days) for SL compuations

TimeVector = StartTime:TimeIntSec/(3600*24):EndTime;

%Total number of Events
NumSpec = length(TimeVector);

ReferenceTime = StartTime;

%% Compute spectra (and track acoustically)

% Received Spectra with 1 Hz bins
S =  cell(NumRec,1);
for RecInd=1:NumRec
    S{RecInd} = zeros(NumSpec,SegLengthSamples(RecInd)/2+1);
end    

%Build Location Matrix from selected recoders
CINMS_B_30SelArr = CINMS_B_30(RefRecInd).xyz(SelRec,:);

%loop to read data piece by piece
for EventIndex = 1:NumSpec
    for RecInd=1:NumRec        
        
        %Load data
        %========================================================          
        %choose correct data file
        FileIndex = ...
            find((xwavReferenceTime{RecInd}(:) - TimeVector(EventIndex))<0);
        switch(CINMS_B_30(RecInd).recType)    
            case('HARP')                
                %skip to right datachunk and read            
                data = xwavread1ch(filename{RecInd}{FileIndex(end)},...
                    TimeVector(EventIndex),DataLengthSamples(RecInd));
                %make colomns=channels
                data=data.';
            case('SoundTrap')
                StartSample = round(Fs(RecInd)*(TimeVector(EventIndex) - ...
                    xwavReferenceTime{RecInd}(FileIndex(end)))*(24*3600));
                EndSample = StartSample+DataLengthSamples(RecInd)-1;
                siz = wavread(filename{RecInd}{FileIndex(end)},'size');
                if EndSample <= siz(1)
                    data = wavread(filename{RecInd}{FileIndex(end)},...
                        [StartSample EndSample]);
                else
                    disp('SoundTrap Data File End. Getting data from next file.')
                    data = wavread(filename{RecInd}{FileIndex(end)},...
                        [StartSample siz(1)]);
                    data = [data wavread(filename{RecInd}{FileIndex(end)+1},...
                        [1 DataLengthSamples(RecInd)-(siz(1)-StartSample)])];
                end    
            case('Bprobe')
                StartSample = round(Fs(RecInd)*(TimeVector(EventIndex) - ...
                    xwavReferenceTime{RecInd}(FileIndex(end)))*(24*3600));
                [data,~,~] = MTRead(filename{RecInd}{FileIndex(end)},...
                    StartSample,StartSample+DataLengthSamples(RecInd)-1,'p');                                
        end
         
        %Compute received spectra
        %==============================
                
        %In case is less than anticipated due to end of file
        NumSeg = floor(length(data)/SegLengthSamples(RecInd));
        %Fourier Transform
        if strcmp('CINMS_B_30_0530',CINMS_B_30(RecInd).label)
            %avoid noisy times
            data = reshape(data(1:NumSeg*SegLengthSamples(RecInd)),SegLengthSamples(RecInd),NumSeg);
            %S{RecInd}(EventIndex,:) = zeros(1,1:SegLengthSamples(RecInd)/2+1);
            for SegInd=1:NumSeg                
                goodindex(SegInd) = length(find(abs(data(:,SegInd)) < 7*mean(abs(data(:,SegInd)))));
                if goodindex(SegInd) == SegLengthSamples(RecInd)
                    X = fft(data(:,SegInd));
                    %compute onesided spectra and calibrate with HARP specific hydrophone transfer functions
                    S{RecInd}(EventIndex,:) = S{RecInd}(EventIndex,:) + ...
                        (((2*abs(X(1:SegLengthSamples(RecInd)/2+1,:)).^2).*G{RecInd}.^2)...
                        /(SegLengthSamples(RecInd)*Fs(RecInd))).';                   
                end
            end
            %normalize by number of Segments used (that had no noise)
            S{RecInd}(EventIndex,:) = S{RecInd}(EventIndex,:)/length(find(goodindex==SegLengthSamples(RecInd)));
        else        
            X = fft(reshape(data(1:NumSeg*SegLengthSamples(RecInd)),...
                SegLengthSamples(RecInd),NumSeg));
            %compute onesided spectra and calibrate with HARP specific hydrophone transfer functions
            S{RecInd}(EventIndex,:) = ...
                (sum(2*abs(X(1:SegLengthSamples(RecInd)/2+1,:)).^2,2).*G{RecInd}.^2)...
                /(NumSeg*SegLengthSamples(RecInd)*Fs(RecInd));            
        end
    end
    
    %display progress
    if mod(EventIndex,round(NumSpec/10))==0
        display([num2str(round(100*EventIndex/NumSpec)) '% done'])   
    end
end

%% PostCorrect Instrument Calibration

%based on roughly the same source levels of ALL HARPs at great distances
S{find(strcmp({CINMS_B_30.label},'CINMS_B_30_0545'))} = ...
    S{find(strcmp({CINMS_B_30.label},'CINMS_B_30_0545'))}.*10^(7/10);%7 dB
    
% %HF soundtrap is ~5dB different then STD(standard) soundtrap
% S{find(strcmp({CINMS_B_30.label},'CINMS_B_30_0275'))} = ...
%     S{find(strcmp({CINMS_B_30.label},'CINMS_B_30_0275'))}.*10^(12/10);%12 dB


%% Spectrograms

FontSizeNumber=14;
%for RecInd=1:NumRec
figure;
    RecInd = RefRecInd;    
    imagesc(TimeVector,fRspectra{RecInd},10*log10(S{RecInd}).');
    title([CINMS_B_30(RecInd).label ': MMSI' MMSI],'Interpreter','None',...
        'FontSize',FontSizeNumber)
    axis xy; 
    h = colorbar;caxis([50 130])
    set(get(h,'ylabel'),'string','RL [dB re 1\muPa//Hz]',...
        'FontSize',FontSizeNumber);
    datetick('x','keeplimits');
    xlabel(['GMT on ' datestr(TimeVector(1),'dd-mmm-yyyy')],...
        'FontSize',FontSizeNumber)
    ylabel('Frequency [Hz]','FontSize',FontSizeNumber)
    ylim([0 1000])
    set(gca, 'FontSize',FontSizeNumber)
%end
%% Average Spectra

figure; 
NumLinesPlot = 3;
NumSubPlots = ceil(NumRec/NumLinesPlot);
colorarray = hsv(NumRec);
for RecInd=1:NumRec
    SubPlotInd = floor(RecInd/(NumLinesPlot+1))+1;
    if RecInd == 7; SubPlotInd = 3; end;
    subplot(NumSubPlots,1,SubPlotInd)
        plot(fRspectra{RecInd},10*log10(sum(S{RecInd})/NumSpec),...
            'Color',colorarray(RecInd,:));
        hold on
        ylim([60 130])
        xlim([5 100])
        grid on        
        if SubPlotInd == 1
            title(['MMSI' MMSI ': Spectra time-averaged over '...
                num2str((TimeVector(end)-TimeVector(1))*24*60) 'min'])
        end            
        ylabel('RL [dB re 1\muPa^2//Hz]')
        if (floor((RecInd+1)/(NumLinesPlot+1))+1) > SubPlotInd || RecInd == 6
            legend({CINMS_B_30(RecInd-2:RecInd).label},'Interpreter','None')
        end   
        if RecInd == NumRec
            legend({CINMS_B_30((SubPlotInd-1)*NumLinesPlot+1:RecInd).label},...
                'Interpreter','None')
            xlabel('f [Hz]')
        end
end

%% Process AIS data

%find ship passage in AIS data (give +/- 3 min extra time for smoothing)
AIS = findShips(Ship,[StartTime-(3/(60*24)) EndTime+(3/(60*24))],[],MMSI);
AIS = AIS{1};
NumAISlocs = length(AIS.time);

%map AIS lats/longs to x-y    
[AIS.x, AIS.y] = latlon2xy(AIS.lat,AIS.long,...
        CINMS_B_30(RefRecInd).lat,CINMS_B_30(RefRecInd).long);

%fix potential errors of time stamp due to AIS computer clock drift
%==================================================================

%determine closest point of approach by finding the maximum RL at propeller blade line
PropFreqLineInd = 6; %8Hz  blade rate fundamental (orginates from propeller)
PropCPAInd = zeros(NumRec,1);
AISCPAInd = zeros(NumRec,1);
for RecInd=1:NumRec
    [~,PropCPAInd(RecInd)]=max(S{RecInd}(:,PropFreqLineInd));
    [~,AISCPAInd(RecInd)]=min((AIS.x-CINMS_B_30SelArr(RecInd,1)).^2 + (AIS.y-CINMS_B_30SelArr(RecInd,2)).^2);
end    
[PropCPAInd AISCPAInd (TimeVector(PropCPAInd)-AIS.time(AISCPAInd)).'*24*60]

%determine correct draft
[~,DraftIndex] = min(abs((StartTime+(EndTime-StartTime)/2) - AIS.DraftDay));

%plot AIS data
h=plotShips({AIS},1);
FontSizeNumber=16;
set(gca, 'FontSize',FontSizeNumber)
xlabel('Eastings [km]','FontSize',FontSizeNumber)
ylabel('Northings [km]','FontSize',FontSizeNumber)    
set(gca,'XTick',[-6:3:6])

%plot time-resolution of AIS data
figure;
    subplot(2,1,1)
    plot(AIS.time(1:end-1),diff(AIS.time*24*3600),'.');
    datetick('x')
    ylabel('\DeltaT [sec]')    
subplot(2,1,2)
    plot(AIS.time,AIS.speed,'.');
    datetick('x')
    ylabel('Ship speed [m/s]')
    xlabel('HH:MM')    
    
%compute ship length and width
AIS.length = AIS.GPS2bow + AIS.GPS2stern;
AIS.width = AIS.GPS2port + AIS.GPS2stbd;    

% ISO 17208-1:2016(E): 1/4 of ship length away from stern
AIS.GPS2RefPoint.ISO = [AIS.width/2-AIS.GPS2port AIS.length/4-AIS.GPS2stern];

% ANSI/ASA S12.64-2000: halfway between propeller and engine room
AIS.GPS2RefPoint.ANSI = [AIS.width/2-AIS.GPS2port Stern2EngRoom/2-AIS.GPS2stern];

%Interpolate ship location & speed at the time stamps of the acous. spectra  
MaxAISdt = Inf;%65/(3600*24); % [decimal days]
AIS.xint = zeros(1,NumSpec);
AIS.yint = zeros(1,NumSpec);
AIS.speedint = zeros(1,NumSpec);
AIS.xintISO = zeros(1,NumSpec);
AIS.yintISO = zeros(1,NumSpec);
AIS.xintANSI = zeros(1,NumSpec);
AIS.yintANSI = zeros(1,NumSpec);

warning('off')
for TimeInd=1:NumSpec
    [dt dtInd] = min(abs(AIS.time-TimeVector(TimeInd)));
        if dt < MaxAISdt
            if dtInd >= 3 && dtInd <= NumAISlocs-2
                p.x = polyfit(AIS.time(dtInd-2:dtInd+2),AIS.x(dtInd-2:dtInd+2),1);
                AIS.xint(TimeInd) = polyval(p.x,TimeVector(TimeInd));                 
                p.y=polyfit(AIS.time(dtInd-2:dtInd+2),AIS.y(dtInd-2:dtInd+2),1);
                AIS.yint(TimeInd) = polyval(p.y,TimeVector(TimeInd)); 
                p.s = polyfit(AIS.time(dtInd-2:dtInd+2),AIS.speed(dtInd-2:dtInd+2),1);
                AIS.speedint(TimeInd) = polyval(p.s,TimeVector(TimeInd));
                
                rotangle.x = atan(p.x(1));
                rotangle.y = atan(p.y(1));
                
                % ISO 17208-1:2016(E): 1/4 of ship length away from stern
                AIS.xintISO(TimeInd) = AIS.xint(TimeInd) ...
                    + cos(rotangle.x).*AIS.GPS2RefPoint.ISO(1) + sin(rotangle.y)*AIS.GPS2RefPoint.ISO(2);
                AIS.yintISO(TimeInd) = AIS.yint(TimeInd) ...
                    -sin(rotangle.x)*AIS.GPS2RefPoint.ISO(1) + cos(rotangle.y)*AIS.GPS2RefPoint.ISO(2);

                % ANSI/ASA S12.64-2000: halfway between propeller and engine room
                AIS.xintANSI(TimeInd) = AIS.xint(TimeInd) ...
                    + cos(rotangle.x).*AIS.GPS2RefPoint.ANSI(1) + sin(rotangle.y)*AIS.GPS2RefPoint.ANSI(2);
                AIS.yintANSI(TimeInd) = AIS.yint(TimeInd) ...
                    -sin(rotangle.x)*AIS.GPS2RefPoint.ANSI(1) + cos(rotangle.y)*AIS.GPS2RefPoint.ANSI(2);                
            end
        else
            %longer gaps in the AIS data
            disp([datestr(TimeVector(TimeInd)) ': AIS data gap > '...
                num2str(MaxAISdt) 'sec.'])  
        end    
end    
warning('on')

switch ShipStandard 
    case 'ANSI'
        AIS.xint = AIS.xintANSI;
        AIS.yint = AIS.yintANSI;
        disp('Ship locations computed for Ship Reference Point according to ANSI/ASA.')
    case 'ISO'
        AIS.xint = AIS.xintISO;
        AIS.yint = AIS.yintISO;
        disp('Ship locations computed for Ship Reference Point according to ISO.')
    case 'GPS'    
        disp('Ship locations computed for Ship Reference Point according to GPS antenna location.')
end        
    
%% Compute Speed vectors from track
%needed for computing the radiation pattern

%speed from interpolated AIS track
AIS.vx = diff(AIS.xint)./(diff(TimeVector)*24*3600); %m/s
AIS.vy = diff(AIS.yint)./(diff(TimeVector)*24*3600); %m/s
AIS.v = sqrt(AIS.vx.^2 + AIS.vy.^2); %m/s

%speed from raw AIS track
AIS.vxraw = diff(AIS.x)./(diff(AIS.time)*24*3600); %m/s
AIS.vyraw = diff(AIS.y)./(diff(AIS.time)*24*3600); %m/s
AIS.vraw = sqrt(AIS.vxraw.^2 + AIS.vyraw.^2); %m/s

figure;
subplot(3,1,1)
    plot(TimeVector,AIS.xint,'.r');
    hold on
    plot(AIS.time,AIS.x,'.b');    
    datetick('x','keeplimits');
    xlabel(['GMT on ' datestr(TimeVector(1),'dd-mmm-yyyy')]);    
subplot(3,1,2)
    plot(TimeVector,AIS.yint,'.r');
    hold on
    plot(AIS.time,AIS.y,'.b');    
    datetick('x','keeplimits');
    xlabel(['GMT on ' datestr(TimeVector(1),'dd-mmm-yyyy')]);    
subplot(3,1,3)
    plot(TimeVector,AIS.speedint,'.c');
    hold on
    plot(AIS.time,AIS.speed,'.b');  
    plot(TimeVector(1:end-1),AIS.v,'.m')
    plot(AIS.time(1:end-1),AIS.vraw,'.r')
    datetick('x','keeplimits');
    xlabel(['GMT on ' datestr(TimeVector(1),'dd-mmm-yyyy')]);    
    legend('AIS int.','AIS raw','AIS int. track','AIS raw track')
    ylim([0 15])
    grid on
    
%% Compute radiation pattern

alpha = zeros(NumRec,NumSpec-1);        %horizontal angle [deg]
phi = zeros(NumRec,NumSpec-1);          %vertical angle [deg]
R = zeros(NumRec,NumSpec-1);            %slant range [m]
RLship = zeros(NumRec,NumSpec-1);       %broad band received level [dB]
SL.broadband = zeros(NumRec,NumSpec-1); %broadband source level [dB]
SL.ALL = cell(NumRec,1);                %source spectra level [dB]

for RecInd = 1:NumRec;    
    SL.ALL{RecInd} = zeros(NumSpec-1,SegLengthSamples(RecInd)/2+1);
    for index=1:NumSpec-1  
            
        %compute beam angles alpha and phi
        thetax=0; % roll angle of Zc =  rotates beampattern around (0deg,0deg)
        thetay = 0;
        thetaz = atan2(AIS.vy(index),AIS.vx(index));
        Rx = [1 0 0; 0 cos(thetax) sin(thetax); 0 -sin(thetax) cos(thetax)];
        Ry = [cos(thetay) 0 sin(thetay); 0 1 0; -sin(thetay) 0 cos(thetay)];
        Rz = [cos(thetaz) sin(thetaz) 0; -sin(thetaz) cos(thetaz) 0; 0 0 1];
        b = CINMS_B_30SelArr(RecInd,:)...
            - [AIS.xint(index) AIS.yint(index) CINMS_B_30(RefRecInd).depthHydro];
        bprime = Rx*Ry*Rz*b.';         
               
        alpha(RecInd,index) =180*atan2(bprime(2),bprime(1))/pi;
        phi(RecInd,index) = 180*atan2(bprime(3),norm(bprime(1:2)))/pi;
                        
        %slant range of each ship location
        R(RecInd,index) = sqrt((AIS.xint(index)-CINMS_B_30SelArr(RecInd,1)).^2 +  (AIS.yint(index)-CINMS_B_30SelArr(RecInd,2)).^2 + CINMS_B_30(RecInd).depthHydro.^2);
                
        %compute source level from spectrum for beam angle           
        RLship(RecInd,index) = 10*log10(trapz(fRspectra{RecInd}(fstartind:fstopind),...
            (S{RecInd}(index,fstartind:fstopind))));%...
        SL.broadband(RecInd,index) = ...
            RLship(RecInd,index) + 20*log10(R(RecInd,index));
        SL.ALL{RecInd}(index,:) = ...
            10*log10(S{RecInd}(index,:))+20*log10(R(RecInd,index));        
    end
end      

%% Indicies of good SL and angles (for plots, further processing)

%bad/no AIS information
PlotIndicies = find(AIS.xint(1:end-1));

%weared angles or source levels
dummy = find(alpha <=3.5 & alpha >=-3.5 &  phi >=-4);
[row col] = ind2sub([NumRec,NumSpec-1],dummy);
PlotIndicies=setdiff(PlotIndicies,col);

%% Plot freq-depend. radiation pattern in GMT

%frequency bin number (=frequency + 1 due to 1 Hz bin spacing)
%fintind = 201:401;  
%fintind = 120;  %119Hz
%fintind = 159;  %158 Hz
%fintind = 11:111;  
%fintind = 201:301;  
%fintind = 801:901;
%fintind = 111:211;

%fintind = 9; %B1 blade line
%fintind = 49; % F4 firing line with no STBD/PORT asymmetry
fintind = 501:601; % OK
%fintind = 401:601; % OK

if length(fintind) == 1
    GMTfilename = [MMSI 'RadPatt' ShipStandard '_' num2str(fRspectra{RefRecInd}(fintind)) 'Hz3m'];
else
    GMTfilename = [MMSI 'RadPatt' ShipStandard '_' num2str(fRspectra{RefRecInd}(fintind(1))) '_'...
                    num2str(fRspectra{RefRecInd}(fintind(end))) 'Hz3m'];
end

% Export  Beampattern  
RadPatGMT=cell(NumRec-1,1);
MaxSL=zeros(NumRec-1,1);
for RecInd=1:NumRec-1
    %find good angles
    dummy = find(abs(alpha(RecInd,:)) <=14 &  phi(RecInd,:) >=-10);    
    RadPattIndicies{RecInd}=setdiff(PlotIndicies,dummy);
    dummy = find(abs(alpha(RecInd,:)) >=160 &  phi(RecInd,:) >=-10);    
    RadPattIndicies{RecInd}=setdiff(RadPattIndicies{RecInd},dummy);    
    MaxSL(RecInd) = max(10*log10((sum(10.^(SL.ALLLloyd{RecInd}(RadPattIndicies{RecInd},fintind)/10),2))));
end
MaxSL = max(MaxSL)
for RecInd=1:NumRec-1
    if strcmp(CINMS_B_30(RecInd).recType,'HARP') %&& RecInd ~= NumRec
        RadPatGMT{RecInd} = [-alpha(RecInd,RadPattIndicies{RecInd});...
                                  phi(RecInd,RadPattIndicies{RecInd}); ...
                                  10*log10(sum(10.^(SL.ALLLloyd{RecInd}(RadPattIndicies{RecInd},fintind)/10),2)).'-MaxSL].';    
        dummy = squeeze(RadPatGMT(RecInd,:,:));
        %save(['D:\Data\gmt\ShipBeamPattern\Ship' MMSI '_' datestr(ReferenceTime,'yyyymmdd_HHMMSS') '_HARP' num2str(RecInd) '.dat'],'dummy','-ascii','-tabs');
    end
end

%plot data stereographically (-Je for azimuthal equial distant)
dos(['psbasemap -P -Xc-2c -Yc -R0/360/-90/0 -Js0/-90/5.5c/0 -Bg10'...
    ' -BcD:\Data\gmt\ShipBeamPattern\SpherLab --MAP_FRAME_TYPE=plain'... 
    ' --FONT_TITLE=18p -B+t"MMSI ' MMSI ' at ' num2str(mean(AIS.speed),'%6.1f')...
    ' m/s on ' datestr(TimeVector(1),'dd/mm/yyyy HH:MM') ' GMT" --ANNOT_FONT_SIZE_PRIMARY=12p -K -V > ' GMTfilename '.ps'])
gmt(['pstext D:\Data\gmt\ShipBeamPattern\LatLab.txt -R0/360/-90/0 -Js0/-90/5.5c/0 -F+f -Gwhite -TO -O -K -V >> ' GMTfilename '.ps'])
%draw colobar via command line to be able to plot greek letter mu
%! psscale -D13/0/9/0.75 -Xa -Ya6c -CSLverylowjet -I0.5 --ANNOT_FONT_SIZE_PRIMARY=12  -B10f5:"SL [dB re 1 @~\155@~Pa \100 1m]":/:: -P -O -K >> shipfreq.ps
dos(['psscale -D13/0/9/0.75 -Xa -Ya6c -CSL0to25 -I0.5 --ANNOT_FONT_SIZE_PRIMARY=12  -B10f5:"Level [dB re 1 @~\155@~Pa]":/:: -P -O -K >> ' GMTfilename '.ps'])

for RecInd=1:NumRec-1
    switch CINMS_B_30(RecInd).recType
        case 'HARP'            
            gmt(['psxy -R -Js -Sc0.3c -O -K -V -CC:\programs\gmt5\share\cpt\SL0to25.cpt >> ' GMTfilename '.ps'],RadPatGMT{RecInd})
        case 'SoundTrap'
            gmt(['psxy -R -Js -Ss0.3c -O -K -V -CC:\programs\gmt5\share\cpt\SL0to25.cpt >> ' GMTfilename '.ps'],squeeze(RadPatGMT(RecInd,:,:)))
        case 'Bprobe'
            gmt(['psxy -R -Js -Ss0.3c -O -K -V -CC:\programs\gmt5\share\cpt\SL0to25.cpt >> ' GMTfilename '.ps'],squeeze(RadPatGMT(RecInd,:,:)))
    end     
end    
        
%% Compute Source Levels according to ANSI/ASA S12.64-2009

%Compute Data Window Length (DWL) in meter from Eq.1
Ind0115 = find(strcmp({CINMS_B_30.label},'CINMS_B_30_0115')); %STBD
Ind0515 = find(strcmp({CINMS_B_30.label},'CINMS_B_30_0515')); %PORT
[RCPA.STBD RCPAind.STBD] = min(R(Ind0115,:));
[RCPA.PORT RCPAind.PORT] = min(R(Ind0515,:));
dCPA.STBD = sqrt(RCPA.STBD^2 + CINMS_B_30(Ind0115).depthHydro.^2);
dCPA.PORT= sqrt(RCPA.PORT^2 + CINMS_B_30(Ind0515).depthHydro.^2);
DWL.STBD = 2*dCPA.STBD*tand(30);
DWL.PORT= 2*dCPA.PORT*tand(30);

%Data Window Period (DWP) in sec from Eq. 2
DWP.STBD.length = DWL.STBD/AIS.v(RCPAind.STBD);
DWP.PORT.length= DWL.PORT/AIS.v(RCPAind.PORT);

[~,DWP.STBD.StartInd]=min(abs(TimeVector-(TimeVector(RCPAind.STBD)...
    -DWP.STBD.length/(2*3600*24))));
[~,DWP.STBD.EndInd]=min(abs(TimeVector-(TimeVector(RCPAind.STBD)...
    +DWP.STBD.length/(2*3600*24))));
[~,DWP.PORT.StartInd]=min(abs(TimeVector-(TimeVector(RCPAind.PORT)...
    -DWP.PORT.length/(2*3600*24))));
[~,DWP.PORT.EndInd]=min(abs(TimeVector-(TimeVector(RCPAind.PORT)...
    +DWP.PORT.length/(2*3600*24))));

% Source Level
SL.STBD = zeros(4,length(fRspectra{Ind0115})); %15,30,45deg and total
SL.PORT = zeros(4,length(fRspectra{Ind0515})); %15,30,45deg and total
SL.TOTAL = zeros(1,length(fRspectra{Ind0515}));
for RecInd=1:3
    TL = R(Ind0115+RecInd-1,DWP.STBD.StartInd:DWP.STBD.EndInd).^2;
    TL = TL(ones(length(fRspectra{Ind0115}),1),:).';
    SL.STBD(RecInd,:) = 10*log10(...
        sum(S{Ind0115+RecInd-1}(DWP.STBD.StartInd:DWP.STBD.EndInd,:).*TL)...
        ./(DWP.STBD.EndInd-DWP.STBD.StartInd+1));
    TL = R(Ind0515+RecInd-1,DWP.PORT.StartInd:DWP.PORT.EndInd).^2;
    TL = TL(ones(length(fRspectra{Ind0515}),1),:).';
    SL.PORT(RecInd,:) = 10*log10(...
        sum(S{Ind0515+RecInd-1}(DWP.PORT.StartInd:DWP.PORT.EndInd,:).*TL)...
        ./(DWP.PORT.EndInd-DWP.PORT.StartInd+1));
end    
% according to Eq. 8 from ANSI/ASA S12.64-2009 
SL.STBD(4,:) = 10*log10(sum(10.^(SL.STBD(1:3,:)/10))/3);
SL.PORT(4,:) = 10*log10(sum(10.^(SL.PORT(1:3,:)/10))/3);
% according to Eq. 9 from ANSI/ASA S12.64-2009 
SL.TOTAL = (SL.STBD(4,:) + SL.PORT(4,:))/2;

%ANSI/ASA/ISO broadband source level
SL.TOTALbroadband = 10*log10(trapz(fRspectra{Ind0115}(fstartind:fstopind),...
            10.^(SL.TOTAL(fstartind:fstopind)/10)))

%% Compute SL with Lloyd's mirror TL

%source depth [m]
zs=10;

SL.ALLLloyd = cell(NumRec,1);
TL = cell(NumRec,1);
% Urick book
% for RecInd=1:NumRec
%     zr = CINMS_B_30(RecInd).depthHydro;    
%     l0 = 4*zs*zr./(c./fRspectra{1});
%     for Rind = 1:length(R(RecInd,:))
%         TL{RecInd}(Rind,:) = (2/R(RecInd,Rind)^2)*(1-cos(pi*l0./R(RecInd,Rind)));
%         SL.ALLLloyd{RecInd}(Rind,:) = 10*log10(S{RecInd}(Rind,:)./TL{RecInd}(Rind,:));
%     end
% end 


%Kuperman book
k = 2*pi./(c./fRspectra{1}); %wavenumber
Rc = ones(length(k),1); %reflection coefficient
% Rc(60:80) = linspace(1,0.9,21); 
% Rc(80:end) = 0.9; 
Rc = Rc.';

for RecInd=1:NumRec
    zr = CINMS_B_30(RecInd).depthHydro;           
    for Rind = 1:length(R(RecInd,:))        
        %compute acoustic pressure p
        R1 = sqrt((R(RecInd,Rind).^2-zr.^2) + (zr-zs).^2);
        R2 = sqrt(R(RecInd,Rind).^2-zr.^2 + (zr+zs).^2);         
        TL{RecInd}(Rind,:) = abs(1*exp(1i*k.*R1)./R1 - Rc.*exp(1i*k.*R2)./R2).^2;
        %curb Lloyds mirrors curve when hit 20*log10R to avoid interference
        %lobes
        findex=1;
        while TL{RecInd}(Rind,findex)<(1/R(RecInd,Rind)^2)
            findex = findex+1;
            if findex == length(fRspectra{1}); break; end;
        end                    
        if findex <= length(fRspectra{1})
            TL{RecInd}(Rind,findex:end) = 1/R(RecInd,Rind)^2;
        end            
        SL.ALLLloyd{RecInd}(Rind,:) = 10*log10(S{RecInd}(Rind,:)./TL{RecInd}(Rind,:));
    end
end 
    
figure;
    imagesc(10*log10(TL{end}).'); axis xy; colorbar
     %semilogx(fRspectra{Ind0115},10*log10(TL{end}(286,:).'));xlim([5 1e3])      

STBDlloyd=0;
for index=1:3    
    STBDlloyd = STBDlloyd+sum(10.^(SL.ALLLloyd{index}(DWP.STBD.StartInd:DWP.STBD.EndInd,:)/10))...
            ./(DWP.STBD.EndInd-DWP.STBD.StartInd+1);

end   
PORTlloyd=0;
for index=1:3
    PORTlloyd = PORTlloyd+sum(10.^(SL.ALLLloyd{4+index}(DWP.PORT.StartInd:DWP.PORT.EndInd,:)/10))...
            ./(DWP.PORT.EndInd-DWP.PORT.StartInd+1);
end        
SL.TOTALlloyd = (10*log10(STBDlloyd/3)+10*log10(PORTlloyd/3))/2;
 
%% PLOT SL(f)

figure;
subplot(2,1,1)    
    semilogx(fRspectra{Ind0115},SL.STBD.'); hold on
    semilogx(fRspectra{Ind0115},SL.TOTAL,'k');    
    title([MMSI ': STBD acc. to ANSI/ASA S12.64-2009 '])
    ylabel('SL [dB re 1 \muPa//Hz @ 1m]')    
    legend([{CINMS_B_30(Ind0115:Ind0115+2).label},'STBD','TOTAL'],...
        'Interpreter','None','Location','SouthWest')
    grid on
subplot(2,1,2)    
    semilogx(fRspectra{Ind0115},SL.PORT.'); hold on
    semilogx(fRspectra{Ind0115},SL.TOTAL,'k');  
    title([MMSI ': PORT acc. to ANSI/ASA S12.64-2009 '])
    xlabel('Frequency [Hz]')    
    ylabel('SL [dB re 1 \muPa//Hz @ 1m]')    
    legend([{CINMS_B_30(Ind0115:Ind0115+2).label},'PORT','TOTAL.'],...
        'Interpreter','None','Location','SouthWest')
    grid on
    
%% Difference

figure;
colorarray = hsv(3);
    for hdrInd=1:3
        semilogx(fRspectra{Ind0115},SL.STBD(hdrInd,:)-SL.PORT(hdrInd,:),...
            'b','LineWidth',1,'Color',colorarray(hdrInd,:));
        hold on
    end        
    legend('15deg','30deg','45deg')
%% Plot STBD and PORT  
        
FontSizeNumber = 10;%
figure;        
    %plot(fRspectra{Ind0115},SL.TOTAL,'k'); xlim([5 100])
    semilogx(fRspectra{Ind0115},SL.STBD(1,:),'b','LineWidth',2);
    hold on
    semilogx(fRspectra{Ind0115},SL.PORT(1,:),'r','LineWidth',2);
    %semilogx(fRspectra{Ind0115},SL.TOTAL,'b','LineWidth',2); 
    set(gca, 'FontSize',FontSizeNumber) 
    xlabel('Frequency [Hz]','FontSize',FontSizeNumber)       
    ylabel('SL [dB re 1 \muPa//Hz @ 1m]','FontSize',FontSizeNumber)       
    title(['MMSI' MMSI ' ' ShipStandardLabel],....
        'FontSize',FontSizeNumber)     
    legend('STBD','PORT')
    grid on          
    ylim([145 190])
    xlim([5 1e3])

    

%% Compute ANSI/ISO 1/3 octave bands

n=7:36;                                                     %number of frequencies  
SL.f_thirdoctave.center = 10.^(n/10);                       %center frequency [Hz]
SL.f_thirdoctave.low = 10^(-1/20)*SL.f_thirdoctave.center;  %lower frequency [Hz]
SL.f_thirdoctave.high = 10^(1/20)*SL.f_thirdoctave.center;  %upper frequency [Hz]

for index=1:length(n)
    startf = round(SL.f_thirdoctave.low(index));
    endf = round(SL.f_thirdoctave.high(index));
%     if endf > 1000 %Hz
%         endf = 1000;
%     end        
    SL.TOTALthirdoctave(index) = ...
        trapz(startf:endf,10.^(SL.TOTAL(startf+1:endf+1)/10));   
end    
SL.TOTALthirdoctave = 10*log10(SL.TOTALthirdoctave);

%% Plot 1Hz and 1/3 octave spectra

FontSizeNumber = 14;%
figure;        
    semilogx(SL.f_thirdoctave.center,SL.TOTALthirdoctave,'k','LineWidth',3);
    hold on
    semilogx(fRspectra{Ind0115},SL.TOTAL,'k','LineWidth',1);
    set(gca, 'FontSize',FontSizeNumber) 
    xlabel('Frequency [Hz]','FontSize',FontSizeNumber)       
    ylabel('Source Level [dB re 1 \muPa @ 1m]','FontSize',FontSizeNumber)       
%     title([ShipStandardLabel ' ' AIS.ShipName ' (MMSI' MMSI...
%         ') at ' num2str(mean(AIS.speed)*1.94384,'%6.1f') ' knots | Draft: '...
%         num2str(AIS.Draft(DraftIndex),'%6.1f') ' m'],'FontSize',FontSizeNumber)     
    legend('1/3 octave','1 Hz')
    grid on          
    ylim([145 190])
    xlim([5 1e3])
    
%% Investigation Differences btw. Megan vs. ANSI/ISO

figure;
HARPindex1 = 1;
HARPindex2 = 8;
semilogx(10*log10(S{HARPindex1}(1:end-1,:)).' + repmat(20*log10(R(HARPindex1,:)).',1,5001).','k');
title([AIS.ShipName ' | MMSI: ' MMSI ' | ' datestr(StartTime)])
hold on; 
semilogx((10*log10(S{HARPindex2}(1:end-1,:)).' + repmat(20*log10(R(HARPindex2,:)).',1,5001).'),'b');
legend()
title([AIS.ShipName ' | MMSI: ' MMSI ' | ' datestr(StartTime)])
legend({CINMS_B_30(HARPindex1).label,CINMS_B_30(HARPindex2).label},'Interpreter','None')   
axis([1 5000 120 230])
grid on
xlabel('Freq [Hz]')
ylabel('Source Level [dB re 1 \muPa @ 1m]')     

figure;
    plot(TimeVector(1:end-1),RLship([HARPindex1 HARPindex2],:).')
    datetick('x')
    ylabel('Received Level [dB re 1 \muPa broadband]')     
    legend({CINMS_B_30(HARPindex1).label,CINMS_B_30(HARPindex2).label},'Interpreter','None')  
    title([AIS.ShipName ' | MMSI: ' MMSI ' | ' datestr(StartTime)])
    grid on
    
%% Investigation of differences in TL between site B (~10deg) and 0515

figure;
semilogx(fRspectra{RecInd},10*log10(sum(S{1}(293:297,:))/5)-10*log10(sum(S{end}(311:315,:))/5),'.-'); 
hold on
semilogx(fRspectra{RecInd},10*log10(sum(TL{1}(293:297,:))/5)-10*log10(sum(TL{end}(311:315,:))/5),'-r'); 
title(['RL 0515-SiteB @ 12sec interval around CPA \newline0515:' ...
    datestr(TimeVector(295)) 'B: ' datestr(TimeVector(313))])

xlabel('Freq [Hz]')
ylabel('RL 0515-site B')
xlim([4 5e3])
grid on
ylim([5 35])

figure;
semilogx(fRspectra{RecInd},10*log10(sum(S{5}(293:297,:))/5)-10*log10(sum(S{7}(293:297,:))/5),'.-'); 
hold on
semilogx(fRspectra{RecInd},10*log10(sum(TL{1}(293:297,:))/5)-10*log10(sum(TL{3}(293:297,:))/5),'-c'); 
xlabel('Freq [Hz]')
ylabel('\Delta TL')
xlim([4 5e3])
grid on
legend('0115-0145','Lloyds z_s=10m','Lloyds z_s=5m','Lloyds z_s=3m','Lloyds z_s=1m')
title(['Source Depth of CSCL South China Sea via Transmission Loss differences\newline12sec interval around CPA:' ...
    datestr(TimeVector(295))])    
%% Compute Megan Spectra

%compute time it takes the ship to travel its own length
IndSiteB = find(strcmp({CINMS_B_30.label},'CINMS_B_30_00')); 
[~,CPAsiteBindex] = min(R(IndSiteB,:));
MeganCPAtraveltimeSec = (AIS.GPS2bow + AIS.GPS2stern)/...
    mean(AIS.v(CPAsiteBindex-2:CPAsiteBindex+2));
NumCPAspectraSiteB = round((MeganCPAtraveltimeSec/2)/TimeIntSec);

%compute mean received 1Hz spectrum
SL.MeganRL = 10*log10(...
        mean((S{IndSiteB}(CPAsiteBindex-NumCPAspectraSiteB:CPAsiteBindex+NumCPAspectraSiteB,:)))).';    
%SL with 20logR
SL.Megan = SL.MeganRL + repmat(20*log10(R(IndSiteB,CPAsiteBindex)).',1,5001).';

%broadband SL
%SL.Meganbroadband = round(10*log10(sum(10.^(SL.Megan(21:1001)/10))));
SL.Meganbroadband = round(10*log10(sum(10.^(SL.Megan(fstartind:fstopind)/10))));

%1/3 octave levels
for index=1:length(SL.f_thirdoctave.center)
    startf = round(SL.f_thirdoctave.low(index));
    endf = round(SL.f_thirdoctave.high(index));
%     if endf > 1000 %Hz
%         endf = 1000;
%     end        
    SL.Megan_thirdoctave(index) = ...
        trapz(startf:endf,10.^(SL.Megan(startf+1:endf+1)/10));   
end    
SL.Megan_thirdoctave = 10*log10(SL.Megan_thirdoctave);
    
%% Compute keel spectra

%compute time it takes the ship to travel its own length
IndSiteKEEL = find(strcmp({CINMS_B_30.label},'CINMS_B_30_03')); 
[~,CPAsiteKEELindex] = min(R(IndSiteKEEL,:));
keelCPAtraveltimeSecBefore = (AIS.GPS2bow+abs(AIS.GPS2RefPoint.ANSI(2)))...
    /mean(AIS.v(CPAsiteKEELindex-2:CPAsiteKEELindex+2));
keelCPAtraveltimeSecAfter = (AIS.GPS2stern-abs(AIS.GPS2RefPoint.ANSI(2)))...
    /mean(AIS.v(CPAsiteKEELindex-2:CPAsiteKEELindex+2));
NumCPAspectraSiteKEELBefore = round((keelCPAtraveltimeSecBefore)/TimeIntSec);
NumCPAspectraSiteKEELAfter = round((keelCPAtraveltimeSecAfter)/TimeIntSec);


%constant range (from CPA of reference point) for the whole ship lenght
SL.KeelRL = 10*log10(...
        mean((S{IndSiteKEEL}(CPAsiteKEELindex-keelCPAtraveltimeSecBefore:CPAsiteKEELindex+keelCPAtraveltimeSecAfter,:)))).';    
SL.Keel = SL.KeelRL + repmat(20*log10(R(IndSiteKEEL,CPAsiteKEELindex)).',1,5001).';

%range varies as reference point transits over KEEL HARP
SL.Keel = 10*log10(...
        mean(S{IndSiteKEEL}(CPAsiteKEELindex-keelCPAtraveltimeSecBefore:CPAsiteKEELindex+keelCPAtraveltimeSecAfter,:)...
         .*repmat(R(IndSiteKEEL,CPAsiteKEELindex-keelCPAtraveltimeSecBefore:CPAsiteKEELindex+keelCPAtraveltimeSecAfter).^2.',1,5001))).';    

   
%broadband SL
SL.Keelbroadband = round(10*log10(sum(10.^(SL.Keel(fstartind:fstopind)/10))));
 
%1/3 octave levels
for index=1:length(SL.f_thirdoctave.center)
    startf = round(SL.f_thirdoctave.low(index));
    endf = round(SL.f_thirdoctave.high(index));
%     if endf > 1000 %Hz
%         endf = 1000;
%     end        
    SL.Keel_thirdoctave(index) = ...
        trapz(startf:endf,10.^(SL.Keel(startf+1:endf+1)/10));   
end    
SL.Keel_thirdoctave = 10*log10(SL.Keel_thirdoctave);

%% Plot RL 1 Hz Spectra: ANSI/ISO, Megan, Keel
 
FontSizeNumber = 10;%
figure;
    semilogx(fRspectra{IndSiteB},SL.KeelRL,'r');
    hold on
    semilogx(fRspectra{IndSiteB},10*log10(sum(S{3}(DWP.STBD.StartInd:DWP.STBD.EndInd,:))...
        ./(DWP.STBD.EndInd-DWP.STBD.StartInd+1)),'k');
    semilogx(fRspectra{IndSiteB},10*log10(sum(S{2}(DWP.STBD.StartInd:DWP.STBD.EndInd,:))...
        ./(DWP.STBD.EndInd-DWP.STBD.StartInd+1)),'-.k');
    semilogx(fRspectra{IndSiteB},10*log10(sum(S{1}(DWP.STBD.StartInd:DWP.STBD.EndInd,:))...
        ./(DWP.STBD.EndInd-DWP.STBD.StartInd+1)),'k');
    semilogx(fRspectra{IndSiteB},SL.MeganRL,'b');        
    set(gca, 'FontSize',FontSizeNumber) 
    xlabel('Frequency [Hz]','FontSize',FontSizeNumber)       
    ylabel('RL [dB re 1 \muPa//Hz]','FontSize',FontSizeNumber)       
     title([AIS.ShipName ' (MMSI' MMSI...
        ') at ' num2str(mean(AIS.speed)*1.94384,'%6.1f') ' knots | Draft: '...
        num2str(AIS.Draft(DraftIndex),'%6.1f') ' m'],'FontSize',FontSizeNumber)  
    legend('KEEL','45°','30°','15°','Site B')
    grid on          
    ylim([70 140])
    xlim([5 1e3])
    

%% PLOT FOR PAPER: 20logR SL 1 Hz Spectra: ANSI/ISO, Megan, Keel     
FontSizeNumber = 14;%
figure;                
    semilogx(fRspectra{IndSiteB},SL.Keel,'r','LineWidth',2);           
    hold on  
    %semilogx(fRspectra{Ind0115},10*log10(sum(10.^(SL.STBD(1:3,:)/10))/3),'k','LineWidth',2);            
    semilogx(fRspectra{Ind0115},SL.TOTAL,'k','LineWidth',2);            
    semilogx(fRspectra{IndSiteB},SL.Megan,'b','LineWidth',2);       
    set(gca, 'FontSize',FontSizeNumber) 
    xlabel('Frequency [Hz]','FontSize',FontSizeNumber)       
    ylabel('Spher. Spread. SL [dB re 1 \muPa//Hz @ 1m]','FontSize',FontSizeNumber)       
    %title([AIS.ShipName ' (MMSI' MMSI...
    %    ') at ' num2str(mean(AIS.speed)*1.94384,'%6.1f') ' knots | Draft: '...
    %    num2str(AIS.Draft(DraftIndex),'%6.1f') ' m'],'FontSize',FontSizeNumber)  
    legend( ['KEEL: ' num2str(round(SL.Keelbroadband)) 'dB'],...
            ['ANSI: ' num2str(round(SL.TOTALbroadband)) 'dB'],...
            ['Site B: ' num2str(round(SL.Meganbroadband)) 'dB'])
    grid on          
    ylim([140 220])
    xlim([5 1e3])        

%% PLOT FOR PAPER: Lloyd-Mirror Source Levels

SL.KeelLloyd = 10*log10(mean(10.^(SL.ALLLloyd{4}...
    (CPAsiteKEELindex-keelCPAtraveltimeSecBefore:CPAsiteKEELindex+keelCPAtraveltimeSecAfter,:)/10)));        
SL.MeganLloyd = 10*log10(mean(10.^(SL.ALLLloyd{end}...
    (CPAsiteBindex-NumCPAspectraSiteB:CPAsiteBindex+NumCPAspectraSiteB,:)/10)));
        

FontSizeNumber = 14;%
figure;
    semilogx(fRspectra{IndSiteB},SL.KeelLloyd,'r','LineWidth',1.5);            
    hold on
    semilogx(fRspectra{IndSiteB},SL.TOTALlloyd,'k','LineWidth',1.5);            
    semilogx(fRspectra{IndSiteB},SL.MeganLloyd,'b','LineWidth',1.5);            
    set(gca, 'FontSize',FontSizeNumber) 
    xlabel('Frequency [Hz]','FontSize',FontSizeNumber)       
    ylabel('Lloyd Mirror SL [dB re 1 \muPa//Hz @ 1m]','FontSize',FontSizeNumber)       
%      title([AIS.ShipName ' (MMSI' MMSI...
%         ') at ' num2str(mean(AIS.speed)*1.94384,'%6.1f') ' knots | Draft: '...
%         num2str(AIS.Draft(DraftIndex),'%6.1f') ' m | SourceDepth: '...
%         num2str(zs,'%6.1f') ' m'],'FontSize',FontSizeNumber)  
    legend( ['KEEL: ' num2str(round(10*log10(sum(10.^(SL.KeelLloyd(fstartind:fstopind)/10))))) 'dB'],...
            ['ANSI: ' num2str(round(10*log10(sum(10.^(SL.TOTALlloyd(fstartind:fstopind)/10))))) 'dB'],...
            ['Site B: ' num2str(round(10*log10(sum(10.^(SL.MeganLloyd(fstartind:fstopind)/10))))) 'dB']);
    grid on          
    ylim([140 220])
    xlim([5 1e3])
    
%% PLOT FOR PAPER: 1/3 octrave STANDARD VS MEGAN

FontSizeNumber = 16;%
figure;        
    semilogx(SL.f_thirdoctave.center,SL.TOTALthirdoctave,'-k','LineWidth',3,'MarkerFaceColor','k');
    hold on
    %semilogx(fRspectra{Ind0115},SL.TOTAL,'k','LineWidth',1);
    %semilogx(fRspectra{IndSiteB},MeganSpectra,'b');
    semilogx(SL.f_thirdoctave.center,SL.Megan_thirdoctave,'--k','LineWidth',3,'MarkerFaceColor','b');
    set(gca, 'FontSize',FontSizeNumber) 
    xlabel('Frequency [Hz]','FontSize',FontSizeNumber)       
    ylabel('1/3 Oct. Source Level [dB re 1 \muPa @ 1m]','FontSize',FontSizeNumber)       
%     title([ShipStandardLabel ' ' AIS.ShipName ' (MMSI' MMSI...
%         ') at ' num2str(mean(AIS.speed)*1.94384,'%6.1f') ' knots | Draft: '...
%         num2str(AIS.Draft(DraftIndex),'%6.1f') ' m'],'FontSize',FontSizeNumber)         
    grid on
    %axis square
    ylim([160 190])
    xlim([5 3e3])
%     legend(['1/3 octave bands'],...
%         ['1 Hz bands'],...
%         ['McKenna et al., 2012: ' num2str(round(10*log10(sum(10.^(MeganSpectra(21:fstopind)/10))))) 'dB'])
     legend(['ASA/ANSI/ISO:' num2str(round(SL.TOTALbroadband)) 'dB'],...
        ['Site B: ' num2str(SL.Meganbroadband) 'dB'])
     
%% Plot Slant Ranges

figure;
    for RecInd=1:NumRec
        scatter(alpha(RecInd,:),phi(RecInd,:),[],R(RecInd,:));
        hold on
    end
    colorbar
    caxis([200 3000])
    %datetick('x')

%% Blade lines and diesel firing lines    

MeanSpeedKnots = mean(AIS.speed)*1.94384;
MaxRPM = 84; %RPM at max speed for 58,100 kW engine with 10 cylinders at 5810kW/cyl.
%from MAN book p.28: EngineRPMfromMAN.pdf 
%MaxRPM = 94; %MSC Monteerey RPM at max speed from MAN book: RPMfromMAN.pdf 
NumBlades = 6;
NumCyl = 10;
NumSpecLines=25;
% frequency,blade line, firing line
SpecLines.computed = [...
            [1:NumSpecLines];
            NumBlades*((MeanSpeedKnots./str2double(AIS.ContShipInfo.MaxSpeedKnt))*MaxRPM/60)*[1:NumSpecLines]; 
            NumCyl*((MeanSpeedKnots./str2double(AIS.ContShipInfo.MaxSpeedKnt))*MaxRPM/60)*[1:NumSpecLines];
            ].'

SpecLines.measured = [...                
                8.0000  188.0183;
               14.0000  178.3970;
               21.0000  183.6209;
               29.0000  189.9931;
               36.0000  186.0054;
               43.0000  183.3621;
               49.0000  180.2905;
               65.0000  178.1960;
               74.0000  176.8122;
               84.0000  176.4894;
               95.0000  174.0171;
               98.0000  172.4258;
              117.0000  172.1908;
              122.0000  170.7690;
              158.0000  167.2722;
              ];
          
%% Cavitation tip index Kt
% the lower Kt, the more likely for caviation to occur

%Draft usually measured amidship unless special ship construction!!!
% ==> draft at propeller is less than the measured one
PropBladeTipBelowSurfM = 1; %[meter] distance btw. tip of prop. blade and surface
PropDiameter = 5.8; % [meter]
%Tip Index Kt (eq. 8.28 on p.271 in Ross)
Kt = (2/3)*(PropBladeTipBelowSurfM+9)...
    ./(((MaxRPM*(MeanSpeedKnots./str2double(AIS.ContShipInfo.MaxSpeedKnt)))*PropDiameter/100).^2)
          
%% Plot Angle versus Freq

%ignore weard angles in beginning and end

GoodAnglesIndex = 1:NumSpec-1;%PlotIndicies;

FreqEnd = [100 500];% 2000];
ColorBarLim = [120 180; 120 180];% 110 170];
for FreqEndIndex=1:length(FreqEnd)
    figure;
    imagesc(fRspectra{RecInd},abs(alpha(RecInd,GoodAnglesIndex)),...
        SL.ALL{RecInd}(:,GoodAnglesIndex))        
    colorbar
    for RecInd=1:NumRec
        subplot(2,4,RecInd)        
        imagesc(fRspectra{RecInd},abs(alpha(RecInd,GoodAnglesIndex)),...
            SL.ALL{RecInd}(GoodAnglesIndex,:))               
        %axis([0 FreqEnd(FreqEndIndex) 0 180])
        xlim([0 500])
        caxis(ColorBarLim(FreqEndIndex,:))        
        %set(gca,'YTick',[0:45/2:180])
        %set(gca,'YTickLabel',{'Bow','','45°','','0°','','45°','','Stern'})
        title(CINMS_B_30(RecInd).label,'Interpreter','None')        
    end  
    B=colorbar;
    set(B, 'Position', [.8314 .09 .0381 .8150]) 
end    

%% PAPER PLOT: Plot Angle versus Freq

%ignore weard angles in beginning and end

GoodAnglesIndex = 1:NumSpec-1;%PlotIndicies;
FreqEnd = [500];% 2000];
ColorBarLim = [20 140];% 110 170];
RecInd = fliplr([1 2 3 5 6 7]); %needs to be flipped for the tight spacing plotting
figure;
nRows = 3 ;
nCols = 2 ;
%set( gcf, 'Units', 'normalized', 'Position', [0.1,0.1,0.5,0.8] ) ;
% - Create grid of axes.
[blx, bly] = meshgrid( 0.06:0.93/nCols:0.95, 0.05:0.95/nRows:0.95 ) ;
hAxes = arrayfun( @(x,y) axes( 'Position', [x, y, 0.95*0.95/nCols, 0.95*0.95/nRows] ), blx, bly, 'UniformOutput', false ) ;
    for k=1:numel(hAxes)
        %axes(hAxes{k}) ;
        figure;
        imagesc(10*log10(S{RecInd(k)}(GoodAnglesIndex,:).'));  
        10*log10(min(S{RecInd(k)}(:)))
        axis xy
        ylim([5 FreqEnd])        
        caxis(ColorBarLim)
        anglevector = abs(alpha(RecInd(k),GoodAnglesIndex));
        Ticks = [10 30 60 90 90+30 90+60 90+80];
        angleind=[];
        for dummyindex=1:length(Ticks)
            [~,angleind(dummyindex)] = min(abs(anglevector-Ticks(dummyindex)));
        end        
        set(gca,'XTick',[1 angleind length(anglevector)])
        set(gca,'TickLength',[0.03 0.035])
%         if k == 1 || k == 4;
%             set(gca,'XTickLabel',{'      Bow','  80° ','30°     ','0°','     30°','  80°','Stern     '})
%             xlabel('Azimuth [°]')
%         else
            set(gca,'XTickLabel',{'','','','','','','','',''})            
%        end
        %if k == 2; ylabel('Frequency [Hz]'); end
        set(gca,'YTick',[0:100:FreqEnd])
%         if k <= 3 
%             set(gca,'YTickLabel',{'0','','200','','400',''});
%         else
            set(gca,'YTickLabel',{'','','','','',''});
%        end        
        %title(CINMS_B_30(RecInd).label,'Interpreter','None')        
    end  
    %B=colorbar;
    %set(B, 'Position', [.8314 .09 .0381 .8150])
    %set(B, get(B,'ylabel'),'string','SL' )
    %ylabel(B,'SL [dB re 1 \muPa//Hz @ 1m]')
    
%% Keel apsect
RecInd=4;
figure;
        imagesc(10*log10(S{RecInd}(GoodAnglesIndex,:).'))        
        B = colorbar;
        colormap
        ylabel(B,'SL [dB re 1 \muPa//Hz @ 1m]')
        axis xy
        ylim([5 FreqEnd])
        caxis(ColorBarLim)
        anglevector = abs(phi(RecInd,GoodAnglesIndex));
        [~,minphiindex]=max(anglevector);
        Ticks = [10 30 60 90];
        angleind=[];
        for dummyindex=1:length(Ticks)
            [~,angleind(dummyindex)] = min(abs(anglevector(1:minphiindex)-Ticks(dummyindex)));
        end
        Ticks = fliplr(Ticks);
        Ticks(1) = [];
        for dummyindex=1:length(Ticks)
            [~,unsinn] = min(abs(anglevector(minphiindex:end)-Ticks(dummyindex)));
            angleind(dummyindex+length(Ticks)+1) = unsinn + minphiindex;
        end             
        set(gca,'XTick',[1 angleind length(anglevector)])
        set(gca,'TickLength',[0.03 0.035])
        set(gca,'XTickLabel',{'','','','','','','','',''})
        set(gca,'YTick',[0:100:FreqEnd])
        set(gca,'YTickLabel',{'','','','','',''});
%        set(gca,'YTickLabel',{'0','','200','','400',''})  
%         xlabel('Vertical Angle [°]')
%         ylabel('Frequency [Hz]');

%% site B
RecInd=8;
figure;
        imagesc(10*log10(S{RecInd}(GoodAnglesIndex,:).'))     
        10*log10(max(S{RecInd}(:)))
        B = colorbar;
        %ylabel(B,'SL [dB re 1 \muPa//Hz @ 1m]')
        axis xy
        ylim([0 FreqEnd])
        caxis([ColorBarLim])
        anglevector = abs(alpha(RecInd,GoodAnglesIndex));
        [~,minphiindex]=max(anglevector);
        Ticks = [30 60 90 90+30 90+60];
        angleind=[];
        for dummyindex=1:length(Ticks)
            [~,angleind(dummyindex)] = min(abs(anglevector-Ticks(dummyindex)));
        end        
        set(gca,'XTick',[1 angleind length(anglevector)])
        set(gca,'TickLength',[0.03 0.035])
        set(gca,'XTickLabel',{'','','','','','','','',''})
        set(gca,'YTick',[0:100:FreqEnd])
        set(gca,'YTickLabel',{'','','','','',''});
%        set(gca,'YTickLabel',{'0','','200','','400',''})  
%         xlabel('Vertical Angle [°]')
%         ylabel('Frequency [Hz]');

%% Plot Lloyds SL 1 Hz Spectra: ANSI/ISO, Megan, Keel     

%ignore weard angles in beginning and end

GoodAnglesIndex = 1:NumSpec-1;%PlotIndicies;
FreqEnd = [500];
ColorBarLim = [150 220; 150 220];% 110 170];
RecInd = fliplr([1 2 3 5 6 7]); %needs to be flipped for the tight spacing plotting
figure;
nRows = 3 ;
nCols = 2 ;
%set( gcf, 'Units', 'normalized', 'Position', [0.1,0.1,0.5,0.8] ) ;
% - Create grid of axes.
[blx, bly] = meshgrid( 0.06:0.93/nCols:0.95, 0.05:0.95/nRows:0.95 );
hAxes = arrayfun( @(x,y) axes( 'Position', [x, y, 0.95*0.95/nCols, 0.95*0.95/nRows] ), blx, bly, 'UniformOutput', false ) ;
    for k=1:numel(hAxes)
        axes(hAxes{k}) ;
        imagesc((SL.ALLLloyd{RecInd(k)}(GoodAnglesIndex,:)).')     ;           
        axis xy
        ylim([0 FreqEnd])
        %caxis(ColorBarLim(FreqEndIndex,:))
        anglevector = abs(alpha(RecInd(k),GoodAnglesIndex));
        Ticks = [10 30 60 90 90+30 90+60 90+80];
        angleind=[];
        for dummyindex=1:length(Ticks)
            [~,angleind(dummyindex)] = min(abs(anglevector-Ticks(dummyindex)));
        end        
        set(gca,'XTick',[1 angleind length(anglevector)])
        set(gca,'TickLength',[0.03 0.035])
%         if k == 1 || k == 4;
%             set(gca,'XTickLabel',{'      Bow','  80° ','30°     ','0°','     30°','  80°','Stern     '})
%             xlabel('Azimuth [°]')
%         else
            set(gca,'XTickLabel',{'','','','','','','','',''})            
%        end
        %if k == 2; ylabel('Frequency [Hz]'); end
        set(gca,'YTick',[0:100:FreqEnd])
%         if k <= 3 
%             set(gca,'YTickLabel',{'0','','200','','400',''});
%         else
            set(gca,'YTickLabel',{'','','','','',''});
%        end        
        %title(CINMS_B_30(RecInd).label,'Interpreter','None')        
    end  
    %B=colorbar;
    %set(B, 'Position', [.8314 .09 .0381 .8150])
    %set(B, get(B,'ylabel'),'string','SL' )
    %ylabel(B,'SL [dB re 1 \muPa//Hz @ 1m]')
    
% Keel apsect
RecInd=4;
figure;
        imagesc((SL.ALLLloyd{RecInd}(GoodAnglesIndex,:)).')        
        B = colorbar;
        ylabel(B,'SL [dB re 1 \muPa//Hz @ 1m]')
        axis xy
        ylim([0 FreqEnd])
        caxis(ColorBarLim(1,:))
        anglevector = abs(phi(RecInd,GoodAnglesIndex));
        [~,minphiindex]=max(anglevector);
        Ticks = [10 30 60 90];
        angleind=[];
        for dummyindex=1:length(Ticks)
            [~,angleind(dummyindex)] = min(abs(anglevector(1:minphiindex)-Ticks(dummyindex)));
        end
        Ticks = fliplr(Ticks);
        Ticks(1) = [];
        for dummyindex=1:length(Ticks)
            [~,unsinn] = min(abs(anglevector(minphiindex:end)-Ticks(dummyindex)));
            angleind(dummyindex+length(Ticks)+1) = unsinn + minphiindex;
        end             
        set(gca,'XTick',[1 angleind length(anglevector)])
        set(gca,'TickLength',[0.03 0.035])
        set(gca,'XTickLabel',{'','','','','','','','',''})
        set(gca,'YTick',[0:100:FreqEnd])
        set(gca,'YTickLabel',{'','','','','',''});
%        set(gca,'YTickLabel',{'0','','200','','400',''})  
%         xlabel('Vertical Angle [°]')
%         ylabel('Frequency [Hz]');

%%
figure;
    semilogx(10*log10(S{1}(290,:)));
    grid on
    xlim([5 3e3])
    
%% Navy style spectrogram: KEEL

RecInd=1;
RecInd2=8;
figure;
        imagesc(10*log10(S{RecInd}(GoodAnglesIndex,:))-10*log10(S{RecInd2}(GoodAnglesIndex,:)))        
        B = colorbar;
        colormap
        ylabel(B,'SL [dB re 1 \muPa//Hz @ 1m]')
        axis xy
        xlim([5 1000])
        %caxis([80 130])
        %caxis([60 110])
        %caxis([-20 20])
        title([CINMS_B_30(RecInd).label '-' CINMS_B_30(RecInd2).label],'Interpreter','None')
        


%%

figure;
    imagesc(10*log10(TL{1})-10*log10(TL{2})); axis xy; colorbar
    xlim([0 1000])
    