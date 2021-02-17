
clear all
close all

%MMSI = 220397000 % Grete
%MMSI = 220379000 % Gudrun
%MMSI = 220413000 % Gunvor
%MMSI = 220414000 % Gjertrud
%MMSI = 220415000 % Gerd
%MMSI = 220416000 % George
%MMSI = 220593000 % Gerner
%MMSI = 220594000 % Gunde
MMSI = 220595000 % Gunhilde
%MMSI = 220596000 % Gustav
%MMSI = 220597000 % Guthorm
%MMSI = 220598000 % Gerda

StartTime = datenum('11/25/2015 13:50:00');
EndTime = datenum('11/25/2015 14:10:00');


load 'G:\Ch.2_MAERSK_Retrofit\AIS\COP_AIS_Data.mat'
load 'G:\MartinGassman\MartinGassman_Ship\SBC\CINMS_B_30_ArrayDefinition.mat'
load 'G:\Ch.2_MAERSK_Retrofit\CINMS_B_DepInfo.mat'

Index = find(contains(CINMS_B_DepInfo.names,'CINMS29-B'));
tffn = char(CINMS_B_DepInfo.tf(1, Index))
draft = 12.08
AISformat = 'COP'

%% Initialize 
Stern2EngRoom = 60;                       %Distance from Stern to engine room [m]
ShipStandard = 'ANSI'                     
Num1secSeg = 5;                             %number of segments for averaging
TimeIntSec = 3;                             %time-difference [sec] btw. beginning of segments
fstartind = 6;                              %10Hz for computing received levels
fstopind = 1001;                            %1000Hz for computing received levels
c = 1488.360552;                             %sound speed [m/s]

if StartTime < 736395
    propDiam = 9  %propeller diameter
else
    propDiam = 9.3;
end

zs = draft - (.85*propDiam)


%% Acoustic Data Files

% d20 file path for each acoustic recorder
filepath  = 'D:\Ch.2_MAERSK_Retrofit\CINMS_B_xwavs\'
%build cell array with filenames
ftype = '*.wav'; % Find all .txt files
listing = dir([filepath, ftype]); 
fn = {listing.name};
fn = char(fn)
filename = {};
for name = 1:size(fn, 1)
    filename{name} = [filepath fn(name, :)];
end


%% Pre-process Input parameters
Fs = 10000;              %df20 data
SegLengthSamples = 1*Fs; % for 1 Hz bins
DataLengthSamples = Num1secSeg*SegLengthSamples; %the length of the data for averaging in samples

switch ShipStandard 
    case 'ANSI'
        ShipStandardLabel = 'ANSI S12.64-2009';
    case 'ISO'
        ShipStandardLabel = 'ISO 17208-1:2016';
end 

%% Transfer Function for HARP/Deployment




tfdir = 'G:\TFs_08282020';          % path to TF files
                                    % manual input for now
fid = fopen(fullfile(tfdir,tffn),'r');
[A,~] = fscanf(fid,'%f %f',[2,inf]);
TFf = A(1,:);
TFdb = A(2,:);
fclose(fid);
TFflag = 1;
[~,ia,ic] = unique(TFf);
if length(ia) == length(ic)
    freq = TFf;
    uppc = TFdb;
else
    freq = TFf(ia);
    uppc = TFdb(ia);
end

%interpolating for everyon 1 Hz & converting transfer function from dB to uPa
F = 0:5000
Ptf = interp1(freq,uppc,F,'linear','extrap');
G = 10.^(Ptf/20)% transfer function in uPa 
G = G' 


%converting transfer function from dB to uPa Martin's way
%fRspectra =  Fs/2*linspace(0,1,SegLengthSamples/2+1);
% 
%         for findex=1:length(fRspectra)
%             [~,minindex] = ...
%                 min(abs(freq-abs(fRspectra(findex))));
%             Gtest(findex, 1) = 10.^(uppc(minindex)/20);
%         end



%% Get File Start Times from File Headers or File Name
    %Get x.wav file Start Times
    %xwavReferenceTime is the start time for each x.wav file
for FileIndex=1:length(filename)         
                %Get Timestamp of first sample from xwavheader
                %==============================================
                % timestamp of first sample = timestamp of first raw file    

                %open file
                fid = fopen(filename{FileIndex},'r');
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
                xwavReferenceTime(FileIndex) = datenum(xwavStartTimeVector);
end
%% Time stamps (decimal days) for SL compuations

TimeVector = StartTime:TimeIntSec/(3600*24):EndTime; %start time to end time in 3 second chunks

%Total number of Events
NumSpec = length(TimeVector); %how many 3 second chunks
ReferenceTime = StartTime;    %start time of the transit
%% Compute spectra (and track acoustically)

% Received Spectra with 1 Hz bins
% S will be the Received Spectra in 1 Hz bins for each 3 second chunk for
% the duration of the transit
S = zeros(NumSpec,SegLengthSamples/2+1); %3 sec chunks for 1 Hz bins for duration of transit
for EventIndex = 1:NumSpec
        
        %Load data
        %========================================================          
        %choose correct data file
        FileIndex = ...
            find((xwavReferenceTime(:) - TimeVector(EventIndex))<0); %finding file that includes start and end time

         
                StartSample = round(Fs*(TimeVector(EventIndex) - ...
                    xwavReferenceTime(FileIndex(end)))*(24*3600)); % the start sample of the transit in the .x.wav that the transit is in
                EndSample = StartSample+DataLengthSamples-1; % 5 seconds away from start sample including 0-10000 Hz in 1 Hz bins
                siz = audioinfo(filename{FileIndex(end)}); %reading the entire .x.wav that the transit is in
                siz = siz.TotalSamples;
                if EndSample <= siz(1)
                    data = double(audioread(filename{FileIndex(end)},...
                        [StartSample EndSample], 'native'));
                else
                    disp('SoundTrap Data File End. Getting data from next file.')
                    data = double(audioread(filename{FileIndex(end)},...
                        [StartSample siz(1)], 'native'));
                    data = [data double(audioread(filename{FileIndex(end)+1},...
                        [1 DataLengthSamples-(siz(1)-StartSample)], 'native'))];
                end    
         
        %Compute received spectra
        %==============================
                
        %In case is less than anticipated due to end of file
        NumSeg = floor(length(data)/SegLengthSamples); % should be 5 if everything good
        %Fourier Transform    
            X = fft(reshape(data(1:NumSeg*SegLengthSamples),...
                SegLengthSamples,NumSeg)); %reshapes to 5 second length 1-1000 Hz and takes fft

            %compute onesided spectra and calibrate with HARP specific hydrophone transfer functions
            S(EventIndex,:) = ...
                (sum(2*abs(X(1:SegLengthSamples/2+1,:)).^2,2).*G.^2)...
                /(NumSeg*SegLengthSamples*Fs);            
end


%% Spectrograms

FontSizeNumber=14;
%for RecInd=1:NumRec
figure;
    imagesc(TimeVector,F,10*log10(S).');
    axis xy; 
    h = colorbar;caxis([50 110])
    set(get(h,'ylabel'),'string','RL [dB re 1\muPa//Hz]',...
        'FontSize',FontSizeNumber);
    datetick('x','keeplimits');
    xlabel(['GMT on ' datestr(TimeVector(1),'dd-mmm-yyyy')],...
        'FontSize',FontSizeNumber)
    ylabel('Frequency [Hz]','FontSize',FontSizeNumber)
    ylim([0 1000])
    %vline(AIS.dnums(AISCPAInd), 'black', '-')
    %vline(TimeVector(PropCPAInd), 'red', '-')
    set(gca, 'FontSize',FontSizeNumber)


%% Average Spectra

figure;
semilogx(F,10*log10(sum(S)/NumSpec))
ylim([60 100])
xlim([5 1000])
grid on
ylabel('RL [dB re 1\muPa^2//Hz]')
xlabel('f [Hz]')


%%
%% Process AIS data

%find ship passage in AIS data (give +/- 3 min extra time for smoothing)
AIS_MMSI = [shipTracks.MMSI];
AIS_MMSI_idx = ismember([shipTracks.MMSI], MMSI);
AIS = (shipTracks(AIS_MMSI_idx));
StartTimeHour = datenum(dateshift(datetime(datestr(StartTime)), 'start', 'hour'));

clear indT
indT = {};
for i = 1:length(AIS)
    indT = AIS(i).dnums > StartTime & AIS(i).dnums < EndTime;
    sumindT(i,1) = sum(indT)
    if sumindT(i) > 0 
        AISTime = AIS(i)
    end
    
end

AIS = AISTime

NumAISlocs = length(AIS.dnums);

%map AIS lats/longs to x-y    
[AIS.x, AIS.y] = latlon2xy(AIS.lats,AIS.lons,...
        CINMS_B_DepInfo.lats(Index),CINMS_B_DepInfo.lons(Index));
%[CINMSB.x, CINMSB.y] = latlon2xy(CINMS_B_DepInfo.lats(Index), CINMS_B_30(12).long,...
    %CINMS_B_30(12).lat, CINMS_B_30(12).long);

%fix potential errors of time stamp due to AIS computer clock drift
%==================================================================

%determine closest point of approach by finding the maximum RL at propeller blade line
PropFreqLineInd = 7; %8Hz  blade rate fundamental (orginates from propeller)
PropCPAInd = zeros(1,1);
AISCPAInd = zeros(1,1);

[~,PropCPAInd]=max(S(:,PropFreqLineInd));
[AISCPA,AISCPAInd]=min(sqrt((AIS.x).^2 + (AIS.y).^2));


 
[PropCPAInd AISCPAInd (TimeVector(PropCPAInd)-AIS.dnums(AISCPAInd)).'*24*60]% difference in CPA in minutes

%determine correct draft
[~,DraftIndex] = min(abs((StartTime+(EndTime-StartTime)/2) - AIS.vData(1,3)));

%plot AIS data
% h=plotShips({AIS},1);
% FontSizeNumber=16;
% set(gca, 'FontSize',FontSizeNumber)
% xlabel('Eastings [km]','FontSize',FontSizeNumber)
% ylabel('Northings [km]','FontSize',FontSizeNumber)    
% set(gca,'XTick',[-6:3:6])

%plot time-resolution of AIS data
figure;
    subplot(2,1,1)
    plot(AIS.dnums(1:end-1),diff(AIS.dnums*24*3600),'.');
    datetick('x')
    ylabel('\DeltaT [sec]')    
subplot(2,1,2)
    plot(AIS.dnums,AIS.SOG,'.');
    datetick('x')
    ylabel('Ship speed [m/s]')
    xlabel('HH:MM')    
    
%compute ship length and width 
% COP total length in 4 and total width in 5
% COP GPS2port is 7 COP GPS2stern is 6
% SBARC total length is 4 + 5 and total width is 6 + 7

AIS.length = 366.89;
AIS.width = 42.94;  
AIS.GPS2port = AIS.vData(1,7)

if strcmp(AISformat, 'COP')
    AIS.GPS2stern = AIS.length - AIS.vData(1,6);
elseif strcmp(AISformat, 'SBARC')
    AIS.GPS2stern = AIS.vData(1,5);
elseif strcmp(AISformat, 'SCI')
    AIS.GPS2stern = AIS.length - AIS.vData(1,6);
end


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
    [dt dtInd] = min(abs(AIS.dnums-TimeVector(TimeInd)));
        if dt < MaxAISdt
            if dtInd >= 3 && dtInd <= NumAISlocs-2
                p.x = polyfit(AIS.dnums(dtInd-2:dtInd+2),AIS.x(dtInd-2:dtInd+2),1);
                AIS.xint(TimeInd) = polyval(p.x,TimeVector(TimeInd));                 
                p.y=polyfit(AIS.dnums(dtInd-2:dtInd+2),AIS.y(dtInd-2:dtInd+2),1);
                AIS.yint(TimeInd) = polyval(p.y,TimeVector(TimeInd)); 
                p.lats=polyfit(AIS.dnums(dtInd-2:dtInd+2),AIS.lats(dtInd-2:dtInd+2),1);
                AIS.latsint(TimeInd) = polyval(p.lats,TimeVector(TimeInd));
                p.lons=polyfit(AIS.dnums(dtInd-2:dtInd+2),AIS.lons(dtInd-2:dtInd+2),1);
                AIS.lonsint(TimeInd) = polyval(p.lons,TimeVector(TimeInd));
                p.s = polyfit(AIS.dnums(dtInd-2:dtInd+2),AIS.SOG(dtInd-2:dtInd+2),1);
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
AIS.vxraw = diff(AIS.x)./(diff(AIS.dnums)*24*3600); %m/s
AIS.vyraw = diff(AIS.y)./(diff(AIS.dnums)*24*3600); %m/s
AIS.vraw = sqrt(AIS.vxraw.^2 + AIS.vyraw.^2); %m/s

figure;
subplot(3,1,1)
    plot(TimeVector,AIS.xint,'.r');
    hold on
    plot(AIS.dnums,AIS.x,'.b');    
    datetick('x','keeplimits');
    xlim([TimeVector(1) TimeVector(end)])
    xlabel(['GMT on ' datestr(TimeVector(1),'dd-mmm-yyyy')]);    
    vline(AIS.dnums(AISCPAInd), 'black', '-')
    vline(TimeVector(PropCPAInd), 'red', '-')
subplot(3,1,2)
    plot(TimeVector,AIS.yint,'.r');
    hold on
    plot(AIS.dnums,AIS.y,'.b');    
    datetick('x','keeplimits');
    xlim([TimeVector(1) TimeVector(end)])
    xlabel(['GMT on ' datestr(TimeVector(1),'dd-mmm-yyyy')]);    
    vline(AIS.dnums(AISCPAInd), 'black', '-')
    vline(TimeVector(PropCPAInd), 'red', '-')
subplot(3,1,3)
    plot(TimeVector,AIS.speedint,'.c');
    hold on
    plot(AIS.dnums,AIS.SOG,'.b');  
    plot(TimeVector(1:end-1),AIS.v,'.m')
    plot(AIS.dnums(1:end-1),AIS.vraw,'.r')
    datetick('x','keeplimits');
    xlim([TimeVector(1) TimeVector(end)])
    xlabel(['GMT on ' datestr(TimeVector(1),'dd-mmm-yyyy')]);    
    legend('AIS int.','AIS raw','AIS int. track','AIS raw track')
    vline(AIS.dnums(AISCPAInd), 'black', '-')
    vline(TimeVector(PropCPAInd), 'red', '-')
    ylim([0 15])
    grid on

%%
alpha = zeros(1,NumSpec-1);        %horizontal angle [deg]
phi = zeros(1,NumSpec-1);          %vertical angle [deg]
R = zeros(1,NumSpec-1);            %slant range [m]
RLship = zeros(1,NumSpec-1);       %broad band received level [dB]
SL.broadband = zeros(1,NumSpec-1); %broadband source level [dB]
SL.ALL = cell(1,1);                %source spectra level [dB]

 
for index=1:NumSpec-1
    
    %compute beam angles alpha and phi
%     thetax=0; % roll angle of Zc =  rotates beampattern around (0deg,0deg)
%     thetay = 0;
%     thetaz = atan2(AIS.vy(index),AIS.vx(index));
%     Rx = [1 0 0; 0 cos(thetax) sin(thetax); 0 -sin(thetax) cos(thetax)];
%     Ry = [cos(thetay) 0 sin(thetay); 0 1 0; -sin(thetay) 0 cos(thetay)];
%     Rz = [cos(thetaz) sin(thetaz) 0; -sin(thetaz) cos(thetaz) 0; 0 0 1];
%     b = 0 - [AIS.xint(index) AIS.yint(index) CINMS_B_30(12).depthHydro];
%     bprime = Rx*Ry*Rz*b.';
%     
%     alpha(1,index) =180*atan2(bprime(2),bprime(1))/pi;
%     phi(1,index) = 180*atan2(bprime(3),norm(bprime(1:2)))/pi;
    
    %slant range of each ship location
    R(1,index) = sqrt((AIS.xint(index)^2)+  (AIS.yint(index)^2) + CINMS_B_30(12).depthHydro^2);
    
    %compute source level from spectrum for beam angle
    RLship(1,index) = 10*log10(trapz(F(fstartind:fstopind),...
        (S(index,fstartind:fstopind))));%...
    SL.broadband(1,index) = ...
        RLship(1,index) + 20*log10(R(1,index));
    SL.ALL{index,:} = ...
        10*log10(S(index,:))+20*log10(R(1,index));
end

%% Indicies of good SL and angles (for plots, further processing)

% %bad/no AIS information
% PlotIndicies = find(AIS.xint(1:end-1));
% 
% %weared angles or source levels
% dummy = find(alpha <=3.5 & alpha >=-3.5 &  phi >=-4);
% [row col] = ind2sub([1,NumSpec-1],dummy);
% PlotIndicies=setdiff(PlotIndicies,col);


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
% fintind = 501:601; % OK
% %fintind = 401:601; % OK
% 
% if length(fintind) == 1
%     GMTfilename = [MMSI 'RadPatt' ShipStandard '_' num2str(F(fintind)) 'Hz3m'];
% else
%     GMTfilename = [MMSI 'RadPatt' ShipStandard '_' num2str(F(fintind(1))) '_'...
%                     num2str(F(fintind(end))) 'Hz3m'];
% end

% Export  Beampattern  
% RadPatGMT=cell(1,1);
% MaxSL=zeros(1,1);
% 
%     %find good angles
%     dummy = find(abs(alpha(1,:)) <=14 &  phi(1,:) >=-10);    
%     RadPattIndicies{1}=setdiff(PlotIndicies,dummy);
%     dummy = find(abs(alpha(1,:)) >=160 &  phi(1,:) >=-10);    
%     RadPattIndicies=setdiff(RadPattIndicies,dummy);    
%     MaxSL = max(10*log10((sum(10.^(SL.ALLLloyd{1}(RadPattIndicies{1},fintind)/10),2))));
% 
 %MaxSL = max(MaxSL)
%% Compute Source Levels according to ANSI/ASA S12.64-2009

%Compute Data Window Length (DWL) in meter from Eq.1
% Ind0115 = find(strcmp({CINMS_B_30.label},'CINMS_B_30_0115')); %STBD
% Ind0515 = find(strcmp({CINMS_B_30.label},'CINMS_B_30_0515')); %PORT
[RCPA RCPAind] = min(R);
%RCPAind = RCPAind - 380
dCPA = sqrt(RCPA^2 + CINMS_B_30(12).depthHydro.^2);
DWL = 2*dCPA*tand(30); % using the distance to CPA
DWL2 = 2*366.89*tand(30); % using the length of the ship

%Data Window Period (DWP) in sec from Eq. 2
DWP.length = DWL/AIS.v(RCPAind);
DWP.length2= DWL2/AIS.v(RCPAind);

[~,DWP.StartInd]=min(abs(TimeVector-(TimeVector(RCPAind)...
    -DWP.length/(2*3600*24))));
[~,DWP.EndInd]=min(abs(TimeVector-(TimeVector(RCPAind)...
    +DWP.length/(2*3600*24))));

[~,DWP.StartInd2]=min(abs(TimeVector-(TimeVector(RCPAind)...
    -DWP.length2/(2*3600*24))));
[~,DWP.EndInd2]=min(abs(TimeVector-(TimeVector(RCPAind)...
    +DWP.length2/(2*3600*24))));


% Source Level with spherical spreading
    TL = R(1, DWP.StartInd2:DWP.EndInd2).^2; %spherical spreading TL in linear space
    TL = TL(ones(length(F),1),:).'; % same TL for every frequency is the same (could change to frequency dependent)
    %TL = 20*log10(TL)
    SL = 10*log10(...
        sum(S(DWP.StartInd2:DWP.EndInd2,:).*TL)...
        ./(DWP.EndInd2-DWP.StartInd2+1));% applying TL to RL and taking the average over the DWP, converting to dB


    
%ANSI/ASA/ISO broadband source level
SL_TotalBroadband = 10*log10(trapz(F(fstartind:fstopind),...
            10.^(SL(fstartind:fstopind)/10)))
figure(6)       
semilogx(F, SL)
grid on
xlim([5 1000])
ylim([120 190])
ylabel('Spectrum Level (dB re 1uPa^2/Hz')
xlabel('Frequency Hz')

FontSizeNumber=14;
%for RecInd=1:NumRec
figure(10);
    imagesc(TimeVector,F,10*log10(S).');
    axis xy; 
    h = colorbar;caxis([50 110])
    set(get(h,'ylabel'),'string','RL [dB re 1\muPa//Hz]',...
        'FontSize',FontSizeNumber);
    datetick('x','keeplimits');
    xlabel(['GMT on ' datestr(TimeVector(1),'dd-mmm-yyyy')],...
        'FontSize',FontSizeNumber)
    ylabel('Frequency [Hz]','FontSize',FontSizeNumber)
    ylim([0 1000])
    vline(TimeVector(RCPAind), 'black', '-')
    vline(TimeVector(PropCPAInd), 'red', '-')
    set(gca, 'FontSize',FontSizeNumber)
        
    %% Compute SL with Lloyd's mirror TL


%source depth [m]

%SL.ALLLloyd = cell(1,1);
% TL = cell(1,1);
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
k = 2*pi./(c./F); %wavenumber
Rc = ones(length(k),1); %reflection coefficient
% Rc(60:80) = linspace(1,0.9,21); 
% Rc(80:end) = 0.9; 
Rc = Rc.';

Zs = zs
Zr = 565
for n = 1:length(R)
    for k = 1:length(F)
        TLLM(k, n) = 10*log10((R(n).^2)/(2*(1-(cos((4*pi*F(k)*Zs*Zr)/(c*R(n)))))));
    end
end

Sal = 34.2481 %salinity
pH = 7.183 %pH at hydrophone depth from CalCofi
T = 6.6477
f0 = F(1)
f1 = F(end)
r = sqrt((AIS.xint.^2)+  (AIS.yint.^2));
zr = Zr

for i = 1:length(R)
    TL_SS(:, i) = range2SS_TL(f0, f1, R(i), zs, zr, pH, T, Sal, c);
end
TLSS = TL_SS;
%%%finding location of intersection
ModTLLM = abs(TLLM-TLSS);
ModTLLM = 1./ModTLLM;

for star = 1:size(ModTLLM,2)
    [~, indtemp] = findpeaks(ModTLLM(:, star));
    index(star) = min(indtemp);
end


%%USE USE USE USE USE USE
%getting TLLM/TLSS for each transit
TLLMFinal=[];

for iTest=1:size(TLLM,2)
    tlmmTemp = TLLM(1:index(iTest), iTest);
    stS = index(iTest)+1;
    endS = size(TLSS,1);
    tlssTemp = TLSS(stS:endS, iTest);
    TLLMFinal{iTest} = [tlmmTemp;tlssTemp];
end


TLLMFinal = cell2mat(TLLMFinal);


%     zr = 565;           
%     for Rind = 1:length(R)        
%         %compute acoustic pressure p
%         R1 = sqrt((R(1,Rind).^2-zr.^2) + (zr-zs).^2);
%         R2 = sqrt(R(1,Rind).^2-zr.^2 + (zr+zs).^2);         
%         TL(Rind,:) = abs(1*exp(1i*k.*R1)./R1 - Rc.*exp(1i*k.*R2)./R2).^2;
%         %curb Lloyds mirrors curve when hit 20*log10R to avoid interference
%         %lobes
%         findex=1;
%         while TL(Rind,findex)<(1/R(1,Rind)^2)
%             findex = findex+1;
%             if findex == length(F); break; end;
%         end                    
%         if findex <= length(F)
%             TL(Rind,findex:end) = 1/R(1,Rind)^2;
%         end            
%         SL.ALLLloyd(Rind,:) = 10*log10(S(Rind,:)./TL(Rind,:));
%     end

    
figure;
    imagesc(TLLM.'); axis xy; colorbar
    semilogx(F,TLLM.');xlim([5 1e3])      


    
    
TLLMFinal = (10.^(TLLMFinal(:, DWP.StartInd2:DWP.EndInd2)/10)); %TLLM in linear space
SL_TLLM = 10*log10(...
        sum(S(DWP.StartInd2:DWP.EndInd2,:).*TLLMFinal')...
        ./(DWP.EndInd2-DWP.StartInd2+1))
SL_TLLM_TotalBroadband = 10*log10(trapz(F(fstartind:fstopind),...
            10.^(SL_TLLM(fstartind:fstopind)/10)))

figure(3)     
semilogx(F, SL);
hold on
semilogx(F, SL_TLLM);
grid on
xlim([5 1000])
ylim([120 190])
ylabel('Spectrum Level (dB re 1uPa^2/Hz)')
xlabel('Frequency Hz')

AIS.CPA = min(R)
AIS.latlonCPA = [AIS.latsint(RCPAind) AIS.lonsint(RCPAind)]

SLfinal.SLSS = SL;
SLfinal.SLSSBB = SL_TotalBroadband;
SLfinal.SLLM = SL_TLLM;
SLfinal.SLLMBB = SL_TLLM_TotalBroadband;
SLfinal.RL = RLship;
SLfinal.Spectra = S;
SLfinal.StartDWP = DWP.StartInd2;
SLfinal.EndDWP = DWP.EndInd2;
SLfinal.CPA = min(R)
SLfinal.CPAcoords = [AIS.latsint(RCPAind) AIS.lonsint(RCPAind)]
SLfinal.SOG = AIS.speedint(RCPAind)
SLfinal.sourceDepth = zs
SLfinal.retrofit = 'post'
SLfinal.datenum = AIS.dnums(AISCPAInd)

formatout = 'dd_mmm_yyyy'
ofn = [num2str(MMSI) '_' datestr(AIS.dnums(AISCPAInd), formatout) '.mat']

save(ofn, 'AIS', 'SLfinal')

%%


% 
% %Kuperman book
% k = 2*pi./(c./F); %wavenumber
% Rc = ones(length(k),1); %reflection coefficient
% % Rc(60:80) = linspace(1,0.9,21); 
% % Rc(80:end) = 0.9; 
% Rc = Rc.';
% 
% 
%     zr = CINMS_B_30(12).depthHydro;           
%     for Rind = 1:length(R)        
%         %compute acoustic pressure p
%         R1 = sqrt((R(Rind).^2-zr.^2) + (zr-zs).^2);
%         R2 = sqrt(R(Rind).^2-zr.^2 + (zr+zs).^2);         
%         TL(Rind,:) = abs(1*exp(1i*k.*R1)./R1 - Rc.*exp(1i*k.*R2)./R2).^2;
%         %curb Lloyds mirrors curve when hit 20*log10R to avoid interference
%         %lobes
%         findex=1;
%         while TL(Rind,findex)<(1/R(Rind)^2)
%             findex = findex+1;
%             if findex == length(F); break; end;
%         end                    
%         if findex <= length(F)
%             TL(Rind,findex:end) = 1/R(Rind)^2;
%         end            
%         SL.ALLLloyd(Rind,:) = 10*log10(S(Rind,:)./TL(Rind,:));
%     end
% 
% figure
% semilogx(F, -10*log10(TL), 'black')
% hold on
% semilogx(F, TLLMFinal, 'green')
% grid on
% ylabel('Transmission Loss (dB re 1uPa^2)')
% xlabel('Frequency (Hz)')
% legend('Martins Code', 'Vanessas Code')
% 
% 
