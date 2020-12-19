%% Frequency dependent TL from RAMSGeo model
% Need the output from BathRunRAMSGEO_withBathymetry.m 
% (both .in and p.grid files)

% Vanessa ZoBell 18-Dec-2020
% MATLAB version 2016b
%% Initialize
clear all 
close all 


% dir where output from BathRunRAMSGEO_withBathymetry.m is stored
idir = 'C:\Users\HARP\Desktop\Code\ramsgeoWithBath-main\Test_Bath'; %

% depth of receiver
hd = 565 



% reads in .in files to get range to receiver and frequency for specific
% model
inFile = dir(fullfile(idir, '*_ramsgeo.in')) 
Freq = zeros(1, length(inFile))
for i = 1:length(inFile)
    inFileName = inFile(i).name;
    fullinFileName = fullfile(idir, inFileName);
    fid = fopen(fullinFileName);
    linenumRange = 3;
    range = floor(cell2mat(textscan(fid, '%f', 1, 'delimiter', '\n', 'headerlines', linenumRange-1))); %range of source to receiver
    frewind(fid)
    linenumFreq = 2;
    Freq(i) = cell2mat(textscan(fid, '%f', 1, 'delimiter', '\n', 'headerlines', linenumFreq-1)); %frequency used in model
    frewind(fid)
    fclose(fid)
    clear inFileName
    clear fullinFileName
end

% reads p.grid files and converts to frequency dependent TL
pgridFiles = dir(fullfile(idir, '*.grid'));   
FreqDep_TL = zeros(1, length(pgridFiles));
for k = 1:length(pgridFiles)
    baseFileName = pgridFiles(k).name;
    fullFileName = fullfile(idir, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    Pgrid = ReadRamPGridmod(fullFileName);
    TL = -20*log10(abs(Pgrid));
    FreqDep_TL(k) = TL(hd, range);
    clear baseFileName
    clear fullFileName
end

Freq_TL = [Freq; FreqDep_TL]';


semilogx(Freq_TL(:, 1), Freq_TL(:, 2), '.'); 
set(gca, 'YDir', 'reverse')
xlim([5 1000])
grid on
ylabel('Transmission Loss (dB)')
xlabel('Frequency (Hz)')


% 
% Testing stuff...
% 
% 
% pgridBath = 'C:\Users\HARP\Desktop\Code\ramsgeoWithBath-main\Test_SSP\Test_SSP_400Hz_p.grid'
% pgridFlat = 'G:\MartinGassman\MartinGassman_PropagationModels\RAM\RAMSGeo\SBC5m_VZ_Jan\SBC5m_VZ_Jan400.5Hz_p.grid'
% 
% TLBath = ReadRamPGridmod(pgridBath);
% TLBath = -20*log10(abs(TLBath));
% TLFlat = ReadRamPGridmod(pgridFlat);
% TLFlat = -20*log10(abs(TLFlat));
% 
% TLFlatdot = TLFlat(565, 3365)
% TLBathdot = TLBath(565, 3365)
% 
% figure(1);
%     %raw
%     pcolor(TLBath);
%     axis ij
%     %matlab interpolated
%     pcolor(TLBath);
%     axis ij
%     shading interp;
%     hold on    
%     plot(3300, 565, 'o', 'y')
%     %plot(rr,858*ones(nrr,1),'-k','LineWidth',5);
%     %title([pltitl '\newline f=' num2str(freq) ' Hz']);
%     xlabel('Range [m]')
%     ylabel('Depth [m]')    
%     t=colorbar;
%     test=flipud(colormap('jet'));
%     colormap(test);
%     set(get(t,'ylabel'),'String', ['\fontsize{10} TL [dB]']);        
%     caxis([15 75])
% 
% figure(2);
%     %raw
%     pcolor(TLFlat);
%     axis ij
%     %matlab interpolated
%     pcolor(TLFlat);
%     axis ij
%     shading interp;
%     hold on    
%     plot(3300, 565, 'o')
%     %plot(rr,858*ones(nrr,1),'-k','LineWidth',5);
%     %title([pltitl '\newline f=' num2str(freq) ' Hz']);
%     xlabel('Range [m]')
%     ylabel('Depth [m]')    
%     t=colorbar;
%     test=flipud(colormap('jet'));
%     colormap(test);
%     set(get(t,'ylabel'),'String', ['\fontsize{10} TL [dB]']);        
%     caxis([15 75])
