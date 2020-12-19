%% BatchRun RAM(S)geo with bathymetry
%
% Martin Gassmann, 2-June-2017
% Adapted and rewritten for models with bathymetry by Eric Snyder, 3-Dec-2020
% Fixed output directory bug and file name Vanessa ZoBell, 17-Dec-2020

%% Inilize

clear all
close all

ProjectName = 'Test_Bath'
f = [5:1000];
%f = 200:250;

hydloc = [34.2755, -120.0185, 565]; % hydrophone [lat, lon, depth]
sloc = [34.305, -120.01, 5]; % source [lat, lon, depth]

% [sx, sy] = latlon2xy(sloc(1), sloc(2), hydloc(1), hydloc(2));

zs = sloc(3); %m
zrr = 5; %m (doesn't matter)

RAM.path = 'C:\Users\HARP\Desktop\Code\ramsgeoWithBath-main';
RAM.type = 'RAMSGEO.exe';


%% Post-Initialize

RAM.outputfilename.pgrid = 'p.grid'
switch RAM.type
    case 'RAMSGEO.exe'
        RAM.inputfilename = 'ramsgeo.in';
    case 'RAMGEO.exe'
        RAM.inputfilename = 'ramgeo.in';
end        

fpn = fullfile(RAM.path, RAM.inputfilename);
FileNamePrefix=ProjectName;

%make Project directory
shellcmd=['mkdir ' ProjectName];
dos(shellcmd);

%% Bathymetry and SSP
% Extract bathymetry between source and receiver
[r, bty] = extractBTY(sloc(1), sloc(2), hydloc(1), hydloc(2));

bty = fliplr(bty)

load('sspTest.mat');
z = ssp(:,1);
c = ssp(:,2);
%% Run RAM



for findex=1:length(f)
    
    make_in_file(FileNamePrefix, fpn, f(findex), zs, zrr, r, bty, z, c)
    
    %run RAM
    shellcmd=[RAM.type ' ' RAM.inputfilename];
    dos(shellcmd);
    %rename output files
    shellcmd=['rename ' RAM.outputfilename.pgrid ' ' FileNamePrefix '_' num2str(f(findex)) 'Hz_' RAM.outputfilename.pgrid];
    dos(shellcmd);
    %move output files
    shellcmd=['move ' FileNamePrefix '_' num2str(f(findex)) 'Hz_' RAM.outputfilename.pgrid ' ' RAM.path '\' ProjectName];
    dos(shellcmd);
    %copy ram(s)geo.in
    shellcmd=['copy ' RAM.inputfilename ' ' RAM.path '\' ProjectName];
    dos(shellcmd);            
    shellcmd=['rename ' RAM.path '\' ProjectName '\' RAM.inputfilename ' ' FileNamePrefix num2str(f(findex)) 'Hz_' RAM.inputfilename];
    dos(shellcmd);            
end    


%% Extract TL(r,z) from RAM(S)Geo.exe output file p.grid or AcTup output file xxx.pgrid
% Uses modified ReadRamPGridmod.m from Bruce, not the one that comes with ActUP
% complex pressure is for unity source strength 
% tested&works with RAMSGeo.exe (v0.5C01.01.01) and RAMGeo.exe (v1.5C00.03)
% that were compiled by ActUp people

%Read in Complex Pressure Grid from RAMSgeo output file (p.grid)
PGrid = ReadRamPGridmod( [RAM.path '\' ProjectName '\' FileNamePrefix '_' num2str(f(findex)) 'Hz_' RAM.outputfilename.pgrid] );

% Plot
figure;
    %raw
    pcolor(-20*log10(abs(PGrid)));
    axis ij
    %matlab interpolated
    pcolor(-20*log10(abs(PGrid)));
    axis ij
    shading interp;
    hold on    
    hold on
    plot(3300, 565, 'o')
    %plot(rr,858*ones(nrr,1),'-k','LineWidth',5);
    %title([pltitl '\newline f=' num2str(freq) ' Hz']);
    xlabel('Range [m]')
    ylabel('Depth [m]')    
    t=colorbar;
    test=flipud(colormap('jet'));
    colormap(test);
    set(get(t,'ylabel'),'String', ['\fontsize{10} TL [dB]']);        
    caxis([15 75])