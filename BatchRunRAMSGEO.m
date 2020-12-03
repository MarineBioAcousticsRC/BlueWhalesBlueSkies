%% BatchRun RAM(S)geo 
%
% Martin Gassmann, 2-June-2017
%

%% Inilize

clear all
close all

ProjectName = 'SBC_2016'
f = [20];
%f = 200:250;
zs = 151; %m
zr = 5; %m (doesn't matter)

RAM.path = 'G:\MartinGassman_PropagationModels\RAM\RAMSGeo\';
RAM.type = 'RAMSGEO.exe';


%% Post-Initialize

RAM.outputfilename.pgrid = 'p.grid'
switch RAM.type
    case 'RAMSGEO.exe'
        RAM.inputfilename = 'ramsgeo.in';
    case 'RAMGEO.exe'
        RAM.inputfilename = 'ramgeo.in';
end        

FileNamePrefix=ProjectName;

%make Project directory
shellcmd=['mkdir ' ProjectName];
dos(shellcmd);

%% Run RAM

%Line in ram(s)geo.in that has frequency
replaceLineFreq = 2;

for findex=1:length(f)
    %adjust frequency in input file
    fid = fopen(RAM.inputfilename,'r+');
    for k=1:(replaceLineFreq-1);
        tline=fgetl(fid);
    end;
    fseek(fid, 0, 'cof'); % call fseek between read and write operations
    fprintf(fid, [sprintf( '%5f', f(findex))  sprintf('\t')...
        sprintf( '%5f',zs)  sprintf('\t') sprintf( '%5f',zr)  sprintf('\t')...
         sprintf('f, zs, zr')]); % print the new values
    fclose(fid);
    %run RAM
    shellcmd=[RAM.type ' ' RAM.inputfilename];
    dos(shellcmd);
    %rename output files
    shellcmd=['rename ' RAM.outputfilename.pgrid ' ' FileNamePrefix num2str(f(findex)) 'Hz_' RAM.outputfilename.pgrid];
    dos(shellcmd);
    %move output files
    shellcmd=['move ' FileNamePrefix num2str(f(findex)) 'Hz_' RAM.outputfilename.pgrid ' ' RAM.path ProjectName];
    dos(shellcmd);
    %copy ram(s)geo.in
    shellcmd=['copy ' RAM.inputfilename ' ' RAM.path ProjectName];
    dos(shellcmd);            
    shellcmd=['rename ' RAM.path ProjectName '\' RAM.inputfilename ' ' FileNamePrefix num2str(f(findex)) 'Hz_' RAM.inputfilename];
    dos(shellcmd);            
end    


%% Extract TL(r,z) from RAM(S)Geo.exe output file p.grid or AcTup output file xxx.pgrid
% Uses modified ReadRamPGridmod.m from Bruce, not the one that comes with ActUP
% complex pressure is for unity source strength 
% tested&works with RAMSGeo.exe (v0.5C01.01.01) and RAMGeo.exe (v1.5C00.03)
% that were compiled by ActUp people

%Read in Complex Pressure Grid from RAMSgeo output file (p.grid)
PGrid = ReadRamPGridmod( [RAM.path ProjectName '\' FileNamePrefix num2str(f(findex)) 'Hz_' RAM.outputfilename.pgrid] );

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
    %plot(rr,858*ones(nrr,1),'-k','LineWidth',5);
    %title([pltitl '\newline f=' num2str(freq) ' Hz']);
    xlabel('Range [m]')
    ylabel('Depth [m]')    
    t=colorbar;
    test=flipud(colormap('jet'));
    colormap(test);
    set(get(t,'ylabel'),'String', ['\fontsize{10} TL [dB]']);        
    caxis([15 75])