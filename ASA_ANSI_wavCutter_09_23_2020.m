clear all
fpr = 'G:\Ch.1_ShipProject\CINMS_B_Shipping\COP\2015-05\'; % path to folder containing files
fpw = 'G:\CPA_wav_ASA\fixedTimeOffset_11_23_2020\';% path to folder where you want to save shorter .wav files
addpath(genpath('G:\Ch.1_ShipProject\CINMS_B_Shipping\CombinedPassageData\VSR_2014_2017'));

ftype = '*.txt'; % Find all .txt files
listing = dir([fpr, ftype]); 
fn = {listing.name}; % getting the file names of all the .txt files
wavFilename = strrep(fn,'.txt','.wav');
%getting wav file names from corresponding txt files
n = length(fn);



%This loop reads in the .wav files, CPA Time, CPA Distance, and Mean Speed
for k = 1:n
    
    fprintf(['\n', wavFilename{k}]) %printing what file it's on
    
    filenameChar = wavFilename{k};
    
    finfo(k) = audioinfo([fpr, wavFilename{n}]); %getting the info of each
    
    
    soundFile = fullfile(listing(k).folder, wavFilename(k));
    % check if txt file exists
    txtFile = strrep(soundFile,'.wav','.txt');
    
    if ~exist(char(txtFile),'file')
        warning('Could not find file %s, skipping.','txtFile')
        continue
    end
    
    textData = importdata(char(txtFile));
    
    
    if ~isfield(textData,'textdata')
        textTemp = textData;
        textData = struct;
        textData.textdata = textTemp;
    end
    
    
    if isfield(textData,'data') && size(textData.data,2)>3
        meanSOG(k) = mean(textData.data(:,4));
    end
    
    
    k6Idx = find(~cellfun(@isempty,strfind((textData.textdata(:,1)),'CPATime'))==1);
    if ~isempty(k6Idx)
        [k6,~]= regexp(textData.textdata{k6Idx,1},'(\d*/\d*/\d*\W\d*:\d*:\d*)','tokens','match');
        CPATime(k) = datenum(k6{1,1});
    end
end

A = readtable('TimeOffset_Transits11_23_2020.csv');%getting all the ship lengths
shipLength = A(:,17); %just ship lengths
textFilecolumn = table2array(A(:, 2));
textFilecolumn = char(textFilecolumn); %.txt files to which the ship lengths correspond to
%textFilecolumnStr = textFilecolumn(:, 48:74)
textFilecolumnCell = cellstr(textFilecolumn)
shipLength_File = horzcat(shipLength, textFilecolumnCell);
shipLength_File.Properties.VariableNames{'Var2'} = 'txtFile'


fnTable = array2table(fn')
CPATimeTable = array2table(CPATime')
fnTable.Properties.VariableNames{'Var1'} = 'txtFile'
T = innerjoin(fnTable, shipLength_File); %getting the ship lengths for the files i have here
shipLength_Deployment = T(:, 2);
shipLength_Deployment = table2array(shipLength_Deployment);
shipLength_Deployment = shipLength_Deployment';
shipLength = shipLength_Deployment


%Getting the exact length in time I want for each specific transit
DWL = 2*shipLength*tand(30);         %data window length, entire length (in meters)
DWP = DWL./meanSOG                  % Desired file length in seconds, getting length in time (seconds)
DWchunks = DWP/2                    % divided by two because we want chunk on either side


%%%Changing the DWchunks from seconds into days.
spd = 60*60*24;
DWchunkdays = DWchunks/spd





for fnum = 1:length(fn)
    fprintf(['\n', wavFilename{fnum}])
    
    fs = finfo.SampleRate;
    dfs = DWP*fs;                               % desired file size
    
    nfnbase = wavFilename{fnum};                       % basis for new file names
    
    startCPATimechunk = CPATime(fnum) - DWchunkdays(fnum);
    endCPATimechunk = CPATime(fnum) + DWchunkdays(fnum);
    
    ftime = datenum([2000 + str2num(wavFilename{fnum}(11:12)), str2num(wavFilename{fnum}(13:14)), str2num(wavFilename{fnum}(15:16)), str2num(wavFilename{fnum}(18:19)), str2num(wavFilename{fnum}(20:21)), str2num(wavFilename{fnum}(22:23))]);
    
    nstart = floor((startCPATimechunk - ftime)*spd*fs);
    nend = ceil((endCPATimechunk - ftime)*spd*fs);
    
    info = audioinfo([fpr wavFilename{fnum}]);
    
    if  nstart > 0 && nend < info.TotalSamples

            [y, ~] = audioread([fpr wavFilename{fnum}], [nstart, nend]);
            audiowrite([fpw, wavFilename{fnum}], y, fs);

    end

 
    clear nstart
    clear nend
    clear y


end
