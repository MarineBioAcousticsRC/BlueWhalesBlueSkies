%%reading in datasets
COP = readtable('COP_combinedPassageData.csv');
SBARC = readtable('SBARC_CombinedPassageData.csv');

%%combining them

both = [COP;SBARC];

%%%%%picking out the unique transit date times
[~,idx] = unique(both(:,8));
uniquerows = both(idx,:);

%getting rid of any rows with NaNs
uniquerows2 = uniquerows(~any(ismissing(uniquerows),2),:);