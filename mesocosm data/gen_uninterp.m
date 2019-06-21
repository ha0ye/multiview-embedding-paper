% Generates time series with a specified time gap, dt. Start with 1st
% point, construct as long a chunk of timeseries with dt in [gap].
% record start and end points in lib (library). start with next
% unused point, repeat until all points used. keep only chunks at least
% c days in length.

cd('/Users/hye/Dropbox/projects/multiembed/analysis/mesocosm data/');

gap = 6:8;
c = 12; % minimum length of each chunk in days
basename = 'mesocosm_data';

M = importdata('ELE.full.raw.txt');

% remove duplicated data point on day 2301
t = M.data(:,1);
duplicated_indices = find(t(2:end)-t(1:(end-1)) == 0);
M.data(duplicated_indices, :) = [];

% add nutrient data
q = importdata('ELE.nutrients.raw.txt');
t2 = q.data(:,1);

[lia, loc] = ismember(t2, t);
added_cols = size(M.data,2) + (2:size(q.data, 2)) - 1;
M.data(:, added_cols) = NaN;
M.colheaders = [M.colheaders q.colheaders(2:end)];
M.data(loc(lia), added_cols) = q.data(lia,2:end);

% find closest time index for non-exact date matches
unused = q.data(~lia, :);
for idex = 1:size(unused,1)
    [xi,Ii] = min(abs(unused(idex,1) - t));

    if abs(xi) <= 3 && isnan(M.data(Ii(1), end))
        M.data(Ii(1), added_cols) = unused(idex, 2:end);
    end    
end

% remove times without nutrient data
%valid_index = isfinite(M.data(:, end));
%M.data = M.data(valid_index,:);

% start finding data segments
tpointsleft = M.data(:,1);
npointsused = 0;
tseries = [];
lib = [];

while ~isempty(tpointsleft)

    % start a new chunk
    
    chunk = M.data(M.data(:,1) == tpointsleft(1),:);
    tpointsleft = tpointsleft(2:end); % remove start point of new chunk from tpointsleft
    
    while ~isempty(tpointsleft) && any(ismember(tpointsleft(:) - chunk(end,1),gap))
        
        i_found = find(ismember(tpointsleft(:) - chunk(end,1),gap));
        chunk = [chunk ; M.data(M.data(:,1) == tpointsleft(i_found(1)),:)];
        tpointsleft(i_found) = []; % remove point from tpointsleft if used
        
    end
    
    if size(chunk,1) >= c
    
        lib = [lib ; size(tseries,1)+1 size(tseries,1)+size(chunk,1)];
        tseries = [tseries; chunk];
    
    end
    
end

% save data and library files as tab delimited text
labels = sprintf('%s\t', M.colheaders{:});
dlmwrite([basename '.txt'], labels, 'delimiter', '');
dlmwrite([basename '.txt'], tseries, 'delimiter', '\t', '-append');
dlmwrite([basename '.lib'], lib, 'delimiter', '\t');

% clear tseries*


%%
%#### APPEND NUTRIENTS

labels = [M.colheaders,'N','P'];
M = importdata('ELE.nutrients.raw.txt');

tseries = [tseries nan(size(tseries,1),2)];

for idex = 1:size(tseries,1)
    [xi,Ii] = min(abs(tseries(idex,1)-M.data(:,1)));
    
    if abs(xi) <= 3
        tseries(idex,end-1:end) = M.data(Ii(1),2:3);
    end
    
end

dlmwrite([basename '+nutrients.txt'],tseries,'delimiter','\t')

%% spline_interpolated version

clear('all');

m = importdata('ELE.full.raw.txt');
spline_tp = 3.35;

x = m.data(:,1);
duplicated_indices = find(x(2:end)-x(1:(end-1)) == 0);
m.data(duplicated_indices, :) = [];

x = m.data(:,1);
xx = [1:spline_tp:max(x)];
block = nan(length(xx), size(m.data,2));

for(j = 1:size(m.data,2))
    y = m.data(:,j);
    valid_index = isfinite(y);
    block(:,j) = interp1(x(valid_index), y(valid_index), xx, 'pchip', NaN);
end

labels = sprintf('%s\t',m.colheaders{:});
dlmwrite('ELE.spline.txt', labels, 'delimiter', '');
dlmwrite('ELE.spline.txt', block, 'delimiter', '\t', '-append')










