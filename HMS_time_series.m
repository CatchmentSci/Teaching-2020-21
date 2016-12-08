function [timestamp, normalisedQ, ssc] = HMS_time_series(riverId)

% Usage:  [timestamp, normalisedQ, ssc] = HMS_time_series({'1001','9004'});

% Required input:
% riverId: The riverId of interest. This corresponds to the monitoring
% codes for each HMS within the database. This should be input as a text
% e.g. '6009' (including the quotation marks)

% Outputs:
% timestamp = time of the sample taking place
% normalisedQ = the normalised discharge estimate for each measurement
% ssc = suspended sediment concentrations

% Load the datasets
outfilename = websave('FlowSed.csv','https://github.com/CatchmentSci/Teaching-2016-17/raw/master/Flow_Sed_Trimmed.csv');
websave('readtext.m', 'https://raw.githubusercontent.com/CatchmentSci/Glaisdale-Beck-diversion-scheme/master/readtext.m'); % Download dependancy
websave('replace.m', 'https://raw.githubusercontent.com/CatchmentSci/Glaisdale-Beck-diversion-scheme/master/replace.m'); % Download dependancy
[data_text,~] = readtext(outfilename, ',', '','','textual'); % read in the comma delimeted data
sample_id = str2double(data_text(1:end,1));
det = str2double(data_text(1:end,2));
result = str2double(data_text(1:end,3));

outfilename2 = websave('tblSample.csv','https://github.com/CatchmentSci/Teaching-2016-17/raw/master/tblSample.csv');
[data_text2,~] = readtext(outfilename2, ',', '','','textual'); % read in the comma delimeted data
data_text2(:,3) = cellfun(@(x) x(1:10), data_text2(1:end,3),'UniformOutput', false); % Extract the dates from the third column
data_text2(:,4) = cellfun(@(x) x(12:end), data_text2(1:end,4),'UniformOutput', false); % Extract the times from the forth column
data_text2(:,5) = strcat(data_text2(:,3), {' '}, data_text2(:,4)); % merge the dates and times

for iter = 1:length(riverId) % loop through each of the stations
    
    cell_id = double(find(strcmp(data_text2(:,2),riverId(iter)) == 1)); % Find the measurement associated with the Coquet
    meas_id = str2double(data_text2(cell_id,1)); % Extract the measurement identifier
    [~,idx]=ismember(meas_id,sample_id); % find the locations of each coquet sample id within the main dataet
    timestamp = data_text2(cell_id,5);
    
    Out = cell(1, numel(meas_id));
    X = [];
    for iA = 1:numel(meas_id)
        Out{iA} = find(sample_id == meas_id(iA));
        splitOut = Out{iA};
        for a = 1:length(splitOut)
            tempDet = det(splitOut(a));
            if isempty(X)
                X(1,1:3) = NaN;
            elseif length(X(:,1)) < iA
                X(iA,1:3) = NaN;
            end
            
            if tempDet == 107
                X(length(X(:,1)),1) = result(splitOut(a)); % average flow in first column
            elseif tempDet == 108
                X(length(X(:,1)),2) = result(splitOut(a)); % instantaneous flow in second column
            elseif tempDet == 114
                X(length(X(:,1)),3) = result(splitOut(a)); % suspended sediment concentration in third column
            end
        end
    end
    
    X = replace(X,0,NaN);
    indUse = find(~isnan(X(:,2))); % Find existing instantaneous Q measurements
    tempQ(1:length(X),1) = NaN;
    tempQ(indUse) = X(indUse,2); % add these measurements to the temporary array
    for a = 1:length(X)
        if isnan(tempQ(a)) && ~isnan(X(a,1)) % If there are no instantaneous measurements, only averaged, use the average ones
            tempQ(a) = X(a,1);
        end
    end
    
    indUse2 = find(~isnan(tempQ(:,1))); % Find Q measurements
    geoMeanQ = geomean(tempQ(indUse2,1));
    normalisedQ = tempQ(:,1)./geoMeanQ(:,1);
    ssc = X(:,3);
    Tnum = datenum(timestamp,'dd/mm/yyyy HH:MM:SS');
    [~,ind] = sort(Tnum);
    
    h = figure();
    set(h,'units','normalized','outerposition',[0 0 1 1]) % set full screen mode
    set(h,'DefaultTextFontname', 'cmr12');
    set(h,'DefaultAxesFontName', 'cmr12');
    set(h,'DefaultAxesFontSize',14);
    
    [ax, p1, p2] = plotyy(Tnum(ind),tempQ(ind), Tnum(ind), ssc(ind), @line, @line); hold on
    set(p1, 'LineWidth', 1, ...
        'Marker', '+')
    set(p2, 'LineWidth', 1, ...
        'LineStyle', '-', ...
        'Marker', '+')
    box off
    grid on
    
    axis(ax(:), 'tight')
    set (ax, 'TickDir', 'out')
    set (ax(1), 'YLim', [-100 max(get(ax(1),'YLim'))])
    set (ax(2), 'YLim', [0 max(get(ax(2),'YLim'))])
    legend ('off')
    yLabel1 = sprintf('River Discharge (m^{3} s^{-1})');
    ylabel(ax(1),yLabel1);
    yLabel2 = sprintf('Suspended sediment concentration (mg L^{-1})');
    ylabel(ax(2),yLabel2);
    L = get(gca,'XLim');
    NumTicks = 9;
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    datetick('x','yyyy')
    
    titleStr = [strcat('Station Number=' , riverId(iter))];
    hTitle = title (titleStr);
    set(hTitle,'FontWeight', 'normal');
    
    
end























