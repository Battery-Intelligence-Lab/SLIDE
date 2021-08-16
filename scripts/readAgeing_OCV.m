% Script to read and plot the half-cell OCV curves from the check-ups in 
% the degradation simulations. 
%
% This script should not be executed on its own, but is called by one of
% three higher-level scripts: 
%   readCalendarAgeing
%   readCycleAgeing
%   readProfileAgeing
%
%
% Copyright (c) 2019, The Chancellor, Masters and Scholars of the University 
% of Oxford, VITO nv, and the 'Slide' Developers.
% See the licence file LICENCE.txt for more information.

fileName = 'DegradationData_OCV.csv';               % Name of the file which contains the OCV curves

%% Read the OCV curves
% A loop to read the files of each ageing regime
for i=1:length(IDs)
    fol = strcat(pref,'_',ageingID,'_',IDs{i});     % Folder in which to find the file
    na2 = fullfile(pathvar.results_folder, fol, fileName);                   % Full name of the file
    try
        B = csvread(na2);
    catch    
        B = nan(50*3, 30);                          % If we can't read the file, print a warning to the user, and fill the data with NaNs
        warning(['warning no file ' na2 ' could be found'])
    end

    % Store the OCV curves in a struct
    OCV{i}.curves = B;
end

%% Plot the OCV curves
% Plot the half-cell OCV curves for the different ageing regimes. There is
% one subplot per ageing regime.
% The curves indicate the half cell OCV curves. The 0 of the x-axis is
% where the difference between both OCV curves is 4.2V (i.e. when the cell
% is fully charged).
% There are stars at the points on the OCV curves where the cell was
% operating before the check-up. They allow to track the slipping of the 
% usable voltage window (when cycling) or the exact point on the OCV curves 
% where the cell is resting (when doing calendar ageing). Note that the
% y-point of the start is exact (i.e. the cell is exactly operating at this
% voltage), but the x-value is rounded to the nearest value of the OCV
% curve.

% Find how to organise the subplots
nregimes = length(IDs);
if nregimes <= 15
    nrows = 3;                      % number of rows
else 
    nrows = 4;                      % number of columns
end
ncols = ceil(nregimes/nrows);

figure()
for i=1:length(IDs)
    try
        A = OCV{i}.curves;
            % Every OCV measurement consist of 3 rows:
            %   - the cumulative discharged Ah
            %   - the OCV of the cathode
            %   - the OCV of the anode
            % First first x columns give the OCV curves
            % then there are two empty columns
            % then, on the row of the OCV anode, there are two values. They
            %   represent the points on the OCV curves where the cell was
            %   opoerating before the checkup was done. 
        A(A == 0) = nan;                % replace the zeros with nan so they are not plotted

        % Define the colours
        [a,b] = size(A);                % size of the OCV file
        col = winter(a/3);              % Every OCV measurement consist of 3 rows, so the number of measurements is a/3
        col2 = autumn(a/3);

        % Get the half cell curves
        Ah = A(1:3:end,:);              % discharged Ah, zero when the cell is at 4.2V
        OCVp = A(2:3:end,:);            % cathode OCV
        OCVn = A(3:3:end,:);            % anode OCV

        subplot(nrows,ncols,i)
            for j=1:a/3                 % loop for every OCV measurement
                % Remove the points which were out of range of the OCV curve.
                % (If we reached the end of one OCV curve, the value written in
                % the CSV file was just the last value. Now remove these
                % values)
                ocvp = OCVp(j,:);
                ip = find(ocvp == ocvp(1));     % find where the OCVp values were the same as the first value
                ocvp(ip(1:end-1)) = nan;        % remove these points, except the last one
                ip = find(ocvp == ocvp(end));   % find where the OCVp values were the same as the last value
                ocvp(ip(2:end)) = nan;          % remove these points, except the first one
                ocvn = OCVn(j,:);
                in = find(ocvn == ocvn(1));
                ocvn(in(1:end-1)) = nan;
                in = find(ocvn == ocvn(end));
                ocvn(ip(2:end)) = nan;

                % Plot the OCV curves
                plot(Ah(j,:),ocvp,'color',col(j,:));
                hold on
                plot(Ah(j,:),ocvn,'color',col(j,:));

                % Get the operating point where the cell was before we started the
                % checkUp. 
                OCVpi = nan;
                OCVni = nan;
                for m=1:length(OCVn(j,:))
                    if isnan(OCVn(j,m))         % this was the last column of the OCV measurement
                        OCVpi = OCVn(j,m+2);    % the cathode operating point is 2 columns further
                        OCVni = OCVn(j,m+3);    % the anode operating point is 3 columns further
                        break
                    end
                end
                % Find where on the x-axis these points are
                inp = find(OCVp(j,:)<OCVpi);    % indices of the points to the right on the cathode OCV curves
                xr = Ah(j,inp(1));              % Ah(inp(1)) is the x-value of the data point to the right on the cathode OCV
                xl = Ah(j,inp(1)-1);            % the previous data point on the OCV curve is the data point to the left on the cathode OCV
                yr = ocvp(inp(1));              % the cathode OCV at the data point to the right
                yl = ocvp(inp(1)-1);            % the cathode OCV at the data point to the left
                xp = interp1([yl yr],[xl xr],OCVpi); % interpolate linearly between both to find the x-point of the cathode operating point
                inn = find(OCVn(j,:)>OCVni);    % indices of the points to the right on the anode OCV curves
                xr = Ah(j,inn(1));              % Ah(inn(1)) is the x-value of the data point to the right on the anode OCV
                xl = Ah(j,inn(1)-1);            % the previous data point on the OCV curve is the data point to the left on the anode OCV
                yr = ocvn(inn(1));              % the anode OCV at the data point to the right
                yl = ocvn(inn(1)-1);            % the anode OCV at the data point to the left
                xn = interp1([yl yr],[xl xr],OCVni); % interpolate linearly between both to find the x-point of the anode operating point
                % Plot these operating points
                try
                    plot(xp,OCVpi,'*','color',col2(j,:))
                    plot(xn,OCVni,'*','color',col2(j,:))
                catch % if we can't find them, just skip this step
                end
            end
    catch
        % No OCV data was available for this ageing regime
    end
    % Labels and title
    ylim([0 5])
    xlabel('[Ah]')
    ylabel('[V]')
    title(IDs{i},'Interpreter', 'none')
end