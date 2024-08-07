% Similar to plotHistogram but will plot bars for multiple cells
% The width of the extreme bins is however the same as the width of the
% bins in the middle, since Matlab's bar does not allow to change the widht
% of individual bars

function plotHistogramMultiple(edges,bins, lin)
% IN
% edges     the boundaries between the bins
%           must be ONE ARRAY (boundaries must be the same for all cells)
% bins      counts in each bin
%               how often was edge(i-1) <= data < edge(i)
%               while the first bin has no left edge
%               and the right edge has no right edge
%           must be a MATRIX, one column = one cell; one row = one bin
% lin       boolean, if true, the plot is made with one line per cell
%           if false, a bar plot is used with one bar per cell
% The plot is drawn with the width of the outer bins equal to the
% width of a regular bin

% x-axis
dx = abs(edges(2) - edges(1));
xmin = edges(1) - 1*dx;
xmax = edges(end)+1*dx;
x = [xmin ; edges ; xmax];

% yaxis, make frequency
[~,n] = size(bins);
for i=1:n
    s = sum(bins(:,i));
    bins(:,i) = bins(:,i)/s;
end

if lin
    % Plot as stairs
    stairs(x(1:end-1),bins)
else
    % plot as bars
    xmid = x(1:end-1) + 1/2*diff(x);
    bar(xmid,bins)
end

xlim([xmin xmax])


