% Function to plot one histogram for the script readTestStatsData.m
% it plots on whatever axis is active when the function is called

function plotHistogram(edges,bins)
% IN
% edges     the boundaries between the bins
% bins      counts in each bin
%           how often was edge(i-1) <= data < edge(i)
%           while the first bin has no left edge
%           and the right edge has no right edge
% The plot is drawn with the width of the outer bins equal to 10 times the
% width of a regular bin

dx = abs(edges(2) - edges(1));
xmin = edges(1) - 10*dx;
xmax = edges(end)+10*dx;
x = [xmin ; edges ; xmax];


% Plot as stairs
stairs(x(1:end-1),bins)

% Plot with filled area under stairs
N = length(x);
xf = reshape(repmat(x',2,1),2*N,1);         % double x values -> [0 0 1 1 ... N N]
xf = [xf(2:end) ; xf(1)];                   % go round in circle -> [0 1 1 2 2 ... N N 0]
N = length(bins);
yf = reshape(repmat(bins',2,1),2*N,1);
yf = [yf ; 0 ; 0];                          % -> [y(1) y(1) y(2) y(2) y(3) .... y(N) 0 0]
xf = [xf ; xf(1)];                          % repeat first point to fill the area between 0 and y(1) at x(1)
yf = [yf ; yf(1)];

area(xf,yf)


% plot as bars
% xmid = x(1:end-1) + 1/2*diff(x);
% bar(xmid,bins)


xlim([xmin xmax])


