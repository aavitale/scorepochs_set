function plot_multichan_nonormalize( x, bad_chan_table, freq_thresh )
%function plot_multichan( x, y, interval, normalize, bad_chan_idx )
%function plot_multichan( x, bad_chan_idx )

% Simple MATLAB function for plot multi-channel time-series data
% 
% Usage:
%    plot_multichan( y )    % <- y: signal
%    plot_multichan( x, y ) % <- x: time
% 
% Example: 
%    y = randn([20, 2000]); 
%    plot_multichan(y);
%
% Written by Hio-Been han, hiobeen.han@kaist.ac.kr, 2020-03-07
% modified by andrea.vitale@gmail.com  20210325

% parse arguments
if nargin==1, y=x; x = 1:length(y); end
y=x; x = 1:length(y); 
nChan = size(y,1);
if nChan > size(y,2),  y = y'; nChan = size(y,1); end
%if nargin < 4, normalize = 1; end
%if nargin < 4, normalize = 0; end
normalize = 0
if normalize
    stds = nanstd( y, 0, 2 );
    for chIdx = 1:size(y,1), y(chIdx,:) = nanmean(stds) * (y(chIdx,:) / stds(chIdx)); end
end
% if nargin < 3
    interval = nanmean(range(y, 2)) * nChan / 2.5;
% end
y_center = linspace( -interval, interval, nChan );
y_center = linspace( -interval*1.15, interval*1.15, nChan );

if nargin < 3
    lowfreq_lim = 10; highfreq_lim = 65
else
    lowfreq_lim = freq_thresh(1); 
    highfreq_lim = freq_thresh(2);
end

[bad_chan_idx, ~] = find(sum(bad_chan_table,2) >2 );
[bad_chan_lowfreq_idx, ~] = find(sum(bad_chan_table(:,1:lowfreq_lim),2) > 2);
[bad_chan_highfreq_idx, ~] = find(sum(bad_chan_table(:,highfreq_lim:end),2) > 2);

% set colormap
color_template =...
   [843 088 153;
    992 750 280;
    400 200 030;
    573 716 350;
    055 538 083]*.001;
c_space = repmat( color_template, [ ceil( nChan/size(color_template,1)), 1]);
% c_space=imresize(colormap('lines'),[nChan,3],'nearest');

% main plot
chanlab = {}; chanlab_pos = [];
lw = 0.5 %1;
chan_scale = 50 %micro_volt
for chanIdx = 1:nChan
    shift = y_center(chanIdx) + nanmean( y( chanIdx, : ), 2);
    
    if ismember(chanIdx, bad_chan_highfreq_idx) 
        plot( x, y( chanIdx, : ) - shift, 'r', 'LineWidth', lw*2); 
    %elseif ismember(chanIdx, bad_chan_idx) 
    elseif ismember(chanIdx, bad_chan_lowfreq_idx) 
        plot( x, y( chanIdx, : ) - shift, 'k', 'LineWidth', lw*2); 
    else
        plot( x, y( chanIdx, : ) - shift, 'b', 'LineWidth', lw/2);
        %plot( x, y( chanIdx, : ) - shift, 'Color', c_space( chanIdx,: ) , 'LineWidth', lw);
    end
    %ylim([ 0 chan_scale ])
    
%     if ismember(chanIdx, bad_chan_idx) 
%         chanIdx_reverse = nChan-chanIdx+1;
%         chanlab{chanIdx} = sprintf( 'Ch %02d', chanIdx_reverse);
%         chanlab_pos(chanIdx) =  y_center(chanIdx) ;
%     else
        chanIdx_reverse = nChan-chanIdx+1;
        if ismember(chanIdx_reverse, bad_chan_idx)
            chanlab{chanIdx} = sprintf('Ch %02d', chanIdx_reverse);
        else
            chanlab{chanIdx} = sprintf(' ');
        end
        chanlab_pos(chanIdx) =  y_center(chanIdx) ;
%    end
    
    
    if chanIdx ==1, hold on; end
end
hold off;

% enhance visibility
set(gca, 'YTick', chanlab_pos, 'YTickLabel', chanlab, 'FontSize', 5, ...
    'Clipping', 'on', 'Box', 'off', 'LineWidth', 2);
ylim([-1 1]*interval*1.2);
