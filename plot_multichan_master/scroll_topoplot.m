function scroll_topoplot(epoch_data, bad_chan_table, chanloc, freq_thresh)
    %INPUT example:
    % chanloc = eeg_struct.chanlocs;
    
    lowfreq_thresh = freq_thresh(1);
    highfreq_thresh = freq_thresh(2);
    
    figure; hold on;
    
    subplot(3,3,3); hold on   
    marker_size = 20;
    for i_chan = 1:size(bad_chan_table,1)
        for i_freq = 1:size(bad_chan_table,2)
            if bad_chan_table(i_chan, i_freq) > 0
                if i_freq <= lowfreq_thresh; scatter(i_freq, i_chan, marker_size, 'k.')
                elseif i_freq >= highfreq_thresh; scatter(i_freq, i_chan, marker_size, 'r.')
                else; scatter(i_freq, i_chan, marker_size/5, 'b.')
                end
            end
        end
    end
    xlim([ 1 size(bad_chan_table,2) ]) 
    ylim([ 1 size(bad_chan_table,1) ]) 
    ylabel('channels')
        
    % LOW frequency (< 10Hz)
    bad_lowfreq_xchan = sum(bad_chan_table(:,1:lowfreq_thresh),2);
    subplot(3,3,6); hold on    
    %single vector of channel values
    topoplot(bad_lowfreq_xchan, chanloc, 'electrodes','on');
    colormap(flipud(gray));
    caxis([0 max(bad_lowfreq_xchan)]); %colorbar
    freezeColors;
    title(['chan outlier: LOW freq <' num2str(lowfreq_thresh) 'Hz'])
    
    % HIGH frequency (> 40Hz)
    bad_highfreq_xchan = sum(bad_chan_table(:,highfreq_thresh:end),2);
    subplot(3,3,9); hold on    
    topoplot(bad_highfreq_xchan, chanloc, 'electrodes','on');
    colormap(flipud(hot));
    caxis([0 max(bad_lowfreq_xchan)*2.5]); %colorbar
    freezeColors;
   title(['chan outlier: HIGH freq <' num2str(highfreq_thresh) 'Hz'])
    
    % MULTICHANNEL SCROLL - - - - - - - - - -
    subplot(3,3,[1,2,4,5,7,8])
    plot_multichan_nonormalize(epoch_data, bad_chan_table); %colorbar
    %plot_multichan_nonormalize(epoch_data, bad_chan_table, freq_thresh); %colorbar
    xlabel('sample points')
    ylabel('channels')
end