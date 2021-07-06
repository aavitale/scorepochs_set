% Scorepochs 
%
% Function to select the best (most homogenoous) M/EEG epochs from a
% resting-state recordings. 
%
%     Copyright (C) 2020 Matteo Demuru, Matteo Fraschini
%     last modified by andrea.vitale@gmail.com  20210706
%
% INPUT
%    cfg struct with the following fields
%           freqRange    - array with the frequency range used to compute the power
%                          spectrum (see MATLAB pwelch function)
%           fs           - integer representing sample frequency         
%           windowL      - integer representing the window length (in seconds)  
%           smoothFactor - smoothing factor for the power spectrum
% 
%    set_file            - file .set (as in the eeglab structure with already the channels location) 
%    (eeglab_dir          - directory with eeglab package) !!! not required
%     
%
%
% OUTPUT
%      
%    epoch       -  cell array of the data divided in equal length epochs 
%                   of length windowL (channels X time samples)
%                  
%    idx_best_ep - array of indexes sorted according to the best score
%                  this array should be used for the selection of the best
%                  epochs
%
%    score_Xep   - array of score per epoch
%
%    ADDITIONAL Output:
%           score_chXep:  scorepoch values not averaged across channels
%
%           psd_chXep: psd computed for each epoch and channel
%
%           2 FIGURES
%
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.


%function [idx_best_ep,epoch,score_Xep] = scorEpochs(cfg,data) 
function [idx_best_ep, epoch, score_Xep, score_chXep, psd_chXep, pxx] = scorEpochs_set(cfg, set_file) %, eeglab_dir) 

    % !!! EEGLAB and scorepoch dir (with accessories)
    % is assumed to be already in the path
    
    %eeglab_dir = 'D:\_TOOLBOX\eeglab_20201226'; cd(eeglab_dir); eeglab('nogui');

    % set_file = fullfile('C:\Users\Utente\Desktop\CMI_EEG_PREProcess\HAPPE_study\intermediate2_ICAclean','EEG_concat_ICAcleanedwithMARA.set')
    % set_file = fullfile('C:\Users\Utente\Desktop\CMI_EEG_PREProcess\data_prep','NDARUC980NZ5_RestingState_cleanraw_avgref_nobadICA.set')
    [ filepath , set_file_name , ext ] = fileparts( set_file )
    set_file_name = regexprep(set_file_name, '_', ' ');
    
    try
        eeg_struct = pop_loadset(set_file)
    catch
        error_message
        disp(error_message)
    end
    
    
    % check channel location - - - - - - --  -
    % figure; topoplot([], eeg_struct.chanlocs, 'electrodes','labelpoint','chaninfo', eeg_struct.chaninfo);
    
    
    %% CONFIGURATION file  
    if ~exist('cfg') || isempty('cfg')
        cfg = []; 
        % <<<<<<<<<<<<<<<<<< ENTIRE FREQUENCY RANGE<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        cfg.freqRange = [ 2 : 80 ];
        % <<<<<<<<<<<<<<<<<< ONLY ALPHA BAND <<<<<<<<<<<<<<<<<<<<<<<<<<<<
        %cfg.freqRange = [ 8 : 13 ]; 
        
        cfg.fs = eeg_struct.srate;
        cfg.windowL = 2; % in sec <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        cfg.smoothFactor = 0;
    end
    
    
    %% divide the data in epochs of windowL length (in sec)
    data = eeg_struct.data;
    
    epLen   = cfg.windowL * cfg.fs;
    dataLen = size(data,2);
    nCh     = size(data,1);
    idx_ep  = 1 : epLen : dataLen - epLen + 1;
    nEp     = numel(idx_ep);
    
    epoch   = cell(1,nEp);
     
    pxx     = cell(1,nEp);
    
    
        
    for e = 1 : nEp 
        epoch{e} = data(:,idx_ep(e):idx_ep(e)+epLen-1);
        % compute POWER SPECTRUM
        % pxx dimensions = n_channel x n_freq 
        % (power for each frequency in a specific epoch)
        pxx{e} = pwelch(epoch{e}',[],[],cfg.freqRange,cfg.fs)';
       
        if(cfg.smoothFactor ~= 0)
                pxx{e} = movmean(pxx{e}',cfg.smoothFactor)';
        end
    end
    
    % compute SCORE across channels and across epochs 
    pxxXch      = zeros(nEp,numel(cfg.freqRange));
    score_chXep = zeros(nCh,nEp);
    psd_chXep = [];
    for c = 1 : nCh
        
        for e = 1 : nEp
            
            pxxXch(e,:) = pxx{e}(c,:);  
            
            psd_chXep(c,e,:) =  pxx{e}(c,:);  
            
        end
        score_ch = corr(pxxXch','type','Spearman');
        
        % score_chXep : % correlation across epochs for a single channel
        score_chXep(c,:) = mean(score_ch);  
    end    
    
    % score_Xep : % mean of the score_epoch across all channels
    score_Xep = mean(score_chXep,1);
    
    [~,idx_best_ep] = sort(score_Xep,'descend');
    
    
    %% <<<<<<< PLOT <<<<<<<<<<<<<<<<<<<<<<
    do_plot = 1;
    if do_plot
        
        % FIGURE 1: BAR GRAPH of the score_epoch values (not sorted)
        do_mean = 1; thresh_level = 2; %std deviation
        %do_median = 0; thresh_level = 2; %mad
        line_width = 1.5;

        figure; subplot(121); hold on; %subplot(2,5,1:4); hold on
        pl1 = bar(score_Xep);
        set(pl1,'FaceColor',[ 1 1 1 ]); %[0.5 0.5 0.5]);
        xlim([1 length(score_Xep)])
        xlabel(['epochs (with window length of ' num2str(cfg.windowL) ' sec)'])
        ylabel('score x ALLChannels')
        ylim([0 1.05]);
        
        if do_mean
            yline(mean(score_Xep),'k-', 'LineWidth',line_width)
            yline(mean(score_Xep) - std(score_Xep),'k--', 'LineWidth',line_width)
            yline(mean(score_Xep) - 2*std(score_Xep),'k:', 'LineWidth',line_width)
            title([ set_file_name ' (mean = ' num2str(mean(score_Xep)) ')'])
        
%         elseif do_median
%             score_median = median(score_Xep);
%             yline(score_median,'k-', 'LineWidth',line_width)
%             yline(score_median - thresh_level*mad(score_Xep),'k--', 'LineWidth',line_width)
%             yline(score_median - thresh_level*1.5*mad(score_Xep),'k:', 'LineWidth',line_width)
%             title(['resting state scorepochs (median = ' num2str(score_median)])
            
        end
        
        % = = = = = = =  = = = = = = = = =  = = == =  ==  == = = = = 
        % identify CHANNELS with low score_epoch as compared to other channels:
        subplot 122; hold on; 
        
        if do_mean
            score_xchan_mean = mean(score_chXep,2);
            score_xchan_std = std(score_chXep,1,2);

            score_allchan_mean = mean(score_xchan_mean,1);
            score_allchan_std = std(score_xchan_mean,1,1);
                  
            errorbar(score_xchan_mean, score_xchan_std,'ko')
            
            hold on;
            yline(score_allchan_mean, 'k','LineWidth',line_width)
            yline(score_allchan_mean - score_allchan_std,'k--','LineWidth',line_width); 
            yline(score_allchan_mean - 2*score_allchan_std, 'k:','LineWidth',line_width); 
            xlabel('channels')
            ylabel(['score (mean - ' num2str(thresh_level) ' std)'])
            ylim([0 1.05]);
           
            % find channel outliers:
            bad_chan_idx = find(score_xchan_mean <= (score_allchan_mean - thresh_level*score_allchan_std))';
            for ii = 1:length(bad_chan_idx)
                hold on
                errorbar(bad_chan_idx(ii), score_xchan_mean(bad_chan_idx(ii)), ...
                    score_xchan_std(bad_chan_idx(ii)),'rx', 'LineWidth', 2)
            end
            title(['channels with score < ' num2str(thresh_level) 'STD = '...
                    num2str(reshape(bad_chan_idx,[1,length(bad_chan_idx)]))])
        
%         elseif do_median
%             % MEDIAN - - - - - - - - - - - 
%             score_xchan_median = median(score_chXep,2);
%             score_xchan_mad = mad(score_chXep,1,2);
% 
%             score_allchan_median = median(score_xchan_median,1);
%             score_allchan_mad = mad(score_xchan_median,1,1);
% 
%             hold on; %subplot 122
%             errorbar(score_xchan_median, score_xchan_mad,'ko')
%             hold on;
%             yline(score_allchan_median,'k-','LineWidth',line_width)
%             %yline(score_allchan_median + thresh_level*score_allchan_mad, '--'); 
%             yline(score_allchan_median - thresh_level*score_allchan_mad, 'k--','LineWidth',line_width); 
% 
%             %yline(score_allchan_median + thresh_level*1.5*score_allchan_mad, ':'); 
%             yline(score_allchan_median - thresh_level*1.5*score_allchan_mad, 'k:','LineWidth',line_width); 
%             xlabel('channels')
%             ylabel(['score (median - ' num2str(thresh_level) ' mad)'])
%             %ylim([y_lim(1) 1]);
%             ylim([0 1.1]);
%            
%             % find channel outliers:
%             bad_chan_idx = find(score_xchan_median <= (score_allchan_median - thresh_level*1.5*score_allchan_mad))';
%             for ii = 1:length(bad_chan_idx)
%                 hold on
%                 errorbar(bad_chan_idx(ii), score_xchan_median(bad_chan_idx(ii)), ...
%                     score_xchan_mad(bad_chan_idx(ii)),'rx', 'LineWidth', 2)
%             end
%             title(['channel with score < ' num2str(thresh_level) 'MAD = ' num2str(bad_chan_idx)])
        
        end
        % -- - - - - - - - -  
            
            
        %% MULTICHANNEL plot + TOPOPLOT
        [ tmp, i_epoch ] = min(score_Xep);  %<<<<<<< epoch with the LOWEST score
        i_epoch = i_epoch(1);
        %     % or find epoch with score value closest to the thresh-level:
        %     score_median = median(score_Xep);
        %     score_thresh = score_median - thresh_level*1.5*mad(score_Xep);
        %     [ tmp, i_epoch ] = min(abs(score_Xep - score_thresh));

        % change color of the histogram - - - - - - - - - - - - 
        subplot(121); hold on
        %pl1(i_epoch).FaceColor = 'r';
        bar(i_epoch,score_Xep(i_epoch),'r');

        
        % FIGURE 2 - - - - - - - - - - - - - - - - - - 
        epoch_data = epoch{1,i_epoch};
        %epoch_psd = pwelch(epoch_data,[],[],cfg.freqRange(1):cfg.freqRange(2),cfg.fs)';
        %epoch_psd = pwelch(epoch_data,[],[],cfg.freqRange(1):cfg.freqRange(end),cfg.fs)';
        epoch_psd_xch = squeeze(psd_chXep(:,i_epoch,:));

        % threshold level for selecting NOISY CHANNELS:
        % +/- n  standard deviation
        %<<<<<<<<<< TO SET <<<<<<<<<<<<<<<<<<<<
        %thresh_level = 3;  
        thresh_level = 2;
        %thresh_level = 1;
       
        lowfreq_thresh = 10;
        highfreq_thresh = 40;
        freq_thresh = [ lowfreq_thresh  highfreq_thresh];
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
        chanloc = eeg_struct.chanlocs;
        n_chan = length(chanloc);

        % - - - - -  - - - - - - -
        % compare the psd of of all the other channels in a single epoch)
        % compute the psd mean chan by chan
        n_freq = size(epoch_psd,2);
        psd_mean = mean(epoch_psd_xch,1);
        psd_std = std(epoch_psd_xch,1);
        %psd_mean = []; psd_std = []; %psd_median = [];
        %for i_freq = 1:n_freq
            %psd_median(i_freq) = median(power_spectrum(:,i_freq));
            %psd_mean(i_freq) = mean(epoch_psd(:,i_freq));
            %psd_std(i_freq) = std(epoch_psd(:,i_freq));
        %end

        thresh_xfreq = thresh_level * psd_std; % noise threshold = +/- standard deviation
        
        % identify bad channels - - - - - -  -
        bad_chan_table = zeros(n_chan, n_freq);
        for i_chan = 1:n_chan
            %i_chan = 67
            for i_freq = 1:n_freq
                % if power_spectrum(i_chan, i_freq) > (psd_median(i_freq)+thresh(i_freq)) || ...
                %        power_spectrum(i_chan, i_freq) < (psd_median(i_freq)-thresh(i_freq))
%                 if epoch_psd(i_chan, i_freq) > (psd_mean(i_freq)+thresh_xfreq(i_freq)) || ...
%                         epoch_psd(i_chan, i_freq) < (psd_mean(i_freq)-thresh_xfreq(i_freq))
%                    bad_chan_table(i_chan, i_freq) = 1;
                if epoch_psd_xch(i_chan, i_freq) > (psd_mean(i_freq)+thresh_xfreq(i_freq)) || ...
                        epoch_psd_xch(i_chan, i_freq) < (psd_mean(i_freq)-thresh_xfreq(i_freq))
                   bad_chan_table(i_chan, i_freq) = 1;
                end
            end
        end
        
        scroll_topoplot(epoch_data, bad_chan_table, chanloc, freq_thresh)   
        title([' epoch ' num2str(i_epoch) ])
    end
end



    
    
  
    
    
    
