% Scorepochs 
%
% Function to select the best (most homogenoous) M/EEG epochs from a
% resting-state recordings. 
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
%    eeglab_dir          - directory with eeglab package
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

     

%     Copyright (C) 2020 Matteo Demuru, Matteo Fraschini
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

function [idx_best_ep,epoch,score_Xep,score_chXep, psd_chXep, pxx] = scorEpochs_set(cfg, set_file, eeglab_dir) 

    % ADD EEGLAB DIR and LOAD the .SET file
    if ~exist('set_file') || isempty('set_file')
        project_dir = 'D:\IIT\_PROJECT\SCORE_epoch'
        set_file = fullfile(project_dir, 'data', 'S001R01.set')
    end
    if ~exist('eeglab_dir') || isempty('eeglab_dir')
        eeglab_dir = 'D:\_TOOLBOX\eeglab_20201226'
        
        project_dir = 'D:\IIT\_PROJECT\SCORE_epoch'
        code_dir = fullfile(project_dir, 'code')
        addpath(genpath(code_dir));
    end

    cd(eeglab_dir);
    eeglab('nogui');

    eeg_struct = pop_loadset(set_file)
    % check channel location
    figure; topoplot([], eeg_struct.chanlocs, 'electrodes','labelpoint','chaninfo', eeg_struct.chaninfo);
    
    
    % CONFIGURATION file  
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
    
    
    % divide the data in epochs of windowL length
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
        
        score_chXep(c,:) = mean(score_ch);
       
    end    
    
    score_Xep = mean(score_chXep,1);
    
    [~,idx_best_ep] = sort(score_Xep,'descend');
    
    
    %<<<<<<< PLOT <<<<<<<<<<<<<<<<<<<<<<
    do_plot = 1;
    do_mean = 1; thresh_level = 1.5; %std deviation
    do_median = 0; thresh_level = 2; %mad
    line_width = 1.5;
          
    if do_plot
        % BAR GRAPH with schore epochs 
        figure; subplot(121); hold on; %subplot(2,5,1:4); hold on
        pl1 = bar(score_Xep);
        set(pl1,'FaceColor',[0.5 0.5 0.5]);
        xlim([1 length(score_Xep)])
        xlabel(['epochs (with window length of ' num2str(cfg.windowL) ' sec)'])
        ylabel('score x ALLChannels')
        ylim([0 1.1]);
        
        if do_mean
            score_mean = mean(score_Xep);
            yline(score_mean,'k-', 'LineWidth',line_width)
            yline(score_mean - thresh_level*std(score_Xep),'k--', 'LineWidth',line_width)
            yline(score_mean - thresh_level*1.5*std(score_Xep),'k:', 'LineWidth',line_width)
            title(['resting state scorepochs (mean = ' num2str(score_mean)])
        
        elseif do_median
            score_median = median(score_Xep);
            yline(score_median,'k-', 'LineWidth',line_width)
            yline(score_median - thresh_level*mad(score_Xep),'k--', 'LineWidth',line_width)
            yline(score_median - thresh_level*1.5*mad(score_Xep),'k:', 'LineWidth',line_width)
            title(['resting state scorepochs (median = ' num2str(score_median)])
            
        end
        
        % = = = = = = =  = = = = = = = = =  = = == =  ==  == = = = = 
        % identify BAD CHANNELS (with low score as compared to other channels:
        %figure; 
        subplot 122
        hold on; 
        
        if do_mean
            score_xchan_mean = mean(score_chXep,2);
            score_xchan_std = std(score_chXep,1,2);

            score_allchan_mean = mean(score_xchan_mean,1);
            score_allchan_std = std(score_xchan_mean,1,1);
                  
            errorbar(score_xchan_mean, score_xchan_std,'ko')
            
            hold on;
            yline(score_allchan_mean, 'k','LineWidth',line_width)
            yline(score_allchan_mean - thresh_level*score_allchan_std,'k--','LineWidth',line_width); 
            yline(score_allchan_mean - thresh_level*1.5*score_allchan_std, 'k:','LineWidth',line_width); 
            xlabel('channels')
            ylabel(['score (mean - ' num2str(thresh_level) ' std)'])
            %y_lim = ylim;
            ylim([0 1.1]);
           
            % find channel outliers:
            bad_chan_idx = find(score_xchan_mean <= (score_allchan_mean - thresh_level*score_allchan_std))';
            for ii = 1:length(bad_chan_idx)
                hold on
                errorbar(bad_chan_idx(ii), score_xchan_mean(bad_chan_idx(ii)), ...
                    score_xchan_std(bad_chan_idx(ii)),'rx', 'LineWidth', 2)
            end
            %title(['channel with score < ' num2str(thresh_level) 'STD = ' num2str(bad_chan_idx)])
            title(['channel with score < ' num2str(thresh_level) 'STD = '...
                    num2str(reshape(bad_chan_idx,[1,length(bad_chan_idx)]))])
        
        elseif do_median
            % MEDIAN - - - - - - - - - - - 
            score_xchan_median = median(score_chXep,2);
            score_xchan_mad = mad(score_chXep,1,2);

            score_allchan_median = median(score_xchan_median,1);
            score_allchan_mad = mad(score_xchan_median,1,1);

            hold on; %subplot 122
            errorbar(score_xchan_median, score_xchan_mad,'ko')
            hold on;
            yline(score_allchan_median,'k-','LineWidth',line_width)
            %yline(score_allchan_median + thresh_level*score_allchan_mad, '--'); 
            yline(score_allchan_median - thresh_level*score_allchan_mad, 'k--','LineWidth',line_width); 

            %yline(score_allchan_median + thresh_level*1.5*score_allchan_mad, ':'); 
            yline(score_allchan_median - thresh_level*1.5*score_allchan_mad, 'k:','LineWidth',line_width); 
            xlabel('channels')
            ylabel(['score (median - ' num2str(thresh_level) ' mad)'])
            %ylim([y_lim(1) 1]);
            ylim([0 1.1]);
           
            % find channel outliers:
            bad_chan_idx = find(score_xchan_median <= (score_allchan_median - thresh_level*1.5*score_allchan_mad))';
            for ii = 1:length(bad_chan_idx)
                hold on
                errorbar(bad_chan_idx(ii), score_xchan_median(bad_chan_idx(ii)), ...
                    score_xchan_mad(bad_chan_idx(ii)),'rx', 'LineWidth', 2)
            end
            title(['channel with score < ' num2str(thresh_level) 'MAD = ' num2str(bad_chan_idx)])
        
        end
        % -- - - - - - - - -  
            
            
        %% scroll_topoplot
        [ tmp, i_epoch ] = min(score_Xep);  %<<<<<<< epoch with the lowest score
        i_epoch = i_epoch(1);
        %     % or find epoch with score value closest to the thresh-level:
        %     score_median = median(score_Xep);
        %     score_thresh = score_median - thresh_level*1.5*mad(score_Xep);
        %     [ tmp, i_epoch ] = min(abs(score_Xep - score_thresh));

        % change color of the histogram - - - - - - - - - - - - 
        thresh_level = 2;

        %hist_plot = gcf;
        figure; hold on; %subplot(2,5,1:4); hold on
        line_width = 1.5;
        pl1 = bar(score_Xep);
        set(pl1,'FaceColor',[0.5 0.5 0.5]);
        xlim([1 length(score_Xep)])
        xlabel(['epochs (with window length of ' num2str(cfg.windowL) ' sec)'])
        ylabel('score x ALLChannels')
        ylim([0 1.1]);

        score_mean = mean(score_Xep);
        yline(score_mean,'k-', 'LineWidth',line_width)
        yline(score_mean - thresh_level*std(score_Xep),'k--', 'LineWidth',line_width)
        %yline(score_mean - thresh_level*1.5*std(score_Xep),'k:', 'LineWidth',line_width)
        title(['resting state scorepochs (mean = ' num2str(score_mean)])

        hold on
        %pl1(i_epoch).FaceColor = 'r';
        bar(i_epoch,score_Xep(i_epoch),'r');

        % - - - - - - - - - - - - - - - - - - 
        epoch_data = epoch{1,i_epoch};
        %epoch_psd = squeeze(psd_chXep(:,:,i_epoch));
        %epoch_psd = pwelch(epoch_data,[],[],cfg.freqRange(1):cfg.freqRange(2),cfg.fs)';
        epoch_psd = pwelch(epoch_data,[],[],2:80,cfg.fs)';

        % threshold level for selecting BAD CHANNELS:
        % +/- n  standard deviation
        %<<<<<<<<<< TO SET <<<<<<<<<<<<<<<<<<<<
        %thresh_level = 3;  
        thresh_level = 1.5;
        %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        lowfreq_thresh = 10;
        highfreq_thresh = 40;
        freq_thresh = [ lowfreq_thresh  highfreq_thresh];

        chanloc = eeg_struct.chanlocs;
        n_chan = length(chanloc);

        % - - - - -  - - - - - - -
        % FIGURE 3 (compared to the psd of all the other channels in a single epoch)
        % compute the median chan by chan
        n_freq = size(epoch_psd,2);
        psd_mean = []; psd_std = []; %psd_median = [];
        for i_freq = 1:n_freq
        %for i_chan = 1:size(power_spectrum,1)
            %psd_median(i_freq) = median(power_spectrum(:,i_freq));
            psd_mean(i_freq) = mean(epoch_psd(:,i_freq));
            psd_std(i_freq) = std(epoch_psd(:,i_freq));
        end

        thresh = thresh_level * psd_std; % noise threshold = +/- 1 standard deviation
        %thresh = thresh_level * 1.5*psd_std; % noise threshold = +/- 1 standard deviation

        % identify bad channels - - - - - -  -
        bad_chan_table = zeros(n_chan, n_freq);
        for i_chan = 1:n_chan
            %i_chan = 67
            for i_freq = 1:n_freq
                % if power_spectrum(i_chan, i_freq) > (psd_median(i_freq)+thresh(i_freq)) || ...
                %        power_spectrum(i_chan, i_freq) < (psd_median(i_freq)-thresh(i_freq))
                if epoch_psd(i_chan, i_freq) > (psd_mean(i_freq)+thresh(i_freq)) || ...
                        epoch_psd(i_chan, i_freq) < (psd_mean(i_freq)-thresh(i_freq))
                   bad_chan_table(i_chan, i_freq) = 1;
                end
            end
        end
        %bad_chXfreq_singlepoch = sum(bad_chan_table,2);

        %figure
        scroll_topoplot(epoch_data, bad_chan_table, chanloc, freq_thresh)   
        title([' epoch ' num2str(i_epoch) ])
    end
end



    
    
  
    
    
    
