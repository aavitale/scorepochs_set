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
%    data             - 2d array with the time-series (channels X time samples)
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

%function [idx_best_ep,epoch,score_Xep,score_chXep] = scorEpochs_av(cfg,data) 
function [idx_best_ep,epoch,score_Xep,score_chXep, psd_chXep, pxx] = scorEpochs_av2(cfg,data) 


    % divide the data in epochs of windowL length
    
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
    
    %<<<<<<< TO SET <<<<<<<<<<<<<<<<<<<<<<
    do_plot = 1;
    do_mean = 1; thresh_level = 1.5; %std deviation
    do_median = 0; thresh_level = 2; %mad
    line_width = 1.5;
          
    if do_plot
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
    end
    

    
    
    
  
    
    
    
