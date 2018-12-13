
%% PPI calculation

%%% We used the Phase Preservation Index (PPI) to assess possible occurrence of
%%% phase reset of the ongoing µ-oscillation (8-12 Hz) after TMS
%%% PPI was implemented as described in: 

%%% Mazaheri A, Jensen O. 2006. Posterior alpha activity is not phase-reset by visual stimuli.
%%% Proc Natl Acad Sci U S A. 103:2948-2952.

%%% Implemented by Debora Desideri 
%%% Department of Neurology & Stroke, and Hertie Institute for Clinical Brain Research,
%%% University of Tübingen, Germany

%%% Last tested on 13th Dec 2018 with MATLAB R2018a
    

%% PPI
%%% NB An implicit assumption of PPI is that the pre-stimulus phase is
%%% random across trials

Nsbj = 12; %%% number of subjects 
Neegchan = 64; %% number of EEG channels


%%% INST_PHASE =  instantaneous µ-phase angle from -0.5 s to 1 s relative to TMS
%%% INST_PHASE is a cell array with Nsbj cells.
%%% In each  cell there is a matrix with dimensions: trials x EEG channels x time points

%%% time_ax = time axes of µ-phase 
%%% in our case time_ax = linspace(-0.5,1,size(INST_PHASE{1},3));

for sbj = 1:Nsbj
   
        for ch = 1:Neegchan 
            
            
            data = squeeze(INST_PHASE{sbj}(:,ch,:));
            
            t_ref = -0.2; %%% reference time, 200 ms before TMS
            [~,idx] = min(abs(time_ax-t_ref));
            phase_ref = data(:,idx); %%% phase at t_ref in each trial for the specific EEG channel
            
            ppi_latency = -0.1:0.1:0.8; %%% latencies at which we want to compute PPI (= every 100 ms from t_ref)
            for tt = 1:length(ppi_latency)
                [~,idx] = min(abs(time_ax-ppi_latency(tt)));
                phase_post_stim(:,tt) = data(:,idx);%%% phase at post_stim_latency in each trial for the specific EEG channel
                
            end
            
            %%% Compute PPI
            for tt = 1:length(ppi_latency)
                PPI(sbj,tt,ch) = abs(mean(exp(1i*(phase_ref-phase_post_stim(:,tt)))));
                PPI_Z(sbj,tt,ch) = size(data,1)*power(PPI(sbj,tt,ch),2); %%% z-tranformed PPI
            end
            
            clear   phase_ref   phase_post_stim
            
            
        end %%% channels
        
        Ntrials(sbj) = size(data,1); %%% store numver of trials for computaiton of PPI critical value
 
end %%% subjects

%% PPI for data shuffled in time 
%%% NB Only post-stimulus data are shuffled becasue we want to simulate
%%% phase reset, which occurs upon stimulus application

for sbj = 1:Nsbj
  
     
        for ch = 1:64
           
            data = squeeze(INST_PHASE{sbj}(:,ch,:));
            
             t_stim = 0; %%% time point of TMS 
            [~,idx_stim] = min(abs(time_ax-t_stim));
            
            data_prestim= data(:,1:idx_stim-1);
            data_poststim = data(:,idx_stim:end);
            
            
            data_shuffled = [];
            
            for rr = 1:100 %%% rr = number of randomization
                
                t_ref = -0.2;
                [~,idx] = min(abs(time_ax-t_ref));
                phase_ref = data(:,idx);
                
     
                for tr = 1:size(data,1) %%% tr = number of trial
     
                    data_shuffled(tr,:) = [data_prestim(tr, :) data_poststim(tr,randperm(size(data_poststim,2)))];
                end
                
                
                ppi_latency = -0.1:0.1:0.8;
                for tt = 1:length(ppi_latency)
                    [~,idx] = min(abs(time_ax-ppi_latency(tt)));
                    phase_post_stim(:,tt) = data_shuffled(:,idx);
                    
                end
                
                %%% PPI
                for tt = 1:length(ppi_latency)
                    PPI_shuffled(rr,sbj,tt,ch) = abs(nanmean(exp(1i*(phase_ref-phase_post_stim(:,tt)))));
                    PPI_Z_shuffled(rr,sbj,tt,ch) = size(data,1)*power(PPI_shuffled(rr,sbj,tt,ch),2);
                end
                
                clear   phase_ref   phase_post_stim
                
                
            end %% randomization
      
            
        end %%% channels
   
end %%% subjects



  
%% Example plot

%%% Critical PPI computed with Rayleigh’s test for circular uniformity 
%%% Fisher NI. 1993. Statistical analysis of circular data. Cambridge: Cambridge University Press.
%%% Zar JH. 1999. Biostatistical analysis. Upper Saddle River, N.J.: Prentice Hall.


  pvalue = 0.01; %%% significance threshold
  PPI_critical = sqrt((-log(pvalue)*sqrt(length(Ntrials)))./sum(Ntrials)); 
  
  
ch = 5; %%% select channel to plot

 ppi_latency = -0.1:0.1:0.8;
 t_ref = -0.2;
figure
hold on 

  line_1 =  plot(ppi_latency,squeeze(mean(PPI(:,:,ch))),'Color','g');
  line_2 = plot(ppi_latency,mean(squeeze(mean(squeeze(PPI_shuffled(:,:,:,ch))))),'-','Color','b');
  
   xlim([-0.3 0.9])
  ylim([0 0.5])
  
  line_3 = plot(xlim,[PPI_critical PPI_critical],'Color','r');
  line_4 = plot([t_ref t_ref],ylim,'Color','k');
 
  xlabel('Time (s)')
  ylabel('PPI')
  title(['Channel ' num2str(ch)])
  legend([line_1; line_2; line_3; line_4],{'PPI'; 'PPI shuffled'; 'PPI critical'; 'Reference time'},'Location','northeast')


