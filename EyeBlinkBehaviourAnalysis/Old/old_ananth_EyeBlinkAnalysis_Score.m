% AUTHOR - Kambadur Ananthamurthy
% PURPOSE - Read saved raw data and analyse

clear all
close all

baselineCorrection = 0;
aoc = 0;
saveFigures = 1;


samplingRate = 100; %Hz
samplePoint = samplingRate/1000;
preTime = 5000; %ms
csTime = 350; %ms
usTime = 50; %ms

%mySessionType = 'NPControl';
mySessionType = 'TEC';

nAnimals = 6;
nSessions = 1;

saveDirec = ('/Users/ananth/Desktop/Work/Analysis/eyeBlinkBehaviourAnalysis/');
fontSize = 10;
for animal = 1:nAnimals
    Z = 28;
    
    %Identify the sessionType and related parameters
    if strcmp(mySessionType,'TEC')
        if animal<4
            sessionType = 11;
        else
            sessionType = 9;
        end
        probe = 1;
    elseif strcmp(mySessionType,'NPControl')
        if animal>3 && animal<7
            sessionType = 5;
        else
            sessionType = 3;
        end
        probe = 0;
    else
        warning('Unable to detect sessionType!');
    end
    
    for session = 1:nSessions
        mouse = ['K' num2str(animal+Z)];
        if animal>3 && animal<7
            traceTime = 500; %ms
        else
            traceTime = 250; %ms
        end
        %session = 1;
        dataset = [mouse '_SessionType' num2str(sessionType) '_Session' num2str(session)];
        saveFolder = [saveDirec, 'Mouse' mouse '/' 'Mouse' dataset '/'];
        
        if baselineCorrection == 1
            csPlusTrials = csvread([saveFolder 'Mouse' dataset '_csPlus_dRR.csv']);
            csMinusTrials = csvread([saveFolder 'Mouse' dataset '_csMinus_dRR.csv']);
        else
            csPlusTrials = csvread([saveFolder 'Mouse' dataset '_csPlus.csv']);
            csMinusTrials = csvread([saveFolder 'Mouse' dataset '_csMinus.csv']);
        end
        
        if probe == 1
            probeTrials = csvread([saveFolder 'Mouse' dataset '_probeTrials.csv']);
            csPlusProbeTrials = csPlusTrials(probeTrials,:);
            csPlusTrials(probeTrials,:) = [];
        end
        
        %Time points
        csOnset = preTime*samplePoint;
        usOnset = (preTime+csTime+traceTime)*samplePoint;
        usOffset = usOnset+(usTime*samplePoint);
        
        if aoc == 1
            %Area Under The Curve
            %CS+
            aoc_csPlusTrials = trapz(csPlusTrials(:,csOnset:usOnset)');
            
            %Probe
            if probe ==1
                aoc_csPlusProbeTrials = trapz(csPlusProbeTrials(:,csOnset:usOnset)');
            end
            
            %CS-
            aoc_csMinusTrials = trapz(csMinusTrials(:,csOnset:usOnset)');
            
            figure(3)
            if probe == 1
                plot(aoc_csPlusProbeTrials,'g')
                hold on
            end
            plot(aoc_csPlusTrials,'b')
            hold on
            plot(aoc_csMinusTrials,'r')
            title([mouse ' ST' num2str(sessionType) ' S' num2str(session)],...
                'FontSize',fontSize,...
                'FontWeight','bold');
            xlabel('Trials',...
                'FontSize',fontSize,...
                'FontWeight','bold');
            xt = get(gca, 'XTick');
            set(gca, 'FontSize', fontSize);
            ylabel('AOC',...
                'FontSize',fontSize,...
                'FontWeight','bold');
            yt = get(gca, 'YTick');
            set(gca, 'FontSize', fontSize);
            %axis([200 800])
        end
        
        if strcmp(mySessionType,'TEC')
            %Peak Ratios
            
            %CS+
            for trial = 1:size(csPlusTrials)
                peak_CR_csPlus(session,trial) = median(csPlusTrials(trial,csOnset:usOnset)-median(csPlusTrials(trial,:),2));
                peak_UR = median(min(csPlusTrials(:,(usOnset:usOffset)),[],2)-median(csPlusTrials(trial,:),2));
                crurRatio_csPlus(session, trial) = abs(peak_CR_csPlus(session,trial)/peak_UR)*100;
            end
            %CS-
            for trial = 1:size(csMinusTrials)
                peak_CR_csMinus(session,trial) = median(csMinusTrials(trial,csOnset:usOnset)-median(csMinusTrials(trial,:),2));
                crurRatio_csMinus(session,trial) = abs(peak_CR_csMinus(session,trial)/peak_UR)*100;
            end
            if probe == 1
                %Probe
                for trial = 1:length(probeTrials)
                    peak_CR_csPlusProbe(session,trial) = median(csPlusProbeTrials(trial,csOnset:usOnset)-median(csPlusProbeTrials(trial,:),2));
                    crurRatio_csPlusProbe(session,trial) = abs(peak_CR_csPlusProbe(session,trial)/peak_UR)*100;
                end
            end
        end
        
        if strcmp(mySessionType,'TEC')
            figure(1)
            %CS+ Trials
            if probe == 1
                subplot(3,3,1:2)
            else
                subplot(2,3,1:2)
            end
            imagesc(csPlusTrials);
            title([mouse ' ST' num2str(sessionType) ' S' num2str(session) ' (AU)'],...
                'FontSize',fontSize+2,...
                'FontWeight','bold');
            xlabel('Time/10 ms',...
                'FontSize',fontSize,...
                'FontWeight','bold');
            xt = get(gca, 'XTick');
            set(gca, 'FontSize', fontSize);
            ylabel('CS+ Trials',...
                'FontSize',fontSize,...
                'FontWeight','bold');
            yt = get(gca, 'YTick');
            set(gca, 'FontSize', fontSize);
            caxis([200 800])
            colorbar;
            colormap(jet);
            
            %CS+ (Probe) Trials
            if probe == 1
                subplot(3,3,4:5)
                
                imagesc(csPlusProbeTrials);
                %         title([mouse ' ST' num2str(sessionType) ' S' num2str(session)],...
                %             'FontSize',fontSize,...
                %             'FontWeight','bold');
                xlabel('Time/10 ms',...
                    'FontSize',fontSize,...
                    'FontWeight','bold');
                xt = get(gca, 'XTick');
                set(gca, 'FontSize', fontSize);
                ylabel('CS+ (Probe) Trials',...
                    'FontSize',fontSize,...
                    'FontWeight','bold');
                yt = get(gca, 'YTick');
                set(gca, 'FontSize', fontSize);
                caxis([200 800])
                colorbar;
                colormap(jet);
            end
            
            %CS- Trials
            if probe == 1
                subplot(3,3,7:8)
            else
                subplot(2,3,4:5)
            end
            imagesc(csMinusTrials);
            %         title([mouse ' ST' num2str(sessionType) ' S' num2str(session)],...
            %             'FontSize',fontSize,...
            %             'FontWeight','bold');
            xlabel('Time/10 ms',...
                'FontSize',fontSize,...
                'FontWeight','bold');
            xt = get(gca, 'XTick');
            set(gca, 'FontSize', fontSize);
            ylabel('CS- Trials',...
                'FontSize',fontSize,...
                'FontWeight','bold');
            yt = get(gca, 'YTick');
            set(gca, 'FontSize', fontSize);
            caxis([200 800])
            colorbar;
            colormap(jet);
            
            %CR/UR (%)
            %CS+ Trials
            subplot(3,3,3)
            plot(crurRatio_csPlus(session,:));
            title(['Session' num2str(session)],...
                'FontSize',fontSize+2,...
                'FontWeight','bold');
            %xlabel(['Session' num2str(session)],...
            xlabel('CS+ Trials',...
                'FontSize',fontSize,...
                'FontWeight','bold');
            %xt = get(gca, 'XTick');
            set(gca, 'FontSize', fontSize);
            ylabel('Peak CR/UR (%)',...
                'FontSize',fontSize,...
                'FontWeight','bold');
            %yt = get(gca, 'YTick');
            set(gca, 'FontSize', fontSize);
            axis([0 size(csPlusTrials,1) 0 max(crurRatio_csPlus(session,:))])
            
            %CS+ (Probe) Trials
            if probe == 1
                subplot(3,3,6)
                plot(crurRatio_csPlusProbe(session,:));
                %xlabel(['Session' num2str(session)],...
                xlabel('CS+ (Probe) Trials',...
                    'FontSize',fontSize,...
                    'FontWeight','bold');
                %xt = get(gca, 'XTick');
                set(gca, 'FontSize', fontSize);
                ylabel('Peak CR/UR (%)',...
                    'FontSize',fontSize,...
                    'FontWeight','bold');
                %yt = get(gca, 'YTick');
                set(gca, 'FontSize', fontSize);
                axis([0 size(csPlusProbeTrials,1) 0 max(crurRatio_csPlusProbe(session,:))])
            end
            
            %CS- Trials
            if probe == 1
                subplot(3,3,9)
            else
                subplot(2,3,6)
            end
            plot(crurRatio_csMinus(session,:));
            %xlabel(['Session' num2str(session)],...
            xlabel('CS- Trials',...
                'FontSize',fontSize,...
                'FontWeight','bold');
            %xt = get(gca, 'XTick');
            set(gca, 'FontSize', fontSize);
            ylabel('Peak CR/UR (%)',...
                'FontSize',fontSize,...
                'FontWeight','bold');
            %yt = get(gca, 'YTick');
            set(gca, 'FontSize', fontSize);
            axis([0 size(csMinusTrials,1) 0 max(crurRatio_csMinus(session,:))])
        else
            figure(1)
            %CS+ Trials
            subplot(2,2,1:2)
            imagesc(csPlusTrials);
            title([mouse ' ST' num2str(sessionType) ' S' num2str(session) ' (AU)'],...
                'FontSize',fontSize+2,...
                'FontWeight','bold');
            xlabel('Time/10 ms',...
                'FontSize',fontSize,...
                'FontWeight','bold');
            xt = get(gca, 'XTick');
            set(gca, 'FontSize', fontSize);
            ylabel('CS+ Trials',...
                'FontSize',fontSize,...
                'FontWeight','bold');
            yt = get(gca, 'YTick');
            set(gca, 'FontSize', fontSize);
            caxis([200 800])
            colorbar;
            colormap(jet);
            
            %CS- Trials
            subplot(2,2,3:4)
            imagesc(csMinusTrials);
            %         title([mouse ' ST' num2str(sessionType) ' S' num2str(session)],...
            %             'FontSize',fontSize,...
            %             'FontWeight','bold');
            xlabel('Time/10 ms',...
                'FontSize',fontSize,...
                'FontWeight','bold');
            xt = get(gca, 'XTick');
            set(gca, 'FontSize', fontSize);
            ylabel('CS- Trials',...
                'FontSize',fontSize,...
                'FontWeight','bold');
            yt = get(gca, 'YTick');
            set(gca, 'FontSize', fontSize);
            caxis([200 800])
            colorbar;
            colormap(jet);
        end
        
        if saveFigures == 1
            if ~isdir(['/Users/ananth/Desktop/Figures/Mouse' mouse])
                mkdir(['/Users/ananth/Desktop/Figures/Mouse' mouse]);
            end
            figure(1);
            if strcmp(mySessionType,'TEC')
                print(['/Users/ananth/Desktop/Figures/Mouse' mouse '/crurRatio_' mouse ...
                    '_' mySessionType ...
                    '_Session' num2str(session)],...
                    '-djpeg');
            else
                print(['/Users/ananth/Desktop/Figures/Mouse' mouse '/rawBlinks_' mouse ...
                    '_' mySessionType ...
                    '_Session' num2str(session)],...
                    '-djpeg');
            end
        end
    end
    
    %Boxplots of the CR/UR Ratios
    if probe == 1
        figure(2);
        %CS+ Trials
        if probe == 1
            subplot(1,3,1)
        else
            subplot(1,2,1)
        end
        boxplot(crurRatio_csPlus')
        title(['Mouse' mouse ...
            ' CS+ Trials'],...
            'FontSize',fontSize,...
            'FontWeight','bold');
        xlabel('Session',...
            'FontSize',fontSize,...
            'FontWeight','bold');
        %xt = get(gca, 'XTick');
        set(gca, 'FontSize', fontSize);
        ylabel('Peak CR/UR (%)',...
            'FontSize',fontSize,...
            'FontWeight','bold');
        %yt = get(gca, 'YTick');
        set(gca, 'FontSize', fontSize);
        
        %CS+ (Probe) Trials
        if probe == 1
            subplot(1,3,2)
            boxplot(crurRatio_csPlusProbe')
            title(['Mouse' mouse ...
                ' CS+ (Probe) Trials'],...
                'FontSize',fontSize,...
                'FontWeight','bold');
            xlabel('Session',...
                'FontSize',fontSize,...
                'FontWeight','bold');
            %xt = get(gca, 'XTick');
            set(gca, 'FontSize', fontSize);
            ylabel('Peak CR/UR (%)',...
                'FontSize',fontSize,...
                'FontWeight','bold');
            %yt = get(gca, 'YTick');
            set(gca, 'FontSize', fontSize);
        end
        
        %CS- Trials
        if probe == 1
            subplot(1,3,3)
        else
            subplot(1,2,2)
        end
        boxplot(crurRatio_csMinus')
        title(['Mouse' mouse ...
            ' CS- Trials'],...
            'FontWeight','bold');
        xlabel('Session',...
            'FontSize',fontSize,...
            'FontWeight','bold');
        %xt = get(gca, 'XTick');
        set(gca, 'FontSize', fontSize);
        ylabel('Peak CR/UR (%)',...
            'FontSize',fontSize,...
            'FontWeight','bold');
        %yt = get(gca, 'YTick');
        set(gca, 'FontSize', fontSize);
        
        if saveFigures == 1
            figure(2)
            print(['/Users/ananth/Desktop/Figures/Mouse' mouse '/boxplot_'...
                mouse '_' mySessionType ...
                '_tillSession' num2str(session)],...
                '-djpeg');
        end
    end
end
%pause(5);
beep;
disp('Done!');
%close all