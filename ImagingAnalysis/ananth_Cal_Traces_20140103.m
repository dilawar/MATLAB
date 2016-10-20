%clear all
close all

% always feed in the pixel resolution before manually drawing rois,

% expcellnum, nstims, nsets, nsessiontrials, etc.. Search for "@!"

expcellnum = 150; % @! Change values for every dataset with an expected number of cells
nstims=6; % @! Number of defined stims
nsets=10; % @! Number of sets of trials for a particular stim
nsessiontrials=60; % @! Total number of trials for the session - all stims, all sets..

% @! Always edit for every ROI dataset - the resolution listing is given by (Fiji is just) ImageJ
first=245;
second=192;

triallength = 30;  % @! trial length in seconds - edit for every video untill a standard protocol is established

fid = fopen('/Users/ananth/Desktop/Work/Imaging/analyze.txt');  % in analyze.txt, add url till the ROIindex folder name, ending, without "/"

tau = 0.5;  % time-course of calcium transient in seconds

%align_ovsampl_factor = 16;

save_direc = '/Users/ananth/Desktop/Work/Analysis/ImageAnalysis';

while 1 %fid is valid
    direc = fgetl(fid);
    
    if ~ischar(direc)
        break
    else
    end
    
    %     dataset_date=('20131218');
    date = findstr(direc, '201');
    dataset_date = direc(date:date+7);
    
    %     dataset_mousename=('mouse4');
    uscorei=findstr(direc, '_');
    dataset_mousename=direc((date+9):(uscorei-1));
    
    %     dataset_ROIindex=('cells-spnt2-ROI-1r');
    dataset_ROIindex=direc((uscorei+1):length(direc));
    dataset = [dataset_date '/' dataset_mousename '/' dataset_ROIindex];
    mkdir([save_direc '/' dataset]);
    
    tiftag = imfinfo([direc '/1_Trial-1_ROI-1.tif']); %use the tiff file to be used to make cellmask
    img_desc = tiftag.ImageDescription;
    
    framei = findstr(img_desc, 'Frames');
    nframes = str2num(img_desc((framei+9):(framei+12)));
    %nframes=266;
    
    samplingrate = nframes/triallength; % Sampling rate
    
    %     %Opening 5 frames from first trial to form averaged first image, ROIs
    %     %imagei = double(imread([direc 'cells4-ROI-1.tif'], 1));
    %
    %     imagei = double(imread([direc '/1_Trial-1_ROI-1r.tif'], 1));
    %     for w = 2:5
    %         image = double(imread([direc '/1_Trial-1_ROI-1r.tif'], w));
    %         imagei = imagei + image;
    %     end   % adds intensity values of frames 1:5
    %     referenceimage = imagei./5;
    %     save([save_direc '/' dataset '/' 'ReferenceImage.mat'], 'referenceimage');
    %
    %     %Manually Drawing ROIs for cells
    %     BW=zeros(second,first);
    %     rois=zeros(second,first);
    %
    %     referenceimagecopy=referenceimage; %Only to help make ROIs - This value shall be cleared out after use
    %
    %     cellsmarked=0; % only for initialization
    %     for i=1:((expcellnum)-1)
    %
    %         imagesc(referenceimagecopy)
    %         colormap(gray)
    %
    %         %Make ROIs in the grayscale image popping up
    %         BWi=roipoly;
    %         close all
    %
    %         %Just in case an ROI needs to be deleted
    %         if find(BWi)>0
    %             BW=BW+BWi;
    %         else
    %             i=i-1;
    %         end
    %
    %         cellsmarked=cellsmarked+1 % to know how many cells have been marked - no ";" to display value
    %
    %         %To know which ROIs have been marked
    %         referenceimagecopy(find(BWi))=0; % you can't mark a cell that you can't see
    %         reply = input('Do you want to quit adding more ROIs? (Y to quit; Enter to continue adding ROIs) [Y]: ', 's');
    %         if numel(reply)>0
    %             break
    %         else
    %         end
    %
    %     end  % Loop to mark all the cells
    %
    %     clear referenceimagecopy
    %
    %     %Manually obtaining ROIs for cells
    %
    %     rois(find(BW))=referenceimage(find(BW));
    %     rois(~find(BW))=0;
    %
    %     referenceimagei = referenceimage;
    %
    %     cellmask = bwlabel(rois);
    %     ncells = max(max(cellmask)); % Number of cells
    %
    %     save([save_direc '/' dataset '/' 'cellmask.mat'], 'cellmask');
    
    
    %load(['/Users/ananth/Desktop/Work/Analysis/ImageAnalysis/20131218/mouse4/Tone10knPuff-ROI/Cellmask.mat'])
    %load([save_direc '/' dataset '/Cellmask.mat'])
    %ncells = max(max(cellmask)); % Number of cells
    
    cal = zeros(max(max(cellmask)), nsessiontrials, nframes); % Cell Number, Trial Number, Frame Number
    caldf = cal;
    caldfsig = cal;
    calbdf = cal;
    calbdfsig = cal;
    
    % Looping over all the trials (nsessiontrials)
    sessiontrial=0; % for initialization
    for set = 1 : nsets % number of sets of trials for a particular stim
        for stim=1:nstims
            %for stim=1:2:nstims
            % Read in image sequence
            img = zeros(nframes, tiftag(1).Height, tiftag(1).Width);
            
            for frame = 1 : nframes
                tempimg = imread([direc '/' num2str(stim) '_Trial-' num2str(set) '_ROI-1r.tif'], frame);
                
                %In case Motion Correction is not being used
                %                 img(frame,:,:)=tempimg;
                
                % Motion Correction 1
                % Register img to refimg frame by frame by cross corr
                cc = normxcorr2(tempimg, referenceimage);
                [max_cc, imax] = max(abs(cc(:)));
                [ypeak, xpeak] = ind2sub(size(cc),imax(1));
                corr_offset = [(xpeak-size(referenceimage,2)) (ypeak-size(referenceimage,1))];
                
                tform = maketform('affine',[1 0 0; 0 1 0; corr_offset(1) corr_offset(2) 1]);
                img (frame, :, :) = imtransform(tempimg, tform, 'XData', [1 size(referenceimage, 2)], 'YData', [1 size(referenceimage, 1)]);
                
                imagebundle=[frame '.tif'];
                %save([save_direc '/' dataset '/' 'image/' imagebundle], 'img');
                
                % Motion Correction 2 - Figure out the deal with TurboReg and mij
                %for Mac OS.
                
                %             temp = imread([direc 'cells4-ROI-1.tif'], frame); % Temporary Image
                %
                %
                %                 % Register img to referenceimage using Turboreg in ImageJ
                %                 MIJ.createImage(['Img' num2str(frame)], int16(tempimg), 1);
                %
                %                 MIJ.run('TurboReg', ['-align -window Img', num2str(frame), ' 0 0 ', num2str(tiftag(1).Width-1),...
                %                     ' ', num2str(tiftag(1).Height-1), ' -window Refimg 0 0 ', ...
                %                     num2str(tiftag(1).Width-1), ' ', num2str(tiftag(1).Height-1),...
                %                     ' -translation ', num2str(tiftag(1).Width/2), ' ', ...
                %                     num2str(tiftag(1).Height/2), ' ', num2str(tiftag(1).Width/2), ...
                %                     ' ', num2str(tiftag(1).Height/2), ' -showOutput']);
                %
                %                 temp = MIJ.getImage ('Output');
                %             img (frame, :, :) = squeeze(uint16(temp(:, :, 1)));
                %                 MIJ.run('Close', 'Output');
                %                 MIJ.run('Close', ['Img' num2str(frame)]);
                %
                
            end
            
            % Generate intensity traces for all cells
            sessiontrial=sessiontrial+1 % adds all stims and then considers next set - no ";" to display value
            
            for c = 1 : max(max(cellmask)) % 2D matrix for all cells
                if min(min(min(img(:, cellmask == c)))) >= 0
                    cal(c, sessiontrial, :) = squeeze(mean(mean(img(:, cellmask == c), 3), 2));
                else
                    cal(c, sessiontrial, :) = zeros(nframes, 1) * NaN;
                end
            end
        end  % number of stims for a set
    end
    save([save_direc '/' dataset '/cal.mat'], 'cal');
    
    % NOTE: all session trials should have been run through, by this step
    
    caltr = zeros(size(cal)); % Calcium Traces for all Cells
    calbaseline = squeeze(mean(cal, 3)); % Specific baseline fluorescence value for each cell
    nstd = 2; % Standard Deviation for Excitatory signal
    nstdlow = 0.5; % Standard Deviation for (Hypothetical) Inhibitory signal
    min_consec_frames = 3; % min number of consecutive frames for a calcium transient
    dframe = 0;
    
    for t = 1 : sessiontrial
        for c = 1 : max(max(cellmask)) % or ncells
            caldf(c, t, :) = (cal(c, t, :) - calbaseline(c, t))/calbaseline(c, t); % Change in Fluorescence relative to baseline
            caldfsig(c, t, :) = caldf(c, t, :) > (nstd * std(squeeze(caldf(c, t, :)))); % Significant-only changes in Fluorescence
            
            % Baseline correction
            [maxtab mintab] = peakdetect(squeeze(caldf(c, t, :))); % modified peakdetect.m to supress graphical display
            temp = squeeze(cal(c, t, :));
            for m = 1 : size(maxtab, 1)
                temp(uint32(maxtab(m,1)-5)+1: maxtab(m,1)+30) = NaN;
            end
            
            %* 10 percentile value correction
            cabase = prctile(squeeze(cal(c, t, :)), 10);
            
            castd = nanstd(temp-cabase);
            
            calbdf(c, t, :) = ((cal(c, t, :) - cabase)/cabase); % Change in Fluorescence relative to baseline
            calbdfsig(c, t, :) = abs(calbdf(c, t, :)) > (nstd * castd/cabase); % Significant-only changes in Fluorescence
            
            % Generate signifcant traces only in caltr
            
            frame = 1;
            consec = 0;
            consecframe = 0;
            
            for frame = 1 : nframes
                if consec
                    tempstd = nstdlow;
                else
                    tempstd = nstd;
                end
                if calbdf(c, t, frame) > (tempstd * castd/cabase)
                    if ~consec
                        dframe = dframe + 1;
                    end
                    consec = 1;
                    caltr(c, t, frame) = calbdf(c, t, frame);
                    consecframe = consecframe + 1;
                else
                    % Impose min triallength
                    if consec && consecframe < min_consec_frames
                        caltr(c, t, frame-consecframe:frame-1) = 0;
                    end
                    consecframe = 0;
                    
                    consec = 0;
                    
                end
            end
            
        end  % calbdf
    end  % for all cells, c
end

%To get rid of some cells
% for i = 1:ncells
%   edit_calbdf = calbdf;
%   edit_calbdf([24, 70], :,:) = [];
%
%   edit_cellmask=cellmask;
%   edit_cellmask(find(24))=0;
%   edit_cellmask(find(70))=0;
%   %clear('y');   % distract the JIT
% end
%
% calbdf=edit_calbdf;
% cellmask=edit_cellmask);
% ncells=size(calbdf,1);

close all

% for i=1:60
%     i
%     imagesc(squeeze(squeeze(calbdf(:,i,:))))
%     reply = input('Do you want to quit looking at datasets? (Y to quit; Enter to continue adding ROIs) [Y]: ', 's');
%     if numel(reply)>0
%         break
%     else
%     end
% end
% clear i

skipTrials=[7,11,21,23,25,26,27,34,35,42,51,53,57,59];

%Bundle trials categorically starting with a tone trial with every
%alternate trial being identical
puff_count=1;
tone_count=1;
for i=1:size(calbdf,2)
    %tone_session_mouse4_tonenpuff((1:ncells),tone,(1:nframes))=a.calbdf((1:ncells),tone_trials(tone),(1:nframes));
    if mod(i,2)==0 %checking if i is odd or even
        if ismember(i,skipTrials)
            
        else
            puff_mouse6(:,puff_count,:)=calbdf(:,i,:);
            puff_count=puff_count+1;
        end
    else
        if ismember(i,skipTrials)
        else
            tone_mouse6(:,tone_count,:)=calbdf(:,i,:);
            tone_count=tone_count+1;
        end
    end
end
% for i=1:length(skipTrials)
%     calbdf(:,skipTrials(i),:)=[];
% end

% save([save_direc '/' dataset '/calbdf.mat'], 'calbdf');
% save([save_direc '/' dataset '/caltr.mat'], 'caltr');
% frame=1:size(calbdf,3); %or just use nframes!
% time=frame./samplingrate;
disp('Remember to use "squeeze" before plotting');

% % Plotting every cell for the trial
% for j=1:ncells
%     plot(frame,squeeze(calbdf(j,1,:)))
%     hold on
% end

%HeatMap
%imagesc(squeeze(squeeze(calbdf(:,1,:))))