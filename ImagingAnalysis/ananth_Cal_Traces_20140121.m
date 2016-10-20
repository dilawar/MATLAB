% 1. Generate dF/F calcium traces
% - Kambadur Ananthamurthy

clear all
close all

% always feed in the pixel resolution before manually drawing rois
% expcellnum, nstims, nsets, nsessiontrials, etc.. Search for "@!"

% NOTE: In case of single trials saved using SkullVIEW.vi rather than
% SenseVIEW.vi, remember to manually add the "nframes" as "img_desc =
% tiftag.ImageDescription;" doesn't work.

%Add JARs to the current dynamic class path
javaaddpath /Users/ananth/Documents/MATLAB/java/ij.jar
javaaddpath /Users/ananth/Documents/MATLAB/java/mij.jar

% Execution modules
manualROI=1; % boolean - manually mark out ROIs/cells or load from saved file
generateCal=1; % boolean - calculate "cal" or load from saved file
motionCorrection=1; % 0 - no correction, 1 - using 2D normalized crossed correlation, 2 - Turboreg
plotEachCell=0; % boolean
skipCells=0; % boolean
plotEachTrial=0; % boolean
skipTrials=0; % boolean
bundleTrials=0; % 0 - skip bundling, 1 - alternate trials are identical, 2 - Set = "blocksize" or 5 tones then "blocksize" or 5 puff
saveData=0; % boolean

% Important initializations
expcellnum = 150; % @! Change values for every dataset with an expected number of cells
nstims=1;%10; % @! Number of defined stims
nsets=1%5; % @! Number of sets of trials for a particular stim
blocksize=1%5; % @! Number of tone/puff trials in a set
nsessiontrials=1%50; % @! Total number of trials for the session - all stims, all sets..

% @! Always edit for every ROI dataset - the resolution listing is given by (Fiji is just) ImageJ
first=189;
second=176;

triallength = 10;  % @! trial length in seconds - edit for every video untill a standard protocol is established

tau = 0.5;  % time-course of calcium transient in seconds
%align_ovsampl_factor = 16;
save_direc = ('/Users/ananth/Desktop/Work/Analysis/ImageAnalysis');

fid = fopen('/Users/ananth/Desktop/Work/Imaging/analyze.txt');  % in analyze.txt, add url till the ROIindex folder name, ending, without "/"
while 1 %fid is valid when 1
    direc = fgetl(fid);
    
    if ~ischar(direc)
        break
    else
    end
    
    %     dataset_date=('20131218');
    date = findstr(direc, '201');
    dataset_date = direc(date:date+7);
    
    %     dataset_mousename=('mouse4');
    %uscorei=findstr(direc, '_');
    fslashi=findstr(direc, '/'); %change made; used '/' instead of '_'
    dataset_mousename=direc((fslashi(7)+1:length(direc)-4));
    
    %     dataset_ROIindex=('cells-spnt2-ROI-1r');
    %dataset_ROIindex=direc((fslashi+1):length(direc));
    %dataset = [dataset_date '/' dataset_mousename '/' dataset_ROIindex];
    dataset = [dataset_date '/' dataset_mousename];
    mkdir([save_direc '/' dataset]);
    
    %tiftag = imfinfo([direc '/1_Trial-1_ROI-1.tif']); %use the tiff file to be used to make cellmask
    tiftag = imfinfo([direc '/1_Trial-1_ROI-1r.tif']); %use the tiff file to be used to make cellmask
    img_desc = tiftag.ImageDescription;
    
    framei = findstr(img_desc, 'Frames');
    %nframes = str2num(img_desc((framei+9):(framei+12)));
    nframes=251;
    
    samplingrate = nframes/triallength; % Sampling rate
    
    if manualROI==1
        %Opening 5 frames from first trial to form averaged first image, ROIs
        %imagei = double(imread([direc 'cells4-ROI-1.tif'], 1));
        
        imagei = double(imread([direc '/1_Trial-1_ROI-1r.tif'], 1));
        for w = 2:5
            image = double(imread([direc '/1_Trial-1_ROI-1r.tif'], w));
            imagei = imagei + image;
        end   % adds intensity values of frames 1:5
        referenceimage = imagei./5;
        save([save_direc '/' dataset '/' 'referenceimage.mat'], 'referenceimage');
        
        %Manually Drawing ROIs for cells
        BW=zeros(second,first);
        rois=zeros(second,first);
        
        referenceimagecopy=referenceimage; %Only to help make ROIs - This value shall be cleared out after use
        
        cellsmarked=0; % only for initialization
        for i=1:((expcellnum)-1)
            
            imagesc(referenceimagecopy)
            colormap(gray)
            
            %Make ROIs in the grayscale image popping up
            BWi=roipoly;
            close all
            
            %Just in case an ROI needs to be deleted
            if find(BWi)>0
                BW=BW+BWi;
            else
                i=i-1;
            end
            
            cellsmarked=cellsmarked+1 % to know how many cells have been marked - no ";" to display value
            
            %To know which ROIs have been marked
            referenceimagecopy(find(BWi))=0; % you can't mark a cell that you can't see
            reply = input('Quit adding more ROIs? ', 's'); % Y to quit;'                                          Enter to continue
            if numel(reply)>0
                break
            else
            end
            
        end  % Loop to mark all the cells
        
        clear referenceimagecopy
        
        %Manually obtaining ROIs for cells
        
        rois(find(BW))=referenceimage(find(BW));
        rois(~find(BW))=0;
        
        referenceimagei = referenceimage;
        
        cellmask = bwlabel(rois);
        ncells = max(max(cellmask)); % Number of cells
        
        save([save_direc '/' dataset '/' 'cellmask.mat'], 'cellmask');
    else
        %load(['/Users/ananth/Desktop/Work/Analysis/ImageAnalysis/20131218/mouse4/Tone10knPuff-ROI/Cellmask.mat'])
        
        load([save_direc '/' dataset '/cellmask.mat'])
        ncells = max(max(cellmask)); % Number of cells
        load([save_direc '/' dataset '/' 'referenceimage.mat'], 'referenceimage');
    end
    
    if generateCal==1
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
                    if motionCorrection==0
                        img(frame,:,:)=tempimg;
                        
                    elseif motionCorrection==1
                        % Motion Correction 1
                        % Register img to refimg frame by frame by cross
                        % correlation
                        cc = normxcorr2(tempimg, referenceimage);
                        [max_cc, imax] = max(abs(cc(:)));
                        [ypeak, xpeak] = ind2sub(size(cc),imax(1));
                        corr_offset = [(xpeak-size(referenceimage,2)) (ypeak-size(referenceimage,1))];
                        
                        tform = maketform('affine',[1 0 0; 0 1 0; corr_offset(1) corr_offset(2) 1]);
                        img (frame, :, :) = imtransform(tempimg, tform, 'XData', [1 size(referenceimage, 2)], 'YData', [1 size(referenceimage, 1)]);
                        
                        imagebundle=[frame '.tif'];
                        %save([save_direc '/' dataset '/' 'image/' imagebundle], 'img');
                    else
                        % Motion Correction 2 - Figure out the deal with TurboReg and mij
                        %for Mac OS.
                        
                        temp = imread([direc '1_Trial-1_ROI-1r.tif'], frame); % Temporary Image
                        
                        
                        % Register img to referenceimage using Turboreg in ImageJ
                        MIJ.createImage(['Img' num2str(frame)], int16(tempimg), 1);
                        
                        MIJ.run('TurboReg', ['-align -window Img', num2str(frame), ' 0 0 ', num2str(tiftag(1).Width-1),...
                            ' ', num2str(tiftag(1).Height-1), ' -window Refimg 0 0 ', ...
                            num2str(tiftag(1).Width-1), ' ', num2str(tiftag(1).Height-1),...
                            ' -translation ', num2str(tiftag(1).Width/2), ' ', ...
                            num2str(tiftag(1).Height/2), ' ', num2str(tiftag(1).Width/2), ...
                            ' ', num2str(tiftag(1).Height/2), ' -showOutput']);
                        
                        temp = MIJ.getImage ('Output');
                        img (frame, :, :) = squeeze(uint16(temp(:, :, 1)));
                        MIJ.run('Close', 'Output');
                        MIJ.run('Close', ['Img' num2str(frame)]);
                        
                    end
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
    else
        load([save_direc '/' dataset '/cal.mat'], 'cal');
        sessiontrial=nsessiontrials;
        %sessiontrial=60;
    end
    % NOTE: all session trials should have been run through, by this step
    
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
            
        end  % calbdf
    end  % for all cells, c
end

if plotEachCell==1
    % Plotting every cell for the trial
    for j=1:ncells
        figure(2);
        plot(frame,squeeze(calbdf(j,1,:)), 'LineWidth', 2);
        disp(['cell number:' num2str(j)]);
        %hold on
        reply = input('Quit plotting cells? ', 's'); %  Y to quit; Enter to continue
        if numel(reply)>0
            break
        else
        end
    end
else
end

%To get rid of some cells
if skipCells==1
    for i = 1:ncells
        edit_calbdf = calbdf;
        edit_calbdf([24, 70], :,:) = [];
        
        edit_cellmask=cellmask;
        edit_cellmask(find(24))=0;
        edit_cellmask(find(70))=0;
        %clear('y');   % distract the JIT
    end
    
    calbdf=edit_calbdf;
    cellmask=edit_cellmask;
    ncells=size(calbdf,1);
else
end

if plotEachTrial==1
    for i=1:60
        i
        figure(1);
        imagesc(squeeze(squeeze(calbdf(:,i,:))));
        reply = input('Quit looking at datasets? ', 's'); %  Y to quit; Enter to continue
        if numel(reply)>0
            break
        else
        end
    end
    clear i
end

%Trials2Skip=[3,11,13,18,21,25,26,27,29,30,31,34,36,37,38,44,45,51,53,54,55,57];

%Bundle trials categorically starting with a tone trial with every
if bundleTrials==1
    %alternate trial being identical
    tone_count=1;
    puff_count=1;
    for i=1:size(calbdf,2)
        %tone_session_mouse4_tonenpuff((1:ncells),tone,(1:nframes))=a.calbdf((1:ncells),tone_trials(tone),(1:nframes));
        if mod(i,2)==0 %checking if i is odd or even
            puffTrials(:,tone_count,:)=calbdf(:,i,:);
            tone_count=tone_count+1;
        else
            toneTrials(:,puff_count,:)=calbdf(:,i,:);
            puff_count=puff_count+1;
        end
    end
    clear tone_count
    clear puff_count
    
    if skipTrials==1
        calbdf(:,Trials2Skip,:)=[];
        
        a=[];
        b=[];
        for i=1:length(Trials2Skip)
            if mod(Trials2Skip(i),2)==0
                a=[a (Trials2Skip(i))/2];
            else
                b=[b (((Trials2Skip(i)-1)/2)+1)];
            end
        end
        clear i
        puffTrials(:,a,:)=[];
        toneTrials(:,b,:)=[];
        clear a
        clear b
    else
    end
    
elseif bundleTrials==2
    %"blocksize" or 5 tone trials, then "blocksize" or 5 puff trials repeated in sets
    %tone_session_mouse4_tonenpuff((1:ncells),tone,(1:nframes))=a.calbdf((1:ncells),tone_trials(tone),(1:nframes));
    count=1;
    tone_index=1;
    puff_index=tone_index+blocksize;
    for i=1:nsets
        toneTrials(:,(count:(count+blocksize)-1),:)=calbdf(:,(tone_index:(tone_index+blocksize)-1),:);
        tone_index=tone_index+(2*blocksize);
        
        puffTrials(:,(count:(count+blocksize)-1),:)=calbdf(:,(puff_index:(puff_index+blocksize)-1),:);
        puff_index=puff_index+(2*blocksize);
        
        count=count+blocksize;
    end
    clear count
else
end

if saveData==1
    save([save_direc '/' dataset '/calbdf.mat'], 'calbdf');
    %save([save_direc '/' dataset '/toneTrials.mat'], 'toneTrials');
    %save([save_direc '/' dataset '/puffTrials.mat'], 'puffTrials');
else
end

frame=1:size(calbdf,3); %or just use nframes!
time=frame./samplingrate;
disp('Remember to use "squeeze" before plotting');

%HeatMap
figure(3); imagesc(squeeze(squeeze(100*(mean(calbdf,2))))); colormap(jet) %trial-averaged heat map
figure(4); plot(frame,squeeze(100*(calbdf(42,1,:))), 'LineWidth', 2); %separate cells
figure(5); plot(frame,squeeze(100*(calbdf(42,1,:))), 'LineWidth', 2); %separate cells