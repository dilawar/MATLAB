% 1. Generate dF/F calcium traces
% - Kambadur Ananthamurthy
% Always feed in the pixel resolution before manually drawing rois. Search for "@!"

addpath('/Users/ananth/Documents/MATLAB/imagingAnalysis_MATLAB')
addpath('/Users/ananth/Documents/MATLAB/CustomFunctions_MATLAB')

%Add JARs to the current dynamic class path
javaaddpath /Users/ananth/Documents/MATLAB/java/ij.jar
javaaddpath /Users/ananth/Documents/MATLAB/java/mij.jar

clear all
close all

%Operations (0 = Don't Perform, 1 = Perform)
getGeneralInformation = 1;
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
expcellnum = 150; % @! Change values for every dataset with an expected number of cell
% @! Always edit for every ROI dataset - the resolution listing is given by (Fiji is just) ImageJ
first=137;
second=244;

trialLength = 11;  % @! trial length in seconds - edit for every video untill a standard protocol is established

saveDirec = ('/Users/ananth/Desktop/Work/Analysis/ImageAnalysis/');
listPath = ('/Users/ananth/Desktop/Work/Imaging/analyze.txt');
fid = fopen(listPath);  % in analyze.txt, add url till the ROIindex folder name, ending, without "/"

while 1 %fid is valid when 1
    tline = fgetl(fid);
    
    if ~ischar(tline)
        break
    end
    
    direc = [tline '/'];
    
    if getGeneralInformation == 1
        %General Information:
        [mouse, sessionType, session] = getInfo(direc);
        sessionType_description = {'Control'; 'Trace'; 'Delay'};
        dataset = [mouse, '_SessionType', sessionType, '_Session' session];
        disp(dataset);
    end
    
    %tiftag = imfinfo([direc 'Trial1-ROI-1.tif']); %use the tiff file to be used to make cellmask
    tiftag = imfinfo([direc 'Trial1-ROI-1r.tif']); %use the tiff file to be used to make cellmask
    imageDescription = tiftag.ImageDescription;
    
    framei = findstr(imageDescription, 'images');
    nframes = str2num(imageDescription(framei+7:framei+9));
    %nframes=251;
    
    samplingrate = nframes/trialLength; % Sampling rate
    
    if manualROI==1
        %Opening 5 frames from first trial to form averaged first image, ROIs
        imagei = double(imread([direc 'Trial1-ROI-1r.tif'], 1));
        for i = 2:5
            image = double(imread([direc 'Trial1-ROI-1r.tif'], i));
            imagei = imagei + image;
        end   % adds intensity values of frames 1:5
        referenceimage = imagei./5;
        clear i
        
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

            cellsmarked=cellsmarked+1;
            %Just in case an ROI needs to be deleted
            if find(BWi)>0
                BW=BW+BWi;
            else
                cellsmarked = cellsmarked-1;
            end

            disp(['Cells marked: ' num2str(cellsmarked)]);
            
            %To know which ROIs have been marked
            referenceimagecopy(find(BWi))=0; % you can't mark a cell that you can't see
            reply = input('Quit adding more ROIs? ', 's'); % Y to quit;'                                          Enter to continue
            if numel(reply)>0
                break
            else
            end
            
        end  % Loop to mark all the cells
        clear i
        clear referenceimagecopy
        
        %Manually obtaining ROIs for cells
        rois(find(BW))=referenceimage(find(BW));
        rois(~find(BW))=0;
        
        referenceimagei = referenceimage;
        
        cellmask = bwlabel(rois);
        ncells = max(max(cellmask)); % Number of cells
        
        save([saveDirec dataset '/' 'cellmask.mat'], 'cellmask');
    else
        %load(['/Users/ananth/Desktop/Work/Analysis/ImageAnalysis/20131218/mouse4/Tone10knPuff-ROI/Cellmask.mat'])
        load([cellmaskDirec 'cellmask.mat'], 'cellmask');
        ncells = max(max(cellmask)); % Number of cells
        load([referenceImageDirec 'referenceimage.mat'], 'referenceimage');
    end
    
    if generateCal==1
        cal = zeros(max(max(cellmask)), nSessionTrials, nframes); % Cell Number, Trial Number, Frame Number
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
        save([saveDirec '/' dataset '/cal.mat'], 'cal');
    else
        load([saveDirec '/' dataset '/cal.mat'], 'cal');
        sessiontrial=nSessionTrials;
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
    puff_index=tone_index+blockSize;
    for i=1:nsets
        toneTrials(:,(count:(count+blockSize)-1),:)=calbdf(:,(tone_index:(tone_index+blockSize)-1),:);
        tone_index=tone_index+(2*blockSize);
        
        puffTrials(:,(count:(count+blockSize)-1),:)=calbdf(:,(puff_index:(puff_index+blockSize)-1),:);
        puff_index=puff_index+(2*blockSize);
        
        count=count+blockSize;
    end
    clear count
else
end

if saveData==1
    saveFolder = [saveDirec, '/', mouse '/' dataset '/'];
    %mkdir
    if ~isdir(saveFolder)
        mkdir(saveFolder);
    end
    
    save([saveFolder, dataset '/calbdf.mat'], 'calbdf');
    %save([saveFolder, dataset '/toneTrials.mat'], 'toneTrials');
    %save([saveFolder, dataset '/puffTrials.mat'], 'puffTrials');
    save([saveFolder, dataset '/' 'referenceimage.mat'], 'referenceimage');
else
end

frame=1:size(calbdf,3); %or just use nframes!
time=frame./samplingrate;
disp('Remember to use "squeeze" before plotting');

%HeatMap
figure(3); imagesc(squeeze(squeeze(100*(mean(calbdf,2))))); colormap(jet) %trial-averaged heat map
figure(4); plot(frame,squeeze(100*(calbdf(42,1,:))), 'LineWidth', 2); %separate cells
figure(5); plot(frame,squeeze(100*(calbdf(42,1,:))), 'LineWidth', 2); %separate cells