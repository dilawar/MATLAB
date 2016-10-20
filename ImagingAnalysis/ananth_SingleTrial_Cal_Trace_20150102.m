clear all
close all

% always feed in the pixel resolution before manually drawing rois,
% expcellnum ntrials, etc.. Search for "@!"

% javaaddpath '/Users/kambadurananthamurthy/Documents/MATLAB/java/ij.jar'
% javaaddpath '/Users/kambadurananthamurthy/Documents/MATLAB/java/mij.jar'

tau = 0.5;  % time-course of calcium transient in seconds
triallength = 10;  % trial length in seconds - edit for every video untill a standard protocol is established
%align_ovsampl_factor = 16;

%Pixel Resolution of the .tif file being analysed in the form of first x second
first=189;
second=176;

expcellnum = 150; % @! Change values for every dataset with an expected number of cells

% loadcal = 1; % load cal instead of computing it %Good if cal.mat has
% already been created.

%fid = fopen('/Users/ananthamurthy/Desktop/Work/Imaging/analyze.txt');  % in analyze.txt, add url till the ROIindex folder name, ending, without "/"

%tau = 0.5;  % time-course of calcium transient in seconds

%while 1
%direc = fgetl(fid);

%if ~ischar(direc)
%    break
%else
%end
directory = '/Users/ananth/Desktop/cells5-ROI-1r.tif'; % Full path to the .tif file

tiftag = imfinfo(directory); %use the tiff file to be used to make cellmask
img_desc = tiftag.ImageDescription;
%
%     framei = findstr(img_desc, 'Frames');
nframes = 251; % @! Always edit to reflect the differences in tiftag for SkullVEIW and SenseVIEW
clear framei;

samplingrate = nframes/triallength; % Sampling rate

%Opening 5 frames from first trial to form averaged first image, ROIs
%imagei = double(imread([direc 'cells4-ROI-1.tif'], 1));

imagei = double(imread(directory, 1));
for i = 2:5
    image = double(imread(directory, i));
    imagei = imagei + image;
end
referenceimage = imagei./5;
%save([save_direc '/' dataset '/' 'referenceimage.mat'], 'referenceimage');
save('/Users/ananth/Desktop/referenceimage.mat', 'referenceimage');

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

%save([save_direc '/' dataset '/' 'cellmask.mat'], 'cellmask');
save('/Users/ananth/Desktop/cellmask.mat', 'cellmask');

%savedirec = mkdir(save_direc, dataset_date); %Parent Folder where all analysed data shall be saved according to imaging date
%save([save_direc '/' dataset '/' 'ReferenceImage.mat'], 'referenceimage');

%     %House-keeping!:)
%     %re-aligning ROI matrix to referenceimage, in case there has been some systematic shift
%     c = xcorr2(referenceimage, rois);
%     [maxr, maxcol] = find(c == max(max(c)));
%     col_lag = size(rois, 1) - maxr;
%     row_lag = size(rois, 2) - maxcol;
%
%     %Realigning rois as per lags detected
%     %For rows
%     if sign(row_lag) == -1
%         rois(:, ( (end-abs(row_lag) +1 ):end)) = [];
%         pad = zeros(size(rois, 1), abs(row_lag) );
%         rois = [pad, rois];
%     elseif sign(row_lag) == 1
%         rois(:, 1:(abs(row_lag))) = [];
%         pad = zeros(size(rois, 1), abs(row_lag) );
%         rois = [rois, pad];
%     elseif sign(row_lag) == 0
%     end
%
%     %For columns
%     if sign(col_lag) == -1
%         rois( ( (end-abs(col_lag)+1):end), :) = [];
%         pad = zeros(abs(col_lag), size(rois, 2) );
%         rois = [pad; rois];
%     elseif sign(col_lag) == 1
%         rois(1:(abs(col_lag)), :) = [];
%         pad = zeros(abs(col_lag), size(rois, 2) );
%         rois = [rois; pad];
%     elseif sign(col_lag) == 0
%     end
%
%     %close all
%     pause(0.1)      %this is to prevent MIJI from causing the program to hang - the dialog box needs a gap to close before MIJI can run
%
%     ref_img = referenceimage;
%     ref_ft = fft2(ref_img); %Fourier Transformation

ntrials=1; %change for every dataset
cal = zeros(max(max(cellmask)), ntrials, nframes); % Cell Number, Trial Number, Frame Number
caldf = cal;
caldfsig = cal;
calbdf = cal;
calbdfsig = cal;

%     MIJ.start ('F:\CSHL\Code\Fiji.app\');  % Start ImageJ
%     MIJ.createImage('ReferenceImage', int16(referenceimage), 1);

%     fid2=fopen(direc);
%
%     while 1
%         tiffid=fgetl(fid2);
%         if ~ischar(tiffid)
%             break
%         else
%         end
for trial = 1 : ntrials
    
    %             trial
    %         for r = 1 : ncells % so that all the ROIs are dealt with
    
    %             Read in image sequence.
    img = zeros(nframes, tiftag(1).Height, tiftag(1).Width);
    
    for count = 1 : nframes
        tempimg = imread(directory, count);
        
        %In case Motion Correction is not being used
        %                 img(count,:,:)=tempimg;
        
        % Motion Correction 1
        % Register img to refimg frame by frame by cross corr
        cc = normxcorr2(tempimg, referenceimage);
        [max_cc, imax] = max(abs(cc(:)));
        [ypeak, xpeak] = ind2sub(size(cc),imax(1));
        corr_offset = [(xpeak-size(referenceimage,2)) (ypeak-size(referenceimage,1))];
        
        tform = maketform('affine',[1 0 0; 0 1 0; corr_offset(1) corr_offset(2) 1]);
        img (count, :, :) = imtransform(tempimg, tform, 'XData', [1 size(referenceimage, 2)], 'YData', [1 size(referenceimage, 1)]);
        
        imagebundle=[count '.tif'];
        %                 save([save_direc '/' dataset '/' 'image/' imagebundle], 'img');
        
        % Motion Correction 2 - Figure out the deal with TurboReg and mij
        %for Mac OS.
        
        %             temp = imread([direc 'cells4-ROI-1.tif'], count); % Temporary Image
        %
        %
        %                 % Register img to referenceimage using Turboreg in ImageJ
        %                 MIJ.createImage(['Img' num2str(count)], int16(tempimg), 1);
        %
        %                 MIJ.run('TurboReg', ['-align -window Img', num2str(count), ' 0 0 ', num2str(tiftag(1).Width-1),...
        %                     ' ', num2str(tiftag(1).Height-1), ' -window Refimg 0 0 ', ...
        %                     num2str(tiftag(1).Width-1), ' ', num2str(tiftag(1).Height-1),...
        %                     ' -translation ', num2str(tiftag(1).Width/2), ' ', ...
        %                     num2str(tiftag(1).Height/2), ' ', num2str(tiftag(1).Width/2), ...
        %                     ' ', num2str(tiftag(1).Height/2), ' -showOutput']);
        %
        %                 temp = MIJ.getImage ('Output');
        %             img (count, :, :) = squeeze(uint16(temp(:, :, 1)));
        %                 MIJ.run('Close', 'Output');
        %                 MIJ.run('Close', ['Img' num2str(count)]);
        %
        
    end
    % Generate intensity traces for all cells
    for c = 1 : max(max(cellmask)) % 2D matrix for all cells
        if min(min(min(img(:, cellmask == c)))) >= 0
            cal(c, trial, :) = squeeze(mean(mean(img(:, cellmask == c), 3), 2));
        else
            cal(c, trial, :) = zeros(nframes, 1) * NaN;
        end
    end
    
    %             mkdir([save_direc '/' dataset '/' trial],'trial');
    
end    % Intensity Traces for all cells - generates cal(cells,trials,frames)

caltr = zeros(size(cal)); % Calcium Traces for all Cells
calbaseline = squeeze(mean(cal, 3)); % Specific baseline fluorescence value for each cell
nstd = 2; % Standard Deviation for Excitatory signal
nstdlow = 0.5; % Standard Deviation for (Hypothetical) Inhibitory signal
minoncount = 3; % min number of consecutive frames for a calcium transient
dcount = 0;

for c = 1 : max(max(cellmask))
    for t = 1 : ntrials
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
        %             calbdfsig(c, t, :) = abs(calbdf(c, t, :)) > (nstd * castd/cabase); % Significant-only changes in Fluorescence
        %
        %             % Generate signifcant traces only in caltr
        %
        %             count = 1;
        %             on = 0;
        %             oncount = 0;
        %
        %             for count = 1 : nframes
        %                 if on
        %                     tempstd = nstdlow;
        %                 else
        %                     tempstd = nstd;
        %                 end
        %                 if calbdf(c, t, count) > (tempstd * castd/cabase)
        %                     if ~on
        %                         dcount = dcount + 1;
        %                     end
        %                     on = 1;
        %                     caltr(c, t, count) = calbdf(c, t, count);
        %                     oncount = oncount + 1;
        %                 else
        %                     % Impose min triallength
        %                     if on && oncount < minoncount
        %                         caltr(c, t, count-oncount:count-1) = 0;
        %                     end
        %                     oncount = 0;
        %
        %                     on = 0;
        %
        %                 end
        %             end
        
    end
end

%     end
%     mkdir(save)
%     save([save_direc '/' dataset '/calbdf.mat'], 'calbdf');
%end
frame=1:size(calbdf,3); %or just use nframes!
time=frame./samplingrate;
disp('Remember to use "squeeze" before plotting');

% % Plotting every cell for the trial
% for j=1:ncells
%     plot(frame,squeeze(calbdf(j,1,:)))
%     hold on
% end

%HeatMap
%imagesc(squeeze(squeeze(calbdf(:,1,:))))