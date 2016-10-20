% clear all
% close all

tau = 0.5;  % time-course of calcium transient in seconds
triallength = 20;  % trial length in seconds - edit for every video untill a standard protocol is established
align_ovsampl_factor = 16;

% exlist = [64 70 80 125];
exlist = [];

expcellnum = 100; % Change values for every dataset with an expected number of cells. Choose a high enough number
%
% % while 1
% direc = fgetl(fid);
%
% if ~ischar(direc)
%     break
% else
% end

%Making save directories
save_direc = '/Users/ananthamurthy/Desktop/Work/Analysis/ImageAnalysis';

%     date = findstr(direc, '201');
%     dataset_date = direc(date:date+7);
%     mouse = findstr(direc, 'mouse'); %to find the name of the mouse
%     dataset_mousename = direc(mouse:mouse+5); %edit in case the name of the mouse is longer
%     dataset_ROIindex=direc(mouse+7); %edit in case the name of the mouse is longer
dataset_date=('20131218');
dataset_mousename=('mouse4');
dataset_ROIindex=('Tone10knPuff-ROI');
dataset = [dataset_date '/' dataset_mousename '/' dataset_ROIindex];

mkdir('save_direc','dataset');
direc =['/Users/ananthamurthy/Desktop/Work/Imaging/' dataset_date '/' dataset_mousename '_' dataset_ROIindex '/'];

%get tiff tag info
num_def_trials=10; %always specify the number of LabVIEW defined trials
ntrials=1; %always specify. NOTE: ntrials= (max(zod)+(max(set)-1)*num_def_trials), unless there are additional trials
% for set=1:3
set=1;

while 1
    direc = fgetl(fid);
    
    if ~ischar(direc) 
        break
    else
    end

    %direc = 'C:\Data\data\BlinkView\20120113\field2_rand_final_final-ROI';
    %parsing direc to identify control directory
    date = findstr(direc, '201');
    dataset_date = direc(date:date+7);

    dataset_name = direc(1, date+9:length(direc));
    uscorei = findstr(dataset_name, '_');
    mouse_name = dataset_name(1:uscorei(1,1)-1);
    control_direc = direc(1, 1:(date+9 + uscorei(1,1)-1) );
    control_direc = [control_direc 'control_no-puff-ROI'];

    save_direc = 'C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\new_extraction_program';

    direci = direc;
    for control_direc_counter = 1:2
        if control_direc_counter == 1
            direc = control_direc;
        elseif control_direc_counter == 2
            direc = direci;
        end

        %assembling directories split due to acquisition program crashes
        direc = dataset_assembler(direc);

        %parsing directory to find date and dataset names
        datei = findstr(direc, '201');
        data_date = direc(1, datei:(datei+7) );
        dataset_name = direc(1, (datei+9):(length(direc)-1) );
        mkdir([save_direc '\' data_date '\' dataset_name]);
        imsave_direc = [direc(1, 1:(end - 5) ) '_aligned_dft\'];
        direc
        
        mkdir(imsave_direc);
        
        clear datei

        tif_tag = imfinfo([direc 'Trial - 0-1.tif']);
        img_desc = tif_tag.ImageDescription;


        %assigning dataset variables
        no_trials = floor(( length(dir([direc '\*.tif'])) ) );
        


        %checking if this dataset has already been analysed
        if isdir(imsave_direc) == 1

            saved_images = dir([imsave_direc '*.mat' ]);

             if length(saved_images) == no_trials
                 disp(['Dataset Alread Analysed, Skipping - ' direc]);
                 clear saved_images
                 
                 %saving control_direc so that ROI can be recovered for next dataset
                 if control_direc_counter == 1
                     control_save_direc = [save_direc '\' data_date '\' dataset_name];
                 else
                 end

                 continue
             else
                 clear saved_images
                 
             end
        else
        end
        


        framei = findstr(img_desc, 'Frames');
        no_frames = str2num(img_desc( (framei + 9):(framei + 12)));
        clear framei;

        CS_oni = findstr(img_desc, 'CS Onset');
        CS_duri = findstr(img_desc, 'CS Dur'); 
        CS_US_delayi = findstr(img_desc, 'CS-US');
        US_durationi = findstr(img_desc, 'US Dur');
        Pb_Tr_i = findstr(img_desc, 'Pb Tr');
        CS_minus_i = findstr(img_desc, 'CS- Tr');


        CS_onset_time = str2num(img_desc( (CS_oni+16):(CS_oni+19) ));                   %in ms
        CS_duration = str2num(img_desc( (CS_duri+14):(CS_duri+17) ));                   %in ms
        CS_US_delay = str2num(img_desc( (CS_US_delayi+20):(CS_US_delayi+23) ));         %in ms
        US_duration = str2num(img_desc( (US_durationi+14):(US_durationi+16) ));         %in ms
        isProbe_bool = str2num(img_desc( (Pb_Tr_i+16)) );                               %logical
        isCSminus_bool = str2num(img_desc( (CS_minus_i+17)) );                          %logical

        clear CS_oni
        clear CS_duri
        clear CS_US_delayi
        clear US_durationi
        clear PB_TR_i
        clear CS_minus_i


        %Opening 5 frames from first trial to form averaged first image, ROIs
        imagei = double(imread([direc 'Trial - 0-1.tif'], 1));
        for i = 2:5
            image = double(imread([direc 'Trial - 0-1.tif'], i));
            imagei = imagei + image;    
        end
        im_ave = imagei./5;
       
        %Manually obtaining ROIs for cells
        ROI_mat = load([save_direc '\' data_date '\' dataset_name '\ROIs.txt']);
        im_avei = im_ave;

%         cells = find (ROI_mat > 0);
%         im_avei(cells) = mean(mean(im_avei));
%         figure (1)
%         imagesc(im_avei);
%         %set(gcf, 'Position', get(0,'Screensize')); %maximises figure

        labels = bwlabel(ROI_mat);

        clear BW
        clear button
        clear cells
        clear i
        clear image
        clear imagei
        clear img_desc
        clear next
        clear tif_tag

        %re-aligning ROI matrix to im_ave, in case there has been some systematic shift
        c = xcorr2(im_ave, ROI_mat);
        [maxr, maxcol] = find(c == max(max(c)));
        col_lag = size(ROI_mat, 1) - maxr;
        row_lag = size(ROI_mat, 2) - maxcol;

        %moving ROI_mat as per lags detected
        if sign(row_lag) == -1
            ROI_mat(:, ( (end-abs(row_lag) +1 ):end)) = [];
            pad = zeros(size(ROI_mat, 1), abs(row_lag) );
            ROI_mat = [pad, ROI_mat];
        elseif sign(row_lag) == 1
            ROI_mat(:, 1:(abs(row_lag))) = [];
            pad = zeros(size(ROI_mat, 1), abs(row_lag) );
            ROI_mat = [ROI_mat, pad];
        elseif sign(row_lag) == 0
        end

        if sign(col_lag) == -1
            ROI_mat( ( (end-abs(col_lag)+1):end), :) = [];
            pad = zeros(abs(col_lag), size(ROI_mat, 2) );
            ROI_mat = [pad; ROI_mat];
        elseif sign(col_lag) == 1
            ROI_mat(1:(abs(col_lag)), :) = [];
            pad = zeros(abs(col_lag), size(ROI_mat, 2) );
            ROI_mat = [ROI_mat; pad];
        elseif sign(col_lag) == 0
        end
        
        cells = find (ROI_mat > 0);
        im_avei(cells) = mean(mean(im_avei));
        figure (1)
        imagesc(im_avei);
        clear im_avei
        %set(gcf, 'Position', get(0,'Screensize')); %maximises figure

        %close all
        pause(0.1)      %this is to prevent MIJI from causing the program to hang - the dialog box needs a gap to close before MIJI can run
        
        ref_img = im_ave;
        ref_ft = fft2(ref_img);
        
        raw_data_mat = zeros(no_frames, max(max(labels)), no_trials);
        trial_types = zeros(1, no_trials);
        status_bar = zeros(1, no_trials);
        for trial_no = 0:(no_trials-1)

            %finding and saving nature of trial
            tif_tag = imfinfo([direc 'Trial - ' int2str(trial_no) '-1.tif']);
            img_desc = tif_tag.ImageDescription;

            Pb_Tr_i = findstr(img_desc, 'Pb Tr');
            CS_minus_i = findstr(img_desc, 'CS- Tr');
            isProbe_bool = str2num(img_desc( (Pb_Tr_i+16)) );                               %logical
            %isCSminus_bool = str2num(img_desc( (CS_minus_i+17)) );                          %logical


            clear tif_tag
            clear img_desc
            clear Pb_Tr_i
            clear CS_minus_i

            if isProbe_bool == 1
                trial_types(1, (trial_no+1) ) = 1;
            else
            end

            
            figure(2)
            status_bar(1, 1:trial_no) = 1;
            imagesc(status_bar)
            title(direc)
            drawnow

            for frame_no = 1:no_frames
                curr_frame = imread([direc 'Trial - ' int2str(trial_no) '-1.tif'], frame_no);
                im_ft = fft2(curr_frame);

                [output] = dftregistration(ref_ft, im_ft, align_ovsampl_factor);
                
                figure(1)
                imagesc(im_ave)
                
                figure(2)
                imagesc(curr_frame)
                
                output
                
                if output(1, 1) > .5
                    keyboard
                else
                end
                
                img(frame_number, :, :) = curr_frame;
            

                %reading out mean intensity of each cell
                for cell_number = 1:max(max(labels))
                    celli = find(labels == cell_number);
                    cell_pix = img(frame_nnumber, celli);
                    ave_int = mean(mean(cell_pix));
                    raw_data_mat(frame_number, cell_number, (trial_nnumber+1)) = ave_int;
                end

                    clear cell_no
                    clear celli
                    clear cell_pix
                    clear ave_int


            end
            
            %writing aligned images as .mat files
            
            save([imsave_direc '\trial-' int2str(trial_no) '.mat' ], 'img');
            clear img
            


        end
        
       
      
        close all

        info = [no_frames; no_trials; max(max(labels)); CS_onset_time; CS_duration; CS_US_delay; US_duration; isProbe_bool; isCSminus_bool];

        %saving data to .mat files
        save([save_direc '\' data_date '\' dataset_name '\' 'raw_int_data.mat'], 'raw_data_mat')
        save([save_direc '\' data_date '\' dataset_name '\' 'info.mat'], 'info')
        save([save_direc '\' data_date '\' dataset_name '\' 'trial_type_list.txt'], 'trial_types', '-ASCII')

    end
    
    
     clear CS_US_delay
    clear CS_duration
    clear CS_onset_time
    clear ROI_mat
    clear US_duration
    clear a
    clear ans
    clear c
    clear cells
    clear col_lag
    clear control_direc
    clear curr_frame
    clear data_date
    clear dataset_date
    clear dataset_name
    clear date
    clear frame_no
    clear im_ave
    clear imsave_direc
    clear info
    clear isCSminus_bool
    clear isProbe_bool
    clear labels
    clear maxcol
    clear maxr
    clear mouse_name
    clear no_frames
    clear no_trials
    clear pad
    clear raw_data_mat
    clear row_lag
    clear save_direc
    clear status_bar
    clear temo
    clear test
    clear trial_no
    clear trial_types
    clear uscorei
   
   
end
fclose(fid);