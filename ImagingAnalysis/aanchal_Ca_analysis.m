%%%%%%%%%%Read ascii calcium movie files from the Current folder

clear all
%close all


final=[];
final2=[];
CH=[];
off=0;

%trace = zeros(68, N_files-2);


%int_ser=importdata('Ca_Wed Feb 25 2015_18.18.24_0.asc');
int_ser=sifread('KClwash.sif');

int_ser=int_ser-950;

% Plot mean of the entire image over all frames
Mean_all=[];
for N=1:2400
     MM_all=int_ser(64*(N-1)+1:64*(N),2:65);
     Mean_all(N)=mean(mean(MM_all,1));
end
% plot(Mean_all);hold on;


%figure;
%plot(Mean_all);


%%%Taking mean intensity of the selected ROI

N_frames=size(int_ser,1)/64; %binning 8x8, hence 512/8=64

ROI_start_row=25;
ROI_stop_row=35;
ROI_start_col=25;
ROI_stop_col=35;

%Baseline = 500ms = 33 frames
Mean_base=[];
for N=1:33    
%    MM=int_ser(128*(N-1)+1:N*128,ROI_start_col:ROI_stop_col);
    MM=int_ser(64*(N-1)+ROI_start_row:64*(N-1)+ROI_stop_row,ROI_start_col:ROI_stop_col);
    Mean_base(N)=mean(mean(MM,1));
end

MM_prestim=Mean_base(end);

%Signal - average over 500ms post-stimulus just to average over the same
%number of frames
Mean_resp=[];
ind=1;
for N=34:34+2300
%    NN=int_ser(128*(N-1)+1:N*128,ROI_start_col:ROI_stop_col);
    NN=int_ser(64*(N-1)+ROI_start_row:64*(N-1)+ROI_stop_row,ROI_start_col:ROI_stop_col);
    Mean_resp(ind)=mean(mean(NN,1));
    ind=ind+1;
end

MM_stim=Mean_resp(3); % It is the 36th frame which is the peak

%trace(:, Num-2) = [Mean_base Mean_resp];

%%%Calculate change in baseline delta F/F

F=mean(Mean_base(:,2:33));
% F=Mean_base(end-1); % Take baseline as the frame right before the stimulus
Delta_F=MM_stim-MM_prestim;
Delta_F1=MM_stim-F;

Change=100*Delta_F/F;
Change1 = 100*Delta_F1/F;

CH=[CH Change1];

%ROI mean plot for 68 frames
AA=[Mean_base(:,2:end) Mean_resp];
%BB=AA/F;
% figure;
%plot(AA/F-1)
final=[final; AA/F];
off=off+0.05;

figure(1)
plot(AA/F + off, 'LineWidth',2);  hold on;
final2=[final2; AA];



