close all

plot(squeeze(calbdf(87,12,:)), 'Color','red', 'LineWidth', 3)
hold on
plot(squeeze(calbdf(60,20,:)+0.3), 'Color','blue', 'LineWidth', 3)
hold on
plot(squeeze(calbdf(39,37,:)+0.4), 'Color','green', 'LineWidth', 3)
hold on
plot(squeeze(calbdf(15,50,:)+0.6), 'Color','magenta', 'LineWidth', 3)


imagei = double(imread(['/Users/ananthamurthy/Desktop/Work/Imaging/20131218/mouse4_Tone10knPuff-ROI/1_Trial-1_ROI-1r.tif'], 1));
for w = 2:5
    image = double(imread(['/Users/ananthamurthy/Desktop/Work/Imaging/20131218/mouse4_Tone10knPuff-ROI/1_Trial-1_ROI-1r.tif'], w));
    imagei = imagei + image;
end   % adds intensity values of frames 1:5
referenceimage = imagei./5;


cellmask_temp=cellmask;
tempi2=find(cellmask==15);
cellmask_temp(tempi2)=300;
figure(2)
imagesc(cellmask_temp);