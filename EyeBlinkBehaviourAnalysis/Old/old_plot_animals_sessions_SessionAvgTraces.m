close all

figure(1);
for i = 1:4
    subplot(4,1,i)
    imagesc(squeeze(csPlus(i,:,:)));
    title(['MouseK' num2str(i)]);
    colormap(jet);
    caxis([-10 10])
    xlabel('Time (/10 ms)');
    ylabel('Session');
    axis([0 1000 0 5]);
    %title(['MouseK' num2str(i)]);
end

figure(2);
for i = 1:4
    subplot(4,1,i)
    imagesc(squeeze(csPlus(:,i,:)));
    title(['Session' num2str(i)]);
    colormap(jet);
    caxis([-10 10])
    xlabel('Time (/10 ms)');
    ylabel('Animal');
    axis([0 1000 0 5]);
    %title(['MouseK' num2str(i)]);
end

figure(3);
for i = 1:4
    subplot(4,1,i)
    imagesc(squeeze(csMinus(i,:,:)));
    title(['MouseK' num2str(i)]);
    colormap(jet);
    caxis([-10 10])
    xlabel('Time (/10 ms)');
    ylabel('Session');
    axis([0 1000 0 5]);
    %title(['MouseK' num2str(i)]);
end

figure(4);
for i = 1:4
    subplot(4,1,i)
    imagesc(squeeze(csMinus(:,i,:)));
    title(['Session' num2str(i)]);
    colormap(jet);
    caxis([-10 10])
    xlabel('Time (/10 ms)');
    ylabel('Animal');
    axis([0 1000 0 5]);
    %title(['MouseK' num2str(i)]);
end