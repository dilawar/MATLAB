for i = 1:3
    figure(1);
    %subplot(2,3,i);
    subplot(2,3,i+3);
%     plot(c_probe(i,:),'g');
%     hold on;
    plot(a_plus(i,:),'b');
    hold on;
    plot(b_minus(i,:),'r');
    
    title(['MouseK' num2str(i+8) ' Session7'],...
        'FontSize',20,...
        'FontWeight','bold');
    axis([300 700 -30 30]);
    xlabel('Time/10 ms',...
        'FontSize',20,...
        'FontWeight','bold');
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 20)
    ylabel('Blink (fold change)',...
        'FontSize',20,...
        'FontWeight','bold');
    yt = get(gca, 'YTick');
    set(gca, 'FontSize', 20);
end