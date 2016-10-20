A=GoLickPercent;

for i=1:8    
    [h p] = kstest2(A(:,1),A(:,i)); %KS test with p<0.05
    if h==1
        disp(i)
        %disp(p)
    end
end

figure(1)
errorbar(meanPerformance,semPerformance,'-ok', 'LineWidth', 4);
hold on
errorbar(meanGoLickPercent,semGoLickPercent, '-xg', 'LineWidth', 4)

figure(2)
errorbar(meanPTLoops,semPTLoops,'-ob','LineWidth', 4);
hold on
errorbar(meanNGLoops,semNGLoops,'-xr', 'LineWidth', 4)

