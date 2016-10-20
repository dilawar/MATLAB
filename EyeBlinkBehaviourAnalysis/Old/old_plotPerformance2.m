figure(1);
subplot(1,3,1);
plot(perf_blinks_csPlus');
title('CS+')
xlabel('Sessions')
ylabel('Performance (%)')

subplot(1,3,2);
plot(perf_blinks_csPlusProbe');
title('Probe')
xlabel('Sessions')
ylabel('Performance (%)')

subplot(1,3,3);
plot(perf_blinks_csMinus');
title('CS-')
xlabel('Sessions')
ylabel('Performance (%)')