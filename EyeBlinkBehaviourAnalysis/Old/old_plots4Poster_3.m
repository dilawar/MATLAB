clear all
close all

direc = ('/Users/ananth/Desktop/Poster2016/');
nAnimals = 3;

a = [];
b = [];

for i = 1:nAnimals
    load([direc 'csPlus' num2str(i)]);
    a = [a; csPlus];
    
    load([direc 'csMinus' num2str(i)]);
    b = [b; csMinus];
end
clear csPlus
clear csMinus

subplot(1,2,1);
plot(a(1,:)*100,'b');
hold on;
plot(a(2,:)*100,'b');
hold on;
plot(a(3,:)*100,'b');
hold on;
errorbar(median(a*100,1), std(a*100,1),'-bo',...
    'LineWidth', 3,...
    'MarkerSize', 10);
axis([1 8 -55 70]);
title('CS+ Scores', ...
    'FontSize',20,...
    'FontWeight','bold');
xlabel('Sessions',...
    'FontSize',20,...
    'FontWeight','bold');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 20)
ylabel('Performance/%',...
    'FontSize',20,...
    'FontWeight','bold');
yt = get(gca, 'YTick');
set(gca, 'FontSize', 20)

subplot(1,2,2);
plot(b(1,:)*100,'r');
hold on;
plot(b(2,:)*100,'r');
hold on;
plot(b(3,:)*100,'r');
errorbar(median(b*100,1), std(b*100,1), '-ro',...
    'LineWidth', 3,...
    'MarkerSize', 10);
axis([1 8 -55 70]);
title('CS- Scores', ...
    'FontSize',20,...
    'FontWeight','bold');
xlabel('Sessions',...
    'FontSize',20,...
    'FontWeight','bold');
xt = get(gca, 'XTick');
set(gca, 'FontSize', 20)
ylabel('Performance/%',...
    'FontSize',20,...
    'FontWeight','bold');
yt = get(gca, 'YTick');
set(gca, 'FontSize', 20)