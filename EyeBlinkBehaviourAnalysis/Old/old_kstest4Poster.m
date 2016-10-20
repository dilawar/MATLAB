%kstest - 1 sample vs normal distribution

for i = 1:8
    [h1_plus,p1_plus]=kstest(a);
    [h2_minus,p2_minus]=kstest(b);
    
    h1 = [h1 h1_plus];
    p1 = [p1 p1_plus];
    
    h2 = [h2 h2_minus];
    p2 = [p2 p2_minus];
end
% kstest - 2 sample
hp = nan(8,8);
pp = nan(8,8);

hm = nan(8,8);
pm = nan(8,8);

%all sessions compared to control
for j = 1:8
    for i = 1:8
        disp([j i]);
        [h_plus,p_plus] = kstest2(a(:,j),a(:,i)); % CS+
        [h_minus,p_minus] = kstest2(b(:,j),b(:,i)); % CS-
        
        hp(j,i) = h_plus;
        pp(j,i) = p_plus;
        
        hm(j,i) = h_minus;
        pm(j,i) = p_minus;
    end
end

% figure(1);
% imagesc(hp, 'b');
% hold on
% imagesc(hm, 'r');
% title('h');

% figure(2);
% imagesc(pp, 'b');
% hold on
% plot(pm, 'r');
% title('p');