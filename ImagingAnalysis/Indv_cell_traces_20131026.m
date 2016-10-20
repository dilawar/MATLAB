% Code to look at individual cell traces in plot
% Kambadur Ananthamurthy

%Remember to load "calbdf" in case you haven't

%calbdf(cells,trials,frames)

for i=1:sessiontrial
    i
    imagesc(squeeze(squeeze(caltr(:,i,:))))
%     frames=1:size(calbdf,3);
%     trialduration=10;
%     samplingrate=length(frames)/trialduration;
%     time=frames/samplingrate;
%     plot(time,a,'LineWidth',1.2);
    
    %Next Cell
    reply = input('Proceed to the next Cell? y/n ', 's');
    if isempty(reply)
        reply = 'y';
    end
end