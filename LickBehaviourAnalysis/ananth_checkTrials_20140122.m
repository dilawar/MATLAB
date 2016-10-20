F=toneTrials;
for i=1:size(F,2)
    i
    figure(1);
    imagesc(squeeze(F(:,i,:)))
    reply = input('Quit looking at datasets? ', 's'); %  Y to quit; Enter to continue
    if numel(reply)>0
        break
    else
    end
end