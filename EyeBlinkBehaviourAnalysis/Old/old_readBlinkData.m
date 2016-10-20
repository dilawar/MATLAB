% Use this code only if there is some problem assigning trials as CS+/- with
% the original analysis code.

% AUTHOR: Kambadur Ananthamurthy

clear all
close all

saveDirec = ('/Users/ananth/Desktop/Work/Analysis/eyeBlinkBehaviourAnalysis/');
imageSaveDirec = ('/Users/ananth/Desktop/Analysis4Giulia/');
% 
% animal = 5;
sessionType = 0;
session = 1;

getData = 1;
saveData = 1;
plotData = 1;

for animal = 5:8
    %for session = 1:4
        mouseName = ['MouseK' num2str(animal)];
        dataset = [mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session)];
        
        if getData == 1
            data = load([saveDirec mouseName '/' dataset '/' dataset '_filtered.mat']);
            blinkData = data.blinkData_filtered;
            clear data
            CS = load([saveDirec mouseName '/' dataset '/' dataset '_csType.mat']);
            csType = CS.csType;
            clear CS
            blinkData(1,:) = [];
            %csType (1,:) = [];
            new_csMinus = [];
            new_csPlus = [];
            for i = 1:(length(csType)-1)
                if csType(i) == 1
                    new_csPlus = [new_csPlus; blinkData(i,:)];
                    %disp(i);
                else
                    new_csMinus = [new_csMinus; blinkData(i,:)];
                end
            end
        end
        
        %save new files
        if saveData == 1
            save([saveDirec mouseName '/' dataset '/' dataset '_new_csPlus'], 'new_csPlus');
            save([saveDirec mouseName '/' dataset '/' dataset '_new_csMinus'], 'new_csMinus');
        end
        
        if plotData == 1
            csPlus_struc = load([saveDirec mouseName '/' dataset '/' dataset '_new_csPlus']);
            csMinus_struc = load([saveDirec mouseName '/' dataset '/' dataset '_new_csMinus']);
            csPlus = csPlus_struc.new_csPlus;
            csMinus = csMinus_struc.new_csMinus;
            clear csPlus_struc
            clear csMinus_struc
            
            figure(1);
            imagesc(csPlus);
            colormap(jet);
            colorbar;
            title ([mouseName ' SessionType' num2str(sessionType) ' Session' num2str(session) ' CS+']);
            xlabel('Time (/10 ms)');
            ylabel('Trials');
            print([imageSaveDirec mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session) '_csPlus'], '-djpeg')
            
            figure(2);
            imagesc(csMinus);
            colormap(jet);
            colorbar;
            title ([mouseName ' SessionType' num2str(sessionType) ' Session' num2str(session) ' CS-']);
            xlabel('Time (/10 ms)');
            ylabel('Trials');
            print([imageSaveDirec mouseName '_SessionType' num2str(sessionType) '_Session' num2str(session) '_csMinus'], '-djpeg')
            
            % In case of multiple datasets,
            reply = input('Quit looking at datasets? ', 's'); %  Y to quit; Enter to continue
            if numel(reply)>0
                break
            end
        end
    %end
end
