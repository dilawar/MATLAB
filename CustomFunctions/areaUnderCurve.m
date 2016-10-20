% AUTHOR - Kambadur Ananthamurthy
% Area Under the Curve

function [auc_csPlus, auc_csPlusProbe, auc_csMinus] = areaUnderCurve(csPlusTrials, csPlusProbeTrials, csMinusTrials, baselineCorrection, criticalWindow)

%Preallocate
csPlus = zeros(size(csPlusTrials));
csPlusProbe = zeros(size(csPlusProbeTrials));
csMinus = zeros(size(csMinusTrials));
auc_csPlus = zeros(size(csPlusTrials,1),1);
auc_csPlusProbe = zeros(size(csPlusProbeTrials,1),1);
auc_csMinus = zeros(size(csMinusTrials,1),1);

%Set baseline to zero
if baselineCorrection == 0
    for trial = 1:size(csPlusTrials,1)
        csPlus(trial,:) = csPlusTrials(trial,:)-median(csPlusTrials(trial,:));
    end
    for trial = 1:size(csPlusProbeTrials,1)
        csPlusProbe(trial,:) = csPlusProbeTrials(trial,:)-median(csPlusProbeTrials(trial,:));
    end
    for trial = 1:size(csMinusTrials,1)
        csMinus(trial,:) = csMinusTrials(trial,:)-median(csMinusTrials(trial,:));
    end
else    
end

%Rectify
csPlus = abs(csPlus);
csPlusProbe = abs(csPlusProbe);
csMinus = abs(csMinus);

%Area under the curve for the critical window

for trial = 1:size(csPlus)
    auc_csPlus(trial) = trapz(csPlus(trial,criticalWindow));
end
for trial = 1:size(csPlusProbe)
    auc_csPlusProbe(trial) = trapz(csPlusProbe(trial,criticalWindow));
end
for trial = 1:size(csMinus)
    auc_csMinus(trial) = trapz(csMinus(trial,criticalWindow));
end

end