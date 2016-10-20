% AUTHOR - Kambadur Ananthamurthy
% Area Under the Curve of specific data

function [X, Y, Z] = crProbability(csPlusTrials, csPlusProbeTrials, csMinusTrials, sessionType, spontaneousWindow, criticalWindow, blinkThreshold_csPlus, blinkThreshold_csMinus);

if str2double(sessionType) <=5
    %No-Puff Control
    probe = 0;
else
    %TEC
    probe = 1;
end

%Preallocate
csPlus = zeros(size(csPlusTrials));
csPlusProbe = zeros(size(csPlusProbeTrials));
csMinus = zeros(size(csMinusTrials));
%Set baseline to zero
for trial = 1:size(csPlusTrials,1)
    csPlus(trial,:) = csPlusTrials(trial,:)-median(csPlusTrials(trial,:));
end
for trial = 1:size(csPlusProbeTrials,1)
    csPlusProbe(trial,:) = csPlusProbeTrials(trial,:)-median(csPlusProbeTrials(trial,:));
end
for trial = 1:size(csMinusTrials,1)
    csMinus(trial,:) = csMinusTrials(trial,:)-median(csMinusTrials(trial,:));
end

%Rectify
csPlus = abs(csPlus);
csPlusProbe = abs(csPlusProbe);
csMinus = abs(csMinus);

%Set all negative values to zero
csPlus(csPlus<0) = 0;
csPlusProbe(csPlusProbe<0) = 0;
csMinus(csMinus<0) = 0;

nbins = floor(length(spontaneousWindow)/length(criticalWindow));

%Preallocate
auc_spont_csPlus = zeros(size(csPlusTrials,1),nbins);
auc_spont_csPlusProbe = zeros(size(csPlusProbeTrials,1),nbins);
auc_spont_csMinus = zeros(size(csMinusTrials,1),nbins);

%Reshape
for trial = 1:size(csPlus)
    A = reshape(csPlus(trial,spontaneousWindow),[length(criticalWindow),nbins]);
    auc_spont_csPlus(trial,:) = trapz(A);
end
for trial = 1:size(csPlusProbe)
    B = reshape(csPlusProbe(trial,spontaneousWindow),[length(criticalWindow),nbins]);
    auc_spont_csPlusProbe(trial,:) = trapz(B);
end
for trial = 1:size(csMinus)
    C = reshape(csMinus(trial,spontaneousWindow),[length(criticalWindow),nbins]);
    auc_spont_csMinus(trial,:) = trapz(C);
end

for trial = 1:size(auc_spont_csPlus,1)
    for bin = 1:nbins
        if auc_spont_csPlus(trial,bin) > blinkThreshold_csPlus
            X(trial,bin) = 1;
        else
            X(trial,bin) = 0;
        end
    end
end

if probe == 1
    for trial = 1:size(auc_spont_csPlusProbe,1)
        for bin = 1:nbins
            if auc_spont_csPlusProbe(trial,bin) > blinkThreshold_csPlus
                Y(trial,bin) = 1;
            else
                Y(trial,bin) = 0;
            end
        end
    end
else
    Y = [];
end

for trial = 1:size(auc_spont_csMinus,1)
    for bin = 1:nbins
        if auc_spont_csMinus(trial,bin) > blinkThreshold_csMinus
            Z(trial,bin) = 1;
        else
            Z(trial,bin) = 0;
        end
    end
end
end