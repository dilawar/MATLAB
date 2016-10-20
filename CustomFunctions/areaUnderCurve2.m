% AUTHOR - Kambadur Ananthamurthy
% Area Under the Curve of specific data

function [X, Y, Z] = areaUnderCurve2(csPlusTrials, csPlusProbeTrials, csMinusTrials, sessionType, spontaneousWindow, criticalWindow, blinkThreshold_csPlus, blinkThreshold_csMinus);

if str2double(sessionType) <=5
    %No-Puff Control
    probe = 0;
else
    %TEC
    probe = 1;
end

%Rectify
rectifiedcsPlusTrials = abs(csPlusTrials);
rectifiedcsPlusProbeTrials = abs(csPlusProbeTrials);
rectifiedcsMinusTrials = abs(csMinusTrials);

%Preallocate
csPlus = zeros(size(csPlusTrials));
csPlusProbe = zeros(size(csPlusProbeTrials));
csMinus = zeros(size(csMinusTrials));
%Set baseline to zero
for trial = 1:size(rectifiedcsPlusTrials,1)
    csPlus(trial,:) = rectifiedcsPlusTrials(trial,:)-median(rectifiedcsPlusTrials(trial,:));
end
for trial = 1:size(rectifiedcsPlusProbeTrials,1)
    csPlusProbe(trial,:) = rectifiedcsPlusProbeTrials(trial,:)-median(rectifiedcsPlusProbeTrials(trial,:));
end
for trial = 1:size(rectifiedcsMinusTrials,1)
    csMinus(trial,:) = rectifiedcsMinusTrials(trial,:)-median(rectifiedcsMinusTrials(trial,:));
end

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
prob_csPlus_tec(animal,session) = numel(X(X==1))/numel(X);
prob_csPlusProbe_tec(animal,session) = numel(Y(Y==1))/numel(Y);
prob_csMinus_tec(animal,session) = numel(Z(Z==1))/numel(Z);
end