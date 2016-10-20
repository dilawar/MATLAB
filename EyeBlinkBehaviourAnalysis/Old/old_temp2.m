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
auc_csPlus = zeros(size(csPlusTrials,1),nbins);
auc_csPlusProbe = zeros(size(csPlusProbeTrials,1),nbins);
auc_csMinus = zeros(size(csMinusTrials,1),nbins);

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

X = auc_csPlus;
Y = auc_csPlusProbe;
Z = auc_csMinus;

X(X>blinkThreshold_csPlus) = 1;
X(X<blinkThreshold_csPlus) = 0;

Y(Y>blinkThreshold_csPlus) = 1;
Y(Y<blinkThreshold_csPlus) = 0;

Z(Z>blinkThreshold_csMinus) = 1;
Z(Z<blinkThreshold_csMinus) = 0;
