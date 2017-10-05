function [locs] = computeGFP2(EEG)
    [gfp,~] = eeg_gfp(EEG.data',1);
    [~, locs] = findpeaks(double(gfp));
end
