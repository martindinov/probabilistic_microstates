%%Preprocess data
[imagesToCluster, imagesToClusterReal, imagesToClusterImagined, kmGroupCentroids, kmMicrostatesGroup, EEG_imaginedMotor, EEG_realMotor, gfpReal, gfpImagined] = preprocess_data(8, 1, 30, 60, 5, 'sqeuclidean', 'euclidean');

%%read channel locations, used for plotting the maps with topoplot
chanlocs = readlocs('cleanerData/eloc64.txt', 'filetype', 'loc');

opts = statset('Display','iter', 'UseParallel', true);
SSEs_real = zeros(1,10);
meanWCDs_real = zeros(1,10);
for numClusters=1 
    [groupClusterPointsReal, groupCentroidsReal, SUMD, D] = kmeans(zscore(imagesToClusterReal)', numClusters, 'MaxIter', 1000, 'Replicates', 1, 'Distance', 'sqeuclidean' ,'OnlinePhase', 'on','Options', opts);
    SSEs_real(numClusters) = sum(SUMD); 
end

opts = statset('Display','iter', 'UseParallel', true);
%SSEs_imagined = zeros(1,8);
%meanWCDs_imagined = zeros(1,8);
for numClusters=9:12
    [groupClusterPointsImagined, groupCentroidsImagined, SUMD, D] = kmeans(zscore(imagesToClusterImagined)', numClusters, 'MaxIter', 1000, 'Replicates', 1, 'Distance', 'sqeuclidean' ,'OnlinePhase', 'on','Options', opts);
    SSEs_imagined(numClusters) = sum(SUMD);
end

%%%We derive the number of clusters k from looking at the results of the
%%%SSEs_real and SSEs_imagined, from above:
figure(6);plot(SSEs_real);
figure(6);hold on;plot(SSEs_imagined, 'r')

kmTemplatesReal = zeros(9,64);
for i=1:length(imagesToClusterReal)
kmTemplatesReal(kmClusteredImagesReal(i), :) = kmTemplatesReal(kmClusteredImagesReal(i), :) + imagesToClusterReal(:,i)';
end
[counts,~]=hist(kmClusteredImagesReal,unique(kmClusteredImagesReal));
kmTemplatesReal = kmTemplatesReal./counts';

kmTemplatesImagined = zeros(9,64);
for i=1:length(imagesToClusterImagined)
    kmTemplatesImagined(kmClusteredImagesImagined(i), :) = kmTemplatesImagined(kmClusteredImagesImagined(i), :) + imagesToClusterImagined(:,i)';
end
[counts,~]=hist(kmClusteredImagesImagined,unique(kmClusteredImagesImagined));
kmTemplatesImagined = kmTemplatesImagined./counts';


[groupClusterPointsReal, groupCentroidsReal, SUMD, D] = kmeans(imagesToClusterReal', , 'MaxIter', 1000, 'Replicates', 10, 'Distance', 'sqeuclidean' ,'OnlinePhase', 'off','Options', opts);
[groupClusterPointsImagined, groupCentroidsImagined, SUMD, D] = kmeans(imagesToClusterImagined', 8, 'MaxIter', 1000, 'Replicates', 10, 'Distance', 'sqeuclidean' ,'OnlinePhase', 'off','Options', opts);

%%%Silhoutte plot:

figure(7);silhouette(imagesToClusterReal(foo(1:20000))', fcmClusteredImagesReal(foo(1:20000)));

%%%The below takes a while to run for the whole data set
options = [1.1 1000 NaN 1];
fcmReplicates = 10;

max_objective_found = Inf;
bestReplicateImagined = 0;
fcmCentersImagined_all = {};
fcmUImagined_all = {};
fcmObjImagined_all = {};
for replicate=1:fcmReplicates
    [fcmCentersImagined_all{replicate}, fcmUImagined_all{replicate}, fcmObjImagined_all{replicate}] = fcm(imagesToClusterImagined', 5, options);
    curFcmObj = fcmObjImagined_all{replicate};
    if(curFcmObj(end) < max_objective_found)
        max_objective_found = curFcmObj(end);
        bestReplicateImagined = replicate;
    end
end

options = [1.1 1000 NaN 1];
fcmReplicates = 10;

max_objective_found = Inf;
bestReplicateReal = 0;
fcmCentersReal_all = {};
fcmUReal_all = {};
fcmObjReal_all = {};
for replicate=1:fcmReplicates
    [fcmCentersReal_all{replicate}, fcmUReal_all{replicate}, fcmObjReal_all{replicate}] = fcm(imagesToClusterReal', 5, options);
    curFcmObj = fcmObjReal_all{replicate};
    if(curFcmObj(end) < max_objective_found)
        max_objective_found = curFcmObj(end);
        bestReplicateReal = replicate;
    end
end

fcmCentersImagined = fcmCentersImagined_all{bestReplicateImagined};
fcmCentersReal = fcmCentersReal_all{bestReplicateReal};

fcmClusteredImagesImaginedFull = fcmUImagined_all{bestReplicateImagined};
fcmClusteredImagesRealFull = fcmUReal_all{bestReplicateReal};


[fcmClusteredImagesReal_ml,fcmClusteredImagesReal] = max(fcmClusteredImagesRealFull);
[fcmClusteredImagesImagined_ml,fcmClusteredImagesImagined] = max(fcmClusteredImagesImaginedFull);


%%%Shuffle labels so that FCM-imagined, FCM-real, KM-imagined, KM-real all
%%%have the same label for the same microstate

[matrix,diag,row,col] = sortcorrel(corr(fcmCentersImagined(:,:)',kmTemplatesImagined(:,:)'))
mapping = [row;col];

relabeledKmTemplatesImagined = zeros(4,64);
for i=1:4
    relabeledKmTemplatesImagined(mapping(1,i),:) = kmTemplatesImagined(mapping(2,i),:);
end


%%%map FCM state labels to KM state labels, to be directly comparable:

[matrix,diag,row,col] = sortcorrel(corr(fcmCentersImagined(:,:)',kmTemplatesImagined(:,:)'))
mapping = [row;col];

relabeledKmTemplatesImagined = zeros(4,64);
for i=1:4
    relabeledKmTemplatesImagined(mapping(1,i),:) = kmTemplatesImagined(mapping(2,i),:);
end


[matrix,diag,row,col] = sortcorrel(corr(fcmCentersImagined(:,:)',kmTemplatesReal(:,:)'))
mapping = [row;col];

relabeledKmTemplatesReal = zeros(4,64);
for i=1:4
    relabeledKmTemplatesReal(mapping(1,i),:) = kmTemplatesReal(mapping(2,i),:);
end


[matrix,diag,row,col] = sortcorrel(corr(fcmCentersImagined(:,:)',fcmCentersReal(:,:)'))
mapping = [row;col];

relabeledFcmCentersReal = zeros(4,64);
for i=1:4
    relabeledFcmCentersReal(mapping(1,i),:) = fcmCentersReal(mapping(2,i),:);
end




%%%%%%%%%%%%% re-label the clustered images (to use for MLP training) from the FCM and KM %%%%%%%%%%%%%%

[matrix,diag,row,col] = sortcorrel(corr(fcmCentersImagined(:,:)',kmTemplatesImagined(:,:)'));
mapping = [row;col];

kmClusteredImagesImaginedRelabeled = zeros(1,length(kmClusteredImagesImagined));
for i=1:length(kmClusteredImagesImaginedRelabeled)
    kmClusteredImagesImaginedRelabeled(i) = mapping(1,find(mapping(2,:)==kmClusteredImagesImagined(i)));
end

%%%--- works to here ---

[matrix,diag,row,col] = sortcorrel(corr(fcmCentersImagined(:,:)',kmTemplatesReal(:,:)'))
mapping = [row;col];

kmClusteredImagesRealRelabeled = zeros(1,length(kmClusteredImagesReal));
for i=1:length(kmClusteredImagesRealRelabeled)
    kmClusteredImagesRealRelabeled(i) = mapping(1,find(mapping(2,:)==kmClusteredImagesReal(i)));
end

%%done with KM, now FCM-real and done

[matrix,diag,row,col] = sortcorrel(corr(fcmCentersImagined(:,:)',fcmCentersReal(:,:)'))
mapping = [row;col];

fcmClusteredImagesRealFullRelabeled = zeros(size(fcmClusteredImagesRealFull,1),size(fcmClusteredImagesRealFull,2));
for i=1:length(fcmClusteredImagesRealFullRelabeled)
    %%%for each vector in fcmClusteredImagesRealFull, we need to flip the
    %%%ordering...
    for j=1:9
        fcmClusteredImagesRealFullRelabeled(mapping(1,find(j==mapping(2,:))),i) = fcmClusteredImagesRealFull(mapping(1,find(j==mapping(1,:))),i);
    end
end

[a,fcmClusteredImagesRealRelabeled] = max(fcmClusteredImagesRealFullRelabeled);

km_net_imagined = patternnet(15,'trainscg','crossentropy');
km_net_imagined.trainParam.epochs = 2000;
km_net_real = patternnet(15,'trainscg','crossentropy');
km_net_real.trainParam.epochs = 2000;

fcm_net_imagined = patternnet(9, 'trainscg', 'crossentropy');
fcm_net_imagined.trainParam.epochs = 2000;
fcm_net_real = patternnet(9, 'trainscg', 'crossentropy');
fcm_net_real.trainParam.epochs = 2000;



fcm_net_imagined_fullProbabilistic = fitnet(15);
fcm_net_imagined_fullProbabilistic.trainParam.epochs = 2000;
fcm_net_imagined_fullProbabilistic.trainFcn = 'trainscg';
fcm_net_real_fullProbabilistic = fitnet(15);
fcm_net_real_fullProbabilistic.trainParam.epochs = 2000;
fcm_net_real_fullProbabilistic.trainFcn = 'trainscg';



kmClusteredImagesImaginedFullRelabeled = full(ind2vec(kmClusteredImagesImagined));
kmClusteredImagesRealFullRelabeled = full(ind2vec(kmClusteredImagesReal));

fcmClusteredImagesRealFull = fcmClusteredImagesRealFullRelabeled;

%%%%% reshape the full FCM images into equivalent binary examples %%%%%

fcmClusteredImagesImaginedFull_temp = [fcmClusteredImagesImaginedFull(1,:) fcmClusteredImagesImaginedFull(2,:) fcmClusteredImagesImaginedFull(3,:) fcmClusteredImagesImaginedFull(4,:) fcmClusteredImagesImaginedFull(5,:) fcmClusteredImagesImaginedFull(6,:) fcmClusteredImagesImaginedFull(7,:) fcmClusteredImagesImaginedFull(8,:) fcmClusteredImagesImaginedFull(9,:)];



%%%% experimental - not used in the paper %%%%
layers = [fullyConnectedLayer(15) regressionLayer];

km_net_real_future = patternnet(15,'trainscg','crossentropy');
km_net_real_future.trainParam.epochs = 2000;
km_net_real_future = train(km_net_real_future, imagesToClusterReal(:,1:end-10), kmClusteredImagesRealFull(:,11:end));
%kmClusteredImagesRealFullRelabeled


fcm_net_imagined_fullProbabilistic = train(fcm_net_imagined_fullProbabilistic, imagesToClusterImagined, fcmUImagined);
fcm_net_real_fullProbabilistic = train(fcm_net_real_fullProbabilistic, imagesToClusterReal, fcmUReal);

patternnet_km_imagined = train(km_net_imagined, imagesToClusterImagined, kmClusteredImagesImaginedFull);
patternnet_km_real = train(km_net_real, imagesToClusterReal, kmClusteredImagesRealFull);

patternnet_fcm_imagined = train(fcm_net_imagined, imagesToClusterImagined, fcmClusteredImagesImaginedFull);
patternnet_fcm_real = train(fcm_net_real, imagesToClusterReal, fcmClusteredImagesRealFullRelabeled);


%patternnet_gmm = train(gmm_net, imagesToCluster', gmmClusteredImagesFull);
%patternnet_fcm = train(fcm_net, imagesToCluster', fcmClusteredImagesFull);

%%%feed full real and imagined motor imagery data into the MLPs
realMotor_concat = {};
for i=1:length(EEG_realMotor)
    realMotor_concat{i} = EEG_realMotor{i}.data;
end

imaginedMotor_concat = {};
for i=1:length(EEG_imaginedMotor)
    imaginedMotor_concat{i} = EEG_imaginedMotor{i}.data;
end

imaginedMotor_concat = [imaginedMotor_concat{:}];
realMotor_concat = [realMotor_concat{:}];



%%%%%%%%%
imaginedMotor_microstates_km = sim(patternnet_km_imagined, imaginedMotor_concat);
realMotor_microstates_km = sim(patternnet_km_real, realMotor_concat);

[imaginedMotor_microstates_km_probabilities,imaginedMotor_microstates_km_max] = max(imaginedMotor_microstates_km);
[realMotor_microstates_km_probabilities,realMotor_microstates_km_max] = max(realMotor_microstates_km);



imaginedMotor_microstates_fcm = sim(patternnet_fcm_imagined, imaginedMotor_concat);
realMotor_microstates_fcm = sim(patternnet_fcm_real, realMotor_concat);

[imaginedMotor_microstates_fcm_probabilities,imaginedMotor_microstates_fcm_max] = max(imaginedMotor_microstates_fcm);
[realMotor_microstates_fcm_probabilities,realMotor_microstates_fcm_max] = max(realMotor_microstates_fcm);



imaginedMotor_microstates_fcm_fullProbabilistic = sim(fcm_net_imagined_fullProbabilistic, imaginedMotor_concat);
realMotor_microstates_fcm_fullProbabilistic = sim(fcm_net_real_fullProbabilistic, realMotor_concat);

[imaginedMotor_microstates_fcm_probabilities_fullProbabilistic,imaginedMotor_microstates_fcm_max_fullProbabilistic] = max(imaginedMotor_microstates_fcm_fullProbabilistic);
[realMotor_microstates_fcm_probabilities_fullProbabilistic,realMotor_microstates_fcm_max_fullProbabilistic] = max(realMotor_microstates_fcm_fullProbabilistic);

%%%%%%%%


%imaginedMotor_microstates_som = sim(patternnet_som_imagined, imaginedMotor_concat');
%realMotor_microstates_som = sim(patternnet_som_real, realMotor_concat');

%[~,imaginedMotor_microstates_som_max] = max(imaginedMotor_microstates_som);
%[~,realMotor_microstates_som_max] = max(realMotor_microstates_som);



%%%subject-level microstates
subject_microstates_km_imagined = {};
subject_microstates_km_real = {};
subject_microstates_fcm_imagined = {};
subject_microstates_fcm_real = {};

subject_microstates_fcm_real_probabilities = {};
subject_microstates_km_real_probabilities = {};
subject_microstates_fcm_imagined_probabilities = {};
subject_microstates_km_imagined_probabilities = {};

subject_microstates_fcm_real_max = {};
subject_microstates_km_real_max = {};
subject_microstates_fcm_imagined_max = {};
subject_microstates_km_imagined_max = {};

for i=1:length(EEG_imaginedMotor)
    subject_microstates_km_imagined{i} = sim(patternnet_km_imagined, EEG_imaginedMotor{i}.data);
    [subject_microstates_km_imagined_probabilities{i}, subject_microstates_km_imagined_max{i}] = max(subject_microstates_km_imagined{i});
    [~,~,~,subject_microstates_km_imagined_len{i}] = findseq(subject_microstates_km_imagined_max{i});
    %subject_microstates_km_imagined_len{i} = subject_microstates_km_imagined_max{
    subject_microstates_fcm_imagined{i} = sim(patternnet_fcm_imagined, EEG_imaginedMotor{i}.data);
    [subject_microstates_fcm_imagined_probabilities{i}, subject_microstates_fcm_imagined_max{i}] = max(subject_microstates_fcm_imagined{i});
    [~,~,~,subject_microstates_fcm_imagined_len{i}] = findseq(subject_microstates_fcm_imagined_max{i});
end

for i=1:length(EEG_realMotor)
    subject_microstates_km_real{i} = sim(patternnet_km_real, EEG_realMotor{i}.data);
    [subject_microstates_km_real_probabilities{i}, subject_microstates_km_real_max{i}] = max(subject_microstates_km_real{i});
    [~,~,~,subject_microstates_km_real_len{i}] = findseq(subject_microstates_km_real_max{i});
    subject_microstates_fcm_real{i} = sim(patternnet_fcm_real, EEG_realMotor{i}.data);
    [subject_microstates_fcm_real_probabilities{i}, subject_microstates_fcm_real_max{i}] = max(subject_microstates_fcm_real{i});
    [~,~,~,subject_microstates_fcm_real_len{i}] = findseq(subject_microstates_fcm_real_max{i});
end


%%%%%%%%%%%%%%


imaginedMotor_microstates_fcm = sim(patternnet_fcm_imagined, imaginedMotor_concat');
[imaginedMotor_microstates_fcm_probabilities,imaginedMotor_microstates_fcm_max] = max(imaginedMotor_microstates_fcm);

realMotor_microstates_fcm = sim(patternnet_fcm_real, realMotor_concat');
[realMotor_microstates_fcm_probabilities,realMotor_microstates_fcm_max] = max(realMotor_microstates_fcm);





%[~,realMotor_microstates_gmm_max] = max(realMotor_microstates_gmm);


%[~,realMotor_microstates_fcm_max] = max(realMotor_microstates_fcm);


%%%Hurst exponent computation (not used in paper):

subject_hurst_realMotor_fcm_max = cellfun(@hurst_exponent, subject_microstates_fcm_real_max);
subject_hurst_realMotor_km_max = cellfun(@hurst_exponent, subject_microstates_km_real_max);
subject_hurst_imaginedMotor_fcm_max = cellfun(@hurst_exponent, subject_microstates_fcm_imagined_max);
subject_hurst_imaginedMotor_km_max = cellfun(@hurst_exponent, subject_microstates_km_imagined_max);

subject_hurst_realMotor_fcm_probabilities = cellfun(@hurst_exponent, subject_microstates_fcm_real_probabilities);
subject_hurst_realMotor_km_probabilities = cellfun(@hurst_exponent, subject_microstates_km_real_probabilities);
subject_hurst_imaginedMotor_fcm_probabilities = cellfun(@hurst_exponent, subject_microstates_fcm_imagined_probabilities);
subject_hurst_imaginedMotor_km_probabilities = cellfun(@hurst_exponent, subject_microstates_km_imagined_probabilities);

meanGfpReal = cellfun(@mean,gfpReal);
meanGfpImagined = cellfun(@mean,gfpImagined);


%%%Compute mean, median and variance of subject-level mean GFP power
meanProbabilities_fcm_imaginedMotor = cellfun(@mean,subject_microstates_fcm_imagined_probabilities)
meanProbabilities_fcm_realMotor = cellfun(@mean,subject_microstates_fcm_real_probabilities)
meanProbabilities_km_imaginedMotor = cellfun(@mean,subject_microstates_km_imagined_probabilities)
meanProbabilities_km_realMotor = cellfun(@mean,subject_microstates_km_real_probabilities)

medianProbabilities_fcm_imaginedMotor = cellfun(@median,subject_microstates_fcm_imagined_probabilities)
medianProbabilities_fcm_realMotor = cellfun(@median,subject_microstates_fcm_real_probabilities)
medianProbabilities_km_imaginedMotor = cellfun(@median,subject_microstates_km_imagined_probabilities)
medianProbabilities_km_realMotor = cellfun(@median,subject_microstates_km_real_probabilities)


varProbabilities_fcm_imaginedMotor = cellfun(@var,subject_microstates_fcm_imagined_probabilities)
varProbabilities_fcm_realMotor = cellfun(@var,subject_microstates_fcm_real_probabilities)
varProbabilities_km_imaginedMotor = cellfun(@var,subject_microstates_km_imagined_probabilities)
varProbabilities_km_realMotor = cellfun(@var,subject_microstates_km_real_probabilities)


%group_hurst_realMotor_fcm = hurst_exponent(realMotor_microstates_fcm_max);
%group_hurst_realMotor_km = hurst_exponent(realMotor_microstates_km_max);
%group_hurst_imaginedMotor_fcm = hurst_exponent(imaginedMotor_microstates_fcm_max);
%group_hurst_imaginedMotor_km = hurst_exponent(imaginedMotor_microstates_km_max);

transitionProbabilities_realMotor_km = zeros(4,4);
transitionProbabilities_realMotor_fcm = zeros(4,4);
transitionProbabilities_imaginedMotor_km = zeros(4,4);
transitionProbabilities_imaginedMotor_fcm = zeros(4,4);

for i=1:length(realMotor_microstates_km_max)-1
    state = realMotor_microstates_km_max(i);
    nextState = realMotor_microstates_km_max(i+1);
    transitionProbabilities_realMotor_km(state,nextState) = transitionProbabilities_realMotor_km(state,nextState) + 1;
end

for i=1:length(realMotor_microstates_fcm_max)-1
    state = realMotor_microstates_fcm_max(i);
    nextState = realMotor_microstates_fcm_max(i+1);
    transitionProbabilities_realMotor_fcm(state,nextState) = transitionProbabilities_realMotor_fcm(state,nextState) + 1;
end

for i=1:length(imaginedMotor_microstates_km_max)-1
    state = imaginedMotor_microstates_km_max(i);
    nextState = imaginedMotor_microstates_km_max(i+1);
    transitionProbabilities_imaginedMotor_km(state,nextState) = transitionProbabilities_imaginedMotor_km(state,nextState) + 1;
end

for i=1:length(imaginedMotor_microstates_fcm_max)-1
    state = imaginedMotor_microstates_fcm_max(i);
    nextState = imaginedMotor_microstates_fcm_max(i+1);
    transitionProbabilities_imaginedMotor_fcm(state,nextState) = transitionProbabilities_imaginedMotor_fcm(state,nextState) + 1;
end

%%%now normalize the probabilities:
transitionProbabilities_realMotor_km = transitionProbabilities_realMotor_km./sum(transitionProbabilities_realMotor_km,1);
transitionProbabilities_realMotor_fcm = transitionProbabilities_realMotor_fcm./sum(transitionProbabilities_realMotor_fcm,1);
transitionProbabilities_imaginedMotor_km = transitionProbabilities_imaginedMotor_km./sum(transitionProbabilities_imaginedMotor_km,1);
transitionProbabilities_imaginedMotor_fcm = transitionProbabilities_imaginedMotor_fcm./sum(transitionProbabilities_imaginedMotor_fcm,1);

[matrix,diag,row,col] = sortcorrel(transitionProbabilities_realMotor_km)
fig6 = figure(6);heatmap(string(col), string(row), matrix, 'ColorMap', redbluecmap)

[matrix,diag,row,col] = sortcorrel(transitionProbabilities_imaginedMotor_km)
fig7 = figure(7);heatmap(string(col), string(row), matrix, 'ColorMap', redbluecmap)

[matrix,diag,row,col] = sortcorrel(transitionProbabilities_realMotor_fcm)
fig8 = figure(8);heatmap(string(col), string(row), matrix, 'ColorMap', redbluecmap)

[matrix,diag,row,col] = sortcorrel(transitionProbabilities_imaginedMotor_fcm)
fig9 = figure(9);heatmap(string(col), string(row), matrix, 'ColorMap', redbluecmap)