function [imagesToClusterReal, imagesToClusterImagined, EEG_imaginedMotor, EEG_realMotor, uncleanedEEG_imagined, uncleanedEEG_real] = preprocess_data(lowFreq, highFreq, downsampleToFreq, chansToRemove, mapClusteringMethod, mapComparisonMethod)

    filesReal = dir('../cleanerData/realMotor_both');
    filesImagined = dir('../cleanerData/imaginedMotor_both');
    
    subjectCount = 0;
    EEG_realMotor = [];
    EEG_imaginedMotor = [];
    EEG_group = {};
    imagesToCluster = {};
    imagesToClusterReal = [];
    imagesToClusterImagined = [];
    uncleanedEEG_imagined = {};
    uncleanedEEG_real = {};

    numReal = 0;
    for i = 1:length(filesReal)
        [~,name,ext] = fileparts(filesReal(i).name);
        if(strcmp(ext, '.edf') == 1)
            foo = pop_biosig(['../cleanerData/realMotor_both','/',filesReal(i).name], 'importevent','off','importannot','off');
            if(foo.nbchan == 64)
                if(~isempty(foo.data))
                    try
                        numReal = numReal + 1;
                        subjectCount = subjectCount + 1;
                        datasetName = name(1:3);
                        EEG_group{subjectCount} = foo;
                        setName = datasetName;
                        EEG_group{subjectCount} = pop_reref(EEG_group{subjectCount}, []);
                        EEG_group{subjectCount} = pop_eegfiltnew(EEG_group{subjectCount}, 59.5, 61.5, [], 1, [], 0);
                        EEG_group{subjectCount} = pop_eegfiltnew(EEG_group{subjectCount}, 1, 30, [], 0, [], 0);
                        uncleanedEEG_real{numReal} = EEG_group{subjectCount};
                        EEG_group{subjectCount}.data = autobss( EEG_group{subjectCount}.data);
                        EEG_group{subjectCount} = pop_autobssemg( EEG_group{subjectCount}, [], [], 'efica', {'eigratio', [1000000]}, 'emg_psd', {'ratio', [10],'fs', [160],'femg', [15],'estimator',spectrum.welch,'range', [0  32]});
                        EEG_group{subjectCount} = pop_resample( EEG_group{subjectCount}, downsampleToFreq); %%downsample
                        EEG_realMotor{numReal} = EEG_group{subjectCount};
                        [locsReal{numReal}] = computeGFP2(EEG_realMotor{numReal});
                        
                    catch
                        disp('-----> Skipped a data set due to some error...')
                        numReal = numReal - 1;
                        subjectCount = subjectCount - 1;
                    end
                end
            end
        end
    end
    
    numImagined = 0;
    for i = 1:length(filesImagined)
        [~,name,ext] = fileparts(filesImagined(i).name);
        if(strcmp(ext, '.edf') == 1)
            foo = pop_biosig(['../cleanerData/imaginedMotor_both','/',filesImagined(i).name]);
            if(foo.nbchan == 64)
                if(~isempty(foo.data))
                     try
                        numImagined = numImagined + 1;
                        subjectCount = subjectCount + 1;
                        datasetName = name(1:3);
                        EEG_group{subjectCount} = foo;
                        setName = datasetName;
                        EEG_group{subjectCount} = pop_reref(EEG_group{subjectCount}, []);
                        EEG_group{subjectCount} = pop_eegfiltnew(EEG_group{subjectCount}, 59, 61, [], 1, [], 0);
                        EEG_group{subjectCount} = pop_eegfiltnew(EEG_group{subjectCount}, 1, 30, [], 0, [], 0);
                        uncleanedEEG_imagined{numImagined} = EEG_group{subjectCount};
                        EEG_group{subjectCount}.data = autobss( EEG_group{subjectCount}.data);
                        EEG_group{subjectCount} = pop_autobssemg( EEG_group{subjectCount}, [], [], 'efica', {'eigratio', [1000000]}, 'emg_psd', {'ratio', [10],'fs', [160],'femg', [15],'estimator',spectrum.welch,'range', [0  32]});
                        EEG_group{subjectCount} = pop_resample( EEG_group{subjectCount}, downsampleToFreq); %%downsample
                        EEG_imaginedMotor{numImagined} = EEG_group{subjectCount};
                        [locsImagined{numImagined}] = computeGFP2(EEG_imaginedMotor{numImagined});
                     catch
                         disp('Skipped a data set due to some error')
                         numImagined = numImagined - 1;
                         subjectCount = subjectCount - 1;
                     end
                    
                end
            end
            
        end
    end
   
    
   totalLength = 0;
   for subject=1:length(EEG_realMotor)
        try 
            tempCount = 0;
            for i=1:length(locsReal{subject})
                tempCount = tempCount + 1;
                curLoc = locsReal{subject}(i);
                totalLength = totalLength + 1;
                if(size(EEG_group{subject}.data(:,curLoc),1) == 64)
                    imagesToClusterReal(:, totalLength) = EEG_realMotor{subject}.data(:,curLoc);
                    %nothing
                end
            end
        catch     
        end
   end
   
   totalLength = 0;
   for subject=1:length(EEG_imaginedMotor)
        try 
            tempCount = 0;
            for i=1:length(locsImagined{subject})
                tempCount = tempCount + 1;
                curLoc = locsImagined{subject}(i);
                totalLength = totalLength + 1;
                if(size(EEG_group{subject}.data(:,curLoc),1) == 64)
                    imagesToClusterImagined(:, totalLength) = EEG_imaginedMotor{subject}.data(:,curLoc);
                    %nothing
                end
            end
        catch     
        end
   end

end