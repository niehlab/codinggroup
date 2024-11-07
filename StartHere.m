% Pick mouse of interest (Gria4 or Stargazer):)</3
mouse = 'Gria4'; % specify the mouse model in format mouse = 'Gria'; 
%mouse = 'Stargazer'; 
addpath '/media/elaX/intanData/ela/individualExperimentDataBase/2024_09_25__12_57_00'; % path to database 
load('2024_09_25__12_57_00_allDataBases.mat') % Load the database
logicalIndexStrain = strcmp(allDataBases.Strain, mouse);% keep one mouse model only 

[originalDB,numMice,whatMice] = BigLag_240405(logicalIndexStrain,allDataBases); % load all data bases and lump structures 

homeFolder = '/media/elaX/MyFigures/SpikeFieldCoherence/September';  %  Make folder where figures will be saved 

% make the matrix for storing coherence 
% prepare for the heatmap 
[BigTableOfLag,labelsBigTable,numStructures,numMice] = makeBigTableOfLag_ed(originalDB,numMice); 
DominantFrequency = zeros(numStructures,numMice); % make the matrix filled with zeros to be populated with dominant frequency
PeakCoherence= zeros(numStructures,numMice); % make the matrix filled with zeros to be populated with peak coherence
CoherenceMatrix  = cell(numStructures,numMice); % make the struct that will hold the mean coherence values (a vector)

% Important values
Fs = 1000;
dt = 1/Fs; % 1 divided by sampling frequency dt - time step
fNQ = 1 / (2 * dt);  % Nyquist frequency - half od the sampling frequency
%
%T = size(timeForEEG,1)/1000; % duration of the seizure
%df = 1/T; % Determine the frequency resolution
% faxis = (0:df:fNQ); % Construct the frequency axis.




%% 
for iMouse = 1:numMice
    display(iMouse)
    logicalIndex = strcmp(originalDB.MouseID, whatMice{iMouse}); % keep one mouse
    MouseDatabase = originalDB(logicalIndex, :); % logical index to extract data from one mouse
    whatSeizures = unique(MouseDatabase.SeizureNumber); %how many seizures are here?
    realID = MouseDatabase.MouseID(1); % what is real ID of this mouse
    realID = realID{1,1}; % extract the number only 
    mouseFolderPath = fullfile(homeFolder, ['Mouse_' num2str(realID)]);
    if ~exist(mouseFolderPath, 'dir')
        mkdir(mouseFolderPath); % Create the folder if it doesn't exist
    end

    numSeizures = length(whatSeizures); % how many seizures are here? 
    [whatStructures] = extractStructuresForOneMouse(MouseDatabase);% what structures are here? 
    indicesToKeep = ~strcmp(whatStructures, 'DontAnalyze'); % Find the indices where the element is not 'DontAnalyze'
    filteredStructures = whatStructures(indicesToKeep); % Use logical indexing to keep only the elements that are not 'DontAnalyze'
    numStructures = length(filteredStructures); % how many structures are here? 
    [EEGTable] = extractTimeAndEEG_ed(MouseDatabase);  % extract time and voltage data for each seizure 

    [AllNeuronsAllSeizures,AllLabelsAllSeizures,newStructures] = ProcessSeizureDataForSFC(MouseDatabase,numSeizures,whatSeizures,numStructures,filteredStructures,EEGTable);
    % check the brain structures (number and names) before proceeding
    % further. Code only keeps structures that have non empty spikes 

    % check if AllNeuronsAllSeizures is empty. If YES, skip to the next
    % iMouss
    if isempty(AllNeuronsAllSeizures{1,1}) % check 1st field only 
        continue; % Skip to the next iMouse iteration
    end
    
    flattenedCells = AllLabelsAllSeizures{1,1};
    flattenedCells = [flattenedCells{:}];
    while any(cellfun(@iscell, flattenedCells))
        flattenedCells = [flattenedCells{:}];
    end
    
    % Ensure all elements are character vectors (this step is usually redundant for names)
    stringData = cellfun(@char, flattenedCells, 'UniformOutput', false);
    uniqueNames = unique(flattenedCells); % Find unique brain structures 
   
    numStructures = size(uniqueNames,2); % re-evaluate the number of brain structures before proceeding 
    [sumSpikes] = sumSpikesAndPlot(AllLabelsAllSeizures,AllNeuronsAllSeizures,numSeizures,numStructures,uniqueNames,flattenedCells);
        
    % loop through all seizures and plot spike count, spectrogram, and
    % power for all brain structures 
    %addpath('/home/Matlab/Downloaded_Functions/FMAToolbox-master/');
    addpath('/home/Matlab/scottSeizureDetection/FMAToolbox-master/');
    addpath('/home/Matlab/Downloaded_Functions/ChronuxToolbox/');
   
    [SXX,SYY,SXY,COHR,scaledHeight,f] = SEPTplotSpikesSpectraAndPower(EEGTable,numStructures,numSeizures,sumSpikes,uniqueNames,mouseFolderPath);
    
    visualizeAllSpectra(numSeizures,f,numStructures,scaledHeight,SXX,SXY,SYY,COHR,mouseFolderPath); 

    calculateMeanCoherence(COHR,numStructures,numSeizures,f,uniqueNames,mouseFolderPath)



    % store dominant frequency and coherence value 
    [DominantFrequency,CoherenceMatrix,PeakCoherence] = StoreCoherenceAndDominantFreq(COHR,iMouse,uniqueNames,DominantFrequency,numSeizures,f,labelsBigTable,CoherenceMatrix,PeakCoherence,numStructures);
    display(CoherenceMatrix); % check if looks good 
    display(DominantFrequency); % check if looks good 

   
    %clear(uniqueNames); 
end 

% linePlots(AllIterationsCoherence); 
% plotPolarPlots(AllIterationsCoherence); 
%% Visualize the coherence for all 
[allFields,heatmapLabels] = HeatmapCoherence(CoherenceMatrix, f, labelsBigTable, whatMice); % plot all structrues and mice in one figure 

% sort from the highest coherence to lowest and plot the heatmap 
[newHeatmapLabels, sortedCoherenceValues] = rankBrainStructures(allFields, heatmapLabels); % rank the brain structures by mean coherence
orderedHeatmap(heatmapLabels,newHeatmapLabels,allFields,f); % plot the heatmap 