function analyseMyERPs(participants, conditions, baseWin, analysisWin, ...
  analysisType, analysisValue, filetypeIn, filenameOut, varargin)
% analyseMyERPs(participants, conditions, baselineWin, analysisWin,
%   analysisType, analysisValue, filetypeIn, filenameOut, varargin)
%
% Perform standard ERP measures on single trials, erp averaged data,
%   or grandaveraged data. Measures include amplitude measures and
%   measures for component latencies.
%
% Inputs:
% participants: 1:5 or [1 4 7]
% conditions: 1:4 or [1 2 5 6]
% baseWin: [-0.2 0] 200 ms pre-stimulus baseline
% analysisWin: [0.1 0.2] 100 ms interval between 100 and 200 ms
%
% analysisType:
% 'meanAmp'
% 'maxPeakAmp'
% 'minPeakAmp'
% 'maxPeakLat'
% 'minPeakLat'
% 'meanMaxPeakAmp'
% 'meanMinPeakAmp'
% 'critOnsetSmaller'
% 'critOnsetLarger'
% 'perOnsetSmaller'
% 'perOnsetLarger'
% 'minFracAreaLat'
% 'maxFracAreaLat'
% 'minFracPeakLat'
% 'maxFracPeakLat'
% 'maxTime2peak'
% 'minTime2peak'
% 'peak2peakAmp'
%
% analysisValue: [] or value (analysis type dependent)
%
% filetypeIn: 'avg', 'all', 'pre', 'GA_avg', 'GA_lrp'
% filenameOut: name for the analysis results file (no file extension)
%
% OPTIONAL inputs [varargin] entered as key-value pairs
% 'boundary', ...
% 'indElectBound', ...
% 'numConseqCrit', ...
% 'baseType', ...  ('abs' = subtraction vs. 'per' = percentage change)
% 'outputDir', ...
% 'trialinfoHeader', ...
%
% Outputs:
% .mat file with rows (participants & conditions) & columns (electrodes)
% .txt file with rows (participants & conditions) & columns (electrodes)
%
% VP    Cond    Fp1     AF7 ... ...
% 1     1
% 1     2
% 2     1
% 2     2
% ...
%
% or with *pre.mat or *all.mat files
%
% VP    Cond    Trial   Fp1     AF7 ...
% 1     1       1
% 1     1       2
% 1     1       3
% 1     2       1
% 1     2       2
% 1     2       3
% 2     1       1
% ...
%
% Examples:
% analyseMyERPs(1:10, 1:2, [-0.2 0], [0.1 0.2], 'meanAmp', [], 'avg', 'filename')
% analyseMyERPs(1:10, 1:2, [-0.2 0], [0.1 0.2; 0.2 0.3], 'meanAmp', [], 'avg', 'filename')
% analyseMyERPs(1:10, 1:2, [-0.2 0], [0.1 0.2; 0.2 0.3], 'meanAmp', [], 'pre', 'filename')
% analyseMyERPs(1:40, 1:4, [-0.1 0], [0.3 0.5], 'minFracAreaLat', 50, 'avg', 'fracArea', 'outputDir', curDir, 'boundary', 1.5)

%%
%%%%%%%%%%%%%%%%%%%%%%%%% SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputProblem = checkMyNumInputArgs(mfilename, nargin, 8:2:14);
if inputProblem
  return
end

funcCall = getMyFunctionCall(mfilename);
funcLog  = mylog(mfilename, funcCall);

allowed = {'boundary', 'indElectBound', 'numConseqCrit', 'baseType', 'outputDir', 'trialinfoHeader'};
vararginProblem = checkMyArgs(varargin(1:2:end), allowed);
if vararginProblem
  return
end
boundary        = getfield(pmd_option(varargin, 'boundary', 'vector', [], 1), 'boundary');
indElectBound   = getfield(pmd_option(varargin, 'indElectBound', 'bool', 1, 0), 'indElectBound');
numConseqCrit   = getfield(pmd_option(varargin, 'numConseqCrit', 'int', 1, 0), 'numConseqCrit');
baseType        = getfield(pmd_option(varargin, 'baseType', 'string', 'abs', 0, {'abs', 'per'}), 'baseType');
outputDir       = getfield(pmd_option(varargin, 'outputDir', 'string', pwd, 0), 'outputDir');
trialinfoHeader = getfield(pmd_option(varargin, 'trialinfoHeader', 'cell', {''}, 1), 'trialinfoHeader');

fileProblem = checkMyFilesExist(participants, conditions, filetypeIn);
if fileProblem
  return
end

allowed = {'meanAmp', 'maxPeakAmp', 'minPeakAmp', 'maxPeakLat', 'minPeakLat',...
  'meanMaxPeakAmp', 'meanMaxPeakAmp', ...
  'critOnsetSmaller', 'critOnsetLarger', 'perOnsetSmaller', 'perOnsetLarger',...
  'minFracAreaLat', 'maxFracAreaLat', 'minFracPeakLat', 'maxFracPeakLat', ...
  'minTime2peak', 'maxTime2peak', 'peak2peakAmp'};
analysisTypeProblem = checkMyArgs(analysisType, allowed);
if analysisTypeProblem
  return
end

if ismember(analysisType, {'critOnsetSmaller', 'critOnsetLarger', 'meanMaxPeakAmp', 'meanMinPeakAmp', ...
    'perOnsetSmaller', 'perOnsetLarger', 'minFracAreaLat', 'maxFracAreaLat', ...
    'minFracPeakLat', 'maxFracPeakLat'}) && isempty(analysisValue)
  fprintf(2, 'Error (analyseMyERPs): Analysis type requires an analysis value ...\n');
  return
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%% analyseMyERPs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grandaverage analysis with empty participants?
grandAverage = false;
if isempty(participants) && strncmp(filetypeIn, 'GA', 2)
  grandAverage = true;
  participants = 1;
end

% save analysis parameters to txt file
fid = fopen([outputDir filesep filenameOut '.txt'], 'w');
if fid < 0
  fprintf(2, 'Unable to create output file. Does directory exist? ...\n');
  return
end
fprintf(fid, 'Function Call: %s\n\n', funcCall);
fprintf(fid, 'Analysis Parameters\n\n');
if ~grandAverage
  fprintf(fid, ['Participants: ' num2str(participants) '\n']);
else
  fprintf(fid, 'Participants: NA\n');
end
fprintf(fid, ['Conditions: ' num2str(conditions) '\n']);
fprintf(fid, 'Analysis File: %s\n', filetypeIn);
fprintf(fid, 'Analysis Type: %s\n', analysisType);
if ~isempty(analysisValue)
  fprintf(fid, 'Analysis Value: %1.3f\n', analysisValue);
else
  fprintf(fid, 'Analysis Value: NA\n');
end
if ~isempty(baseWin)
  fprintf(fid, 'Baseline Window: %1.3f %1.3f\n', baseWin(1), baseWin(2));
else
  fprintf(fid, 'Baseline Window: NA\n');
end

for i = 1:size(analysisWin, 1)
  fprintf(fid, [strcat('Analysis Window', num2str(i)) ': %1.3f %1.3f\n'], ...
    analysisWin(i, 1), analysisWin(i, 2));
end
fprintf(fid, 'Output File: %s\n\n', filenameOut);
fclose(fid);

% preallocate data matrix a bit bigger than needed
analysisData(length(participants)*length(conditions)*100, 85) = 0;

count = 1;
for vp = participants
  for cond = conditions
    
    % load data
    if ~grandAverage
      fname = sprintf('%d_%d_%s.mat', vp, cond, filetypeIn);
    else
      fname = sprintf('%d_%s.mat', cond, filetypeIn);
    end
    fprintf('reading %s\n', fname);
    load(fullfile(pwd, fname), 'data');
    
    % extra analysis info from data.trialinfo field for single trial data?
    if ismember(filetypeIn, {'all', 'pre'})
      data.trialinfo = data.trialinfo(:, any(data.trialinfo));
      numExtraCols   = size(data.trialinfo, 2) - 1;
    else  % averaged data so no trial info/extra columns (removed later!)
      data.trialinfo = [0 0];
      numExtraCols   = 1;
    end
    
    % baseline
    if ~isempty(baseWin) && strcmp(baseType, 'abs')  % standard subtraction baseline for EEG
      data = pmd_baselineEpoch(data, [baseWin(1) baseWin(2)]);
    elseif ~isempty(baseWin) && strcmp(baseType, 'per')  % percentage baseline for EMG
      data = pmdEMG_baselineEpoch(data, [baseWin(1) baseWin(2)]);
    end
    
    % create fake trial field for averaged data
    if ~ismember(filetypeIn, {'all', 'pre'})
      data.trial{1} = data.avg;
    end
    
    % labels/electrodes/trials
    labels    = data.label;
    numElects = length(data.label);
    numTrls   = length(data.trial);
    
    % analysis windows
    % TO DO: some logic withih here could be shorteded with function handles!
    for winNum = 1:size(analysisWin, 1)
      
      % analysis index start/end
      anaIdx = findStartEndIdx(data.time, analysisWin(winNum, 1), analysisWin(winNum, 2));
      time   = data.time(anaIdx(1):anaIdx(2));
      
      % analysis type
      if strcmp(analysisType, 'meanAmp')
        
        for trl = 1:numTrls
          
          values = mean(data.trial{trl}(:, anaIdx(1):anaIdx(2)), 2)';
          
          analysisData(count, 1:numElects + 6 + numExtraCols) = ...
            [vp, cond, trl, data.trialinfo(trl, 2:end), winNum, analysisWin(winNum, 1), ...
            analysisWin(winNum, 2), values];
          
          count = count + 1;
          
        end
        
      elseif strcmp(analysisType, 'maxPeakAmp')
        
        for trl = 1:numTrls
          
          values = max(data.trial{trl}(:, anaIdx(1):anaIdx(2)), [], 2)';
          
          analysisData(count, 1:numElects + 6 + numExtraCols) = ...
            [vp, cond, trl, data.trialinfo(trl, 2:end), winNum, analysisWin(winNum, 1), ...
            analysisWin(winNum, 2), values];
          
          count = count + 1;
          
        end
        
      elseif strcmp(analysisType, 'minPeakAmp')
        
        for trl = 1:numTrls
          
          values = min(data.trial{trl}(:, anaIdx(1):anaIdx(2)), [], 2)';
          
          analysisData(count, 1:numElects + 6 + numExtraCols) = ...
            [vp, cond, trl, data.trialinfo(trl, 2:end), winNum, analysisWin(winNum, 1), ...
            analysisWin(winNum, 2), values];
          
          count = count + 1;
          
        end
        
      elseif strcmp(analysisType, 'maxPeakLat')
        
        for trl = 1:numTrls
          
          [~, maxIdx] = max(data.trial{trl}(:, anaIdx(1):anaIdx(2)), [], 2);
          values      = round(time(maxIdx)*1000);
          
          analysisData(count, 1:numElects + 6 + numExtraCols) = ...
            [vp, cond, trl, data.trialinfo(trl, 2:end), winNum, analysisWin(winNum, 1), ...
            analysisWin(winNum, 2), values];
          
          count = count + 1;
          
        end
        
      elseif strcmp(analysisType, 'minPeakLat')
        
        for trl = 1:numTrls
          
          [~, minIdx] = min(data.trial{trl}(:, anaIdx(1):anaIdx(2)), [], 2);
          values      = round(time(minIdx)*1000);
          
          analysisData(count, 1:numElects + 6 + numExtraCols) = ...
            [vp, cond, trl, data.trialinfo(trl, 2:end), winNum, analysisWin(winNum, 1), ...
            analysisWin(winNum, 2), values];
          
          count = count + 1;
          
        end
        
      elseif strcmp(analysisType, 'meanMaxPeakAmp')
        
        for trl = 1:numTrls
          
          values = zeros(1, numElects);
          for elect = 1:numElects
            
            % first find peak then adapt analysis window around this peak
            [~, maxIdx] = max(data.trial{trl}(elect, anaIdx(1):anaIdx(2)), [], 2);
            peakAnalysisWin = [time(maxIdx) - analysisValue, time(maxIdx) + analysisValue];
            
            % calculate analysis index around peak
            peakAnaIdx = findStartEndIdx(data.time, peakAnalysisWin(winNum, 1), peakAnalysisWin(winNum, 2));
            
            % now take mean amplitude
            values(elect) = mean(data.trial{trl}(elect, peakAnaIdx(1):peakAnaIdx(2)), 2)';
            
          end
          
          analysisData(count, 1:numElects + 6 + numExtraCols) = ...
            [vp, cond, trl, data.trialinfo(trl, 2:end), winNum, analysisWin(winNum, 1), ...
            analysisWin(winNum, 2), values];
          
          count = count + 1;
          
        end
        
      elseif strcmp(analysisType, 'meanMinPeakAmp')
        
        for trl = 1:numTrls
          
          values = zeros(1, numElects);
          for elect = 1:numElects
            
            % first find peak then adapt analysis window around this peak
            [~, minIdx] = min(data.trial{trl}(elect, anaIdx(1):anaIdx(2)), [], 2);
            peakAnalysisWin = [time(minIdx) - analysisValue, time(minIdx) + analysisValue];
            
            % calculate analysis index around peak
            peakAnaIdx = findStartEndIdx(data.time, peakAnalysisWin(winNum, 1), peakAnalysisWin(winNum, 2));
            
            % now take mean amplitude
            values(elect) = mean(data.trial{trl}(elect, peakAnaIdx(1):peakAnaIdx(2)), 2)';
            
          end
          
          analysisData(count, 1:numElects + 6 + numExtraCols) = ...
            [vp, cond, trl, data.trialinfo(trl, 2:end), winNum, analysisWin(winNum, 1), ...
            analysisWin(winNum, 2), values];
          
          count = count + 1;
          
        end
        
      elseif ismember(analysisType, {'critOnsetSmaller', 'critOnsetLarger', 'perOnsetSmaller', 'perOnsetLarger'})
        
        for trl = 1:numTrls
          
          onsetLocs = findIndex(data.trial{trl}, anaIdx, analysisValue, analysisType);
          
          values = zeros(1, length(data.label));
          for elect = 1:numElects
            
            onsetLoc = find(cumsum(onsetLocs(elect, :)) == numConseqCrit, 1, 'first');
            if ~isempty(onsetLoc)
              values(1, elect) = round(time(onsetLoc)*1000);
            else % print warning to file
              fid = fopen([pwd filesep mfilename '.txt'], 'a+');
              fprintf(fid, 'Problem in participant %d condition %d trial %d electrode %s\n', ...
                vp, cond, trl, data.label{elect});
              fclose(fid);
            end
            
          end
          
          analysisData(count, 1:numElects + 6 + numExtraCols) = ...
            [vp, cond, trl, data.trialinfo(trl, 2:end), winNum, analysisWin(winNum, 1), ...
            analysisWin(winNum, 2), values];
          
          count = count +1;
          
        end
        
      elseif ismember(analysisType, {'minFracAreaLat', 'maxFracAreaLat'})
        
        for trl = 1:numTrls
          
          values    = zeros(1, length(data.label));
          timeIndex = find(anaIdx(1):anaIdx(2));
          for elect = 1:numElects
            
            if isempty(boundary)
              if strcmp(analysisType, 'minFracAreaLat')
                boundary = max(data.trial{trl}(elect, anaIdx(1):anaIdx(2)));
              elseif strcmp(analysisType, 'maxFracAreaLat')
                boundary = min(data.trial{trl}(elect, anaIdx(1):anaIdx(2)));
              end
            end
            areaWhole = trapz(abs(data.trial{trl}(elect, anaIdx(1):anaIdx(2)) - boundary));
            
            indexStart = timeIndex(1);
            indexEnd   = timeIndex(end);
            while indexStart < indexEnd
              
              indexMid = round((indexStart + indexEnd) / 2);
              areaPart = trapz(abs(data.trial{trl}(elect, timeIndex(1):indexMid) - boundary));
              
              if areaPart > (areaWhole * (analysisValue/100))
                indexEnd = indexMid - 1;
              elseif areaPart < (areaWhole * (analysisValue/100))
                indexStart = indexMid + 1;
              else
                break
              end
            end
            
            values(1, elect) = round(time(indexMid)*1000);
            
          end
          
          analysisData(count, 1:numElects + 6 + numExtraCols) = ...
            [vp, cond, trl, data.trialinfo(trl, 2:end), winNum, analysisWin(winNum, 1), ...
            analysisWin(winNum, 2), values];
          
          count = count + 1;
          
        end
        
      elseif ismember(analysisType, {'minFracPeakLat', 'maxFracPeakLat'})
        
        for trl = 1:numTrls
          
          if indElectBound  % adjust signal
            if strcmp(analysisType, 'minFracPeakLat')
              boundary = max(data.trial{trl}(:, anaIdx(1):anaIdx(2)), [], 2);
            elseif strcmp(analysisType, 'maxFracPeakLat')
              boundary = min(data.trial{trl}(:, anaIdx(1):anaIdx(2)), [], 2);
            end
          end
          data.trial{trl}(:, anaIdx(1):anaIdx(2)) = ...
            data.trial{trl}(:, anaIdx(1):anaIdx(2)) - repmat(boundary, 1, length(time));
          
          % find peak
          if strcmp(analysisType, 'minFracPeakLat')
            onsetLocs = findIndex(data.trial{trl}, anaIdx, analysisValue, 'perOnsetSmaller');
          elseif strcmp(analysisType, 'maxFracPeakLat')
            onsetLocs = findIndex(data.trial{trl}, anaIdx, analysisValue, 'perOnsetLarger');
          end
          
          values    = zeros(1, length(data.label));
          for elect = 1:numElects
            onsetLoc = find(cumsum(onsetLocs(elect, :)) == numConseqCrit, 1, 'first');
            if ~isempty(onsetLoc)
              values(1, elect) = round(time(onsetLoc)*1000);
            else  % print warning to file
              values(1, elect) = 9999;
              fid = fopen([pwd filesep mfilename '.txt'], 'a+');
              fprintf(fid, 'Problem in participant %d condition %d trial %d electrode %s\n', ...
                vp, cond, trl, data.label{elect});
              fclose(fid);
            end
          end
          
          analysisData(count, 1:numElects + 6 + numExtraCols) = ...
            [vp, cond, trl, data.trialinfo(trl, 2:end), winNum, analysisWin(winNum, 1), ...
            analysisWin(winNum, 2), values];
          
          count = count + 1;
          
        end
        
      elseif ismember(analysisType, {'minTime2peak', 'maxTime2peak'})
        
        for trl = 1:numTrls
          
          if strcmp(analysisType, 'minTime2peak')
            onsetLocs = findIndex(data.trial{trl}, anaIdx, analysisValue, 'critOnsetSmaller');
            [~, peakIdx] = min(data.trial{trl}(:, anaIdx(1):anaIdx(2)), [], 2);
          elseif strcmp(analysisType, 'maxTime2peak')
            onsetLocs = findIndex(data.trial{trl}, anaIdx, analysisValue, 'critOnsetLarger');
            [~, peakIdx] = max(data.trial{trl}(:, anaIdx(1):anaIdx(2)), [], 2);
          end
          
          values    = zeros(1, length(data.label));
          for elect = 1:numElects
            onsetLoc = find(cumsum(onsetLocs(elect, :)) == numConseqCrit, 1, 'first');
            if ~isempty(onsetLoc)
              values(1, elect) = round(time(peakIdx(elect))*1000) - round(time(onsetLoc)*1000);
            else % print warning to file
              fid = fopen([pwd filesep mfilename '.txt'], 'a+');
              fprintf(fid, 'Problem in participant %d condition %d trial %d electrode %s\n', ...
                vp, cond, trl, data.label{elect});
              fclose(fid);
            end
          end
          
          analysisData(count, 1:numElects + 6 + numExtraCols) = ...
            [vp, cond, trl, data.trialinfo(trl, 2:end), winNum, ...
            analysisWin(winNum, 1), analysisWin(winNum, 2), values];
          
          count = count + 1;
          
        end
        
      elseif strcmp(analysisType, 'peak2peakAmp')
        
        for trl = 1:numTrls
          
          values = range(data.trial{trl}(:, anaIdx(1):anaIdx(2)), 2)';
          
          analysisData(count, 1:numElects + 6 + numExtraCols) = ...
            [vp, cond, trl, data.trialinfo(trl, 2:end), winNum, ...
            analysisWin(winNum, 1), analysisWin(winNum, 2), values];
          
          count = count + 1;
          
        end
        
      end
      
    end
    
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select appropriate columns
if ismember(filetypeIn, {'avg', 'lrp', 'GA_avg', 'GA_lrp'})
  analysisData(:, [3, 4]) = [];  % remove trial number/trial info
  analysisData = analysisData(1:count-1, 1:numElects + 5);
elseif any(ismember({'all', 'pre'}, filetypeIn))
  analysisData = analysisData(1:count-1, 1:numElects + 6 + numExtraCols);
end

% save to .mat file
save([outputDir filesep filenameOut], 'analysisData');

% save to .txt file with appropriate headers
if ismember(filetypeIn, {'all', 'pre'})
  try
    if isempty(trialinfoHeader{1})
      trialinfo = strcat('trialinfo', strtrim(cellstr(num2str(transpose(1:numExtraCols)))'));
    else
      trialinfo = trialinfoHeader{1:numExtraCols};
    end
  catch
    trialinfo = [];
  end
end

% header and format string
if ismember(filetypeIn, {'avg', 'lrp', 'GA_avg', 'GA_lrp'})
  header = [{'VP'}, {'Cond'}, {'WindowNum'}, {'WindowStart'}, ...
    {'WindowEnd'}, transpose(labels)];
elseif ismember(filetypeIn, {'all', 'pre'})
  header = [{'VP'}, {'Cond'}, {'Trial'}, trialinfo, {'WindowNum'}, ...
    {'WindowStart'}, {'WindowEnd'}, transpose(labels)];
end
formatHeader = repmat('%s\t', 1, length(header));
formatHeader(end) = 'n';

formatData = '';
for i = 1:length(header)
  try
    if ~mod(analysisData(:, i), 1) % is integer
      formatData = strcat(formatData, '%1.0f\t');
    else
      formatData = strcat(formatData, '%1.3f\t');
    end
  catch
    formatData = strcat(formatData, '%1.3f\t');
  end
end
formatData(end) = 'n';

fid = fopen([outputDir filesep strcat(filenameOut,'.txt')], 'a+');
fprintf(fid, formatHeader, header{:});
fprintf(fid, formatData, analysisData');
fclose(fid);

funcLog.success = true;
