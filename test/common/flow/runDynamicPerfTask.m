function taskOut = runDynamicPerfTask(taskInfo, sharedData)
%RUNDYNAMICPERFTASK Run one dynamic perf/checkpoint task.
% The full perf script and the perf smoke regression use this same task
% adapter so checkpoint/resume tests exercise the production task path.

context = sharedData.context;
flowOpt = sharedData.flowOpt;
snrDbUse = sharedData.snrDbList(taskInfo.snrIndex);
repeatData = buildDynamicRepeatData(context, snrDbUse, taskInfo.taskSeed);

bundle = buildDoaDopplerDynamicTransitionBundle( ...
  repeatData.periodicFixture, repeatData.subsetFixtureCell, ...
  context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
  localGetFieldOrDefault(sharedData, 'optVerbose', false), flowOpt);

numCase = sharedData.numCase;
idxMsDynUnknown = sharedData.idxMsDynUnknown;
caseList = bundle.caseResult;
if isfield(bundle, 'caseDynMsUnknownFinal') && numel(caseList) >= idxMsDynUnknown
  caseList(idxMsDynUnknown) = bundle.caseDynMsUnknownFinal;
end

angleVec = nan(1, numCase);
fdVec = nan(1, numCase);
fdRateVec = nan(1, numCase);
resolvedVec = false(1, numCase);
for iCase = 1:min(numel(caseList), numCase)
  [angleVec(iCase), fdVec(iCase), fdRateVec(iCase), resolvedVec(iCase)] = ...
    extractDynamicCaseMetric(caseList(iCase), repeatData.truth);
end

if numCase >= idxMsDynUnknown
  finalSummary = localGetFieldOrDefault(bundle, 'selectedFinalSummary', struct());
  angleOverride = localGetFieldOrDefault(finalSummary, 'angleErrDeg', NaN);
  fdOverride = abs(localGetFieldOrDefault(finalSummary, 'fdRefErrHz', NaN));
  fdRateOverride = abs(localGetFieldOrDefault(finalSummary, 'fdRateErrHzPerSec', NaN));
  resolvedOverride = logical(localGetFieldOrDefault(finalSummary, 'isResolved', resolvedVec(idxMsDynUnknown)));
  if isfinite(angleOverride)
    angleVec(idxMsDynUnknown) = angleOverride;
  end
  if isfinite(fdOverride)
    fdVec(idxMsDynUnknown) = fdOverride;
  end
  if isfinite(fdRateOverride)
    fdRateVec(idxMsDynUnknown) = fdRateOverride;
  end
  resolvedVec(idxMsDynUnknown) = resolvedOverride;
end

taskOut = struct();
taskOut.angleVec = angleVec;
taskOut.fdVec = fdVec;
taskOut.fdRateVec = fdRateVec;
taskOut.resolvedVec = resolvedVec;
taskOut.selectedSubsetLabel = string(localGetFieldOrDefault(bundle, 'selectedSubsetLabel', "none"));
taskOut.selectedFinalTag = string(localGetFieldOrDefault(bundle, 'selectedFinalTag', "none"));
taskOut.selectedFinalSolveVariant = string(localGetFieldOrDefault( ...
  localGetFieldOrDefault(bundle, 'selectedFinalSummary', struct()), 'solveVariant', ""));
end

function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one nonempty field with default fallback.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName) && ~isempty(dataStruct.(fieldName))
  value = dataStruct.(fieldName);
end
end
