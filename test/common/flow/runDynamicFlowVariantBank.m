function [variantTable, bundleCell] = runDynamicFlowVariantBank(variantNameList, flowOptCell, ...
  periodicFixture, subsetFixtureCell, pilotWave, carrierFreq, sampleRate, optVerbose)
%RUNDYNAMICFLOWVARIANTBANK Compare multiple dynamic flow profiles on one fixture set.
% This helper reuses the same periodic fixture and subset-fixture bank to
% evaluate several flow-option variants in one shot. It is intended for
% representative/fullDiag work where route comparison is valuable, but the
% caller does not want to rerun snapshot generation or fixture construction
% for every route hypothesis.

arguments
  variantNameList (:, 1) string
  flowOptCell (:, 1) cell
  periodicFixture (1, 1) struct
  subsetFixtureCell cell
  pilotWave
  carrierFreq (1, 1) double
  sampleRate (1, 1) double
  optVerbose (1, 1) logical = false
end

numVariant = numel(variantNameList);
if numel(flowOptCell) ~= numVariant
  error('runDynamicFlowVariantBank:VariantCountMismatch', ...
    'variantNameList and flowOptCell must have the same length.');
end

bundleCell = cell(numVariant, 1);
useParfor = localShouldUseVariantParfor(numVariant);
if useParfor
  parfor iVariant = 1:numVariant
    bundleCell{iVariant} = buildDoaDopplerDynamicTransitionBundle( ...
      periodicFixture, subsetFixtureCell, pilotWave, carrierFreq, sampleRate, ...
      optVerbose, flowOptCell{iVariant});
  end
else
  for iVariant = 1:numVariant
    bundleCell{iVariant} = buildDoaDopplerDynamicTransitionBundle( ...
      periodicFixture, subsetFixtureCell, pilotWave, carrierFreq, sampleRate, ...
      optVerbose, flowOptCell{iVariant});
  end
end

selectedSubsetLabel = strings(numVariant, 1);
selectedFinalTag = strings(numVariant, 1);
finalAngleErrDeg = nan(numVariant, 1);
finalFdRefErrHz = nan(numVariant, 1);
finalFdRateErrHzPerSec = nan(numVariant, 1);
totalBundleSec = nan(numVariant, 1);
subsetSelectTotalSec = nan(numVariant, 1);
inToothRefineSec = nan(numVariant, 1);

for iVariant = 1:numVariant
  bundleUse = bundleCell{iVariant};
  finalSummary = localGetFieldOrDefault(bundleUse, 'selectedFinalSummary', struct());
  selectedSubsetLabel(iVariant) = string(localGetFieldOrDefault(bundleUse, 'selectedSubsetLabel', "none"));
  selectedFinalTag(iVariant) = string(localGetFieldOrDefault(bundleUse, 'selectedFinalTag', "none"));
  finalAngleErrDeg(iVariant) = localGetFieldOrDefault(finalSummary, 'angleErrDeg', NaN);
  finalFdRefErrHz(iVariant) = localGetFieldOrDefault(finalSummary, 'fdRefErrHz', NaN);
  finalFdRateErrHzPerSec(iVariant) = localGetFieldOrDefault(finalSummary, 'fdRateErrHzPerSec', NaN);
  stageTiming = localGetFieldOrDefault(bundleUse, 'stageTiming', struct());
  totalBundleSec(iVariant) = localGetFieldOrDefault(stageTiming, 'totalBundleSec', NaN);
  subsetSelectTotalSec(iVariant) = localGetFieldOrDefault(stageTiming, 'subsetSelectTotalSec', NaN);
  inToothRefineSec(iVariant) = localGetFieldOrDefault(stageTiming, 'inToothRefineSec', NaN);
end

variantTable = table(variantNameList(:), selectedSubsetLabel, selectedFinalTag, ...
  finalAngleErrDeg, finalFdRefErrHz, finalFdRateErrHzPerSec, ...
  totalBundleSec, subsetSelectTotalSec, inToothRefineSec, ...
  'VariableNames', {'variantName', 'selectedSubsetLabel', 'selectedFinalTag', ...
  'angleErrDeg', 'fdRefErrHz', 'fdRateErrHzPerSec', ...
  'totalBundleSec', 'subsetSelectTotalSec', 'inToothRefineSec'});
end


function tf = localShouldUseVariantParfor(numVariant)
%LOCALSHOULDUSEVARIANTPARFOR Use parfor only for real variant banks.

tf = false;
if numVariant < 2
  return;
end
if ~localCanUseParfor()
  return;
end
tf = true;
end


function tf = localCanUseParfor()
%LOCALCANUSEPARFOR Check whether parfor can be used safely here.

tf = false;
if ~isempty(getCurrentTask())
  return;
end
if ~license('test', 'Distrib_Computing_Toolbox')
  return;
end
try
  tf = ~isempty(ver('parallel'));
catch
  tf = false;
end
end


function value = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Return one struct field or a default value.

value = defaultValue;
if isstruct(dataStruct) && isfield(dataStruct, fieldName)
  rawValue = dataStruct.(fieldName);
  if ~isempty(rawValue)
    value = rawValue;
  end
end
end
