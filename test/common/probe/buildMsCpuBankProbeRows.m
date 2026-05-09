function [solveRows, doaBasinEntryRows] = buildMsCpuBankProbeRows(method, staticSeedCase, ssMfSeedCase, truth, task, ...
  viewUse, context, fdRangeUse, fdRateRangeUse, dynOpt, debugTruthUse, optVerbose, config, baselineObj)
%BUILDMSCPUBANKPROBEROWS Run replay-side MS CP-U bank probes.
% The bank is diagnostic only: CP-K candidates are used for cheap
% preselection, then the top candidates are released with CP-U. This helper
% does not change estimator winners or define a formal adoption rule.

solveRows = repmat(emptyMfSolveProbeRow(), 0, 1);
doaBasinEntryRows = repmat(emptyMfDoaBasinEntryRow(), 0, 1);
if string(method.satMode) ~= "multi" || method.isKnownRate
  return;
end

bankList = localBuildMsAxisBankCandidateList(staticSeedCase, ssMfSeedCase, config, task);
if isempty(bankList)
  return;
end
truthFdRefHz = resolveMfProbeTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = resolveMfProbeTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
if ~(isfinite(truthFdRefHz) && isfinite(truthFdRateHzPerSec))
  return;
end

methodKnown = method;
methodKnown.isKnownRate = true;
methodKnown.fdRateMode = "known";
dynOptKnown = dynOpt;
dynOptKnown.fdRateMode = 'known';
dynOptKnown.fdRateKnown = truth.fdRateFit;
preRows = repmat(emptyMfSolveProbeRow(), 0, 1);
preBasinRows = repmat(emptyMfDoaBasinEntryRow(), 0, 1);
for iBank = 1:numel(bankList)
  tag = "ms-bank-cpk-" + bankList(iBank).tag;
  [row, basinRows] = localRunOneMsBankCandidate(methodKnown, tag, bankList(iBank).centerDoa, ...
    bankList(iBank).halfWidth, truthFdRefHz, truthFdRateHzPerSec, truth, task, ...
    viewUse, context, fdRangeUse, [], dynOptKnown, debugTruthUse, optVerbose, baselineObj, "ms-bank-cpk-preselect");
  preRows = [preRows; row]; %#ok<AGROW>
  preBasinRows = [preBasinRows; basinRows(:)]; %#ok<AGROW>
end
solveRows = [solveRows; preRows(:)];
doaBasinEntryRows = [doaBasinEntryRows; preBasinRows(:)];
if isempty(preRows)
  return;
end

preTable = struct2table(preRows(:));
validMask = preTable.evalOk & isfinite(preTable.finalObj) & isfinite(preTable.finalLatDeg) & isfinite(preTable.finalLonDeg);
if ~any(validMask)
  return;
end
validIdx = find(validMask);
[~, order] = sort(preTable.finalObj(validMask), 'ascend');
releaseTopK = max(0, round(localGetFieldOrDefault(config, 'msBankReleaseTopK', 1)));
releaseIdx = validIdx(order(1:min(releaseTopK, numel(order))));
for iRel = 1:numel(releaseIdx)
  preRow = preTable(releaseIdx(iRel), :);
  releaseCenter = [preRow.finalLatDeg(1); preRow.finalLonDeg(1)];
  releaseFdRefHz = truthFdRefHz + preRow.finalFdRefErrHz(1);
  if ~isfinite(releaseFdRefHz)
    releaseFdRefHz = truthFdRefHz;
  end
  releaseTag = "ms-bank-cpu-release-" + erase(string(preRow.probeTag(1)), "ms-bank-cpk-");
  [row, basinRows] = localRunOneMsBankCandidate(method, releaseTag, releaseCenter, ...
    localExpandDoaHalfWidth(localGetFieldOrDefault(config, 'msBankCompactHalfWidthDeg', [0.002; 0.002])), ...
    releaseFdRefHz, truthFdRateHzPerSec, truth, task, viewUse, context, fdRangeUse, ...
    fdRateRangeUse, dynOpt, debugTruthUse, optVerbose, baselineObj, "ms-bank-cpu-release");
  solveRows = [solveRows; row]; %#ok<AGROW>
  doaBasinEntryRows = [doaBasinEntryRows; basinRows(:)]; %#ok<AGROW>
end
end


function [row, doaBasinEntryRows] = localRunOneMsBankCandidate(method, tag, centerDoa, halfWidth, fdRefHz, fdRateHzPerSec, truth, task, ...
  viewUse, context, fdRangeUse, fdRateRangeUse, dynOpt, debugTruthUse, optVerbose, baselineObj, probeGroup)
%LOCALRUNONEMSBANKCANDIDATE Run one compact replay-only bank candidate.

row = emptyMfSolveProbeRow();
doaBasinEntryRows = repmat(emptyMfDoaBasinEntryRow(), 0, 1);
initParam = makeMfProbeOptVar(method, centerDoa, fdRefHz, fdRateHzPerSec);
initCandidate = struct('initParam', initParam, 'initDoaParam', centerDoa(:), ...
  'initDoaHalfWidth', halfWidth(:), 'startTag', string(tag));
try
  probeCase = runDynamicDoaDopplerCase(method.displayName, method.satMode, ...
    viewUse, truth, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
    fdRangeUse, fdRateRangeUse, optVerbose, dynOpt, method.isKnownRate, ...
    debugTruthUse, initCandidate);
  row = buildMfSolveProbeRowFromCase(method, probeCase, truth, task, ...
    string(tag), halfWidth, baselineObj, string(probeGroup));
  doaBasinEntryRows = buildMfDoaBasinEntryRowsFromCase( ...
    method, probeCase, truth, task, string(tag), string(probeGroup), baselineObj);
catch ME
  row.displayName = method.displayName;
  row.snrDb = task.snrDb;
  row.taskSeed = task.taskSeed;
  row.probeTag = string(tag);
  row.probeGroup = string(probeGroup);
  row.doaHalfWidthLatDeg = halfWidth(1);
  row.doaHalfWidthLonDeg = halfWidth(min(2, numel(halfWidth)));
  row.evalOk = false;
  row.message = string(ME.message);
end
end


function bankList = localBuildMsAxisBankCandidateList(staticSeedCase, ssMfSeedCase, config, task)
%LOCALBUILDMSAXISBANKCANDIDATELIST Build truth-free static/SS-MF axis-cross DoA centers.

bankList = repmat(struct('tag', "", 'centerDoa', NaN(2, 1), 'halfWidth', NaN(2, 1)), 0, 1);
compactHalfWidth = localExpandDoaHalfWidth(localGetFieldOrDefault(config, 'msBankCompactHalfWidthDeg', [0.002; 0.002]));
stepList = reshape(double(localGetFieldOrDefault(config, 'msBankAxisStepDegList', [0.006; 0.012; 0.024])), [], 1);
wideStep = double(localGetFieldOrDefault(config, 'msBankWideAxisStepDeg', NaN));
wideMaxSnr = double(localGetFieldOrDefault(config, 'msBankWideStepMaxSnrDb', -Inf));
if isfinite(wideStep) && wideStep > 0 && isfield(task, 'snrDb') && task.snrDb <= wideMaxSnr
  stepList = unique([stepList(:); wideStep], 'stable');
end
staticEst = localGetFieldOrDefault(staticSeedCase, 'estResult', struct());
staticDoa = takeFirstTwoValues(localGetFieldOrDefault(staticEst, 'doaParamEst', []));
bankList = [bankList; localBuildOneAxisBankSource("static", staticDoa, stepList, compactHalfWidth)]; %#ok<AGROW>
if ~isempty(fieldnames(ssMfSeedCase))
  ssMfEst = localGetFieldOrDefault(ssMfSeedCase, 'estResult', struct());
  ssMfDoa = takeFirstTwoValues(localGetFieldOrDefault(ssMfEst, 'doaParamEst', []));
  bankList = [bankList; localBuildOneAxisBankSource("ssmf", ssMfDoa, stepList, compactHalfWidth)]; %#ok<AGROW>
end
end


function itemList = localBuildOneAxisBankSource(sourceName, centerDoa, stepList, compactHalfWidth)
%LOCALBUILDONEAXISBANKSOURCE Build center and axis-cross offsets for one source.

itemList = repmat(struct('tag', "", 'centerDoa', NaN(2, 1), 'halfWidth', NaN(2, 1)), 0, 1);
if numel(centerDoa) < 2 || any(~isfinite(centerDoa(1:2)))
  return;
end
centerDoa = centerDoa(1:2);
item = struct('tag', sourceName + "-center", 'centerDoa', centerDoa(:), 'halfWidth', compactHalfWidth(:));
itemList = [itemList; item]; %#ok<AGROW>
for iStep = 1:numel(stepList)
  step = stepList(iStep);
  if ~(isfinite(step) && step > 0)
    continue;
  end
  token = localBankStepToken(step);
  offsetMat = [step, 0; -step, 0; 0, step; 0, -step];
  tagList = ["latp"; "latm"; "lonp"; "lonm"];
  for iOff = 1:size(offsetMat, 1)
    item = struct('tag', sourceName + "-" + tagList(iOff) + token, ...
      'centerDoa', centerDoa(:) + offsetMat(iOff, :).', 'halfWidth', compactHalfWidth(:));
    itemList = [itemList; item]; %#ok<AGROW>
  end
end
end


function token = localBankStepToken(step)
%LOCALBANKSTEPTOKEN Build a compact tag token for a DoA step.

token = replace(string(sprintf('%.4g', step)), '.', 'p');
token = replace(token, '-', 'm');
end


function halfWidth = localExpandDoaHalfWidth(value)
%LOCALEXPANDDOAHALFWIDTH Convert scalar/two-vector half-width to 2x1.

value = reshape(double(value), [], 1);
if isempty(value)
  halfWidth = [NaN; NaN];
elseif numel(value) == 1
  halfWidth = repmat(value, 2, 1);
else
  halfWidth = value(1:2);
end
end


function value = localScalarOrDefault(valueIn, defaultValue)
%LOCALSCALARORDEFAULT Return a finite scalar or default.

value = defaultValue;
valueIn = reshape(double(valueIn), [], 1);
if ~isempty(valueIn) && isfinite(valueIn(1))
  value = valueIn(1);
end
end


function val = localGetFieldOrDefault(s, fieldName, defaultVal)
%LOCALGETFIELDORDEFAULT Get a field with a default fallback.

if nargin < 3
  defaultVal = [];
end
val = defaultVal;
if isstruct(s) && isfield(s, fieldName)
  val = s.(fieldName);
end
end
