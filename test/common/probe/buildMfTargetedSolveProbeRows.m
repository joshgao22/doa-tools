function [solveRows, doaBasinEntryRows] = buildMfTargetedSolveProbeRows(method, baselineCase, staticSeedCase, ssMfSeedCase, truth, task, ...
  viewUse, context, fdRangeUse, fdRateRangeUse, dynOpt, debugTruthUse, optVerbose, config, objectiveProbeRows)
%BUILDMFTARGETEDSOLVEPROBEROWS Rerun controlled MF solver probes.
% This helper keeps the replay's fixed targeted route: every method records a
% baseline row, expensive reruns are gated to MS failures, and the replay-only
% MS CP-U bank is delegated to buildMsCpuBankProbeRows. The ms-bank-only route
% skips deep solve reruns and keeps only baseline plus MS bank rows. It does
% not alter the estimator result or final winner.

solveRows = repmat(emptyMfSolveProbeRow(), 0, 1);
doaBasinEntryRows = repmat(emptyMfDoaBasinEntryRow(), 0, 1);
solveRows = [solveRows; buildMfSolveProbeRowFromCase(method, baselineCase, truth, task, ...
  "baseline", config.staticLocalDoaHalfWidthDeg, NaN, "baseline")];
runDeepSolveProbe = localShouldRunDeepSolveProbe(method, baselineCase, truth, objectiveProbeRows, config);
runMsBankProbe = localShouldRunMsBankProbe(method, baselineCase, objectiveProbeRows, config);
if ~(runDeepSolveProbe || runMsBankProbe)
  return;
end

truthDoa = takeFirstTwoValues(localGetFieldOrDefault(truth, 'latlonTrueDeg', []));
staticEst = localGetFieldOrDefault(staticSeedCase, 'estResult', struct());
staticDoa = takeFirstTwoValues(localGetFieldOrDefault(staticEst, 'doaParamEst', []));
truthFdRefHz = resolveMfProbeTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = resolveMfProbeTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});

probeList = repmat(struct('tag', "", 'centerDoa', NaN(2, 1), 'halfWidth', NaN(2, 1)), 0, 1);
probeItem = struct('tag', "truth-doa-truth-fd", 'centerDoa', truthDoa, ...
  'halfWidth', localExpandDoaHalfWidth(config.probeTruthDoaHalfWidthDeg));
probeList = [probeList; probeItem]; %#ok<AGROW>
runFullWideProbe = ~(string(method.satMode) == "multi" && ~method.isKnownRate);
wideList = reshape(double(config.probeWideDoaHalfWidthDegList), [], 1);
if runFullWideProbe
  for iWide = 1:numel(wideList)
    probeItem = struct('tag', string(sprintf('static-doa-wide%.4g', wideList(iWide))), ...
      'centerDoa', staticDoa, 'halfWidth', [wideList(iWide); wideList(iWide)]);
    probeList = [probeList; probeItem]; %#ok<AGROW>
  end
end

if string(method.satMode) == "multi" && ~isempty(fieldnames(ssMfSeedCase))
  ssMfEst = localGetFieldOrDefault(ssMfSeedCase, 'estResult', struct());
  ssMfDoa = takeFirstTwoValues(localGetFieldOrDefault(ssMfEst, 'doaParamEst', []));
  if all(isfinite(ssMfDoa))
    probeItem = struct('tag', "ss-mf-seed-truth-fd", 'centerDoa', ssMfDoa, ...
      'halfWidth', localExpandDoaHalfWidth(config.probeTruthDoaHalfWidthDeg));
    probeList = [probeList; probeItem]; %#ok<AGROW>
  end
end

if runDeepSolveProbe
  for iProbe = 1:numel(probeList)
    initParam = makeMfProbeOptVar(method, probeList(iProbe).centerDoa, truthFdRefHz, truthFdRateHzPerSec);
    initCandidate = struct('initParam', initParam, 'initDoaParam', probeList(iProbe).centerDoa, ...
      'initDoaHalfWidth', probeList(iProbe).halfWidth, 'startTag', probeList(iProbe).tag);
    try
      probeCase = runDynamicDoaDopplerCase(method.displayName, method.satMode, ...
        viewUse, truth, context.pilotWave, context.carrierFreq, context.waveInfo.sampleRate, ...
        fdRangeUse, fdRateRangeUse, optVerbose, dynOpt, method.isKnownRate, ...
        debugTruthUse, initCandidate);
      solveRows = [solveRows; buildMfSolveProbeRowFromCase(method, probeCase, truth, task, ...
        probeList(iProbe).tag, probeList(iProbe).halfWidth, solveRows(1).finalObj, "probe")]; %#ok<AGROW>
      doaBasinEntryRows = [doaBasinEntryRows; buildMfDoaBasinEntryRowsFromCase( ...
        method, probeCase, truth, task, probeList(iProbe).tag, "probe", solveRows(1).finalObj)]; %#ok<AGROW>
    catch ME
      row = emptyMfSolveProbeRow();
      row.displayName = method.displayName;
      row.snrDb = task.snrDb;
      row.taskSeed = task.taskSeed;
      row.probeTag = probeList(iProbe).tag;
      row.probeGroup = "probe";
      row.doaHalfWidthLatDeg = probeList(iProbe).halfWidth(1);
      row.doaHalfWidthLonDeg = probeList(iProbe).halfWidth(2);
      row.evalOk = false;
      row.message = string(ME.message);
      solveRows = [solveRows; row]; %#ok<AGROW>
    end
  end
end

if runMsBankProbe
  [bankRows, bankBasinRows] = buildMsCpuBankProbeRows(method, staticSeedCase, ssMfSeedCase, truth, task, ...
    viewUse, context, fdRangeUse, fdRateRangeUse, dynOpt, debugTruthUse, optVerbose, config, solveRows(1).finalObj);
  solveRows = [solveRows; bankRows(:)]; %#ok<AGROW>
  doaBasinEntryRows = [doaBasinEntryRows; bankBasinRows(:)]; %#ok<AGROW>
end
end


function shouldRun = localShouldRunDeepSolveProbe(method, baselineCase, truth, objectiveProbeRows, config)
%LOCALSHOULDRUNDEEPSOLVEPROBE Gate expensive controlled solver reruns.

probeRoute = string(localGetFieldOrDefault(config, 'solveProbeRoute', "ms-targeted"));
if probeRoute == "disabled" || probeRoute == "ms-bank-only"
  shouldRun = false;
  return;
end
if string(method.satMode) ~= "multi"
  shouldRun = probeRoute == "ss-parity";
  return;
end

estResult = localGetFieldOrDefault(baselineCase, 'estResult', struct());
estAux = localGetFieldOrDefault(estResult, 'aux', struct());
estDebug = localGetFieldOrDefault(estAux, 'debug', struct());
finalEval = localGetFieldOrDefault(estDebug, 'finalEval', struct());
optimInfo = localGetFieldOrDefault(estResult, 'optimInfo', struct());
[~, nonRefCoherenceFloor] = extractMfProbeCoherenceFloor(finalEval);
firstOrderOpt = localGetFieldOrDefault(optimInfo, 'firstorderopt', NaN);
fdRefErrHz = abs(localGetFieldOrDefault(estResult, 'fdRefEst', NaN) - ...
  resolveMfProbeTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'}));

finalObj = localProbeObj(objectiveProbeRows, "finalPoint");
truthObj = localProbeObj(objectiveProbeRows, "truthPoint");
truthDoaFinalFdObj = localProbeObj(objectiveProbeRows, "truthDoaFinalFdPoint");
objTol = max(1e-6, 1e-9 * max(abs(finalObj), 1));
truthPointBetter = isfinite(truthObj) && isfinite(finalObj) && truthObj < finalObj - objTol;
truthDoaBetter = isfinite(truthDoaFinalFdObj) && isfinite(finalObj) && truthDoaFinalFdObj < finalObj - objTol;
nonRefCollapse = isfinite(nonRefCoherenceFloor) && ...
  nonRefCoherenceFloor < config.deepProbeNonRefCoherenceThreshold;
conditioningBad = isfinite(firstOrderOpt) && firstOrderOpt > config.deepProbeFirstOrderOptThreshold;

if isfinite(localProbeFdAbsErr(objectiveProbeRows, "finalPoint"))
  fdRefErrHz = localProbeFdAbsErr(objectiveProbeRows, "finalPoint");
end
fdRefBranch = isfinite(fdRefErrHz) && fdRefErrHz > config.deepProbeFdRefAbsHzThreshold;
shouldRun = truthPointBetter || truthDoaBetter || nonRefCollapse || conditioningBad || fdRefBranch;
end


function shouldRun = localShouldRunMsBankProbe(method, baselineCase, objectiveProbeRows, config)
%LOCALSHOULDRUNMSBANKPROBE Decide whether to run replay-only MS CP-U bank probes.

shouldRun = false;
if string(localGetFieldOrDefault(config, 'msBankRoute', "cpu-release-only")) == "disabled"
  return;
end
if string(method.satMode) ~= "multi" || method.isKnownRate
  return;
end
probeRoute = string(localGetFieldOrDefault(config, 'solveProbeRoute', "ms-targeted"));
if probeRoute == "disabled"
  return;
end
if probeRoute == "ms-bank-only"
  shouldRun = true;
  return;
end
truthPointBetter = false;
finalObj = localProbeObj(objectiveProbeRows, "finalPoint");
truthObj = localProbeObj(objectiveProbeRows, "truthPoint");
if isfinite(finalObj) && isfinite(truthObj)
  truthPointBetter = truthObj < finalObj - max(1e-6, 1e-9 * max(abs(finalObj), 1));
end
estResult = localGetFieldOrDefault(baselineCase, 'estResult', struct());
estAux = localGetFieldOrDefault(estResult, 'aux', struct());
estDebug = localGetFieldOrDefault(estAux, 'debug', struct());
finalEval = localGetFieldOrDefault(estDebug, 'finalEval', struct());
[~, nonRefCoherenceFloor] = extractMfProbeCoherenceFloor(finalEval);
optimInfo = localGetFieldOrDefault(estResult, 'optimInfo', struct());
firstOrderOpt = localGetFieldOrDefault(optimInfo, 'firstorderopt', NaN);
nonRefCollapse = isfinite(nonRefCoherenceFloor) && ...
  nonRefCoherenceFloor < config.deepProbeNonRefCoherenceThreshold;
conditioningBad = isfinite(firstOrderOpt) && firstOrderOpt > config.deepProbeFirstOrderOptThreshold;
shouldRun = truthPointBetter || nonRefCollapse || conditioningBad;
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


function obj = localProbeObj(probeRows, probeTag)
%LOCALPROBEOBJ Extract one objective value from objective probe rows.

obj = NaN;
if isempty(probeRows)
  return;
end
if istable(probeRows)
  if height(probeRows) == 0
    return;
  end
  idx = find(probeRows.probeTag == string(probeTag), 1, 'first');
  if ~isempty(idx)
    obj = probeRows.objective(idx);
  end
elseif isstruct(probeRows)
  tagList = string({probeRows.probeTag});
  idx = find(tagList == string(probeTag), 1, 'first');
  if ~isempty(idx)
    obj = probeRows(idx).objective;
  end
end
end


function fdAbsErrHz = localProbeFdAbsErr(probeRows, probeTag)
%LOCALPROBEFDABSERR Return absolute fdRef error from an objective probe row.

fdAbsErrHz = NaN;
if isempty(probeRows)
  return;
end
if istable(probeRows)
  if height(probeRows) == 0
    return;
  end
  idx = find(probeRows.probeTag == string(probeTag), 1, 'first');
  if ~isempty(idx)
    fdAbsErrHz = abs(probeRows.fdRefErrHz(idx));
  end
elseif isstruct(probeRows)
  tagList = string({probeRows.probeTag});
  idx = find(tagList == string(probeTag), 1, 'first');
  if ~isempty(idx)
    fdAbsErrHz = abs(probeRows(idx).fdRefErrHz);
  end
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
