function objectiveProbeRows = buildMfObjectiveProbeRows(method, caseUse, staticSeedCase, truth, task, ...
  viewUse, pilotWave, carrierFreq, sampleRate, fdRangeUse, fdRateRangeUse, dynOpt, initParamOverride)
%BUILDMFOBJECTIVEPROBEROWS Evaluate fixed MF objective probe points.
% This helper builds the same MF model as the replay estimator path, evaluates
% static-seed, final, truth, and mixed DoA/frequency points, and returns a
% compact typed row array for replay / scan diagnostics. It does not run a
% solver, choose a winner, or change estimator adoption.

objectiveProbeRows = repmat(emptyMfObjectiveProbeRow(), 0, 1);
try
  [model, ~, ~] = buildDoaDopplerMfModel(viewUse.sceneSeq, viewUse.rxSigMf, ...
    pilotWave, carrierFreq, sampleRate, viewUse.doaGrid, fdRangeUse, fdRateRangeUse, dynOpt);
  [initParamUse, ~, model] = buildDoaDopplerMfInit(model, initParamOverride);
catch ME
  row = emptyMfObjectiveProbeRow();
  row.displayName = method.displayName;
  row.snrDb = task.snrDb;
  row.taskSeed = task.taskSeed;
  row.probeTag = "model-build-failed";
  row.evalOk = false;
  row.evalMessage = string(ME.message);
  objectiveProbeRows = row;
  return;
end

pointMap = localBuildObjectiveProbePointMap(method, caseUse, staticSeedCase, truth, initParamUse);
tagList = string(fieldnames(pointMap));
tagList = tagList(:);
rowList = repmat(emptyMfObjectiveProbeRow(), numel(tagList), 1);
finalObj = localEvaluateProbeObjective(model, localGetStructField(pointMap, 'finalPoint', []));
for iTag = 1:numel(tagList)
  optVar = pointMap.(char(tagList(iTag)));
  rowList(iTag) = localBuildObjectiveProbeRow(model, method, truth, task, tagList(iTag), optVar, finalObj);
end
objectiveProbeRows = rowList(:);
end


function pointMap = localBuildObjectiveProbePointMap(method, caseUse, staticSeedCase, truth, initParamUse)
%LOCALBUILDOBJECTIVEPROBEPOINTMAP Build static/final/truth/mixed objective points.

truthDoa = reshape(localGetStructField(truth, 'latlonTrueDeg', []), [], 1);
truthDoa = truthDoa(1:min(2, numel(truthDoa)));
truthFdRefHz = resolveMfProbeTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = resolveMfProbeTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});

staticEst = localGetStructField(staticSeedCase, 'estResult', struct());
staticDoa = reshape(localGetStructField(staticEst, 'doaParamEst', []), [], 1);
staticDoa = staticDoa(1:min(2, numel(staticDoa)));

estResult = localGetStructField(caseUse, 'estResult', struct());
finalDoa = reshape(localGetStructField(estResult, 'doaParamEst', []), [], 1);
finalDoa = finalDoa(1:min(2, numel(finalDoa)));
finalFdRefHz = localGetStructField(estResult, 'fdRefEst', NaN);
finalFdRateHzPerSec = localGetStructField(estResult, 'fdRateEst', truthFdRateHzPerSec);

if method.isKnownRate
  pointMap = struct( ...
    'initPoint', makeMfProbeOptVar(method, takeFirstTwoValues(initParamUse), initParamUse(3), truthFdRateHzPerSec), ...
    'staticSeedPoint', makeMfProbeOptVar(method, staticDoa, truthFdRefHz, truthFdRateHzPerSec), ...
    'finalPoint', makeMfProbeOptVar(method, finalDoa, finalFdRefHz, finalFdRateHzPerSec), ...
    'truthPoint', makeMfProbeOptVar(method, truthDoa, truthFdRefHz, truthFdRateHzPerSec), ...
    'truthDoaFinalFdPoint', makeMfProbeOptVar(method, truthDoa, finalFdRefHz, finalFdRateHzPerSec), ...
    'finalDoaTruthFdPoint', makeMfProbeOptVar(method, finalDoa, truthFdRefHz, truthFdRateHzPerSec));
else
  pointMap = struct( ...
    'initPoint', makeMfProbeOptVar(method, takeFirstTwoValues(initParamUse), initParamUse(3), initParamUse(4)), ...
    'staticSeedPoint', makeMfProbeOptVar(method, staticDoa, truthFdRefHz, truthFdRateHzPerSec), ...
    'finalPoint', makeMfProbeOptVar(method, finalDoa, finalFdRefHz, finalFdRateHzPerSec), ...
    'truthPoint', makeMfProbeOptVar(method, truthDoa, truthFdRefHz, truthFdRateHzPerSec), ...
    'truthDoaFinalFdPoint', makeMfProbeOptVar(method, truthDoa, finalFdRefHz, finalFdRateHzPerSec), ...
    'finalDoaTruthFdPoint', makeMfProbeOptVar(method, finalDoa, truthFdRefHz, truthFdRateHzPerSec));
end
end


function row = localBuildObjectiveProbeRow(model, method, truth, task, probeTag, optVar, finalObj)
%LOCALBUILDOBJECTIVEPROBEROW Build one point-evaluation row.

row = emptyMfObjectiveProbeRow();
row.displayName = method.displayName;
row.snrDb = task.snrDb;
row.taskSeed = task.taskSeed;
row.probeTag = string(probeTag);
row.objective = localEvaluateProbeObjective(model, optVar);
row.finalObjective = finalObj;
row.objMinusFinal = row.objective - finalObj;
row.evalOk = isfinite(row.objective);
if ~row.evalOk
  row.evalMessage = "objective-eval-failed";
end
row.angleErrDeg = localProbeAngleErr(optVar, truth);
[row.fdRefErrHz, row.fdRateErrHzPerSec] = localProbeFdErr(method, optVar, truth);
end


function obj = localEvaluateProbeObjective(model, optVar)
%LOCALEVALUATEPROBEOBJECTIVE Safely evaluate one MF objective point.

obj = NaN;
try
  optVar = reshape(double(optVar), [], 1);
  if numel(optVar) ~= localGetProbeNumVar(model) || any(~isfinite(optVar))
    return;
  end
  obj = evaluateDoaDopplerMfObjective(model, optVar);
catch
  obj = NaN;
end
end


function numVar = localGetProbeNumVar(model)
%LOCALGETPROBENUMVAR Number of optimizer variables for an MF model.

numVar = 3;
if isfield(model, 'fdRateMode') && strcmp(model.fdRateMode, 'unknown')
  numVar = 4;
end
end


function angleErrDeg = localProbeAngleErr(optVar, truth)
%LOCALPROBEANGLEERR Compute angular error for one probe point.

angleErrDeg = NaN;
truthDoa = reshape(localGetStructField(truth, 'latlonTrueDeg', []), [], 1);
optVar = reshape(double(optVar), [], 1);
if numel(optVar) >= 2 && numel(truthDoa) >= 2 && all(isfinite(optVar(1:2)))
  angleErrDeg = calcLatlonAngleError(optVar(1:2), truthDoa(1:2));
end
end


function [fdRefErrHz, fdRateErrHzPerSec] = localProbeFdErr(method, optVar, truth)
%LOCALPROBEFDERR Compute fdRef/fdRate errors for one probe point.

fdRefErrHz = NaN;
fdRateErrHzPerSec = NaN;
truthFdRefHz = resolveMfProbeTruthScalar(truth, {'fdRefFit', 'fdRefTrueHz'});
truthFdRateHzPerSec = resolveMfProbeTruthScalar(truth, {'fdRateFit', 'fdRateTrueHzPerSec'});
optVar = reshape(double(optVar), [], 1);
if numel(optVar) >= 3
  fdRefErrHz = optVar(3) - truthFdRefHz;
end
if method.isKnownRate
  fdRateErrHzPerSec = 0;
elseif numel(optVar) >= 4
  fdRateErrHzPerSec = optVar(4) - truthFdRateHzPerSec;
end
end


function val = localGetStructField(s, fieldName, defaultVal)
%LOCALGETSTRUCTFIELD Get one struct field with a default fallback.

if nargin < 3
  defaultVal = [];
end

val = defaultVal;
if isstruct(s) && isfield(s, fieldName)
  val = s.(fieldName);
end
end
