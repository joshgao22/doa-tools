function row = buildMfSolveProbeRowFromCase(method, caseUse, truth, task, probeTag, halfWidthDeg, baselineObj, probeGroup)
%BUILDMFSOLVEPROBEROWFROMCASE Convert one probe solve to a compact row.

row = emptyMfSolveProbeRow();
row.displayName = method.displayName;
row.snrDb = task.snrDb;
row.taskSeed = task.taskSeed;
row.probeTag = string(probeTag);
row.probeGroup = string(probeGroup);
row.doaHalfWidthLatDeg = halfWidthDeg(1);
row.doaHalfWidthLonDeg = halfWidthDeg(min(2, numel(halfWidthDeg)));
info = summarizeDoaDopplerCase(caseUse, truth);
estResult = localGetFieldOrDefault(caseUse, 'estResult', struct());
optimInfo = localGetFieldOrDefault(estResult, 'optimInfo', struct());
estAux = localGetFieldOrDefault(estResult, 'aux', struct());
estDebug = localGetFieldOrDefault(estAux, 'debug', struct());
finalEval = localGetFieldOrDefault(estDebug, 'finalEval', struct());
row.finalLatDeg = info.latEstDeg;
row.finalLonDeg = info.lonEstDeg;
row.finalAngleErrDeg = info.angleErrDeg;
row.finalFdRefErrHz = info.fdRefErrHz;
row.finalFdRateErrHzPerSec = info.fdRateErrHzPerSec;
row.finalObj = localGetFieldOrDefault(finalEval, 'obj', localGetFieldOrDefault(estResult, 'fval', NaN));
row.objMinusBaseline = row.finalObj - baselineObj;
row.iterations = info.iterations;
row.exitflag = info.exitflag;
row.solveVariant = string(localGetFieldOrDefault(estResult, 'solveVariant', localGetFieldOrDefault(optimInfo, 'solveVariant', "")));
row.candidateCount = numel(localGetFieldOrDefault(optimInfo, 'candidateVariant', strings(0, 1)));
row.boundaryHit = false;
row.firstOrderOpt = localGetFieldOrDefault(optimInfo, 'firstorderopt', NaN);
[~, row.nonRefCoherenceFloor] = extractMfProbeCoherenceFloor(finalEval);
row.warningFlag = localGetFieldOrDefault(info, 'warningFlag', false);
row.evalOk = isfinite(row.finalObj);
row.message = "";
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
