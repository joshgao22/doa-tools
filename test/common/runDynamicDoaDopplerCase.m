function caseInfo = runDynamicDoaDopplerCase(displayName, satMode, view, ...
  truth, pilotWave, carrierFreq, sampleRate, fdRange, fdRateRange, ...
  verbose, dynBaseOpt, isKnownRate, debugTruth, initParamOverride)
%RUNDYNAMICDOADOPPLERCASE Run one multi-frame continuous-phase case.

if nargin < 13 || isempty(debugTruth)
  debugTruth = [];
end
if nargin < 14
  initParamOverride = [];
end

dynOpt = dynBaseOpt;
if isKnownRate
  dynOpt.fdRateMode = 'known';
  dynOpt.fdRateKnown = truth.fdRateFit;
  fdRateRangeUse = [];
  dynamicMode = "cp-known";
else
  dynOpt.fdRateMode = 'unknown';
  fdRateRangeUse = fdRateRange;
  dynamicMode = "cp-unknown";
end

if ~isempty(debugTruth)
  dynOpt.debugTruth = debugTruth;
end

[estResult, ~, ~] = estimatorDoaDopplerMlePilotMfOpt( ...
  view.sceneSeq, view.rxSigMf, pilotWave, carrierFreq, sampleRate, ...
  view.doaGrid, fdRange, fdRateRangeUse, initParamOverride, verbose, dynOpt);

caseInfo = buildDoaDopplerCaseResult(displayName, satMode, "multi", ...
  "doa-doppler", dynamicMode, estResult);
end
