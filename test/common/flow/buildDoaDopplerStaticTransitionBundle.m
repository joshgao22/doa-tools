function bundle = buildDoaDopplerStaticTransitionBundle(viewRefOnly, viewOtherOnly, viewMs, ...
  wavelen, pilotWave, carrierFreq, sampleRate, fdRange, truth, ...
  otherSatIdxGlobal, optVerbose, doaOnlyOpt, staticBaseOpt, weightSweepAlpha, staticMsHalfWidth, bundleMode)
%BUILDDOADOPPLERSTATICTRANSITIONBUNDLE Build the shared SF transition cases.
% This helper keeps the single-frame SS/MS DoA and static DoA-Doppler
% ladder identical between the static and dynamic dev scripts by default.
% Replay/scan callers may request bundleMode="ms-seed-only" to build only
% the SS/MS static seeds needed by MS dynamic init diagnostics. The static
% sat-weight sweep stays serial by default because it only evaluates a few
% small cases and parfor overhead usually dominates.

arguments
  viewRefOnly (1, 1) struct
  viewOtherOnly (1, 1) struct
  viewMs (1, 1) struct
  wavelen (1, 1) double
  pilotWave
  carrierFreq (1, 1) double
  sampleRate (1, 1) double
  fdRange (1, 2) double
  truth (1, 1) struct
  otherSatIdxGlobal (1, 1) double
  optVerbose (1, 1) logical = false
  doaOnlyOpt (1, 1) struct = struct()
  staticBaseOpt (1, 1) struct = struct()
  weightSweepAlpha (:, 1) double = [0; 0.25; 0.5; 1]
  staticMsHalfWidth (:, 1) double = [0.002; 0.002]
  bundleMode (1, 1) string = "full"
end

bundleMode = string(bundleMode);
if ~(bundleMode == "full" || bundleMode == "ms-seed-only")
  error('buildDoaDopplerStaticTransitionBundle:UnsupportedBundleMode', ...
    'Unsupported bundleMode "%s".', char(bundleMode));
end

bundle = struct();
bundle.doaOnlyOpt = doaOnlyOpt;
bundle.staticBaseOpt = staticBaseOpt;
bundle.weightSweepAlpha = weightSweepAlpha;
bundle.bundleMode = bundleMode;

bundle.caseRefDoa = localRunDoaOnlyCase("SS-SF-DoA", "single", ...
  viewRefOnly, wavelen, pilotWave, optVerbose, doaOnlyOpt);
bundle.caseMsDoa = localRunDoaOnlyCase("MS-SF-DoA", "multi", ...
  viewMs, wavelen, pilotWave, optVerbose, doaOnlyOpt);

bundle.staticRefOpt = staticBaseOpt;
bundle.staticRefOpt.initDoaParam = bundle.caseRefDoa.estResult.doaParamEst(:);

bundle.staticMsOpt = staticBaseOpt;
bundle.staticMsOpt.initDoaParam = bundle.caseMsDoa.estResult.doaParamEst(:);
bundle.staticMsOpt.initDoaHalfWidth = reshape(staticMsHalfWidth, [], 1);

bundle.caseStaticRefOnly = localRunStaticDoaDopplerCase("SS-SF-Static", "single", ...
  viewRefOnly, pilotWave, carrierFreq, sampleRate, fdRange, optVerbose, bundle.staticRefOpt);
bundle.caseStaticMs = localRunStaticDoaDopplerCase("MS-SF-Static", "multi", ...
  viewMs, pilotWave, carrierFreq, sampleRate, fdRange, optVerbose, bundle.staticMsOpt);

if bundleMode == "ms-seed-only"
  bundle.caseOtherDoa = buildDoaDopplerCaseResult( ...
    sprintf("SS-SF-DoA-otherSat(g%d)", otherSatIdxGlobal), "single", "single", ...
    "doa", "none", struct());
  bundle.staticOtherOpt = staticBaseOpt;
  bundle.staticOtherOpt.initDoaParam = NaN(2, 1);
  bundle.caseStaticRefAbl = bundle.caseStaticRefOnly;
  bundle.caseStaticRefAbl.displayName = sprintf("MS-SF-Static-refSat(g%d)", ...
    getDoaDopplerFieldOrDefault(truth, 'refSatIdxGlobal', NaN));
  bundle.caseStaticRefAbl.satMode = "multi";
  bundle.caseStaticOtherOnly = buildDoaDopplerCaseResult( ...
    sprintf("MS-SF-Static-otherSat(g%d)", otherSatIdxGlobal), "multi", "single", ...
    "doa-doppler", "static", struct());
  bundle.weightCase = repmat(buildDoaDopplerCaseResult("", "multi", "single", ...
    "doa-doppler", "static", struct()), 1, 0);
  bundle.bestStaticMsCase = bundle.caseStaticMs;
  return;
end

bundle.caseOtherDoa = localRunDoaOnlyCase(sprintf("SS-SF-DoA-otherSat(g%d)", otherSatIdxGlobal), ...
  "single", viewOtherOnly, wavelen, pilotWave, optVerbose, doaOnlyOpt);
bundle.staticOtherOpt = staticBaseOpt;
bundle.staticOtherOpt.initDoaParam = bundle.caseOtherDoa.estResult.doaParamEst(:);
bundle.caseStaticRefAbl = bundle.caseStaticRefOnly;
bundle.caseStaticRefAbl.displayName = sprintf("MS-SF-Static-refSat(g%d)", ...
  getDoaDopplerFieldOrDefault(truth, 'refSatIdxGlobal', NaN));
bundle.caseStaticRefAbl.satMode = "multi";

bundle.caseStaticOtherOnly = localRunStaticDoaDopplerCase( ...
  sprintf("MS-SF-Static-otherSat(g%d)", otherSatIdxGlobal), "multi", ...
  viewOtherOnly, pilotWave, carrierFreq, sampleRate, fdRange, optVerbose, bundle.staticOtherOpt);

bundle.weightCase = repmat(buildDoaDopplerCaseResult("", "multi", "single", ...
  "doa-doppler", "static", struct()), 1, numel(weightSweepAlpha));
useParforWeight = localShouldUseParfor(staticBaseOpt, numel(weightSweepAlpha), 'minWeightSweepForParfor');
if useParforWeight
  weightCaseCell = cell(1, numel(weightSweepAlpha));
  parfor iCase = 1:numel(weightSweepAlpha)
    alpha = weightSweepAlpha(iCase);
    currentOpt = bundle.staticMsOpt;
    currentOpt.satWeight = [1; alpha];
    weightCaseCell{iCase} = localRunStaticDoaDopplerCase( ...
      sprintf('MS-SF-Static-W%.2f', alpha), "multi", viewMs, ...
      pilotWave, carrierFreq, sampleRate, fdRange, optVerbose, currentOpt);
  end
  bundle.weightCase = [weightCaseCell{:}];
else
  for iCase = 1:numel(weightSweepAlpha)
    alpha = weightSweepAlpha(iCase);
    currentOpt = bundle.staticMsOpt;
    currentOpt.satWeight = [1; alpha];
    bundle.weightCase(iCase) = localRunStaticDoaDopplerCase( ...
      sprintf('MS-SF-Static-W%.2f', alpha), "multi", viewMs, ...
      pilotWave, carrierFreq, sampleRate, fdRange, optVerbose, currentOpt);
  end
end

bundle.bestStaticMsCase = localSelectBestStaticSeedCase(bundle.caseStaticMs, bundle.weightCase, truth);
end


function caseInfo = localRunDoaOnlyCase(displayName, satMode, view, wavelen, ...
  pilotWave, verbose, modelOpt)
%LOCALRUNDOAONLYCASE Run one DoA-only estimator case using a private wrapper.

[estRaw, ~, ~] = estimatorDoaMlePilotOpt( ...
  view.sceneRef.array, wavelen, view.rxSigSf, pilotWave, ...
  view.doaGrid, [], verbose, modelOpt);

estResult = localWrapDoaOnlyResult(estRaw, view);
caseInfo = buildDoaDopplerCaseResult(displayName, satMode, "single", ...
  "doa", "none", estResult);
end

function estOut = localWrapDoaOnlyResult(estIn, view)
%LOCALWRAPDOAONLYRESULT Convert DoA-only output to the common summary format.

estOut = estIn;
estOut.modelType = 'doa-only';
estOut.fdRefEst = NaN;
estOut.fdRateEst = NaN;
if ~isfield(estOut, 'aux') || ~isstruct(estOut.aux)
  estOut.aux = struct();
end
if isfield(estIn, 'latlonEst') && ~isfield(estOut.aux, 'latlonEst')
  estOut.aux.latlonEst = estIn.latlonEst;
end

timeOffsetSec = 0;
if isfield(view, 'sceneSeq') && ~isempty(view.sceneSeq)
  timeOffsetSec = view.sceneSeq.timeOffsetSec(view.sceneSeq.refFrameIdx);
end
estOut.timeOffsetSec = timeOffsetSec;
end

function caseInfo = localRunStaticDoaDopplerCase(displayName, satMode, view, ...
  pilotWave, carrierFreq, sampleRate, fdRange, verbose, modelOpt)
%LOCALRUNSTATICDOADOPPLERCASE Run one single-frame static DoA-Doppler case.

[estResult, ~, ~] = estimatorDoaDopplerMlePilotSfOpt( ...
  view.sceneRef, view.rxSigSf, pilotWave, carrierFreq, sampleRate, ...
  view.doaGrid, fdRange, 1, [], verbose, modelOpt);

caseInfo = buildDoaDopplerCaseResult(displayName, satMode, "single", ...
  "doa-doppler", "static", estResult);
end

function bestCase = localSelectBestStaticSeedCase(baseCase, weightCase, truth)
%LOCALSELECTBESTSTATICSEEDCASE Choose the most reliable static multi-sat seed.
% Dynamic MF refinement should start from the best resolved static MS case
% rather than always from the equal-weight branch.

candidateCase = [baseCase, weightCase];
bestCase = baseCase;
bestScore = [inf, inf];

truthLatlon = reshape(getDoaDopplerFieldOrDefault(truth, 'latlonTrueDeg', nan(2, 1)), [], 1);
truthFdRef = localResolveTruthFdRef(truth);

for iCase = 1:numel(candidateCase)
  currentCase = candidateCase(iCase);
  if ~localCaseResolved(currentCase)
    continue;
  end

  estResult = currentCase.estResult;
  angleErrDeg = calcLatlonAngleError(estResult.doaParamEst(:), truthLatlon);
  fdErrHz = abs(estResult.fdRefEst - truthFdRef);
  scoreNow = [angleErrDeg, fdErrHz];
  if localLexicoLess(scoreNow, bestScore)
    bestScore = scoreNow;
    bestCase = currentCase;
  end
end
end

function truthFdRef = localResolveTruthFdRef(truth)
%LOCALRESOLVETRUTHFDREF Resolve reference Doppler truth with fallback.

truthFdRef = getDoaDopplerFieldOrDefault(truth, 'fdRefTrueHz', NaN);
if ~isfinite(truthFdRef)
  truthFdRef = getDoaDopplerFieldOrDefault(truth, 'fdRefFit', NaN);
end
end

function isResolved = localCaseResolved(caseInfo)
%LOCALCASERESOLVED Return true when one case contains a usable estimate.

isResolved = false;
if ~isstruct(caseInfo) || ~isfield(caseInfo, 'estResult') || isempty(caseInfo.estResult)
  return;
end

estResult = caseInfo.estResult;
if isfield(estResult, 'isResolved') && ~isempty(estResult.isResolved)
  isResolved = logical(estResult.isResolved);
else
  isResolved = true;
end
if ~isResolved
  return;
end

if ~isfield(estResult, 'doaParamEst') || numel(estResult.doaParamEst) ~= 2 || ...
    any(~isfinite(estResult.doaParamEst(:)))
  isResolved = false;
  return;
end
if ~isfield(estResult, 'fdRefEst') || ~isscalar(estResult.fdRefEst) || ~isfinite(estResult.fdRefEst)
  isResolved = false;
end
end

function isLess = localLexicoLess(scoreNow, scoreBest)
%LOCALLEXICOLESS Lexicographic compare for [angleErr, fdErr].

if isempty(scoreBest) || any(~isfinite(scoreBest))
  isLess = all(isfinite(scoreNow));
  return;
end
if any(~isfinite(scoreNow))
  isLess = false;
  return;
end
if scoreNow(1) < scoreBest(1) - 1e-12
  isLess = true;
  return;
end
if scoreNow(1) > scoreBest(1) + 1e-12
  isLess = false;
  return;
end
isLess = scoreNow(2) < scoreBest(2);
end

function useParfor = localShouldUseParfor(optStruct, numCase, minFieldName)
%LOCALSHOULDUSEPARFOR Decide whether the current loop should use parfor.

useParfor = logical(localGetFieldOrDefault(optStruct, 'enableWeightSweepParfor', ...
  localGetFieldOrDefault(optStruct, 'enableParfor', false)));
if ~useParfor
  return;
end
minCase = localGetFieldOrDefault(optStruct, minFieldName, 8);
useParfor = useParfor && (numCase >= minCase) && localCanUseParfor();
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

function fieldValue = localGetFieldOrDefault(dataStruct, fieldName, defaultValue)
%LOCALGETFIELDORDEFAULT Read one field or property with a default value.

fieldValue = defaultValue;
if nargin < 3
  defaultValue = [];
  fieldValue = defaultValue;
end

if isempty(dataStruct)
  return;
end

if isstruct(dataStruct)
  if isfield(dataStruct, fieldName)
    fieldValue = dataStruct.(fieldName);
  end
  return;
end

if isobject(dataStruct)
  if isprop(dataStruct, fieldName)
    fieldValue = dataStruct.(fieldName);
  end
  return;
end

try
  fieldValue = dataStruct.(fieldName);
catch
end
end
