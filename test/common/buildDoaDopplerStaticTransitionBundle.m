function bundle = buildDoaDopplerStaticTransitionBundle(viewRefOnly, viewOtherOnly, viewMs, ...
  wavelen, pilotWave, carrierFreq, sampleRate, fdRange, truth, ...
  otherSatIdxGlobal, optVerbose, doaOnlyOpt, staticBaseOpt, weightSweepAlpha, staticMsHalfWidth)
%BUILDDOADOPPLERSTATICTRANSITIONBUNDLE Build the shared SF transition cases.
% This helper keeps the single-frame SS/MS DoA and static DoA-Doppler
% ladder identical between the static and dynamic dev scripts.

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
  staticMsHalfWidth (:, 1) double = [0.01; 0.01]
end

bundle = struct();
bundle.doaOnlyOpt = doaOnlyOpt;
bundle.staticBaseOpt = staticBaseOpt;
bundle.weightSweepAlpha = weightSweepAlpha;

bundle.caseRefDoa = runDoaOnlyCase("SS-SF-DoA", "single", ...
  viewRefOnly, wavelen, pilotWave, optVerbose, doaOnlyOpt);
bundle.caseMsDoa = runDoaOnlyCase("MS-SF-DoA", "multi", ...
  viewMs, wavelen, pilotWave, optVerbose, doaOnlyOpt);
bundle.caseOtherDoa = runDoaOnlyCase(sprintf("SS-SF-DoA-otherSat(g%d)", otherSatIdxGlobal), ...
  "single", viewOtherOnly, wavelen, pilotWave, optVerbose, doaOnlyOpt);

bundle.staticRefOpt = staticBaseOpt;
bundle.staticRefOpt.initDoaParam = bundle.caseRefDoa.estResult.doaParamEst(:);

bundle.staticOtherOpt = staticBaseOpt;
bundle.staticOtherOpt.initDoaParam = bundle.caseOtherDoa.estResult.doaParamEst(:);

bundle.staticMsOpt = staticBaseOpt;
bundle.staticMsOpt.initDoaParam = bundle.caseMsDoa.estResult.doaParamEst(:);
bundle.staticMsOpt.initDoaHalfWidth = reshape(staticMsHalfWidth, [], 1);

bundle.caseStaticRefOnly = runStaticDoaDopplerCase("SS-SF-Static", "single", ...
  viewRefOnly, pilotWave, carrierFreq, sampleRate, fdRange, optVerbose, bundle.staticRefOpt);
bundle.caseStaticRefAbl = bundle.caseStaticRefOnly;
bundle.caseStaticRefAbl.displayName = sprintf("MS-SF-Static-refSat(g%d)", ...
  getDoaDopplerFieldOrDefault(truth, 'refSatIdxGlobal', NaN));
bundle.caseStaticRefAbl.satMode = "multi";

bundle.caseStaticOtherOnly = runStaticDoaDopplerCase( ...
  sprintf("MS-SF-Static-otherSat(g%d)", otherSatIdxGlobal), "multi", ...
  viewOtherOnly, pilotWave, carrierFreq, sampleRate, fdRange, optVerbose, bundle.staticOtherOpt);

bundle.caseStaticMs = runStaticDoaDopplerCase("MS-SF-Static", "multi", ...
  viewMs, pilotWave, carrierFreq, sampleRate, fdRange, optVerbose, bundle.staticMsOpt);

bundle.weightCase = repmat(buildDoaDopplerCaseResult("", "multi", "single", ...
  "doa-doppler", "static", struct()), 1, numel(weightSweepAlpha));
for iCase = 1:numel(weightSweepAlpha)
  alpha = weightSweepAlpha(iCase);
  currentOpt = bundle.staticMsOpt;
  currentOpt.satWeight = [1; alpha];
  bundle.weightCase(iCase) = runStaticDoaDopplerCase( ...
    sprintf('MS-SF-Static-W%.2f', alpha), "multi", viewMs, ...
    pilotWave, carrierFreq, sampleRate, fdRange, optVerbose, currentOpt);
end

bundle.bestStaticMsCase = selectBestStaticSeedCase(bundle.caseStaticMs, bundle.weightCase, truth);
end
