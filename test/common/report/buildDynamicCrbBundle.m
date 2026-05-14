function crbBundle = buildDynamicCrbBundle(periodicFixture, pilotWave, carrierFreq, sampleRate, noiseVar, crbOptOverride)
%BUILDDYNAMICCRBBUNDLE Build one compact CRB bundle for a periodic fixture.
% The dev/perf scripts share the same CRB anchor construction so the entry
% points can stay focused on orchestration and display.

arguments
  periodicFixture (1,1) struct
  pilotWave
  carrierFreq (1,1) double
  sampleRate (1,1) double
  noiseVar (1,1) double
  crbOptOverride (1,1) struct = struct()
end

truth = periodicFixture.truth;
sceneRef = periodicFixture.sceneRef;
sceneRefOnly = periodicFixture.sceneSeqRefOnly.sceneCell{periodicFixture.sceneSeqRefOnly.refFrameIdx};

crbSfOpt = struct();
crbSfOpt.doaType = 'latlon';
[crbSfDoaRef, auxCrbSfDoaRef] = localTryBuildDoaOnlyCrb( ...
  sceneRefOnly, pilotWave, carrierFreq, sampleRate, truth.latlonTrueDeg, truth.fdRefTrueHz, 1, noiseVar, crbSfOpt);
[crbSfDoaMs, auxCrbSfDoaMs] = localTryBuildDoaOnlyCrb( ...
  sceneRef, pilotWave, carrierFreq, sampleRate, truth.latlonTrueDeg, truth.fdRefTrueHz, 1, noiseVar, crbSfOpt);
[crbSfRef, auxCrbSfRef] = localTryBuildStaticCrb( ...
  sceneRefOnly, pilotWave, carrierFreq, sampleRate, truth.latlonTrueDeg, truth.fdRefTrueHz, 1, noiseVar, crbSfOpt);
[crbSfMs, auxCrbSfMs] = localTryBuildStaticCrb( ...
  sceneRef, pilotWave, carrierFreq, sampleRate, truth.latlonTrueDeg, truth.fdRefTrueHz, 1, noiseVar, crbSfOpt);

pathGainRef = ones(periodicFixture.sceneSeqRefOnly.numSat, periodicFixture.sceneSeqRefOnly.numFrame);
noiseVarRef = noiseVar * ones(periodicFixture.sceneSeqRefOnly.numSat, periodicFixture.sceneSeqRefOnly.numFrame);
pathGainMs = ones(periodicFixture.sceneSeq.numSat, periodicFixture.sceneSeq.numFrame);
noiseVarMs = noiseVar * ones(periodicFixture.sceneSeq.numSat, periodicFixture.sceneSeq.numFrame);

crbMfKnownOpt = struct();
crbMfKnownOpt.doaType = 'latlon';
crbMfKnownOpt.phaseMode = 'continuous';
crbMfKnownOpt.fdRateMode = 'known';
crbMfKnownOpt.steeringMode = 'framewise';
crbMfKnownOpt.steeringRefFrameIdx = periodicFixture.sceneSeq.refFrameIdx;
crbMfKnownOpt = localMergeStruct(crbMfKnownOpt, crbOptOverride);
[crbMfRefKnown, auxCrbMfRefKnown] = localTryBuildDynamicCrb( ...
  periodicFixture.sceneSeqRefOnly, pilotWave, carrierFreq, sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, truth.fdRateFit, pathGainRef, noiseVarRef, crbMfKnownOpt);
[crbMfMsKnown, auxCrbMfMsKnown] = localTryBuildDynamicCrb( ...
  periodicFixture.sceneSeq, pilotWave, carrierFreq, sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, truth.fdRateFit, pathGainMs, noiseVarMs, crbMfKnownOpt);

crbMfStaticOpt = crbMfKnownOpt;
crbMfStaticOpt.fdRateMode = 'zero';
[crbMfRefStatic, auxCrbMfRefStatic] = localTryBuildDynamicCrb( ...
  periodicFixture.sceneSeqRefOnly, pilotWave, carrierFreq, sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, 0, pathGainRef, noiseVarRef, crbMfStaticOpt);
[crbMfMsStatic, auxCrbMfMsStatic] = localTryBuildDynamicCrb( ...
  periodicFixture.sceneSeq, pilotWave, carrierFreq, sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, 0, pathGainMs, noiseVarMs, crbMfStaticOpt);

crbMfUnknownOpt = crbMfKnownOpt;
crbMfUnknownOpt.fdRateMode = 'unknown';
[crbMfRefUnknown, auxCrbMfRefUnknown] = localTryBuildDynamicCrb( ...
  periodicFixture.sceneSeqRefOnly, pilotWave, carrierFreq, sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, truth.fdRateFit, pathGainRef, noiseVarRef, crbMfUnknownOpt);
[crbMfMsUnknown, auxCrbMfMsUnknown] = localTryBuildDynamicCrb( ...
  periodicFixture.sceneSeq, pilotWave, carrierFreq, sampleRate, ...
  truth.latlonTrueDeg, truth.fdRefFit, truth.fdRateFit, pathGainMs, noiseVarMs, crbMfUnknownOpt);

crbBundle = struct();
crbBundle.truth = truth;
crbBundle.crbSfDoaRef = crbSfDoaRef;
crbBundle.auxCrbSfDoaRef = auxCrbSfDoaRef;
crbBundle.crbSfDoaMs = crbSfDoaMs;
crbBundle.auxCrbSfDoaMs = auxCrbSfDoaMs;
crbBundle.crbSfRef = crbSfRef;
crbBundle.auxCrbSfRef = auxCrbSfRef;
crbBundle.crbSfMs = crbSfMs;
crbBundle.auxCrbSfMs = auxCrbSfMs;
crbBundle.crbMfRefStatic = crbMfRefStatic;
crbBundle.auxCrbMfRefStatic = auxCrbMfRefStatic;
crbBundle.crbMfMsStatic = crbMfMsStatic;
crbBundle.auxCrbMfMsStatic = auxCrbMfMsStatic;
crbBundle.crbMfRefKnown = crbMfRefKnown;
crbBundle.auxCrbMfRefKnown = auxCrbMfRefKnown;
crbBundle.crbMfMsKnown = crbMfMsKnown;
crbBundle.auxCrbMfMsKnown = auxCrbMfMsKnown;
crbBundle.crbMfRefUnknown = crbMfRefUnknown;
crbBundle.auxCrbMfRefUnknown = auxCrbMfRefUnknown;
crbBundle.crbMfMsUnknown = crbMfMsUnknown;
crbBundle.auxCrbMfMsUnknown = auxCrbMfMsUnknown;
end

function [crb, aux] = localTryBuildDoaOnlyCrb(scene, pilotWave, carrierFreq, sampleRate, doaParam, fdRefHz, pathGain, noiseVar, crbOpt)
try
  [crb, aux] = crbPilotSfDoaOnlyEffective(scene, pilotWave, carrierFreq, sampleRate, ...
    doaParam, fdRefHz, pathGain, noiseVar, crbOpt);
catch ME
  crb = nan(2, 2);
  aux = struct();
  aux.modelType = "sfDoaOnlyEffective";
  aux.errorMessage = string(ME.message);
end
end

function [crb, aux] = localTryBuildStaticCrb(scene, pilotWave, carrierFreq, sampleRate, doaParam, fdRefHz, numSource, noiseVar, crbOpt)
try
  [crb, aux] = crbPilotSfDoaDoppler(scene, pilotWave, carrierFreq, sampleRate, ...
    doaParam, fdRefHz, numSource, noiseVar, crbOpt);
catch ME
  crb = nan(3, 3);
  aux = struct();
  aux.fdRateMode = "zero";
  aux.phaseMode = "single-frame";
  aux.errorMessage = string(ME.message);
end
end

function [crb, aux] = localTryBuildDynamicCrb(sceneSeq, pilotWave, carrierFreq, sampleRate, doaParam, fdRefHz, fdRateHzPerSec, pathGain, noiseVar, crbOpt)
try
  [crb, aux] = crbPilotMfDoaDoppler(sceneSeq, pilotWave, carrierFreq, sampleRate, ...
    doaParam, fdRefHz, fdRateHzPerSec, pathGain, noiseVar, crbOpt);
catch ME
  crb = nan(3, 3);
  aux = struct();
  aux.fdRateMode = getDoaDopplerFieldOrDefault(crbOpt, 'fdRateMode', "unknown");
  aux.phaseMode = getDoaDopplerFieldOrDefault(crbOpt, 'phaseMode', "continuous");
  aux.errorMessage = string(ME.message);
end
end

function outStruct = localMergeStruct(baseStruct, overrideStruct)
%LOCALMERGESTRUCT Apply optional CRB overrides without changing defaults.

outStruct = baseStruct;
if nargin < 2 || isempty(overrideStruct)
  return;
end
fieldList = fieldnames(overrideStruct);
for iField = 1:numel(fieldList)
  outStruct.(fieldList{iField}) = overrideStruct.(fieldList{iField});
end
end
