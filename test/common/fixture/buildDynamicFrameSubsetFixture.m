function fixture = buildDynamicFrameSubsetFixture(sceneSeqMaster, linkParamCellMaster, rxSigCellMaster, ...
  masterOffsetIdx, subsetOffsetIdx, gridSize, searchRange, E, wavelen, sampleRate, ...
  fdRangeDefault, fdRateRangeDefault)
%BUILDDYNAMICFRAMESUBSETFIXTURE Build one reusable dynamic frame-subset fixture.

frameIdx = localResolveFrameIdxFromOffset(masterOffsetIdx, subsetOffsetIdx);
sceneSeqUse = selectFrameSceneSeq(sceneSeqMaster, frameIdx);
linkParamCellUse = reshape(linkParamCellMaster(frameIdx), 1, []);
rxSigCellUse = selectFrameRxSig(rxSigCellMaster, frameIdx);
sceneRefUse = sceneSeqUse.sceneCell{sceneSeqUse.refFrameIdx};
[refStateUse, refSatIdxLocalUse] = resolveReferenceSatState( ...
  sceneRefUse, sceneRefUse.satPosEci, sceneRefUse.satVelEci);

if isfield(sceneRefUse, 'satIdx') && ~isempty(sceneRefUse.satIdx)
  satIdxGlobalUse = reshape(sceneRefUse.satIdx, 1, []);
else
  error('buildDynamicFrameSubsetFixture:MissingSceneSatIdx', ...
    'sceneRef.satIdx is required in the frame subset fixture.');
end
nonRefSatIdxLocalUse = setdiff(1:sceneSeqUse.numSat, refSatIdxLocalUse, 'stable');
nonRefSatIdxGlobalUse = satIdxGlobalUse(nonRefSatIdxLocalUse);
otherSatIdxLocalUse = [];
otherSatIdxGlobalUse = [];
if sceneSeqUse.numSat == 2
  otherSatIdxLocalUse = nonRefSatIdxLocalUse;
  otherSatIdxGlobalUse = nonRefSatIdxGlobalUse;
end

truthUse = buildDynTruthFromLinkParam(linkParamCellUse, sceneSeqUse.timeOffsetSec, sampleRate, 1);
truthUse.utcRef = sceneSeqUse.utcRef;
truthUse.latlonTrueDeg = reshape(sceneRefUse.usrLla(1:2, 1), [], 1);
truthUse.refSatIdxGlobal = satIdxGlobalUse(refSatIdxLocalUse);
truthUse.refSatIdxLocal = refSatIdxLocalUse;
truthUse.selectedSatIdxGlobal = satIdxGlobalUse(:).';
truthUse.usrElevationDeg = reshape(sceneRefUse.accessInfo.usrElevationDeg(:, 1), 1, []);
truthUse.refWeight = sceneRefUse.ref.weight(:);
truthUse.refStateSource = string(refStateUse.source);
truthUse.fdRefTrueHz = truthUse.fdRefSeries(sceneSeqUse.refFrameIdx);
truthUse.fdRateTrueHzPerSec = truthUse.fdRateFit;
truthUse.fdSatTrueHz = reshape(truthUse.fdSatSeries(:, sceneSeqUse.refFrameIdx), [], 1);
truthUse.deltaFdTrueHz = reshape(truthUse.deltaFdSeries(:, sceneSeqUse.refFrameIdx), [], 1);
truthUse.localDoaRef = reshape(sceneRefUse.localDoa(:, refSatIdxLocalUse), 2, 1);

fdRangeUse = expandRangeToTruth(fdRangeDefault, [truthUse.fdRefFit; truthUse.fdSatTrueHz(:)], 0.1, 2e4);
fdRateTruthCand = [truthUse.fdRateFit; truthUse.fdRateFit + ...
  reshape(getDoaDopplerFieldOrDefault(truthUse, 'deltaFdRate', []), [], 1)];
fdRateRangeUse = expandRangeToTruth(fdRateRangeDefault, fdRateTruthCand, 0.1, 5e2);

rxSigRefUse = rxSigCellUse{sceneSeqUse.refFrameIdx};
sceneSeqRefOnly = selectSatSceneSeq(sceneSeqUse, refSatIdxLocalUse);
sceneRefOnly = sceneSeqRefOnly.sceneCell{sceneSeqRefOnly.refFrameIdx};
rxSigMfRefOnly = selectRxSigBySat(rxSigCellUse, refSatIdxLocalUse, 'multiFrame');
viewRefOnly = buildDoaDopplerEstView(sceneRefOnly, rxSigMfRefOnly, ...
  gridSize, searchRange, E, struct('sceneSeq', sceneSeqRefOnly));

sceneOtherOnly = struct();
viewOtherOnly = struct();
if sceneSeqUse.numSat == 2
  sceneOtherOnly = selectSatScene(sceneRefUse, otherSatIdxLocalUse);
  rxSigOtherOnly = selectRxSigBySat(rxSigRefUse, otherSatIdxLocalUse, 'singleFrame');
  viewOtherOnly = buildDoaDopplerEstView(sceneOtherOnly, rxSigOtherOnly, ...
    gridSize, searchRange, E);
end

viewMsUse = buildDoaDopplerEstView(sceneRefUse, rxSigCellUse, ...
  gridSize, searchRange, E, struct('sceneSeq', sceneSeqUse));

fixture = struct();
fixture.sceneSeq = sceneSeqUse;
fixture.sceneRef = sceneRefUse;
fixture.linkParamCell = linkParamCellUse;
fixture.rxSigCell = rxSigCellUse;
fixture.rxSigRef = rxSigRefUse;
fixture.truth = truthUse;
fixture.fdRange = fdRangeUse;
fixture.fdRateRange = fdRateRangeUse;
fixture.viewRefOnly = viewRefOnly;
fixture.viewOtherOnly = viewOtherOnly;
fixture.viewMs = viewMsUse;
fixture.sceneSeqRefOnly = sceneSeqRefOnly;
fixture.sceneRefOnly = sceneRefOnly;
fixture.sceneOtherOnly = sceneOtherOnly;
fixture.refSatIdxLocal = refSatIdxLocalUse;
fixture.nonRefSatIdxLocal = reshape(nonRefSatIdxLocalUse, 1, []);
fixture.nonRefSatIdxGlobal = reshape(nonRefSatIdxGlobalUse, 1, []);
fixture.otherSatIdxLocal = otherSatIdxLocalUse;
fixture.otherSatIdxGlobal = otherSatIdxGlobalUse;
fixture.refStateSource = string(refStateUse.source);
fixture.wavelen = wavelen;
fixture.debugTruthMs = struct();
localAssertReferenceDopplerInvariant(truthUse);
end

function frameIdx = localResolveFrameIdxFromOffset(masterOffsetIdx, subsetOffsetIdx)
%LOCALRESOLVEFRAMEIDXFROMOFFSET Map frame-offset labels to frame indices.

[isFound, frameIdx] = ismember(subsetOffsetIdx(:).', masterOffsetIdx(:).');
if ~all(isFound)
  missingOffset = subsetOffsetIdx(~isFound);
  error('buildDynamicFrameSubsetFixture:MissingSubsetOffset', ...
    'The requested subset offsets are not present in the master window: %s', ...
    mat2str(missingOffset));
end
frameIdx = reshape(frameIdx, 1, []);
end

function localAssertReferenceDopplerInvariant(truthUse)
%LOCALASSERTREFERENCEDOPPLERINVARIANT Check reference-Doppler bookkeeping.

refSatIdx = double(truthUse.refSatIdxLocal);
fdRef = double(truthUse.fdRefTrueHz);
fdSat = reshape(double(truthUse.fdSatTrueHz), [], 1);
deltaFd = reshape(double(truthUse.deltaFdTrueHz), [], 1);
if ~(isscalar(refSatIdx) && isfinite(refSatIdx) && refSatIdx >= 1 && refSatIdx <= numel(fdSat))
  error('buildDynamicFrameSubsetFixture:InvalidReferenceSat', ...
    'truth.refSatIdxLocal is inconsistent with fdSatTrueHz.');
end
if abs(deltaFd(refSatIdx)) > 1e-6 || abs(fdSat(refSatIdx) - fdRef) > 1e-6 || ...
    max(abs(fdSat - (fdRef + deltaFd))) > 1e-6
  error('buildDynamicFrameSubsetFixture:ReferenceDopplerInvariantFailed', ...
    ['Reference-Doppler truth must satisfy deltaFd(ref)=0, ', ...
    'fdSat(ref)=fdRef, and fdSat=fdRef+deltaFd.']);
end
end
