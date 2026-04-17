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

if sceneSeqUse.numSat ~= 2
  error('buildDynamicFrameSubsetFixture:SubsetNumSatMismatch', ...
    'The subset fixture expects exactly two satellites.');
end

otherSatIdxLocalUse = 3 - refSatIdxLocalUse;
if isfield(sceneRefUse, 'satIdx') && ~isempty(sceneRefUse.satIdx)
  satIdxGlobalUse = reshape(sceneRefUse.satIdx, 1, []);
else
  error('buildDynamicFrameSubsetFixture:MissingSceneSatIdx', ...
    'sceneRef.satIdx is required in the frame subset fixture.');
end
otherSatIdxGlobalUse = satIdxGlobalUse(otherSatIdxLocalUse);

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

sceneOtherOnly = selectSatScene(sceneRefUse, otherSatIdxLocalUse);
rxSigOtherOnly = selectRxSigBySat(rxSigRefUse, otherSatIdxLocalUse, 'singleFrame');
viewOtherOnly = buildDoaDopplerEstView(sceneOtherOnly, rxSigOtherOnly, ...
  gridSize, searchRange, E);

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
fixture.otherSatIdxLocal = otherSatIdxLocalUse;
fixture.otherSatIdxGlobal = otherSatIdxGlobalUse;
fixture.refStateSource = string(refStateUse.source);
fixture.wavelen = wavelen;
fixture.debugTruthMs = struct();
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
