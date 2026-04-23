function context = buildDynamicDualSatEciContext(varargin)
%BUILDDYNAMICDUALSATECICONTEXT Build one reusable dual-satellite dynamic context.
% This helper centralizes the fixed scene / waveform setup used by the
% simplified dynamic dev script and the DoA-profile regression so both
% entry points share the same scenario semantics.

opt = localParseContextOpt(varargin{:});

if opt.numUsr ~= 1
  error('buildDynamicDualSatEciContext:OnlySingleUserSupported', ...
    'The current simplified dynamic context only supports one user.');
end
if ~any(opt.masterOffsetIdx == 0) || ~any(opt.periodicOffsetIdx == 0)
  error('buildDynamicDualSatEciContext:MissingReferenceFrame', ...
    'Both masterOffsetIdx and periodicOffsetIdx must contain one zero offset.');
end

rng(opt.baseSeed);

wavelen = 299792458 / opt.carrierFreq;
elemSpace = wavelen / 2;
arrUpa = createUpa([4, 4], elemSpace);
E = referenceEllipsoid('sphere');
usrLla = [37.78; 36.59; 0];
truthLatlon = usrLla(1:2, 1);
searchRange = [truthLatlon(1) - opt.searchMarginDeg, truthLatlon(1) + opt.searchMarginDeg; ...
               truthLatlon(2) - opt.searchMarginDeg, truthLatlon(2) + opt.searchMarginDeg];

utcRef = datetime([2026, 03, 18, 17, 08, 0], ...
  'TimeZone', 'UTC', ...
  'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
utcVecMaster = utcRef + seconds(opt.masterOffsetIdx * opt.frameIntvlSec);

tle = tleread(resolveTestDataPath("starlink_pair_4154_1165_20260318_170800.tle", "tle"));
[~, satAccessRef] = findVisibleSatFromTle(utcRef, tle, usrLla);
[selectedSatIdxGlobal, satPickAux] = pickVisibleSatByElevation(satAccessRef, 2, 1, "available");
refSatIdxGlobal = selectedSatIdxGlobal(1);

sceneSeqMaster = genMultiFrameScene(utcVecMaster, tle, usrLla, selectedSatIdxGlobal, [], arrUpa, ...
  15, 55, "satellite", refSatIdxGlobal, find(opt.masterOffsetIdx == 0, 1, 'first'));
linkParamCellMaster = cell(1, sceneSeqMaster.numFrame);
for iFrame = 1:sceneSeqMaster.numFrame
  linkParamCellMaster{iFrame} = getLinkParam(sceneSeqMaster.sceneCell{iFrame}, wavelen);
end

refFrameIdxMaster = sceneSeqMaster.refFrameIdx;
sceneRefMaster = sceneSeqMaster.sceneCell{refFrameIdxMaster};
[~, refSatIdxLocal] = resolveReferenceSatState(sceneRefMaster, sceneRefMaster.satPosEci, sceneRefMaster.satVelEci);
if sceneSeqMaster.numSat ~= 2
  error('buildDynamicDualSatEciContext:InvalidNumSat', ...
    'The simplified dual-satellite context expects exactly two satellites.');
end
otherSatIdxLocal = 3 - refSatIdxLocal;
otherSatIdxGlobal = selectedSatIdxGlobal(otherSatIdxLocal);

[pilotSym, ~] = genPilotSymbol(opt.numUsr, opt.numSym, 'pn', 1);
pulseOpt = struct();
pulseOpt.rolloff = 0.25;
pulseOpt.span = 8;
osf = opt.sampleRate / opt.symbolRate;
[pilotWave, waveInfo] = genPilotWaveform(pilotSym, opt.symbolRate, osf, 'rrc', pulseOpt);

simOpt = struct();
simOpt.spatial.model = 'dynamic';
simOpt.spatial.refFrameIdx = sceneSeqMaster.refFrameIdx;
simOpt.phase.timeModel = 'global';
simOpt.phase.frameModel = 'shared';
simOpt.phase.sharedPhase = 2 * pi * rand(sceneSeqMaster.numSat, sceneSeqMaster.numUser);
simOpt.wave.delayModel = 'phaseOnly';
simOpt.wave.timeRef = 'zero';
simOpt.wave.carrierPhaseModel = 'none';
simOpt.precomp.linkParamCell = linkParamCellMaster;

[primarySubsetOffsetCell, primarySubsetLabelList] = getDynamicCuratedSubsetBank();
subsetOffsetCell = primarySubsetOffsetCell;
subsetLabelList = primarySubsetLabelList;

context = struct();
context.frameIntvlSec = opt.frameIntvlSec;
context.periodicOffsetIdx = reshape(opt.periodicOffsetIdx, 1, []);
context.masterOffsetIdx = reshape(opt.masterOffsetIdx, 1, []);
context.sampleRate = opt.sampleRate;
context.symbolRate = opt.symbolRate;
context.numSym = opt.numSym;
context.carrierFreq = opt.carrierFreq;
context.baseSeed = opt.baseSeed;
context.gridSize = opt.gridSize;
context.searchRange = searchRange;
context.fdRangeDefault = opt.fdRangeDefault;
context.fdRateRangeDefault = opt.fdRateRangeDefault;
context.numUsr = opt.numUsr;
context.numSubsetRandomTrial = opt.numSubsetRandomTrial;
context.parallelOpt = opt.parallelOpt;
context.wavelen = wavelen;
context.elemSpace = elemSpace;
context.arrUpa = arrUpa;
context.E = E;
context.usrLla = usrLla;
context.truthLatlon = truthLatlon;
context.utcRef = utcRef;
context.utcVecMaster = utcVecMaster;
context.tle = tle;
context.selectedSatIdxGlobal = reshape(selectedSatIdxGlobal, 1, []);
context.refSatIdxGlobal = refSatIdxGlobal;
context.refSatIdxLocal = refSatIdxLocal;
context.otherSatIdxLocal = otherSatIdxLocal;
context.otherSatIdxGlobal = otherSatIdxGlobal;
context.satPickAux = satPickAux;
context.sceneSeqMaster = sceneSeqMaster;
context.linkParamCellMaster = linkParamCellMaster;
context.sceneRefMaster = sceneRefMaster;
context.pilotWave = pilotWave;
context.waveInfo = waveInfo;
context.simOpt = simOpt;
context.subsetOffsetCell = subsetOffsetCell;
context.subsetLabelList = reshape(string(subsetLabelList), [], 1);
end

function opt = localParseContextOpt(varargin)
opt = struct();
opt.frameIntvlSec = 1 / 750;
opt.periodicOffsetIdx = -4:5;
opt.masterOffsetIdx = -9:10;
opt.sampleRate = 512e6;
opt.symbolRate = 128e6;
opt.numSym = 512;
opt.carrierFreq = 11.7e9;
opt.baseSeed = 253;
opt.gridSize = [50, 50];
opt.searchMarginDeg = 5;
opt.fdRangeDefault = [-2e5, 2e5];
opt.fdRateRangeDefault = [-1e4, 0];
opt.numUsr = 1;
opt.numSubsetRandomTrial = 0;
opt.parallelOpt = struct('enableSubsetEvalParfor', true, 'minSubsetEvalParfor', 4);

if nargin == 1 && isstruct(varargin{1})
  override = varargin{1};
  opt = localMergeStruct(opt, override);
elseif mod(nargin, 2) == 0
  for iArg = 1:2:nargin
    name = varargin{iArg};
    if ~(ischar(name) || (isstring(name) && isscalar(name)))
      error('buildDynamicDualSatEciContext:InvalidNameValue', ...
        'Name-value inputs must use character or string names.');
    end
    opt.(char(name)) = varargin{iArg + 1};
  end
elseif nargin ~= 0
  error('buildDynamicDualSatEciContext:InvalidInput', ...
    'Use either one override struct or name-value pairs.');
end
end

function outStruct = localMergeStruct(baseStruct, overrideStruct)
outStruct = baseStruct;
if nargin < 2 || isempty(overrideStruct)
  return;
end
fieldList = fieldnames(overrideStruct);
for iField = 1:numel(fieldList)
  outStruct.(fieldList{iField}) = overrideStruct.(fieldList{iField});
end
end
