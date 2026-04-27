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
elemSpace = opt.elemSpacingWavelength * wavelen;
arrUpa = createUpa(opt.arraySize, elemSpace);
E = referenceEllipsoid('sphere');
usrLla = reshape(opt.usrLla, [], 1);
truthLatlon = usrLla(1:2, 1);
searchRange = [truthLatlon(1) - opt.searchMarginDeg, truthLatlon(1) + opt.searchMarginDeg; ...
               truthLatlon(2) - opt.searchMarginDeg, truthLatlon(2) + opt.searchMarginDeg];

utcRef = localNormalizeUtcRef(opt.utcRef);
utcVecMaster = utcRef + seconds(opt.masterOffsetIdx * opt.frameIntvlSec);

tle = tleread(resolveTestDataPath(char(opt.tleFileName), "tle"));
[~, satAccessRef] = findVisibleSatFromTle(utcRef, tle, usrLla);
if isempty(opt.selectedSatIdxGlobal)
  [selectedSatIdxGlobal, satPickAux] = pickVisibleSatByElevation( ...
    satAccessRef, opt.satPickCount, opt.satPickAnchorIdx, opt.satPickMode);
else
  selectedSatIdxGlobal = reshape(double(opt.selectedSatIdxGlobal), 1, []);
  satPickAux = struct('source', "explicit-selectedSatIdxGlobal");
end
if isempty(opt.refSatIdxGlobal)
  refSatIdxGlobal = selectedSatIdxGlobal(1);
else
  refSatIdxGlobal = double(opt.refSatIdxGlobal);
end
if ~any(selectedSatIdxGlobal == refSatIdxGlobal)
  error('buildDynamicDualSatEciContext:InvalidReferenceSat', ...
    'refSatIdxGlobal must be one of selectedSatIdxGlobal.');
end

sceneSeqMaster = genMultiFrameScene(utcVecMaster, tle, usrLla, selectedSatIdxGlobal, [], arrUpa, ...
  opt.satElevationMaskDeg(1), opt.satElevationMaskDeg(2), "satellite", refSatIdxGlobal, find(opt.masterOffsetIdx == 0, 1, 'first'));
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
if isfield(opt, 'subsetOffsetCell') && ~isempty(opt.subsetOffsetCell)
  subsetOffsetCell = opt.subsetOffsetCell;
  subsetLabelList = reshape(string(opt.subsetLabelList), [], 1);
  if isempty(subsetLabelList)
    subsetLabelList = "subset" + string((1:numel(subsetOffsetCell)).');
  end
else
  subsetOffsetCell = primarySubsetOffsetCell;
  subsetLabelList = primarySubsetLabelList;
end

context = struct();
context.frameIntvlSec = opt.frameIntvlSec;
context.periodicOffsetIdx = reshape(opt.periodicOffsetIdx, 1, []);
context.masterOffsetIdx = reshape(opt.masterOffsetIdx, 1, []);
context.sampleRate = opt.sampleRate;
context.symbolRate = opt.symbolRate;
context.numSym = opt.numSym;
context.carrierFreq = opt.carrierFreq;
context.baseSeed = opt.baseSeed;
context.tleFileName = string(opt.tleFileName);
context.arraySize = opt.arraySize;
context.elemSpacingWavelength = opt.elemSpacingWavelength;
context.satElevationMaskDeg = opt.satElevationMaskDeg;
context.satPickCount = opt.satPickCount;
context.satPickAnchorIdx = opt.satPickAnchorIdx;
context.satPickMode = string(opt.satPickMode);
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
opt.usrLla = [37.78; 36.59; 0];
opt.utcRef = datetime([2026, 03, 18, 17, 08, 0], ...
  'TimeZone', 'UTC', ...
  'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
opt.tleFileName = "starlink_pair_4154_1165_20260318_170800.tle";
opt.selectedSatIdxGlobal = [];
opt.refSatIdxGlobal = [];
opt.satPickCount = 2;
opt.satPickAnchorIdx = 1;
opt.satPickMode = "available";
opt.satElevationMaskDeg = [15, 55];
opt.arraySize = [4, 4];
opt.elemSpacingWavelength = 0.5;
opt.gridSize = [50, 50];
opt.searchMarginDeg = 5;
opt.fdRangeDefault = [-2e5, 2e5];
opt.fdRateRangeDefault = [-1e4, 0];
opt.numUsr = 1;
opt.numSubsetRandomTrial = 0;
opt.subsetOffsetCell = {};
opt.subsetLabelList = strings(0, 1);
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

function utcRef = localNormalizeUtcRef(utcInput)
if isa(utcInput, 'datetime')
  utcRef = utcInput;
elseif isnumeric(utcInput) && numel(utcInput) >= 6
  utcRef = datetime(reshape(utcInput(1:6), 1, []), ...
    'TimeZone', 'UTC', ...
    'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
elseif ischar(utcInput) || (isstring(utcInput) && isscalar(utcInput))
  utcRef = datetime(char(utcInput), ...
    'TimeZone', 'UTC', ...
    'Format', 'yyyy-MM-dd HH:mm:ss.SSSSSS');
else
  error('buildDynamicDualSatEciContext:InvalidUtcRef', ...
    'utcRef must be a datetime, date vector, or datetime string.');
end
if isempty(utcRef.TimeZone)
  utcRef.TimeZone = 'UTC';
end
utcRef.Format = 'yyyy-MM-dd HH:mm:ss.SSSSSS';
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
