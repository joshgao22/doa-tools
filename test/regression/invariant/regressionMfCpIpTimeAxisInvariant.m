function regressionMfCpIpTimeAxisInvariant(varargin)
% Regression check for MF CP/IP phase-mode time-axis semantics.
% Continuous phase must tie one satellite phase across frames using the
% absolute time axis, while independent phase must use frame-local time and
% profile one phase per valid frame.
opt = parseRegressionCaseOpt(varargin{:});
verbose = opt.verbose;

fprintf('Running regressionMfCpIpTimeAxisInvariant ...\n');

[modelCp, hyp] = localBuildToyModel('continuous');
[~, profCp] = evalDoaDopplerDynProfileLike(modelCp, hyp);

modelIp = modelCp;
modelIp.phaseMode = 'independent';
[~, profIp] = evalDoaDopplerDynProfileLike(modelIp, hyp);

cpFramePhase = profCp.framePhaseEst(1, :);
ipFramePhase = profIp.framePhaseEst(1, :);
cpSpread = localCircularSpread(cpFramePhase);
ipSpread = localCircularSpread(ipFramePhase);

if ~isfinite(profCp.phaseSatEst(1))
  error('regressionMfCpIpTimeAxisInvariant:MissingCommonPhase', ...
    'Continuous phase mode must estimate one finite common satellite phase.');
end
if cpSpread > 1e-8
  error('regressionMfCpIpTimeAxisInvariant:CpFramePhaseNotTied', ...
    'Continuous phase mode must copy the same profiled phase to every frame.');
end
if isfinite(profIp.phaseSatEst(1))
  error('regressionMfCpIpTimeAxisInvariant:IndependentCommonPhase', ...
    'Independent phase mode must not report one shared satellite phase.');
end
if ipSpread < 0.05
  error('regressionMfCpIpTimeAxisInvariant:IndependentFramePhaseCollapsed', ...
    'Independent phase mode should expose different frame-local phase offsets in this toy case.');
end

fprintf('  CP frame phase spread : %.3e rad\n', cpSpread);
fprintf('  IP frame phase spread : %.3e rad\n', ipSpread);
fprintf('PASS: regressionMfCpIpTimeAxisInvariant\n');
end

function [model, hyp] = localBuildToyModel(phaseMode)
%LOCALBUILDTOYMODEL Build one tiny noiseless one-satellite MF case.

array = createUpa([1, 1], 0.5);
wavelength = 1;
numFrame = 3;
numSample = 8;
fs = 1000;
fdHz = 31;
fdRateHzPerSec = 240;
timeOffsetSec = [0, 0.051, 0.093];
pilot = ones(numSample, 1);
localDoa = [0; 0];
aVec = steeringMatrix(array, wavelength, localDoa);
aVec = aVec(:);

rxSig = cell(1, numFrame);
pilotPadCell = cell(1, numFrame);
timeSecCell = cell(1, numFrame);
for iFrame = 1:numFrame
  timeLoc = (0:numSample-1).' / fs;
  timeAbs = timeOffsetSec(iFrame) + timeLoc;
  phase = 2 * pi * (fdHz * timeAbs + 0.5 * fdRateHzPerSec * timeAbs.^2);
  qCol = pilot .* exp(1j * phase);
  rxSig{iFrame} = {aVec * qCol.'};
  pilotPadCell{iFrame} = pilot;
  timeSecCell{iFrame} = timeLoc;
end

model = struct();
model.array = {array};
model.rxSig = rxSig;
model.pilotPadCell = pilotPadCell;
model.timeSecCell = timeSecCell;
model.timeOffsetSec = timeOffsetSec;
model.frameMask = true(1, numFrame);
model.wavelength = wavelength;
model.phaseMode = phaseMode;
model.useLogObjective = false;
model.fdRateMode = 'unknown';
model.phaseGridCount = 72;
model.phaseRefine = true;

hyp = struct();
hyp.localDoaArr = localDoa;
hyp.fdSat = fdHz;
hyp.fdRateSat = fdRateHzPerSec;
end

function spread = localCircularSpread(phaseVec)
%LOCALCIRCULARSPREAD Compute a small circular phase spread.

phaseVec = reshape(phaseVec, 1, []);
phaseVec = phaseVec(isfinite(phaseVec));
if numel(phaseVec) <= 1
  spread = 0;
  return;
end
refPhase = phaseVec(1);
diffVec = angle(exp(1j * (phaseVec - refPhase)));
spread = max(abs(diffVec));
end
