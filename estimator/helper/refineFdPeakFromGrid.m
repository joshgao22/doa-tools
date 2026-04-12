function fdPeak = refineFdPeakFromGrid(scoreVec, fdGrid, bestGridIdx)
%REFINEFDPEAKFROMGRID Refine one 1-D Doppler peak by local parabola fitting.
% The helper first tries a log-parabola fit when the local scores are
% positive, and otherwise falls back to a standard quadratic interpolation.

arguments
  scoreVec (:,1) double
  fdGrid (:,1) double
  bestGridIdx (1,1) {mustBeInteger, mustBePositive}
end

fdPeak = fdGrid(bestGridIdx);
if numel(fdGrid) < 3 || bestGridIdx <= 1 || bestGridIdx >= numel(fdGrid)
  return;
end

scoreTriplet = real(scoreVec(bestGridIdx + (-1:1)));
if ~all(isfinite(scoreTriplet))
  return;
end

stepLeft = fdGrid(bestGridIdx) - fdGrid(bestGridIdx - 1);
stepRight = fdGrid(bestGridIdx + 1) - fdGrid(bestGridIdx);
if ~(isfinite(stepLeft) && isfinite(stepRight))
  return;
end
if abs(stepLeft - stepRight) > 1e-9 * max(abs([stepLeft, stepRight, 1]))
  return;
end
fdStep = 0.5 * (stepLeft + stepRight);

if all(scoreTriplet > 0)
  fitVal = log(max(scoreTriplet, realmin));
else
  fitVal = scoreTriplet;
end

denom = fitVal(1) - 2 * fitVal(2) + fitVal(3);
if ~isfinite(denom) || abs(denom) <= eps(max(abs(fitVal)))
  return;
end

delta = 0.5 * (fitVal(1) - fitVal(3)) / denom;
delta = min(max(delta, -1), 1);
fdPeak = fdGrid(bestGridIdx) + delta * fdStep;
end
