function summaryTable = buildDynamicCaseSummaryTable(caseName, snrDb, angleErrDeg, fdErrHz, fdRateErrHzPerSec, isResolved)
%BUILDDYNAMICCASESUMMARYTABLE Build one compact single-SNR dynamic case summary.

arguments
  caseName (:, 1) string
  snrDb (1, 1) double
  angleErrDeg (:, :) double
  fdErrHz (:, :) double
  fdRateErrHzPerSec (:, :) double
  isResolved (:, :) logical
end

numCase = numel(caseName);
angleRmseDeg = nan(numCase, 1);
angleP95Deg = nan(numCase, 1);
fdRmseHz = nan(numCase, 1);
fdP95Hz = nan(numCase, 1);
fdRateRmseHzPerSec = nan(numCase, 1);
fdRateP95HzPerSec = nan(numCase, 1);
resolveRate = nan(numCase, 1);

for iCase = 1:numCase
  angleFail = ~isResolved(:, iCase) | ~isfinite(angleErrDeg(:, iCase));
  angleStat = summarizeMonteCarloStat(angleErrDeg(:, iCase), angleFail);
  angleRmseDeg(iCase) = angleStat.rmse;
  angleP95Deg(iCase) = angleStat.p95;
  resolveRate(iCase) = 1 - angleStat.failRate;

  if any(isfinite(fdErrHz(:, iCase)))
    fdFail = ~isResolved(:, iCase) | ~isfinite(fdErrHz(:, iCase));
    fdStat = summarizeMonteCarloStat(fdErrHz(:, iCase), fdFail);
    fdRmseHz(iCase) = fdStat.rmse;
    fdP95Hz(iCase) = fdStat.p95;
  end

  if any(isfinite(fdRateErrHzPerSec(:, iCase)))
    fdRateFail = ~isResolved(:, iCase) | ~isfinite(fdRateErrHzPerSec(:, iCase));
    fdRateStat = summarizeMonteCarloStat(fdRateErrHzPerSec(:, iCase), fdRateFail);
    fdRateRmseHzPerSec(iCase) = fdRateStat.rmse;
    fdRateP95HzPerSec(iCase) = fdRateStat.p95;
  end
end

summaryTable = table(caseName(:), repmat(snrDb, numCase, 1), angleRmseDeg, angleP95Deg, ...
  fdRmseHz, fdP95Hz, fdRateRmseHzPerSec, fdRateP95HzPerSec, resolveRate, ...
  'VariableNames', {'displayName', 'snrDb', 'angleRmseDeg', 'angleP95Deg', ...
  'fdRmseHz', 'fdP95Hz', 'fdRateRmseHzPerSec', 'fdRateP95HzPerSec', 'resolveRate'});
end
