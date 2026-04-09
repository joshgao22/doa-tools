function plotDoaDopplerGeometryComparison(latlonTrue, caseResult)
%PLOTDOADOPPLERGEOMETRYCOMPARISON Plot all lat-lon estimates on one figure.

markerList = {'o', 's', '^', 'v', 'x', 'd', '>', '<'};
figure();
plot(latlonTrue(2), latlonTrue(1), 'kp', 'MarkerSize', 12, 'LineWidth', 1.3); hold on;

legendCell = {'truth'};
for iCase = 1:numel(caseResult)
  latlonEst = getDoaDopplerLatlonEst(caseResult(iCase).estResult);
  markerNow = markerList{1 + mod(iCase - 1, numel(markerList))};
  plot(latlonEst(2), latlonEst(1), markerNow, 'MarkerSize', 8, 'LineWidth', 1.2);
  legendCell{end + 1} = char(caseResult(iCase).displayName); %#ok<AGROW>
end

grid on;
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
legend(legendCell, 'Location', 'best');
title('Lat-lon estimates');
end
