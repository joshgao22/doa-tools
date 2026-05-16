function legendList = enableLegendToggle(axList)
%ENABLELEGENDTOGGLE Enable line visibility toggle by clicking legend entries.
%
% legendList = enableLegendToggle(axList) attaches a legend ItemHitFcn to
% each axes in axList. Clicking a legend entry toggles the corresponding
% plotted object visibility while keeping the legend entry available.

if nargin < 1 || isempty(axList)
  axList = gca;
end

axList = axList(:);
legendList = gobjects(0, 1);
for iAx = 1:numel(axList)
  axNow = axList(iAx);
  if ~isgraphics(axNow, 'axes')
    continue;
  end

  lgd = legend(axNow, 'show', 'Location', 'best');
  lgd.AutoUpdate = 'off';
  lgd.ItemHitFcn = @localToggleLegendItem;
  legendList(end + 1, 1) = lgd; %#ok<AGROW>
end
end

function localToggleLegendItem(~, eventData)
%LOCALTOGGLELEGENDITEM Toggle the plotted object tied to a legend entry.

targetObj = eventData.Peer;
if ~isgraphics(targetObj) || ~isprop(targetObj, 'Visible')
  return;
end

if strcmp(targetObj.Visible, 'on')
  targetObj.Visible = 'off';
else
  targetObj.Visible = 'on';
end
end
