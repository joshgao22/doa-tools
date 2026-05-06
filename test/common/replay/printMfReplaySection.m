function printMfReplaySection(sectionTitle, sectionData)
%PRINTMFREPLAYSECTION Print a standard replay summary section banner.

arguments
  sectionTitle (1,:) char
  sectionData = []
end

fprintf('\n========== %s ==========\n', sectionTitle);
if nargin >= 2 && ~isempty(sectionData)
  disp(sectionData);
end
end
