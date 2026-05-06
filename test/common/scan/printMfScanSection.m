function printMfScanSection(sectionTitle, sectionData)
%PRINTMFSCANSECTION Print a standard scan summary section banner.

arguments
  sectionTitle (1,:) char
  sectionData = []
end

fprintf('\n========== %s ==========\n', sectionTitle);
if nargin >= 2 && ~isempty(sectionData)
  disp(sectionData);
end
end
