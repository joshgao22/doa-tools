function perSatProbeTable = buildMfPerSatProbeTable(repeatCell)
%BUILDMFPERSATPROBETABLE Collect compact per-satellite final-eval diagnostics.

rowList = repmat(emptyMfPerSatProbeRow(), 0, 1);
for iRepeat = 1:numel(repeatCell)
  repeatOut = repeatCell{iRepeat};
  if isfield(repeatOut, 'perSatProbeRowList') && ~isempty(repeatOut.perSatProbeRowList)
    rowList = [rowList; repeatOut.perSatProbeRowList(:)]; %#ok<AGROW>
  end
end
perSatProbeTable = struct2table(rowList(:));
end
