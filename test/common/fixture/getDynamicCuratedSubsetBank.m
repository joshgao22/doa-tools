function [primaryOffsetCell, primaryLabelList, rescueOffsetCell, rescueLabelList] = getDynamicCuratedSubsetBank(bankName)
%GETDYNAMICCURATEDSUBSETBANK Return fixed curated subset schedules for dynamic MF runs.
% The fast statistics path should avoid true-random subset fallback because it
% hurts reproducibility and makes regression coverage noisy. This helper keeps
% one small deterministic bank: two primary schedules for the main path and
% two nonperiodic rescue schedules that break the 1/T_f comb symmetry when the
% primary pair disagrees or lands on a bad tooth.

arguments
  bankName (1, 1) string = "default"
end

switch lower(strtrim(bankName))
  case {"default", "dual-sat", "dynamic-mf"}
    primaryOffsetCell = { ...
      [-8, -7, -5, -4, -2, -1, 0, 4, 7, 9], ...
      [-7, -4, -1, 0, 3, 5, 7, 8, 9, 10]};
    primaryLabelList = ["curated1"; "curated2"];

    rescueOffsetCell = { ...
      [-9, -7, -5, -2, -1, 0, 1, 7, 9, 10], ...
      [-6, -3, -2, -1, 0, 2, 3, 6, 8, 10]};
    rescueLabelList = ["curated3"; "curated4"];

  otherwise
    error('getDynamicCuratedSubsetBank:UnknownBank', ...
      'Unknown dynamic curated subset bank "%s".', bankName);
end
end
