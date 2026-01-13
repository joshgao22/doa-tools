function [est, est_idx, resolved] = find_doa_from_spectrum(dg, y, n, varargin)
%FIND_DOA_FROM_SPECTRUM Estimates DOAs by locating spectral peaks.
% Supports both 1D and 2D DOA grids.
%
% Inputs:
%   dg        - DOA grid (1D: [1×N], 2D: [2×N] with azimuth; elevation)
%   y         - Spectrum values (1×N)
%   n         - Number of expected DOAs
%   varargin  - Optional: grid_size = [num_azimuth_points, num_elevation_points]
%
% Outputs:
%   est       - Estimated DOAs (1×n or 2×n)
%   est_idx   - Estimated indices (1×n for 1D, 2×n for 2D)
%   resolved  - Boolean flag

if ~isempty(varargin)
    grid_size = varargin{1};
else
    grid_size = [];
end

dim = size(dg, 1);

if size(dg, 2) ~= numel(y)
    error('Size mismatch between dg and y.');
end

if dim == 1
    % --- 1D case: use findpeaks directly ---
    [peaks, idx_est] = findpeaks(y);
    
    if numel(idx_est) < n
        est = [];
        est_idx = [];
        resolved = false;
        return;
    end

    [~, sort_idx] = sort(peaks, 'descend');
    idx_selected = idx_est(sort_idx(1:n));
    est = dg(idx_selected);
    est_idx = idx_selected;
    resolved = true;

elseif dim == 2
    % --- 2D case: use imregionalmax ---
    if isempty(grid_size)
        [nx, ny] = grid_size_from_dg(dg);  % infer from dg
    else
        nx = grid_size(1);
        ny = grid_size(2);
    end

    y2d = reshape(y, ny, nx);  % [el, az]
    peak_mask = imregionalmax(y2d);
    peak_vals = y2d(peak_mask);

    if numel(peak_vals) < n
        est = [];
        est_idx = [];
        resolved = false;
        return;
    end

    [~, sort_idx] = sort(peak_vals, 'descend');
    [el_idx, az_idx] = find(peak_mask);
    sel_el_idx = el_idx(sort_idx(1:n));
    sel_az_idx = az_idx(sort_idx(1:n));

    lin_idx = sub2ind([ny, nx], sel_el_idx, sel_az_idx);
    est = dg(:, lin_idx);  % 2 × n
    est_idx = [sel_az_idx'; sel_el_idx'];  % az_idx, el_idx
    resolved = true;
else
    error('Only 1D and 2D direction grids are supported.');
end
end


function [nx, ny] = grid_size_from_dg(dg)
%GRID_SIZE_FROM_DG Estimates the original grid size [nx, ny] from flattened dg.
%Assumptions:
%   - The 2D grid is constructed by grid2().
%   - In grid2(), azimuth varies slowly, elevation varies quickly.
%     (i.e., elevation is the faster-changing dimension.)
%Method:
%   - The number of unique elevation values determines ny.
%   - The total number of points divided by ny gives nx.

elev = dg(2,:); 
unique_el = numel(unique(round(elev, 8))); % safe rounding
ny = unique_el;
nx = size(dg,2) / ny;
end
