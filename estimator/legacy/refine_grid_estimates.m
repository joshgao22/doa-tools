function est = refine_grid_estimates(f_obj, grid_list, est_idx)
%REFINE_GRID_ESTIMATES Refines DOA estimates using a subgrid search.
%
%   est = refine_grid_estimates(f_obj, grid_list, est_idx)
%
%   Inputs:
%     f_obj     - Objective function handle. Should accept a vector or matrix
%                 of candidate DOAs and return scalar values to minimize.
%     grid_list - DOA grid list from default_doa_grid:
%                   • If 1D: a vector of DOA values (length N)
%                   • If 2D: a cell {az_list, el_list}, each a 1D vector
%     est_idx   - Estimated indices corresponding to initial grid:
%                   • If 1D: 1×n vector of indices
%                   • If 2D: n×2 array of [az_idx, el_idx] pairs
%
%   Output:
%     est       - Refined DOA estimates:
%                   • 1×n vector for 1D
%                   • n×2 matrix for 2D

    % === 1D case ===
    if ~iscell(grid_list)
        grid = grid_list(:);  % ensure column
        est = zeros(size(est_idx));
        n_iter = 10;
        subgrid_size = 10;

        for kk = 1:length(est_idx)
            idx = est_idx(kk);
            % Determine bounds
            lb = grid(max(idx - 1, 1));
            ub = grid(min(idx + 1, length(grid)));

            % Refine by iterative search
            for ii = 1:n_iter
                subgrid = linspace(lb, ub, subgrid_size);
                obj_vals = f_obj(subgrid);
                [~, min_idx] = min(obj_vals);

                % Update bounds
                if min_idx > 1
                    lb = subgrid(min_idx - 1);
                else
                    lb = subgrid(1);
                end
                if min_idx < subgrid_size
                    ub = subgrid(min_idx + 1);
                else
                    ub = subgrid(end);
                end
            end
            est(kk) = subgrid(min_idx);
        end

    % === 2D case ===
    elseif iscell(grid_list) && numel(grid_list) == 2
        az_list = grid_list{1};
        el_list = grid_list{2};

        est = zeros(size(est_idx));
        n_iter = 10;
        subgrid_size = [4 4];

        for kk = 1:size(est_idx, 2)
            az_idx = est_idx(1, kk);
            el_idx = est_idx(2, kk);

            % Initial bounds
            lb_az = az_list(max(az_idx - 1, 1));
            ub_az = az_list(min(az_idx + 1, length(az_list)));
            lb_el = el_list(max(el_idx - 1, 1));
            ub_el = el_list(min(el_idx + 1, length(el_list)));

            % Refinement loop
            for ii = 1:n_iter
                az_sub = linspace(lb_az, ub_az, subgrid_size(1));
                el_sub = linspace(lb_el, ub_el, subgrid_size(2));
                [az_grid, el_grid] = meshgrid(az_sub, el_sub);
                candidates = [az_grid(:), el_grid(:)]';
                obj_vals = f_obj(candidates);

                [~, min_idx] = min(obj_vals);
                [min_el_idx, min_az_idx] = ind2sub(subgrid_size, min_idx);

                % Update bounds
                lb_az = az_sub(max(min_az_idx - 1, 1));
                ub_az = az_sub(min(min_az_idx + 1, end));
                lb_el = el_sub(max(min_el_idx - 1, 1));
                ub_el = el_sub(min(min_el_idx + 1, end));
            end

            est(:, kk) = [az_sub(min_az_idx), el_sub(min_el_idx)];
        end
    else
        error('grid_list must be a vector (1D) or a 2-element cell array (2D).');
    end
end
