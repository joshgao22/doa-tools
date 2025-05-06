function [dg_rad, dg_display, dg_range, dg_list] = default_doa_grid(n, unit, dim)
%DEFAULT_DOA_GRID Creates the search grid for DOA estimation.
%
%   [dg_rad, dg_display, dg_range, dg_list] = default_doa_grid(n, unit, dim)
%
%   Inputs:
%     n     - Number of grid points. For 2D DOA, n can be a 2-element vector
%             specifying the number of points for azimuth and elevation respectively.
%     unit  - Unit of the grid values: 'radian', 'degree', or 'sin'.
%     dim   - Dimension of DOA: 1 (1D DOA) or 2 (2D DOA).
%
%   Outputs:
%     dg_rad     - Grid of candidate DOAs in radians.
%     dg_display - Grid of candidate DOAs in specified unit.
%     dg_range   - Range of the DOA candidates.
%                  For 1D: [1×2] vector; for 2D: [2×2] matrix.
%     dg_list    - Grid values per axis:
%                  • If dim == 1: dg_list = dg_display (1D array)
%                  • If dim == 2: dg_list = {az_list, el_list}
%                    where az_list and el_list are 1D vectors used in meshgrid.

if isempty(dim)
    if design.dim > 1
        dim = 2;
    else
        dim = 1;
    end
else
    if dim ~= 1 && dim ~= 2
        error('Incorrect DOA dimension.');
    end
end
if dim == 2
    if isscalar(n)
        n = [n n];
    end
end
switch lower(unit)
    case 'radian'
        if dim == 1
            dg_range = [-pi/2 pi/2];
            dg_rad = -pi/2:pi/n:(pi/2-pi/n);
            dg_display = dg_rad;

            dg_list = dg_display;
        else
            dg_range = [0 2*pi;0 pi/2];
            dg_rad = grid2(0, 2*pi, 0, pi/2, n(1), n(2));
            dg_display = dg_rad;

            az_list = unique(dg_display(1, :));
            el_list = unique(dg_display(2, :));
            dg_list = {az_list, el_list};
        end
    case 'degree'
        if dim == 1
            dg_range = [-pi/2 pi/2];
            dg_rad = -pi/2:pi/n:(pi/2-pi/n);
            dg_display = rad2deg(dg_rad);

            dg_list = dg_display;
        else
            dg_range = [0 2*pi;0 pi/2];
            dg_rad = grid2(0, 2*pi, 0, pi/2, n(1), n(2));
            dg_display = rad2deg(dg_rad);

            az_list = unique(dg_display(1, :));
            el_list = unique(dg_display(2, :));
            dg_list = {az_list, el_list};
        end
    case 'sin'
        if dim == 1
            dg_range = [-1 1];
            dg_display = -1:2/n:(1-2/n);
            dg_rad = asin(dg_display);

            dg_list = dg_display;
        else
            dg_range = [-1 1;0 1];
            dg_display = grid2(-1, 1, 0, 1, n(1), n(2));
            dg_rad = asin(dg_display);
            dg_rad(1,:) = dg_rad(1,:) + pi/2; % azimuth -> [0, pi]

            az_list = unique(dg_display(1, :));
            el_list = unique(dg_display(2, :));
            dg_list = {az_list, el_list};
        end
    otherwise
        error('Unknown unit ''%s''.', unit);
end
end


