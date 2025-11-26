function compute_trigreen_matrix(mesh_file, ncore, varargin)
%COMPUTE_TRIGREEN_MATRIX Compute triangular Green's function stiffness matrix
%
% Computes the elastic Green's function (stiffness matrix) for a triangular
% mesh using the Nikkhoo & Walter (2015) triangular dislocation method.
% Outputs results in binary format compatible with TriBIE Fortran code.
%
% Usage:
%   compute_trigreen_matrix(mesh_file, ncore)
%   compute_trigreen_matrix(mesh_file, ncore, 'param', value, ...)
%
% Inputs:
%   mesh_file - Path to .gts mesh file (GTS format)
%   ncore     - Number of processors (determines load distribution)
%
% Optional Parameters:
%   'mu'      - Shear modulus (default: 30000 MPa)
%   'nu'      - Poisson's ratio (default: 0.25)
%   'slip_ss' - Strike-slip component (default: -1.0)
%   'slip_ds' - Dip-slip component (default: 0.0)
%   'slip_ts' - Tensile-slip component (default: 0.0)
%   'output_dir' - Output directory (default: current directory)
%
% Outputs:
%   trigreen_<id>.bin - Binary files (one per processor)
%                       Each contains Nt*Ncell double precision values
%                       where Nt = local elements for processor id
%                             Ncell = total number of elements
%   position.bin      - Element centroid positions (created by processor 0)
%
% Example:
%   compute_trigreen_matrix('triangular_mesh.gts', 4)
%   compute_trigreen_matrix('mesh.gts', 8, 'mu', 32000, 'nu', 0.28)
%
% Reference:
%   Nikkhoo M. and Walter T.R., 2015. Triangular dislocation: An analytical,
%   artefact-free solution. Geophysical Journal International, 201(2), 1119-1141.
%
% Author: TriBIE Development Team
% Date: 2025-11-18

    %% Parse input arguments
    p = inputParser;
    addRequired(p, 'mesh_file', @ischar);
    addRequired(p, 'ncore', @(x) isnumeric(x) && x > 0);
    addParameter(p, 'mu', 30000, @isnumeric);           % Shear modulus (MPa)
    addParameter(p, 'nu', 0.25, @isnumeric);            % Poisson's ratio
    addParameter(p, 'slip_ss', -1.0, @isnumeric);       % Strike-slip
    addParameter(p, 'slip_ds', 0.0, @isnumeric);        % Dip-slip
    addParameter(p, 'slip_ts', 0.0, @isnumeric);        % Tensile-slip
    addParameter(p, 'output_dir', '.', @ischar);        % Output directory

    parse(p, mesh_file, ncore, varargin{:});

    mu = p.Results.mu;
    nu = p.Results.nu;
    ss = p.Results.slip_ss;
    ds = p.Results.slip_ds;
    ts = p.Results.slip_ts;
    output_dir = p.Results.output_dir;

    % Calculate lambda from Poisson's ratio and shear modulus
    lambda = (2 * nu * mu) / (1 - 2 * nu);

    fprintf('======================================================================\n');
    fprintf('TriGreen Stiffness Matrix Computation (MATLAB)\n');
    fprintf('======================================================================\n');
    fprintf('Mesh file: %s\n', mesh_file);
    fprintf('Number of processors: %d\n', ncore);
    fprintf('Material properties:\n');
    fprintf('  Shear modulus (mu):    %.2e MPa\n', mu);
    fprintf('  Poisson ratio (nu):    %.4f\n', nu);
    fprintf('  Lame parameter (lambda): %.2e MPa\n', lambda);
    fprintf('Slip components:\n');
    fprintf('  Strike-slip (ss): %.2f\n', ss);
    fprintf('  Dip-slip (ds):    %.2f\n', ds);
    fprintf('  Tensile-slip (ts): %.2f\n', ts);
    fprintf('Output directory: %s\n', output_dir);
    fprintf('======================================================================\n\n');

    %% Load mesh
    fprintf('Loading mesh from %s...\n', mesh_file);
    [vertices, cells] = load_gts_mesh(mesh_file);
    n_vertex = size(vertices, 1);
    n_cell = size(cells, 1);

    fprintf('  Vertices: %d\n', n_vertex);
    fprintf('  Elements: %d\n', n_cell);
    fprintf('Mesh loaded successfully.\n\n');

    %% Pre-compute cell centroids and triangle data
    fprintf('Pre-computing element centroids...\n');
    centroids = zeros(n_cell, 3);
    triangles = cell(n_cell, 1);

    for i = 1:n_cell
        v1 = vertices(cells(i,1), :);
        v2 = vertices(cells(i,2), :);
        v3 = vertices(cells(i,3), :);

        % Centroid
        centroids(i,:) = (v1 + v2 + v3) / 3;

        % Store triangle vertices
        triangles{i} = [v1; v2; v3];
    end
    fprintf('Centroids computed.\n\n');

    %% Distribute work across processors
    fprintf('Distributing work across %d processors...\n', ncore);
    [proc_cells, start_idx, end_idx] = distribute_cells(n_cell, ncore);

    for i = 0:ncore-1
        fprintf('  Processor %d: cells %d to %d (%d elements)\n', ...
                i, start_idx(i+1), end_idx(i+1), proc_cells(i+1));
    end
    fprintf('\n');

    %% Compute Green's functions for each processor
    fprintf('Computing Green''s functions...\n');
    fprintf('Progress: ');

    for proc = 0:ncore-1
        fprintf('\n  Processor %d/%d: ', proc, ncore-1);

        % Get cells for this processor
        cells_local = start_idx(proc+1):end_idx(proc+1);
        Nt = length(cells_local);

        % Allocate output matrix: Nt x Ncell
        trigreen_matrix = zeros(Nt, n_cell);

        % Compute Green's function for each source-observation pair
        for j = 1:Nt
            % Global cell index (source element)
            idx_source = cells_local(j);

            % Observation point (centroid of source element)
            obs_point = centroids(idx_source, :);

            % Progress indicator
            if mod(j, max(1, floor(Nt/10))) == 0
                fprintf('.');
            end

            % Compute influence from all elements
            for i = 1:n_cell
                % Get triangle vertices
                tri = triangles{i};
                P1 = tri(1,:);
                P2 = tri(2,:);
                P3 = tri(3,:);

                % Compute stress using Nikkhoo-Walter method
                [Stress, ~] = TDstressHS(obs_point(1), obs_point(2), obs_point(3), ...
                                         P1, P2, P3, ss, ds, ts, mu, lambda);

                % Extract shear stress component
                % trigreen = -mu/100 * dot(local_z, Stress * local_x)
                % For now, use simple scalar extraction (matches Fortran)
                trigreen_matrix(j, i) = -mu/100 * Stress(3,1);
            end
        end

        % Write output file
        output_file = fullfile(output_dir, sprintf('trigreen_%d.bin', proc));
        write_binary_matrix(output_file, trigreen_matrix);
        fprintf(' Done. Wrote %s\n', output_file);
    end

    fprintf('\nAll processors completed.\n\n');

    %% Write position file (processor 0 only)
    fprintf('Writing position file...\n');
    position_file = fullfile(output_dir, 'position.bin');
    fid = fopen(position_file, 'w');
    fwrite(fid, centroids', 'double');  % Write as column-major (Fortran compatible)
    fclose(fid);
    fprintf('  Wrote %s\n', position_file);

    fprintf('\n======================================================================\n');
    fprintf('Computation completed successfully!\n');
    fprintf('======================================================================\n');
end

function [vertices, cells] = load_gts_mesh(filename)
    % Load GTS format mesh file
    % Format:
    %   Line 1: n_vertex n_edge n_cell
    %   Lines 2 to n_vertex+1: x y z (vertex coordinates)
    %   Lines n_vertex+2 to end: v1 v2 v3 (cell connectivity)

    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open mesh file: %s', filename);
    end

    % Read header
    header = fscanf(fid, '%d %d %d', 3);
    n_vertex = header(1);
    n_edge = header(2);
    n_cell = header(3);

    % Read vertices
    vertices = fscanf(fid, '%f %f %f', [3, n_vertex])';

    % Read cells
    cells = fscanf(fid, '%d %d %d', [3, n_cell])';

    fclose(fid);
end

function [proc_cells, start_idx, end_idx] = distribute_cells(n_cell, ncore)
    % Distribute cells across processors with load balancing

    base_cells = floor(n_cell / ncore);
    extra_cells = mod(n_cell, ncore);

    proc_cells = zeros(ncore, 1);
    start_idx = zeros(ncore, 1);
    end_idx = zeros(ncore, 1);

    current_idx = 1;
    for i = 1:ncore
        if i <= extra_cells
            proc_cells(i) = base_cells + 1;
        else
            proc_cells(i) = base_cells;
        end

        start_idx(i) = current_idx;
        end_idx(i) = current_idx + proc_cells(i) - 1;
        current_idx = current_idx + proc_cells(i);
    end
end

function write_binary_matrix(filename, matrix)
    % Write matrix to binary file in Fortran-compatible format
    % Writes row-by-row (each row is Ncell values)

    fid = fopen(filename, 'w');
    if fid == -1
        error('Cannot open output file: %s', filename);
    end

    % Write each row
    for i = 1:size(matrix, 1)
        fwrite(fid, matrix(i,:), 'double');
    end

    fclose(fid);
end
