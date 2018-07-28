%% Distance functional surface
xleft = -6;

cleft_start = 20; cleft_end = 30;
cright_start = 20; cright_end = 30;
cleft_span = [-4 4]; cright_span = [-4 4]; cstep = 0.1;

cleft_grid = cleft_start:cstep:cleft_end;
cright_grid = cright_start:cstep:cright_end;

distance = zeros(length(cleft_grid), length(cright_grid));

for i = 1:length(cleft_grid)
    i
    for j = 1:length(cright_grid)
        distance(i, j) = shots_distance(params, xleft, cleft_grid(i), cright_grid(j));
    end
end

figure('Position', [100 100 500 400]); hold on

[cleft_grid_mesh, cright_grid_mesh] = meshgrid(cleft_grid, cright_grid);
surf(cleft_grid_mesh, cright_grid_mesh, distance);
view(2)