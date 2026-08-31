% Regenerate the scalar-paper datasets and all numerical figures.

code_dir=fileparts(mfilename('fullpath'));
start_dir=pwd;
restore_dir=onCleanup(@() cd(start_dir));
cd(code_dir);
fig_dir=fullfile(code_dir,'figures');
if ~isfolder(fig_dir), mkdir(fig_dir); end

check_barrier_evaluation;

run_reg_sweep;
run_reg_sweep2;
run_exact_all;
run_numerics_extra;
run_numerics_gamma_mesh;

make_sweep_figs;
make_scalar_figs;
make_extra_figs;
make_solution_fig;
make_solution_exact_figs;

fprintf('Scalar-paper data and figures regenerated.\n');
