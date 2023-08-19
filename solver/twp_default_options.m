
function options = twp_default_options()
options = [];
options.select_t = 0;
options.minimum_length = 1e-3;
options.maximum_mismatch = 1e-3;
options.threshold_alignment = 1e-6;
options.minimum_displacement = 1e-3;
options.threshold_planarity = 1e-6;
options.threshold_axis = 1e-6;
options.refine_initial = [0, 0];
options.refine_final = [0, 0];
options.threshold_refinement = [1e-6, 1e-6];
end
