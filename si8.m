% Fast end-to-end smoke test for the portable JSON case launcher.
% All calculation parameters and paths are maintained in the adjacent JSON.

config_file = '/public/home/mxy/work/spinor_gw/speed/KSSOLV-Spinor-GW/si8.json');
[eps, sig] = gw_run_json(config_file);

