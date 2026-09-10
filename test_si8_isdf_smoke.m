% Fast end-to-end smoke test for the portable JSON case launcher.
% All calculation parameters and paths are maintained in the adjacent JSON.

test_root = fileparts(mfilename('fullpath'));
config_file = fullfile(test_root, 'test', 'cases', 'si8_isdf_smoke.json');
[eps, sig] = gw_run_json(config_file);

assert(all(isfinite(sig.sig(:))), 'Si8 smoke test produced non-finite sigma.');
fprintf('SI8_ISDF_SMOKE_OK: %d diagonal bands.\n', size(sig.sig, 1));
