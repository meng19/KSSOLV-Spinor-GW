% Fast end-to-end smoke test for the portable JSON case launcher.
% All calculation parameters and paths are maintained in the adjacent JSON.

config_file = '/public/home/mxy/work/spinor_gw/speed/KSSOLV-Spinor-GW/si8.json';
profile_dir = fullfile(fileparts(config_file), 'si8_profile');

% MATLAB profiler: real (wall-clock) time is appropriate for FFT, I/O,
% GPU, and future parallel sections.  The existing gw_timer report remains
% the lightweight summary for the main GW stages.
profile clear;
profile on -timer real;
try
    [eps, sig] = gw_run_json(config_file);
catch err
    profile off;
    rethrow(err);
end
profile off;

if ~isfolder(profile_dir)
    mkdir(profile_dir);
end
profile_data = profile('info');
save(fullfile(profile_dir, 'si8_profile.mat'), 'profile_data');
profsave(profile_data, profile_dir);
fprintf('MATLAB profile written to: %s\n', profile_dir);
