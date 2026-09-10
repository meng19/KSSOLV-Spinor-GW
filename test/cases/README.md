# JSON case configuration

Run a case from any directory with:

```matlab
addpath('C:/Users/USTC/Documents/GitHub/KSSOLV-Spinor-GW');
gw_run_json('C:/path/to/case.json');
```

For a scheduler submission script, use the equivalent non-interactive command:

```sh
matlab -batch "addpath('C:/Users/USTC/Documents/GitHub/KSSOLV-Spinor-GW'); gw_run_json('C:/path/to/case.json')"
```

`qe_path` and the optional `save_file` are resolved relative to the JSON file,
not the submission script.  The JSON must define `qe_path`, `epsilon`, and
`sigma`.  `epsilon.nv` and `epsilon.nc` are optional: when absent, the runner
uses `options.nv` and `epsilon.nbnd - epsilon.nv` respectively.

`si8_isdf_smoke.json` is a complete small Si8 ISDF example.
