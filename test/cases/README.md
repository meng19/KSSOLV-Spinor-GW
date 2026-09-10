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

`qe_path`, `qp_file`, and the optional `save_file` are resolved relative to the JSON file,
not the submission script.  The JSON must define `qe_path`, `epsilon`, and
`sigma`.  `epsilon.nv` and `epsilon.nc` are optional: when absent, the runner
uses `options.nv` and `epsilon.nbnd - epsilon.nv` respectively.

After sigma, the runner writes `qp.dat` by default.  It contains BerkeleyGW-
style eV columns through `Eqp0`, with an additional `Eqp1` column.  Set
`qp_file` to choose its name or location.

`si8_isdf_smoke.json` is a complete small Si8 ISDF example.
