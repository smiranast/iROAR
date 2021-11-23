# iROAR
Tool to adjust clonal count of VDJtool output table

## Recommended system configuration
* Linux, 2 CPU, 8GB RAM

* python=3.7.3
* matplotlib=3.0.3
* numpy=1.16.2
* pandas=0.24.2
* requests=2.21.0
* tqdm=4.43.0
* scipy==1.3.1

## Installation
1. Install the proper python3 version and all listed packages
2. Add iROAR directory to PATH in .bash_profile
3. Add the runfile execution permission
```
chmod +x .../iROAR/iroar
```
4. On the first run germline V/D/J segments of will be downloaded if aux/germline_* files are absent
