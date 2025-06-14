# Setup
## Run docker
(You should install docker and turn it on first)
```
(under the directory you decompressed our project)
$ docker run --platform linux/amd64 -it -v $(pwd):/workspace --name openroad_container openroad/flow-ubuntu22.04-dev
$ docker start openroad_container
$ docker attach openroad_container
$ cd /workspace/OpenROADForCAD/
```
##  Build
```
$ ./etc/Build.sh
```
# Usage
```
$ openroad eval_def.tcl
```