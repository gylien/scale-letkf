### External link
[Oakforest-PACS webpage](https://www.cc.u-tokyo.ac.jp/supercomputer/ofp/service/)  
[User portal](https://ofp-www.jcahpc.jp/cgi-bin/hpcportal.en/index.cgi)

### Version 

Use the branch [https://gitlab.com/scale-met/scale/tree/5.3.x_LETKF/](https://gitlab.com/scale-met/scale/tree/5.3.x_LETKF/)

### Compile SCALE

* **Set environment variables**
 ```ShellSession
$ export SCALE_SYS=OFP
$ export SCALE_ENABLE_OPENMP=T ### necessary for SCALE-5.3
 ```

* **Compile SCALE-RM**
 ```ShellSession
$ cd scale/scale-rm/src
$ make  # 'make -j 8' for parallel compiling using 8 threads
 ```

It uses configuration files in `scale/sysdep/`.

See [SCALE users Guide](http://r-ccs-climate.riken.jp/scale/doc/scale_users_guide.v5.2.5.pdf) for more details.

### To compile SCALE-LETKF

* **Link or copy the configuration file**
 ```ShellSession
$ cd scale/letkf/scale
$ ln -s arch/configure.user.ofp configure.user  # Note the file name change
 ```

* **Compile all LETKF programs**
 ```ShellSession
$ make
 ```

Example configuration files are in `scale/letkf/scale/arch/`.