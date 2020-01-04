### Simultaneous multiple runs

The scripts for ensemble forecast `run_d1_ext`,`run_d1-2`,`run_d3` and `run_d4_init` support multiple runs at the same time.

### DISK_MODE setting

File staging is not used in current OFP system (possible but not effective). 
Therefore `DISK_MODE = 1` is preferable because it saves time by making symbolic links of input/output paths instead of 
unnecessarily copying large files to a temporal runtime directory. 

### Automatic adjustment of a node size 

