# Scripts and data to test the program

Your can try:

* `bash standard.sh`
* `sbatch slurm.sh` to test it with [Slurm](https://slurm.schedmd.com/sbatch.html)
* `qsub pbs.sh` to test it with [PBS/Torque](https://support.adaptivecomputing.com/support/documentation-index/torque-resource-manager-documentation/)
* `bash getSlurmNodesInfo.sh` to check what [Slurm](https://slurm.schedmd.com/sbatch.html) partitions and node are present on your cluster
* `python informaticResources.py -h` to help you choosing computer resources for sorting
* `python getSlurmOpt -c {TotalCPUPerNode} -m {TotalMemoryPerNode} -s {SizeOfSAM} -t {SAMType}`

getSlurmOpt:
	input: 
	   	-c: from informaticResources.py 
		-m: from informaticResources.py in GB
		-s: Sam size in GB
	        -t: uniq or all (contain 1 chromosom or all)

        output: a string line in forms of 
		-N {} -n {} -c 1 --tasks-per-node={} --mem-per-cpu={}GB 
		to pass it directly to sbatch
