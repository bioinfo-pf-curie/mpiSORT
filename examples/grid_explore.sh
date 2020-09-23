#!/bin/bash


function usage()
{
	echo "Usage: "
        echo "bash test_bash: 0 argument (explore the grid and find all partitions)"
        echo "bash test_bash partition_name: 1 argument (explore the partition called partition_name)"
        echo "bash test_bash partition_name node_name: 2 arguments (explore the partition and node called partition_name and node_name )"
        exit 0
}

function compute_required_cpus()
{
	TotalCPUperNode=$1
	CPURequired=1
        while  [[ ${CPURequired} -lt ${TotalCPUperNode} ]]; do
                CPURequired=$(echo $((${CPURequired}*2)))
        done
        CPURequired=$((CPURequired/2))
	echo ${CPURequired}
}

function total_ram_1_chr(){
	TMP=$((5*$1/2))
	echo ${TMP}
}

function total_ram_all_chr(){
	echo $((3*$1/2))
}

function ceil(){
	RES=`echo $1 $2 | awk 'function ceil(x){return (x == int(x) ? x : int(x) + 1)}{ x=$2/$1; print ceil(x) }'`
	echo ${RES}
}

function explore_partition()
{
	########################
	#Explore partition Level#
	#########################

	PARTITION=$1

	TotalNodesTMP=$(scontrol show partition ${PARTITION} | grep 'TotalNodes' | cut -d ' ' -f 1,6 )
	TotalNodes=${TotalNodesTMP#*=}
	TotalCPUSTMP=$(scontrol show partition ${PARTITION} | grep 'TotalCPUs' | cut -d ' ' -f 1,5 )
	TotalCPUS=${TotalCPUSTMP#*=}
	MaxMemoryPerNodeTMP=$(scontrol show partition ${PARTITION} | grep 'MaxMemPerNode=' | cut -d ' ' -f 1,5 )
	MaxMemoryPerNode=${MaxMemoryPerNodeTMP#*=}
	NodesNameTMP=$(scontrol show partition ${PARTITION} | grep ' Nodes=')
	NodesName=${NodesNameTMP#*=}
	
	echo "After exploration of the partition " ${PARTITION}
	echo "Total number of nodes = " ${TotalNodes}
	echo "Total Number of CPUS  = " ${TotalCPUS}
	echo "Maximum Memory per Nodes = " ${MaxMemoryPerNode}
	echo "NODE names list = " ${NodesName}
	echo "########"
	
}

function explore_node()
{
	####################
	#Explore node level#
	####################
 
	NODE_NAME=$1

	TotalCPUperNodeTMP=$(scontrol show node ${NODE_NAME} | grep 'CPUTot=' | cut -d ' ' -f5)
	TotalCPUperNode=${TotalCPUperNodeTMP#*=}
	TotalRAMperNodeTMP=$(scontrol show node ${NODE_NAME} | grep 'RealMemory=' | cut -d ' ' -f4)
	TotalRAMperNode=$((${TotalRAMperNodeTMP#*=}/1000))
	TotalSocketPerNodeTMP=$(scontrol show node ${NODE_NAME} | grep 'Sockets=' | cut -d ' ' -f7)
	TotalSocketPerNode=${TotalSocketPerNodeTMP#*=}

	CPUPerSocket=$(echo $(( TotalCPUperNode / TotalSocketPerNode )))
	MemPerCPU=$(echo $(( TotalRAMperNode / TotalCPUperNode )))
		
	echo "After exploration of the node " ${NODE_NAME}
	echo "Number of CPUS per node  = " ${TotalCPUperNode}
	echo "Memory per Nodes = " ${TotalRAMperNode}
	echo "Number of sockets per node = " ${TotalSocketPerNode}
	echo "Number of CPUS per socket = " ${CPUPerSocket}
	echo "Average memory per CPU = " ${MemPerCPU} "GB"
	echo "########"
	echo " " 
	echo "Here are some examples of commande lines mpiSORT:"
	echo ""
	
	TotalRAMNeeded=$(total_ram_1_chr 50)
	NodesRequired=$(ceil ${TotalRAMperNode} ${TotalRAMNeeded})
	CPURequiredPerNode=$(compute_required_cpus ${TotalCPUperNode})
	
	echo "For instance if your SAM is 50gb and contains 1 chromosoms you will need around 2.5*50 = " ${TotalRAMNeeded}  "GB"  
	echo "So you need at least ${NodesRequired} nodes and from 2 to ${CPURequiredPerNode} CPUS per node"
	echo ""	
	echo "your PBS script could look like this"
	
	echo "#SBATCH -J JOB_NAME"
	echo "#SBATCH -N" ${NodesRequired}
	echo "#SBATCH -n" $((${CPURequiredPerNode}*${NodesRequired}))
	echo "#SBATCH -c 1"
	echo "#SBATCH --tasks-per-node="${CPURequiredPerNode}
	echo "#SBATCH --mem-per-cpu="$(((${NodesRequired}*${TotalRAMperNode})/(${CPURequiredPerNode}*${NodesRequired})))
	echo ""
	
	TotalRAMNeeded=$(total_ram_all_chr 500)
        NodesRequired=$(ceil ${TotalRAMperNode} ${TotalRAMNeeded})
        CPURequiredPerNode=$(compute_required_cpus ${TotalCPUperNode})


	echo "For instance if your SAM is 500gb and contains all chromosoms from HG19 you will need around 1.5*50 = " ${TotalRAMNeeded}  "GB"		
	echo "So you need at least ${NodesRequired} nodes and from 2 to ${CPURequiredPerNode} CPUS per node"
        echo "" 
        echo "your PBS script could look like this"

        echo "#SBATCH -J JOB_NAME"
        echo "#SBATCH -N" ${NodesRequired}
        echo "#SBATCH -n" $((${CPURequiredPerNode}*${NodesRequired}))
        echo "#SBATCH -c 1"
        echo "#SBATCH --tasks-per-node="${CPURequiredPerNode}
        echo "#SBATCH --mem-per-cpu="$(((${NodesRequired}*${TotalRAMperNode})/(${CPURequiredPerNode}*${NodesRequired})))
        echo ""
	
}

function explore_all_partition()
{
	#######################
        #Explore cluster level#
        #######################

	PartitionList=$(scontrol show partition | grep 'PartitionName=')
	PartitionNumber=$(scontrol show partition | grep 'PartitionName=' | wc -l)
	
	declare -A PartitionNameList
	declare i=1

	while [ ${i} -lt ${PartitionNumber} ]; do
		
		PartitionNameTMP=$(echo ${PartitionList} | cut -d ' ' -f${i})
		PartitionName=${PartitionNameTMP#*=}
		echo "########"
		echo "Partition "$i" name : "${PartitionName}
		PartitionNameList[$i]=${PartitionName}
		explore_partition ${PartitionName}
		((i++))
	done
	echo " "
	echo "Number of partition = "${PartitionNumber}
}

################################
## Test if slurm is installed ##
################################

type scontrol

if [[ ! $? -eq 0 ]]
then
    echo "ERROR: Slurm not found in your PATH. please check its installation."
    exit 1
fi

while getopts ":h" option; do 
	case "$option" in
		h) usage;;
	esac
done

if [[ $# -eq 0  ]]; then
	echo "Slurm is installed we are going to check all the partitions present on your cluster."
	explore_all_partition
	exit 0
elif [[ ( $# -eq 1 ) && ( ( $# != "--help" ) || ( $# != "-h") ) ]]; then
	echo "Slurm is installed we are going to check the partitions $1 on your cluster."
	explore_partition $1
	exit 0
elif [[ $# -eq 2 ]]; then
	echo "########"
	echo "Slurm is installed we are going to check the partitions $1 and node $2 on your cluster."
	explore_partition $1
	explore_node $2
	exit 0
else
	echo "Problem with your arguments see usage (-h)"
	exit 1
fi

exit 0
