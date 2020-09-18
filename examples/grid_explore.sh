#!/bin/bash



#PARTITION=$(echo $1)
#NODE_NAME=$(echo $2)
#FILE=$(echo $3)

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

	echo "########"
	echo " "
	echo "After exploration of the partition " ${PARTITION}
	echo "Total number of nodes = " ${TotalNodes}
	echo "Total Number of CPUS  = " ${TotalCPUS}
	echo "Maximum Memory per Nodes = " ${MaxMemoryPerNode}
	echo "NODE names list = " ${NodesName}
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
	TotalRAMperNode=${TotalRAMperNodeTMP#*=}
	TotalSocketPerNodeTMP=$(scontrol show node ${NODE_NAME} | grep 'Sockets=' | cut -d ' ' -f7)
	TotalSocketPerNode=${TotalSocketPerNodeTMP#*=}

	CPUPerSocket=$(echo $(( TotalCPUperNode / TotalSocketPerNode )))
	MemPerCPU=$(echo $(( TotalRAMperNode / TotalCPUperNode )))
	
	echo "########"
	echo " "
	echo "After exploration of the node " ${NODE_NAME}
	echo "Number of CPUS per node  = " ${TotalCPUperNode}
	echo "Memory per Nodes = " ${TotalRAMperNode}
	echo "Number of sockets per node = " ${TotalSocketPerNode}
	echo "Number of CPUS per socket = " ${CPUPerSocket}
	echo "Memory per CPU = " ${MemPerCPU}
}

function explore_all_partition()
{
	#######################
        #Explore cluster level#
        #######################

	PartitionList=$(scontrol show partition | grep 'PartitionName=')
	PartitionNumber=$(scontrol show partition | grep 'PartitionName=' | wc -l)
	
	echo "Number of partition = "${PartitionNumber}

	declare -A PartitionNameList
	declare i=1

	while [ ${i} -lt ${PartitionNumber} ]; do
		
		PartitionNameTMP=$(echo ${PartitionList} | cut -d ' ' -f${i})
		PartitionName=${PartitionNameTMP#*=}
		echo "############"
		echo " " 
		echo "Partition "$i" name : "${PartitionName}
		PartitionNameList[$i]=${PartitionName}
		explore_partition ${PartitionName}
		((i++))
	done
}

if [[ $# -eq 0  ]]; then
	explore_all_partition
elif [[ $# -eq 1 ]]; then
	explore_partition $1
elif [[ $# -eq 2 ]]; then
	explore_partition $1
	explore_node $2
else
	echo "Problem with your arguments"
	echo "Usage: "
	echo "bash test_bash: 0 argument (explore the grid and find all partitions)"
	echo "bash test_bash partition_name: 1 argument (explore the partition called partition_name)"
	echo "bash test_bash partition_name node_name: 2 arguments (explore the partition and node called partition_name and node_name )"
	
