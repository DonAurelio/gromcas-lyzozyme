# Computers on which gromcas will run on a distributed fashion
HOSTS=localhost,slave_0

# MPI main process PID
MPI_PROCESS=$$(ps -fu mpiuser | grep mpirun | awk 'NR==1{printf "%s", $$2}' )

# Creates a binary file 'nvt_interfaz_SDS.tpr' given the gromcas topology and simulation variables
compile:
	gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr

# Run gromcas distributing the work among computers, but only use CPUs
run: compile
	mpirun --host ${HOSTS}  mdrun_mpi -v -deffnm md_0_1 > cpu_out.log 2> cpu_err.log &

# Show the logs of the simulation running in cpu
log:
	tail -f cpu_err.log

monitor:
	watch -n 1 ./monitor.sh

# Stop simulation
stop:
	kill  ${MPI_PROCESS}

# Remove the backup files generated during the simulation. WARNING: You don't want to do this really !!
clean_backups:
	rm -f \#*\#
