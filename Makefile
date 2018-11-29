# MPI main process PID
MPI_PROCESS=$$(ps -fu mpiuser | grep mdrun | awk 'NR==1{printf "%s", $$2}' )

MPI_FLAGS= --mca btl_tcp_if_include eno1  -np 16 --hostfile hostfile

# ENERGY MINIMIZATION
em_compile:
	gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr

em_seq: clean_backups em_clean em_compile
	gmx mdrun -v -deffnm em

em_mpi: clean_backups em_clean em_compile
	mpirun ${MPI_FLAGS} gmx_mpi mdrun -deffnm em

em_clean:
	rm -f em.* mdout.mdp

# EQUILIBRATION NVT
nvt_compile:
	gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr

nvt_seq: clean_backups nvt_clean nvt_compile
	gmx mdrun -v -deffnm nvt

nvt_mpi: clean_backups nvt_clean nvt_compile
	mpirun ${MPI_FLAGS} mdrun_mpi -deffnm nvt

nvt_clean:
	rm -f nvt.cpt nvt.edr nvt.gro nvt.log nvt.tpr nvt.trr nvt_prev.cpt mdout.mdp

# EQUILIBRATION NPT
npt_compile: 
	gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr

npt_seq: clean_backups npt_clean npt_compile
	gmx mdrun -deffnm npt

npt_mpi: clean_backups npt_clean npt_compile
	mpirun ${MPI_FLAGS} mdrun_mpi -deffnm npt

npt_clean:
	rm -f npt.cpt npt.edr npt.gro npt.log npt.tpr npt.trr npt_prev.cpt

# PRODUTION

# Creates a binary file 'nvt_interfaz_SDS.tpr' given the gromcas topology and simulation variables
prod_compile:
	gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr

# Run gromcas distributing the work among computers, but only use CPUs
prod_mpi: clean_logs prod_compile
	mpirun ${MPI_FLAGS} mdrun_mpi -ntmpi 16 -v -deffnm md_0_1 > prod_out.log 2> prod_err.log &

prod_seq: clean_logs prod_compile
	gmx mdrun -deffnm md_0_1  > prod_out.log 2> prod_err.log &

# Show the logs of the simulation running in cpu
log:
	tail -f prod_err.log

# Stop simulation
stop:
	kill  ${MPI_PROCESS}

# Remove the backup files generated during the simulation. WARNING: You don't want to do this really !!
clean_backups:
	rm -f \#*\#

clean_logs:
	rm -f prod_out.log prod_err.log
