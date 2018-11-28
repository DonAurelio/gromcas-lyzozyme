# gromcas-lyzozyme

Files for each step:

* [1aki.pdb](https://files.rcsb.org/download/1AKI.pdb)
* [minim.mdp](http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp)
* [ions.mdp](http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp)

* [nvt.mdp](http://www.mdtutorials.com/gmx/lysozyme/Files/nvt.mdp)

For equilibarion step 

* [npt.mdp](http://www.mdtutorials.com/gmx/lysozyme/Files/npt.mdp)

For production molecular dynamics simulation step 

* [md.mdp](http://www.mdtutorials.com/gmx/lysozyme/Files/npt.mdp)

## Energy Minimization Step

Instuctions [here](http://www.mdtutorials.com/gmx/lysozyme/05_EM.html)

```sh
gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr
```

```sh
gmx mdrun -v -deffnm em
```

## Equilibration Step

NVT equilibration, stabilized the temperature of the system, [details](http://www.mdtutorials.com/gmx/lysozyme/06_equil.html). **NOTE: This step could takes 30 minutes more or less.**

```sh
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt
```

Equilibration of pressure is conducted under an NPT, [details](http://www.mdtutorials.com/gmx/lysozyme/07_equil2.html). **NOTE: This step could takes 51 minutes more or less.**

```sh
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt
```

## Production Molecular Dynamics Simulation 

Upon completion of the two equilibration phases, the system is now well-equilibrated at the desired temperature and pressure. We are now ready to release the position restraints and run production MD for data collection


```sh
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
gmx mdrun -deffnm md_0_1
```

Assuming you have one GPU available, the mdrun command to make use of it is as simple as:

```sh
gmx mdrun -deffnm md_0_1 -nb gpu
```

# References

[Gromacs Performance Measuring and Considerations](http://www.gromacs.org/Documentation/Performance_checklist)