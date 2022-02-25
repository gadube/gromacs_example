#!/bin/bash -e
#THREADS=$1
THREADS=24
ALIGNMENTS=0 32 64

cd /home/gad4877/proj/bigfoot
rm -f ./output/STATUS.out

echo "Cleaning Apps for new test" >> ./output/STATUS.out 2>&1
make clean-apps
make GROMACS >> ./output/STATUS.out 2>&1

mkdir -p /home/gad4877/proj/bigfoot/apps/GROMACS/{bigfootalign0,bigfootalign32,bigfootalign64}/results
# cp /home/gad4877/proj/gmxtest/lysozyme_in_water/gromacs_run.sh /home/gad4877/proj/bigfoot/apps/GROMACS/{bigfootalign0,bigfootalign32,bigfootalign64}/results/gromacs_run.sh

. /home/gad4877/proj/bigfoot/enable
for $a in $alignments
do
    echo "Running GROMACS test for alignment $a" >> ./output/STATUS.out 2>&1
    cd /home/gad4877/proj/bigfoot/apps/GROMACS/bigfootalign$a/results
    . /home/gad4877/proj/bigfoot/apps/GROMACS/bigfootalign$a/bin/GMXRC
    # Tutorial from (http://www.mdtutorials.com/gmx/lysozyme) , Lysozyme in Water

    # grab 1AKI from RSCB PDB https://www.rcsb.org/structure/1AKI, this is our structure
    wget https://files.rcsb.org/download/1AKI.pdb

    # remove water molecules from the pdb file, don't need them for this tutorial
    grep -v HOH 1AKI.pdb > 1AKI_clean.pdb

    # run pdb2gmx and choose a force field, we choose 15, OPLS all atom force field
    echo "15" | gmx pdb2gmx -f 1AKI_clean.pdb -o 1AKI_processed.gro -water spce

    # define a box and fill with a solvent, for this case a cubic box
    # protein is placed in center of box, at least 1 nm from the edge and defines box as a cube
    gmx editconf -f 1AKI_processed.gro -o 1AKI_newbox.gro -c -d 1.0 -bt cubic

    # fill box with solvent (water)
    echo "Running solvate command for alignment $a" >> ./output/STATUS.out 2>&1
    gmx solvate -cp 1AKI_newbox.gro -cs spc216.gro -o 1AKI_solv.gro -p topol.top

    # grab ions file and use it to generate an atomic description of system
    # Assemble .tpr file using grompp
    wget http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp
    gmx grompp -f ions.mdp -c 1AKI_solv.gro -p topol.top -o ions.tpr

    # use genion to add ions to system, replacing some water molecules, choose group 13 for embedding ions in the solvent
    echo "13" | gmx genion -s ions.tpr -o 1AKI_solv_ions.gro -p topol.top -pname NA -nname CL -neutral

    # ensure appropriate system by relaxing structure in energy minimization,
    # first assemble structure topology using grompp and the provided input parameter file
    echo "Running energy minimization for alignment $a" >> ./output/STATUS.out 2>&1
    wget http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp
    gmx grompp -f minim.mdp -c 1AKI_solv_ions.gro -p topol.top -o em.tpr

    # run on 12 threads and 1 MPI proc, name files em.*
    OMP_PROC_BIND=spread OMP_PLACES=cores OMP_NUM_THREADS=$THREADS gmx mdrun -v -ntmpi 1 -deffnm em

    # create energy term plot
    echo "10 0" | gmx energy -f em.edr -o potential.xvg

    # create structure using provided .mpd input and run equilibration
    echo "Running pressure and density equilibration for alignment $a" >> ./output/STATUS.out 2>&1
    wget http://www.mdtutorials.com/gmx/lysozyme/Files/nvt.mdp
    gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
    OMP_PROC_BIND=spread OMP_PLACES=cores OMP_NUM_THREADS=$THREADS gmx mdrun -deffnm nvt

    # Analyze temperature progression with energy, use 16 0 to select temperature of system
    echo "16 0" | gmx energy -f nvt.edr -o temperature.xvg

    # Previous step stabalized temperature of tye system, now we stabalize pressure and density
    echo "Running temperature equilibration for alignment $a" >> ./output/STATUS.out 2>&1
    wget http://www.mdtutorials.com/gmx/lysozyme/Files/npt.mdp
    gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
    OMP_PROC_BIND=spread OMP_PLACES=cores OMP_NUM_THREADS=$THREADS gmx mdrun -deffnm npt

    # Analyze pressure progression and density
    echo "18 0" | gmx energy -f npt.edr -o pressure.xvg
    echo "24 0" | gmx energy -f npt.edr -o density.xvg

    # Run a 1 ns MD simulation, releasing restraints, using md input config
    echo "Running 1 ns MD sim for alignment $a" >> ./output/STATUS.out 2>&1
    wget http://www.mdtutorials.com/gmx/lysozyme/Files/md.mdp
    gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
    OMP_PROC_BIND=spread OMP_PLACES=cores OMP_NUM_THREADS=$THREADS gmx mdrun -deffnm md_0_1

    # Analyze Results, correct for periodicity
    echo "1 0" | gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center

    # Calculate RMSD on corrected trajectory for minimized, equilibrated system
    echo "4 4" | gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns

    # Calculate RMSD on corrected trajectory for crystal structure
    echo "4 4" | gmx rms -s em.tpr -f md_0_1_noPBC.xtc -o rmsd_xtal.xvg -tu ns

    # plot the radius of gyration fo the protein to measure it's compactness
    echo "1" | gmx gyrate -s md_0_1.tpr -f md_0_1_noPBC.xtc -o gyrate.xvg
done
