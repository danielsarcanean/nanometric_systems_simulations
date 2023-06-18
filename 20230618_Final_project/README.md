# Nanoparticle self-assembly & conjugation modeller

**Summary**

This code aims to produce a Cerium metallic nanoparticle .pdb, by introducing a Montecarlo 'Metropolis-Hastings' algorithm limiting it for spherical dimensions. A better insight can be reviewed at the .pdf proportioned within this folder. 

To run the code, simply arrange your radius and run `python nanoparticle_ensembler.py`. An approximation of 1000 steps ~ 5 minutes can be expected.

The namd scripts can be run by opening a terminal within the `namd_input` folder and writting off:

    namd2 +p$ input_file.inp > output_file.out
    
    # Where +p$ has to be changed to the number of $ logical cores you have

Remember to add the `toppar` folder into the input folder so there is no need to change locations at the inputs. The equilibrium script mantains the topological structure of the box, while the product script runs following a constant piston pressure (NpT calculation). Both are choosen to run at 1000 K. 
