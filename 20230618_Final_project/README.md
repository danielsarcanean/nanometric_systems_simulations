![title]([https://github.com/[username]/[reponame]/blob/[branch]/image.jpg?raw=true](https://github.com/danielsarcanean/nanometric_systems_simulations/blob/628ca907267843b6a270ae1e339f29867114c158/20230618_Final_project/Visual_outputs/title.png))

**Summary**

This code aims to produce a Cerium metallic nanoparticle .pdb, by introducing a Montecarlo 'Metropolis-Hastings' algorithm limiting it to spherical dimensions. Better insight can be reviewed at the .pdf proportioned within this folder. 

To run the code, simply arrange your radius and run `python nanoparticle_ensembler.py`. An approximation of 1000 steps ~ 5 minutes can be expected.

The namd scripts can be run by opening a terminal within the `namd_input` folder and writing off:

    namd2 +p$ input_file.inp > output_file.out
    
    # Where +p$ has to be changed to the number of $ logical cores you have

Remember to add the `toppar` folder into the input folder so there is no need to change locations at the inputs. The equilibrium script maintains the topological structure of the box, while the product script runs following a constant piston pressure (NpT calculation). Both are chosen to run at 1000 K. 

**Basics**

The fundamentals of this code cover a Lennard-Jones approach for atomic arrangement:

$$ V (r) = 4 \cdot \varepsilon \left( \left( \frac {\sigma}{r} \right) ^{12} - \left(  \frac {\sigma}{r} \right) ^{6} \right)$$

Where the values $\varepsilon$ and $\sigma$ come from the combinatory calculations described by [Elliott and Akerson, 2015](https://openkim.org/id/MO_959249795837_003). The discriminator for the Montecarlo simulation can be simplified to:

$$ \sum_{ijk} ^N V (r) = \sum_{ijk} V_{Lennard-Jones}(r) - \sum_{ijk} V_{cristal} (r) $$

And then the resulting $r$ value is confirmed following a Boltzmann probability factor:

$$ f (E) = \frac{1}{e^{\frac{E}{k_B T}}} $$ 

Hence the distribution image comes from both the physical meaning of the distances and a realistic displacement produced due to evironmental conditions of the system.
