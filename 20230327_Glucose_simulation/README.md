1. Launch the first program, arranging your molecule.
2. Create a 'conda' environment in a terminal at the same folder you have your program using the following sequence:
```
conda create -n NAME
conda activate NAME
conda install siesta
```
3. Run your ``.fdf`` file by this sequence: ``siesta RUN.fdf``
4. Now you can launch any of the programs that form the second part. Use an external software like VMD to visualize the molecular orbitals.

