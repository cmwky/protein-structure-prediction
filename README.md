# protein-structure-prediction
A branch-and-bound A.I. approach to predicting tertiary protein structure from simulated sparse pairwise NMR distance data.

Determining the 3-dimensional structure of a protein is a significant problem that still poses a major challenge. 
The ability to determine the spatial orientation of a protein is highly desired, as it offers great insight into its functionality. 
Current mainstream experimental methods include X-Ray Crystallography and NMR Spectroscopy. 
Using NMR results, the primary objective is to obtain a solution such that the distances between atoms within the predicted, computed structure, are as close as possible to the experimental distances; this is referred to as the "Molecular Distance Geometry Problem". 
Once a true structure is acquired, it is then uploaded to an online protein database (PDB).

I first began by collecting statistical data from the PDB in order to more accurately characterize bond lengths between various atoms within the target protein. 
From there, a "Branch-and-Prune" A.I. algorithm is implemented that considers various physical and chemical assumptions when constructing a potential solution. 
Without these considerations (i.e., by generating a purely mathematical solution) the computational time is exponential, however, by necessitating the validation of these assumptions, the algorithm can decide when to stop constructing a potential solution and discard it. 
Once a family of potential structures is obtained, the task then becomes identifying the "true" solution. 
This is to say that, even though a structure may be mathematically correct, it still may not be the true structure of the target protein (i.e., when considering a biological context).

As of now, the algorithm only parses the backbone of the target protein, and will output a list of x,y,z-coordinates for roughly 30-35 atoms in PDB file format. 
Currently I'm working on optimizing some of the data structures used to store the potential positions for each atom, which will dramatically improve computational time and should allow an entire backbone to be determined.

See https://onlinelibrary.wiley.com/doi/full/10.1111/j.1475-3995.2007.00622.x for supprting documentation on the physics and mathematics behind this algorithm.
