This is a program that determines hydrogen-bonding affinity from amide-groups across a peptide backbone between all other amide groups on the peptide. This is done by calculating the distance 
between the hydrogen associated with the nitrogen in the amide group and all other amide group oxygens. This matrix is boltzmann weighted by relative potential energy accross 500-700 conformers
of each peptide, and is used as an input for a random forest machine learning model. This model attemps to categorize a peptide's likelihood of cyclization based upon amide-amide distances 
(potential hydrogen bonds). We hypothesize that the Boltzmann-averaged secondary structure of the linearized peptide could reveal important information about its cyclization propensity.
