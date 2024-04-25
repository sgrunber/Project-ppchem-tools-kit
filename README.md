# ppchem

#rdkit.Chem.rdMolDescriptors.CalcFractionCSP3((Mol)mol) → float :
returns the fraction of C atoms that are SP3 hybridized
#The fraction of sp3 hybridised carbon atoms is thought to give a good indicator of the “planarity” of a molecule, with those having low fractions being fairly “planar” (2-dimensional) and so, possibly, less good candidates for drugs.
#Ertl, P. and Schuffenhauer A. “Estimation of Synthetic Accessibility Score of Drug-like Molecules based on Molecular
Complexity and Fragment Contributions” Journal of Cheminformatics 1:8 (2009)
#“Natural Product Likeness Score and Its Application for Prioritization of Compound Libraries” Peter Ertl, Silvio Roggo, and Ansgar Schuffenhauer Journal of Chemical Information and Modeling 48:68-74 (2008) http://pubs.acs.org/doi/abs/10.1021/ci700286x