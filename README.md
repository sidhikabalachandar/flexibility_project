# flexibility_project

Docking:
Contains code to conduct all ligand to all structure docking for a group of proteins

  Calculations:
  dock_rmsd_delete.py -> code to conduct all ligand to all structure docking
                         if run_dock flag is provided, docking is conducted (see lxpowers33/docking)
                          produces a folder each protein and a subfolder for each ligand_to_structure pair
                          each subfolder contains a .in file, .log file, .rept file, and a .csv file of the rmsds
                          the pose viewer file is deleted
                         if check flag is provided, code checks if all ligand_to_structure pair subfolders have been created
                         if results flag is provided, the rmsds of each pair is condensed into a pickled dictionary
                          format of dictionary is {PROTEIN : { LIGAND : { STRUCTURE : [RMSD_LIST] } } }
                          
  aggregate_ligand_poses_MAPK14.py -> code to aggregate all poses into a single .mae file
  
  Analysis:
  docking_evaluation.ipynb -> contains graphs for glide benchmark by docking type, protein, and protein group
  
  best_pose_without_4DLI_with_threshold -> contains 2D color graph of RMSDs clustered by ligands
  
  docking_analysis.ipynb -> contains basic accuracy metrics on MAPK14 glide all ligand to all structure docking
  
Flexibility prediction:
Contains code to crate X and Y data for the train and test sets
Also contains code to run regression based prediction algorithm 

  Calculations:
  residue_feature_vector_creator.py -> obtains features dependent on only the protein structure
                                       features collected: residue name, residue num, raw bfactor, normalized bfactor, 
                                        previous residue's normalized bfactor, previous previous residue's normalized bfactor, 
                                        next residue's normalized bfactor, next next residue's normalized bfactor, molecular 
                                        weight, residue type's general number of rotamers, residue type's general average rmsd 
                                        between different rotamers, residue's specific number of non conflicting rotamers, 
                                        residue's specific average rmsd between different nonconflicting rotamers, solvent 
                                        accessibility, secondary structure
                                        
  rmsd_calculator.py -> obtains features dependent on both the protein structure and ligand
                        calculates the complete rmsd, backbone rmsd, and sidechain rmsd (sidechain is to be improved) of each 
                          individual residue in the binding pocket
                        additional features collected: ligand similarity (maximum common substructure), ligand similarity 
                          ratio, ligand size difference (difference in number of heavy atoms), ligand size ratio, 
                          
  Analysis:
  all_protein_flexibility_prediction_analysis.ipynb -> graph of correlation between bfactor and residue flexibility
                                                       graph of average number of total residues and average number of      
                                                        flexible residues in the binding pocket
                                                        
  cycling_test_set.ipynb -> most up-to-date prediction algorithm (code decomposed in prediction.py)
                            cycles through each protein as a test set and all of the rest of the proteins as the training set
                            uses linear regression, polynomial regression, and various decision tree prediction models
                            calculates the accuracy, recall, and precision using a cutoff of 2 Angstroms as flexible
                            also contains code to plot color bars to compare prediction outputs of a certain model for a 
                              certain protein, ligand_to_structure pair
                              
  flexibility_prediction_analysis.ipynb -> graph of MAPK14 bfactor distributions before and after normalization
  
Mutations:
Contains code to mutate residues to alanine (residues in the binding pocket), create .zip files, and then rerun docking to 
obtain RMSD results

  Calculations:
  flexible_mutation.py -> mutates all residues that have a flexibility RMSD above a certain cutoff
  
  conflict_mutation.py -> mutates all residues that conflict with the target ligand
  
  MAPK14_prediction_mutation.py -> obtains the RMSDs as predicted by the cycling_test_set.ipynb code
                                   mutates at most 6 residues with RMSDs above a cutoff
                                   
  zipper.py -> produces the .zip files from the .mae structure files
               .zip file is necessary for docking
               
  dock_rmsd_delete.py -> same as above (reruns docking)
  
  Analysis:
  mutated_rmsds_analysis.ipynb -> creates graph of cumulative frequency to portray docking accuracy
  
  mutation_analysis.ipynb -> calculates the average number of residues mutated
  
Protein flexibility:
Contains code to calculate the rmsd of all of the residues in the binding pocket between each pair of structures for a given 
protein

  Calculations:
  sequence_finder.py -> obtains the amino acid sequence of a given protein structure
  
  pairwise_alignment.py	-> obtains two aligned strings for each pair of starting structure and the structure corresponding to 
                            the target ligand
                            
  rmsd_calculator.py -> calculates the rmsd of all of the residues in the binding pocket between each pair of structures for a 
                          given protein
                          
  Analysis:
  rmsd_analysis.ipynb -> contains graph of rmsd of all of the residues in the binding pocket between each pair vs docking 
                          performance
                          
Similarity:
Contains code to calculate the maximum common substructure and tanimoto coefficient between every pair of ligands for a given 
protein

  Calculations:
  mcss_similarity.py -> creates a csv file containing the number of atoms in each ligand and the mcss of the pair of ligands
  
  similarity_vs_performance.py -> calculates the tanimoto coefficient between every pair of ligands
  
  Analysis:
  similarity_analysis.ipynb -> graph of tanimoto coefficient vs docking performance
