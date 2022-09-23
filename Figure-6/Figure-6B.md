# The scripts for Figure-6B
We used Orthofinder to identify ortholog genes among six species and constructed phylogenetic tree.

## Step 1  OrtherFinder installation

conda install -y orthofinder

## Step 2 Run ortherfinder
save protein file of each species in a directory. Be sure to use phased genome, or the single-copy-ortholog-gene is zero
```
orthofinder -f orthofinder -M msa -S diamond -t 16 -a 16
```
The defaul e value for allignment is loose (0.001), please increase it to at least e-05
There is no bootstrp value in the tree constructed by orthofinder.
To update two issues above, please updat config.json
```
vim /PATH/miniconda/envs/orthofinder/bin/scripts_of/config.json
```
line 30 iqtree add -bb 1000
line 43 diamond blastp update -evalue to 0.00001
line 46 blast_gz update -evalue to 0.00001

The ouput director is Results_XXX  XXX represents the running date

## Step 3 Plot phylogenetic tree
view the /Results_XXX/Species_Tree/SpeciesTree_rooted_node_labels.txt in figtree 
For more information, please see http://tree.bio.ed.ac.uk/software/figtree/

The numbers on the tree represented the increased/reduced gene number between two closely related species. 
