# The scripts for Figure-6B
We used Orthofinder to identify ortholog genes among six species and constructed phylogenetic tree.

## Step 1  OrtherFinder installation

conda install -y orthofinder

## Step 2 Run ortherfinder
save protein files of each species in a directory. Be sure to use phased genome, otherwise the single-copy-ortholog-gene is zero
```
orthofinder -f orthofinder -M msa -S diamond -t 16 -a 16
```
The defaut e value for allignment is loose (0.001), please increase it to at least e-05
There is no bootstrp value in the tree constructed by orthofinder.
To update above two issues, please updat config.json
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

The numbers on the tree represent the increased/reduced gene number between two closely related species. 

```
Orthogroup	apsi.p	ppz.p	10822/-11468	pca203.p	7954/-10170	pst104.p	8055/-6524	pgt210.p	7227/-6160
OG0000062	0	7	7	16	9	15	-1	7	-8
OG0000272	0	10	10	29	19	4	-25	18	14
OG0001672	0	0	0	6	6	4	-2	5	1
OG0000817	0	0	0	5	5	4	-1	2	-2
OG0000945	0	1	1	24	23	3	-21	4	1
OG0001219	0	0	0	6	6	3	-3	5	2
OG0001519	0	0	0	5	5	3	-2	1	-2
OG0001200	0	0	0	14	14	2	-12	6	4
OG0001997	0	0	0	12	12	2	-10	2	0
OG0000549	0	0	0	10	10	2	-8	3	1
OG0000814	1	2	1	8	6	2	-6	5	3
OG0000075	0	0	0	8	8	2	-6	1	-1
OG0000173	0	0	0	7	7	2	-5	1	-1
OG0000392	0	0	0	7	7	2	-5	4	2
OG0000480	1	2	1	6	4	2	-4	4	2
OG0000538	3	0	-3	6	6	2	-4	2	0
OG0000548	0	0	0	6	6	2	-4	2	0
OG0000554	0	0	0	5	5	2	-3	1	-1
OG0000950	0	0	0	5	5	2	-3	2	0
OG0001454	0	0	0	5	5	2	-3	3	1
OG0000706	3	1	-2	4	3	2	-2	2	0
OG0001230	0	0	0	4	4	2	-2	2	0
OG0001260	0	0	0	4	4	2	-2	2	0
OG0001240	1	2	1	3	1	2	-1	2	0
OG0001504	1	2	1	3	1	2	-1	2	0

```

-The 4th column represents the increased/reduced gene number of ppz evolved from apsi. The 4th colum = the 3rd colum - the 2nd colum.

-The 6th column represents the increased/reduced gene number of pca203 evovled from ppz. The 6th colum = the 5th colum - the 3rd colum.

-The 8th column represents the increased/reduced gene number of pst104 evovled from pca203. The 8th colum = the 7th colum - the 5th colum.

-The 10th column represents the increased/reduced gene number of pgt210 and ptt76 evovled from pst104. The 10th colum = the 9th colum - the 7th colum.
