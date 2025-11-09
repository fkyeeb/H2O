#################################################################################################
############################################# ERIC ##############################################
#################################################################################################

cd ...scripts/phylomod2
python prune_tips.py

h2o infer_ortho -t ../../0_homolog_trees/ERIC -of ERIC_outgroup.txt -e .subtree > infer_ortho.log
# 10325 homologs in 1m51s to 8768 unpruned orthologs
h2o map_dupl -t processed_trees -s ERIC_ASTRAL_rooted_unpruned.tre
# 8768 unpruned orthologs, 11 sec
h2o extract_wgd_trees -t processed_trees -n 3,4

cat processed_trees/unpruned/*ortho*.tre > other_output/ASTRAL_in_unpruned.tre
cat processed_trees/pruned/*ortho*.tre > other_output/ASTRAL_in_pruned.tre

cd other_output
astral4 -i ASTRAL_in_unpruned.tre -o ASTRAL_out_unpruned.tre -t 40
# 8768
astral4 -i ASTRAL_in_pruned.tre -o ASTRAL_out_pruned.tre -t 40
# 5269
astral4 -i ASTRAL_in_unpruned_wgd.tre -o ASTRAL_out_unpruned_wgd.tre -t 40
# 3346
astral4 -i ASTRAL_in_pruned_wgd.tre -o ASTRAL_out_pruned_wgd.tre -t 40
# 1686

pxrr -t ASTRAL_out_unpruned.tre -f ../ERIC_outgroup.txt -o ASTRAL_rooted_unpruned.tre
pxrr -t ASTRAL_out_pruned.tre -f ../ERIC_outgroup.txt -o ASTRAL_rooted_pruned.tre
pxrr -t ASTRAL_out_unpruned_wgd.tre -f ../ERIC_outgroup.txt -o ASTRAL_rooted_unpruned_wgd.tre
pxrr -t ASTRAL_out_pruned_wgd.tre -f ../ERIC_outgroup.txt -o ASTRAL_rooted_pruned_wgd.tre

nohup bp -c ASTRAL_rooted_unpruned.tre -t ASTRAL_in_unpruned.tre -v -tv -w 30 > bp_output_unpruned.txt 2> bp.log &
# 1h34m trees read: 8768 edges read: 346255
nohup bp -c ASTRAL_rooted_pruned.tre -t ASTRAL_in_pruned.tre -v -tv -w 30 > bp_output_pruned.txt 2> bp_pruned.log &
# 28m trees read: 5269 edges read: 235020
bp -c ASTRAL_rooted_unpruned_wgd.tre -t ASTRAL_in_unpruned_wgd.tre -v -tv -w 40 > bp_output_unpruned_wgd.txt
# 6m trees read: 3346 edges read: 99259
bp -c ASTRAL_rooted_pruned_wgd.tre -t ASTRAL_in_pruned_wgd.tre -v -tv -w 40 > bp_output_pruned_wgd.txt
# 24s trees read: 1686 edges read: 30767

h2o bp2pie -f bp_output_unpruned.txt
h2o bp2pie -f bp_output_pruned.txt
h2o bp2pie -f bp_output_unpruned_wgd.txt
h2o bp2pie -f bp_output_pruned_wgd.txt

# /Users/fky/Desktop/PhyloMod/h2o_results/ERIC
gokstad -s -d -b -pie gokstad_pie_unpruned.tre -o unpruned.svg
gokstad -s -d -b -pie gokstad_pie_pruned.tre -o pruned.svg
gokstad -s -d -b -pie gokstad_pie_unpruned_wgd.tre -o unpruned_wgd.svg
gokstad -s -d -b -pie gokstad_pie_pruned_wgd.tre -o pruned_wgd.svg

h2o map_dupl -t processed_trees -s other_output/ASTRAL_rooted_pruned_wgd.tre
h2o constraint -s summary_tree_numbered.tre -n 142,2,Ericales_Sapotaceae_Sarcosperma_laurinum,Ericales_Sapotaceae_Diploknema_yunnanensis,122,119,105,88,12
astral4 -o ASTRAL_out_constraint.tre -c constraint_tree.tre -i ASTRAL_in_unpruned.tre -t 40
pxrr -t ASTRAL_out_constraint.tre -f ../ERIC_outgroup.txt -o ASTRAL_rooted_constraint.tre
h2o map_dupl -t processed_trees -s other_output/ASTRAL_rooted_constraint.tre
h2o gene_loss -t processed_trees -n 3,4

# duplication counts -  seems like the same event just sometimes both lost balsamoniods
node3 Counter({0: 2389, 1: 140, 2: 35, 3: 13, 4: 9, 5: 3, 6: 2, 14: 1, 8: 1, 7: 1, 9: 1})
node4 Counter({0: 2235, 1: 329, 2: 17, 3: 5, 4: 2, 5: 2, 13: 1, 11: 1, 8: 1, 18: 1, 10: 1})
overlap 5

############### ##############################astral-pro3 ##############################################
cd /home/microway/keyi/homolog2ortholog/4_ASTRAL-Pro_runs/ERIC
cat ../../0_homolog_trees/ERIC/*.subtree > homologs.tre
astral-pro3 -R -i homologs.tre -o ASTRAL-pro_out.tre -u 2 -t 50 2>test.log
# 6 min -R more rounds
cd /Users/fky/Desktop/PhyloMod/h2o_results
pxrr -t ERIC/ASTRAL-pro_out.tre -f ERIC/ERIC_outgroup.txt -o ERIC/ASTRAL-pro_rooted.tre
python read_astral_pro.py
# change colors in ~/gokstad/src/gokstad_conf.py
gokstad -s -d -b -pie ERIC/ASTRAL-pro_pie.tre -o ERIC/ASTRAL-pro.svg
# q2 #E69F00 q3 #F0E442 q1 #648fff

############################################# yang&smith ##############################################
cd /home/microway/keyi/homolog2ortholog/5_yang_smith
python3 prune_paralogs_MO.py ../0_homolog_trees/ERIC .subtree 5 ERIC/ERIC_outgroup.txt ERIC/MO
# about 4 min
python3 prune_paralogs_RT.py ../0_homolog_trees/ERIC .subtree ERIC/RT 5 ERIC/ERIC_taxa.txt
# about 19 min
cat ERIC/MO/*.tre > ERIC/MO_astral_in.tre
# 1769
cat ERIC/RT/*.tre > ERIC/RT_astral_in.tre
# 28609

astral4 -i ERIC/MO_astral_in.tre -o ERIC/MO_astral_out.tre -t 50    
astral4 -i ERIC/RT_astral_in.tre -o ERIC/RT_astral_out.tre -t 50

pxrr -t ERIC/MO_astral_out.tre -f ERIC/ERIC_outgroup.txt -o ERIC/MO_astral_rooted.tre
pxrr -t ERIC/RT_astral_out.tre -f ERIC/ERIC_outgroup.txt -o ERIC/RT_astral_rooted.tre

cd ERIC
bp -c MO_astral_rooted.tre -t MO_astral_in.tre -v -tv -w 30 > MO_bp_output.txt 2> MO_bp.log
# 11m25s trees read: 1769 edges read: 159657
nohup bp -c RT_astral_rooted.tre -t RT_astral_in.tre -rng 0-14300 -v -tv -w 30 > RT_bp_output1.txt 2> RT_bp1.log &
# about 1h13m trees read: 14300 edges read: 369768
nohup bp -c RT_astral_rooted.tre -t RT_astral_in.tre -rng 14301-28609 -v -tv -w 30 > RT_bp_output2.txt 2> RT_bp2.log &
# about 1hr25m trees read: 14308 edges read: 383101
h2o bp2pie -f MO_bp_output.txt -n MO
h2o bp2pie -f RT_bp_output1.txt,RT_bp_output2.txt -n RT

# /Users/fky/Desktop/PhyloMod/h2o_results/ERIC
gokstad -s -d -b -pie gokstad_pie_MO.tre -o MO.svg
gokstad -s -d -b -pie gokstad_pie_RT.tre -o RT.svg

#################################################################################################
############################################# AMAR ##############################################
#################################################################################################

# ~/keyi/homolog2ortholog/scripts/phylomod2
python rename_tips.py 
# ~/keyi/homolog2ortholog/3_h2o_runs/AMAR
# h2o infer_ortho -t ../../0_homolog_trees/AMAR_test -of AMAR_outgroup.txt -e .tm
h2o infer_ortho -t ../../0_homolog_trees/AMAR -of AMAR_outgroup.txt -e .tm > infer_ortho.log
# 14583 homologs in 2m, 11934 processed, 23682 unpruned orthologs, 18298 pruned orthologs

cat processed_trees/unpruned/*ortho*.tre > other_output/ASTRAL_in_unpruned.tre
cat processed_trees/pruned/*ortho*.tre > other_output/ASTRAL_in_pruned.tre

cd other_output
astral4 -i ASTRAL_in_unpruned.tre -o ASTRAL_out_unpruned.tre -t 50
astral4 -i ASTRAL_in_pruned.tre -o ASTRAL_out_pruned.tre -t 50
pxrr -t ASTRAL_out_unpruned.tre -f ../AMAR_outgroup.txt -o ASTRAL_rooted_unpruned.tre
pxrr -t ASTRAL_out_pruned.tre -f ../AMAR_outgroup.txt -o ASTRAL_rooted_pruned.tre

cd ..
h2o map_dupl -t processed_trees -s other_output/ASTRAL_rooted_unpruned.tre
cd other_output
nohup bp -c ASTRAL_rooted_unpruned.tre -t ASTRAL_in_unpruned.tre -v -tv -w 30 > bp_output_unpruned.txt 2> bp.log &
# 3h13m trees read: 23682 edges read: 913512
nohup bp -c ASTRAL_rooted_pruned.tre -t ASTRAL_in_pruned.tre -v -tv -w 30 > bp_output_pruned.txt 2> bp_pruned.log &
# 2h32m trees read: 18298 edges read: 861055

# only node 3 is WGD with high conflict
cd ..
h2o extract_wgd_trees -t processed_trees -n 3
cd other_output
astral4 -i ASTRAL_in_unpruned_wgd.tre -o ASTRAL_out_unpruned_wgd.tre -t 30
# 11412
astral4 -i ASTRAL_in_pruned_wgd.tre -o ASTRAL_out_pruned_wgd.tre -t 30
# 8107
pxrr -t ASTRAL_out_unpruned_wgd.tre -f ../AMAR_outgroup.txt -o ASTRAL_rooted_unpruned_wgd.tre
pxrr -t ASTRAL_out_pruned_wgd.tre -f ../AMAR_outgroup.txt -o ASTRAL_rooted_pruned_wgd.tre

nohup bp -c ASTRAL_rooted_unpruned_wgd.tre -t ASTRAL_in_unpruned_wgd.tre -v -tv -w 10 > bp_output_unpruned_wgd.txt 2> bp_unpruned_wgd.log &
# 17m trees read: 11412 edges read: 326246
nohup bp -c ASTRAL_rooted_pruned_wgd.tre -t ASTRAL_in_pruned_wgd.tre -v -tv -w 10 > bp_output_pruned_wgd.txt 2> bp_pruned_wgd.log &
# 12m trees read: 8107 edges read: 284007

h2o bp2pie -f bp_output_unpruned.txt
h2o bp2pie -f bp_output_pruned.txt -n pruned
h2o bp2pie -f bp_output_unpruned_wgd.txt -n wgd
h2o bp2pie -f bp_output_pruned_wgd.txt -n pruned_wgd

# ~/Desktop/PhyloMod/h2o_results/AMAR
gokstad -s -d -b -pie gokstad_pie.tre -o unpruned.svg
gokstad -s -d -b -pie gokstad_pie_pruned.tre -o pruned.svg
gokstad -s -d -b -pie gokstad_pie_wgd.tre -o unpruned_wgd.svg
gokstad -s -d -b -pie gokstad_pie_pruned_wgd.tre -o pruned_wgd.svg

# duplication counts, seems like two separate WGDs
node3 Counter({0: 8231, 1: 2357, 2: 1058, 3: 254, 4: 24, 5: 5, 7: 3, 6: 1, 11: 1}) - could be a triplication
node24 Counter({0: 9950, 1: 1872, 2: 95, 3: 14, 4: 2, 11: 1})
overlap 1257

# topology the same, more support

############### ##############################astral-pro3 ##############################################
cd /home/microway/keyi/homolog2ortholog/4_ASTRAL-Pro_runs/AMAR
cat ../../0_homolog_trees/AMAR/*.tm > homologs.tre
nohup astral-pro3 -R -i homologs.tre -o ASTRAL-pro_out.tre -u 2 -t 50 2> astral-pro.log &
# about 6 min -R
cd /Users/fky/Desktop/PhyloMod/h2o_results
pxrr -t AMAR/ASTRAL-pro_out.tre -f AMAR/AMAR_outgroup.txt -o AMAR/ASTRAL-pro_rooted.tre
# change folder name
python read_astral_pro.py
# change colors in ~/gokstad/src/gokstad_conf.py
gokstad -s -d -b -pie AMAR/ASTRAL-pro_pie.tre -o AMAR/ASTRAL-pro.svg
# q2 #E69F00 q3 #F0E442 q1 #648fff

############################################# yang&smith ##############################################
cd /home/microway/keyi/homolog2ortholog/5_yang_smith
python3 prune_paralogs_MO.py ../0_homolog_trees/AMAR .tm 5 AMAR/AMAR_outgroup.txt AMAR/MO
# about 2 min
python3 prune_paralogs_RT.py ../0_homolog_trees/AMAR .tm AMAR/RT 5 AMAR/AMAR_taxa.txt
# about 5 min
cat AMAR/MO/*.tre > AMAR/MO_astral_in.tre
# 9071
cat AMAR/RT/*.tre > AMAR/RT_astral_in.tre
# 28489

astral4 -i AMAR/MO_astral_in.tre -o AMAR/MO_astral_out.tre -t 50
astral4 -i AMAR/RT_astral_in.tre -o AMAR/RT_astral_out.tre -t 50

pxrr -t AMAR/MO_astral_out.tre -g Polygonaceae_Fagopyrum_tataricum -o AMAR/MO_astral_rooted.tre
pxrr -t AMAR/RT_astral_out.tre -g Polygonaceae_Fagopyrum_tataricum -o AMAR/RT_astral_rooted.tre

cd AMAR
nohup bp -c MO_astral_rooted.tre -t MO_astral_in.tre -v -tv -w 30 > MO_bp_output.txt 2> MO_bp.log &
# about 1.5hr trees read: 9071 edges read: 646993
nohup bp -c RT_astral_rooted.tre -t RT_astral_in.tre -rng 0-14245 -v -tv -w 30 > RT_bp_output1.txt 2> RT_bp1.log &
# about 32m trees read: 14245 edges read: 419415
nohup bp -c RT_astral_rooted.tre -t RT_astral_in.tre -rng 14246-28489 -v -tv -w 30 > RT_bp_output2.txt 2> RT_bp2.log &
# about 1hr8m trees read: 14243 edges read: 554872

h2o bp2pie -f MO_bp_output.txt -n MO
h2o bp2pie -f RT_bp_output1.txt,RT_bp_output2.txt -n RT

# /Users/fky/Desktop/PhyloMod/h2o_results/AMAR
gokstad -s -d -b -pie gokstad_pie_MO.tre -o MO.svg
gokstad -s -d -b -pie gokstad_pie_RT.tre -o RT.svg