PRE_SEQ=GCTTCATATAATCCGAATGATATGGTTTCGGAGCTTCCACCAAGAGCCTTAAACTCTTGATTATGAAGTCTGTCGCTTTATCCGAAATTTTATAAAGAGAAGACT
LENGTH=24
TAG=MGSSHH
POST_SEQ=CATCATCATCACAGCAGCGGCCTGGTGCCGCGCGGCAGCCAT
OUTPUT_FILE=results.xls
sudo pip install synbiochem-py
python riboswitch_predicter.py $PRE_SEQ $LENGTH $TAG $POST_SEQ $OUTPUT_FILE