# 2016-04-08
# batch nucleosome calls with DANPOS

DPOS_EXEC=/Users/benayoun/Softwares/danpos-2.2.2/

LIVER=/Volumes/MyBook_3/BD_aging_project/ChIP-seq/New_H3_lanes/Liver_H3_reseq/Nucleosomes/
CEREB=/Volumes/MyBook_3/BD_aging_project/ChIP-seq/New_H3_lanes/Cerebellum_H3_reseq/Nucleosomes
HEART=/Volumes/MyBook_3/BD_aging_project/ChIP-seq/New_H3_lanes/Heart_H3_reseq/Nucleosomes 
NPC=/Volumes/MyBook_3/BD_aging_project/ChIP-seq/New_H3_lanes/NPC_5-6_H3_reseq/Nucleosomes 
OB=/Volumes/MyBook_3/BD_aging_project/ChIP-seq/New_H3_lanes/OB_H3_reseq/Nucleosomes 

ln -s $LIVER Liver
ln -s $CEREB Cerebellum
ln -s $HEART Heart
ln -s $NPC NPCs
ln -s $OB OlfactoryBulb


python2.7 $DPOS_EXEC/danpos.py dpos Cerebellum/Cerebellum_29m/:Cerebellum/Cerebellum_3m,Cerebellum/Cerebellum_12m/:Cerebellum/Cerebellum_3m/,Cerebellum/Cerebellum_29m/:Cerebellum/Cerebellum_12m/ -o Cerebellum_H3_positionning_aging_1e-15_2016-04-14 -t 1e-15
python2.7 $DPOS_EXEC/danpos.py dpos OlfactoryBulb/OB_29m/:OlfactoryBulb/OB_3m,OlfactoryBulb/OB_12m/:OlfactoryBulb/OB_3m/,OlfactoryBulb/OB_29m/:OlfactoryBulb/OB_12m/ -o OB_H3_positionning_aging_1e-15_2016-04-14 -t 1e-15
python2.7 $DPOS_EXEC/danpos.py dpos Heart/Heart_29m/:Heart/Heart_3m,Heart/Heart_12m/:Heart/Heart_3m/,Heart/Heart_29m/:Heart/Heart_12m/ -o Heart_H3_positionning_aging_1e-15_2016-04-14 -t 1e-15
python2.7 $DPOS_EXEC/danpos.py dpos Liver/Liver_29m/:Liver/Liver_3m,Liver/Liver_12m/:Liver/Liver_3m/,Liver/Liver_29m/:Liver/Liver_12m/ -o Liver_H3_positionning_aging_1e-15_2016-04-14 -t 1e-15
python2.7 $DPOS_EXEC/danpos.py dpos NPCs/NPCs_29m/:NPCs/NPCs_3m,NPCs/NPCs_12m/:NPCs/NPCs_3m/,NPCs/NPCs_29m/:NPCs/NPCs_12m/ -o NPCs_H3_positionning_aging_1e-15_2016-04-14 -t 1e-15
