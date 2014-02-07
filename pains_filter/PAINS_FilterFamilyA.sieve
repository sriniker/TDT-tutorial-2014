#*************************************************************************************#
#                                                                                     #
# Sieve definitions derived from the PAINS                                            #
#    Created by www.macinchem.org based on http://blog.rguha.net/?p=850               #
#  New Substructure Filters for Removal of Pan Assay Interference Compounds (PAINS)   #
# from Screening Libraries and for Their Exclusion in Bioassays                       #
# FILTER FAMILY A (Freq_Hit_5_morethan150.hits)  									  #
#    DOI: 10.1021/jm901137j                     									  #
#  Tested on sieve v3.0.0                                                             #
#*************************************************************************************#

FRAGMENT	regId=ene_six_het_A(483)	[#6]-1(-[#6](~[!#6&!#1]~[#6]-[!#6&!#1]-[#6]-1=[!#6&!#1])~[!#6&!#1])=[#6;!R]	0	0
FRAGMENT	regId=hzone_phenol_A(479)	c:1:c:c(:c(:c:c:1)-[#6]=[#7]-[#7])-[O;H1]	0	0
FRAGMENT	regId=anil_di_alk_A(478)	[C;H2]N([C;H2])c1cc([$([H]),$([C;H2]),$([O][C;H2][C;H2])])c(N)c([H])c1	0	0
FRAGMENT	regId=indol_3yl_alk(461)	n:1(c(c(c:2:c:1:c:c:c:c:2-[H])-[C;D4]-[H])-[$([C;H2]),$([C]=,:[!C]),$([C;H1][N]),$([C;H1]([C;H2])[N;H1][C;H2]),$([C;H1]([C;H2])[C;H2][N;H1][C;H2])])-[$([H]),$([C;H2])]	0	0
FRAGMENT	regId=quinone_A(370)	[!#6&!#1]=[#6]-1-[#6]=,:[#6]-[#6](=[!#6&!#1])-[#6]=,:[#6]-1	0	0
FRAGMENT	regId=azo_A(324)	[#7;!R]=[#7]	0	0
FRAGMENT	regId=imine_one_A(321)	[#6]-[#6](=[!#6&!#1;!R])-[#6](=[!#6&!#1;!R])-[$([#6]),$([#16](=[#8])=[#8])]	0	0
FRAGMENT	regId=mannich_A(296)	[#7]-[C;X4]-c1ccccc1-[O;H1]	0	0
FRAGMENT	regId=anil_di_alk_B(251)	c:1:c:c(:c:c:c:1-[#7](-[#6;X4])-[#6;X4])-[#6]=[#6]	0	0
FRAGMENT	regId=anil_di_alk_C(246)	c:1:c:c(:c:c:c:1-[#8]-[#6;X4])-[#7](-[#6;X4])-[$([#1]),$([#6;X4])]	0	0
FRAGMENT	regId=ene_rhod_A(235)	[#7]-1-[#6](=[#16])-[#16]-[#6](=[#6])-[#6]-1=[#8]	0	0
FRAGMENT	regId=hzone_phenol_B(215)	c:1(:c:c:c(:c:c:1)-[#6]=[#7]-[#7])-[#8]-[#1]	0	0
FRAGMENT	regId=ene_five_hetA1(201A)	[#6]-1(=[#6])-[#6]=[#7]-[#7,#8,#16]-[#6]-1=[#8]	0	0
FRAGMENT	regId=ene_five_het_A(201)	[#6]-1(=[#6])-[#6]=[#7]-[!#6&!#1]-[#6]-1=[#8]	0	0
FRAGMENT	regId=anil_di_alk_D(198)	c:1:c:c(:c:c:c:1-[#7](-[#6;X4])-[#6;X4])-[#6;X4]-[$([#8]-[#1]),$([#6]=[#6]-[#1]),$([#7]-[#6;X4])]	0	0
FRAGMENT	regId=imine_one_isatin(189)	[#8]=[#6]-2-[#6](=!@[#7]-[#7])-c:1:c:c:c:c:c:1-[#7]-2	0	0
FRAGMENT	regId=anil_di_alk_E(186)	[#6](-[#1])-[#7](-[#6](-[#1])-[#1])-c:1:c(:c(:c(:c(:c:1-[#1])-[$([#1]),$([#6](-[#1])-[#1])])-[#6](-[#1])-[$([#1]),$([#6]-[#1])])-[#1])-[#1]	0	0
