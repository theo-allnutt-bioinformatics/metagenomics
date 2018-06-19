heatmap.R

Theo Allnutt
School of Medicine
Deakin University
2016

An R script to draw a heatmap with row and column dendrograms and a column grouping label. Data format is
tab delimited table of abundance data, e.g.

otu	F	G	H	I	J	Y	T	U	V	W	X	O
group	1	1	1	1	1	1	2	2	2	2	2	2
f:Fusobacteriaceae	26	39	98	70	224	26	32	10642	11200	2005	18602	41
g:Leuconostoc	1653	638	2306	2405	2090	3286	164	27	114	110	22	32
f:Enterobacteriaceae	96	23	121	195	162	585	5533	5188	7	9	4	338
o:Bacillales	678	171	4387	2900	3177	0	0	2	0	90	4	0
d:Bacteria	56	6	26	72	50	119	2	8	289	2413	39	1109
f:Planococcaceae	985	413	2289	12	37	0	0	0	0	51	1	0
g:Lactococcus	387	136	389	839	416	1507	31	13	13	27	5	24
g:Clostridium_sensu_stricto	16	0	16	21	1	11	52	45	508	429	2322	25
g:Weissella	42	15	56	89	57	2704	115	21	28	36	6	25

First row is sample ids.
Second row is a group identifier for groups of samples.
Following rows are otu labels followed by abundances for each sample.

Edit line 26 to add groups if required, example is only two groups.

Output is a .pdf 
Edit line 1 if Rscript is not in /usr/bin/Rscript
Uncomment lines 2 and 3 if you need to install packages.

Useage:

heatmap.R input_data.txt output.pdf 25

where '25' is the number of otus from the top of the table to draw.
Sort your data first to obtain, e.g. 25 most abundant.
