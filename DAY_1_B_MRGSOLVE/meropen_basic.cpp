[ prob ]
1: Li C, Kuti JL, Nightingale CH, Nicolau DP. Population pharmacokinetic analysis
and dosing regimen optimization of meropenem in adult patients. J Clin Pharmacol.
2006 Oct;46(10):1171-8. PubMed PMID: 16988206.

https://www.ncbi.nlm.nih.gov/pubmed/16988206

[ set ] delta=0.1, end=8

[ pkmodel ] cmt = "CENT, PERIPH"  

[ param ] @annotated

CL : 1.50E+01 : Typical value of clearance (L/h)
V1 : 1.27E+01 : Typical value of volume 1 (L)
Q  : 1.52E+01 : Intercompartmental clearance (L/h) 
V2 : 1.24E+01 : Typical value of volume 2 (L) 

[ table ] 
capture CP = (CENT/V1);
