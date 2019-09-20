[ prob ]
1: Li C, Kuti JL, Nightingale CH, Nicolau DP. Population pharmacokinetic analysis
and dosing regimen optimization of meropenem in adult patients. J Clin Pharmacol.
2006 Oct;46(10):1171-8. PubMed PMID: 16988206.

https://www.ncbi.nlm.nih.gov/pubmed/16988206

[ set ] delta=0.1, end=8, req=""

[ pkmodel ] cmt = "CENT, PERIPH"  

[ param ] @annotated
WT   : 70 : Weight (kg)
CLCR : 83 : Creatinine clearance (ml/min)
AGE  : 35 : Age (years)
 
[ theta ] @annotated
 1.50E+01 : Typical value of clearance (L/h)
 1.27E+01 : Typical value of volume 1 (L)
 1.52E+01 : Intercompartmental clearance (L/h) 
 1.24E+01 : Typical value of volume 2 (L) 
-4.47E-01 : AGE on CL
 8.20E-01 : WT on V1
 1.88E-01 : Proportional error standard deviation
 4.76E-01 : Additive error standard deviation
 6.20E-01 : CLCR on CL 

[ main ]
double RUV_PROP = THETA7;
double RUV_ADD  = THETA8;

double LOGTWT = log((WT/70.0)); 
  
double LOGTAGE = log((AGE/35.0));
  
double LOGTCLCR = log((CLCR/83.0));
  
double TVCL = log(THETA1) + THETA5 * LOGTAGE + THETA9 * LOGTCLCR;

double CL =  exp(TVCL +  ETA(1)) ;

double TVV1 = log(THETA2) + THETA6 * LOGTWT;
double V1 =  exp(TVV1 +  ETA(2)) ;

double TVVQ = log(THETA3);
double Q =  exp(TVVQ +  ETA(3)) ;

double TVV2 = log(THETA4);
double V2 =  exp(TVV2 +  ETA(4));

[ omega ] @annotated 
ECL : 8.84E-02 : ETA on CL
EV1 : 9.76E-02 : ETA on V1
EQ  : 1.03E-01 : ETA on Q
EV2 : 7.26E-02 : ETA on V2

[sigma ] @annotated
EP : 1 : Not used

[ table ] 
double IPRED = (CENT/V1);
double W = sqrt((RUV_ADD*RUV_ADD)+ (RUV_PROP*RUV_PROP*IPRED*IPRED));
capture Y = IPRED+W*EPS(1);
capture CP = IPRED;

