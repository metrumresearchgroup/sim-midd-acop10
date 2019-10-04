[ prob ]
1: Zierhut ML, Gastonguay MR, Martin SW, Vicini P, Bekker PJ, Holloway D, Leese
PT, Peterson MC. Population PK-PD model for Fc-osteoprotegerin in healthy
postmenopausal women. J Pharmacokinet Pharmacodyn. 2008 Aug;35(4):379-99. 
doi: 10.1007/s10928-008-9093-5. Epub 2008 Jul 17. PubMed PMID: 18633695.
  
[ set ] req = "CP", end = 14*24, rtol=1e-5
  
[GLOBAL]
#define CP (CENT/(VC/1000000.0))
  
[ cmt ] @annotated
SC   : Subcutaneous dosing compartment (mg)
CENT : Central compartment (mg)
P1   : First peripheral compartment (mg)
P2   : Second peripheral compartment (mg)
AUC  : Accumulation compartment (ng*hour/mL)
    
[ param ]  @annotated
IV : 0 : IV dose indicator 
WT : 70 : Body weight (kg)
    
[ param ] @annotated
THETA1  : 168    : Clearance (ml/h)
THETA2  : 2800   : Central volume (ml)
THETA3  : 443    : Volume of first peripheral cmt (ml)
THETA4  : 269    : Volume of second peripheral cmt (ml)
THETA5  : 15.5   : Distribution clearance (ml/h)
THETA6  : 3.02   : Distribution clearance (ml/h)
THETA7  : 0.0131 : Absorption rate constant (1/h)
THETA8  : 13300  : Maximum velocity (ng/h)
THETA9  : 6.74   : Michaelis constant (ng/ml)
THETA10 : 0.0719 : Bioavailability of SC dose (.)
    
[ main ]
double CL   = exp(log(THETA1) + ECL)*pow(WT/70,0.75);
double VC   = exp(log(THETA2) + EVC)*(WT/70);
double VP1  = exp(log(THETA3) + EVP1)*(WT/70);
double VP2  = exp(log(THETA4) + EVP2)*(WT/70);
double Q1   = exp(log(THETA5) + EQ1)*pow(WT/70,0.75);
double Q2   = THETA6*pow(WT/70,0.75);
double KA   = exp(log(THETA7) + EKA);
double VMAX = THETA8*pow(WT/70,0.75);
double KM   = THETA9;
double FSC  = exp(log(THETA10) + EFSC);
  
F_SC = FSC/(1.0+FSC);
  
[ omega ] @annotated
ECL  : 0.0391 : IIV on CL
EVC  : 0.0102 : IIV on VC
EVP1 : 0.0144 : IIV on VP1
EVP2 : 0.0333 : IIV on VP2
EQ1  : 0.0379 : IIV on Q1
EKA  : 0.0457 : IIV on KA
EFSC : 0.263  : IIV on FSC
    
[ sigma ] @annotated 
ADDIV : 0.0193 : Additive error IV dose
ADDSC : 0.7330 : Additive error SC dose 
    
[ ode ]
double CLNL = VMAX/(CP+KM);
dxdt_SC     = -KA*SC;
dxdt_CENT   =  KA*SC - (CL+Q1+Q2+CLNL)*CENT/VC + Q1*P1/VP1 + Q2*P2/VP2;
dxdt_P1     =  CENT*Q1/VC - P1*Q1/VP1;
dxdt_P2     =  CENT*Q2/VC - P2*Q2/VP2;
dxdt_AUC    =  CP;
  
[ table ]
double IPRED = CP;
double PKEPS = IV==1 ? ADDIV : ADDSC;
double PKDV = exp(log(IPRED)+PKEPS);
  
[ capture ] @annotated
IPRED : Individual-predicted concenrtation (ng/mL)
PKDV  : Fc-OPG serum concentration (ng/mL)
AUC   : Accumulated concentration (ng*hour/mL)