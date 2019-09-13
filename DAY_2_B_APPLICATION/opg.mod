[ prob ]
1: Zierhut ML, Gastonguay MR, Martin SW, Vicini P, Bekker PJ, Holloway D, Leese
PT, Peterson MC. Population PK-PD model for Fc-osteoprotegerin in healthy
postmenopausal women. J Pharmacokinet Pharmacodyn. 2008 Aug;35(4):379-99. 
doi: 10.1007/s10928-008-9093-5. Epub 2008 Jul 17. PubMed PMID: 18633695.

[ set ] req = "NTX", end = 14*24, rtol=1e-5

[GLOBAL]
#define CP (CENT/(VC/1000000.0))

[ cmt ] @annotated
SC   : Subcutaneous dosing compartment (mg)
CENT : Central compartment (mg)
P1   : First peripheral compartment (mg)
P2   : Second peripheral compartment (mg)
NTX  : Urinary N-telopeptide

[ param ]  @annotated
IV : 0 : IV dose indicator 
 
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
THETA11 : 0.864  : Biomarker synthesis rate (.)
THETA12 : 0.0204 : Biomarker elimination rate constant (1/h)
THETA13 : 5.38   : Half-maximal inhibitory conc (ng/ml)

[ main ]
double CL   = exp(log(THETA1) + ECL);
double VC   = exp(log(THETA2) + EVC);
double VP1  = exp(log(THETA3) + EVP1);
double VP2  = exp(log(THETA4) + EVP2);
double Q1   = exp(log(THETA5) + EQ1);
double Q2   = THETA6;
double KA   = exp(log(THETA7) + EKA);
double VMAX = THETA8;
double KM   = THETA9;
double FSC  = exp(log(THETA10) + EFSC);
double KSYN = exp(log(THETA11) + EKSYN);
double KDEG = exp(log(THETA12) + EKDEG);
double IC50 = exp(log(THETA13) + EIC50);

NTX_0 = KSYN/KDEG;

F_SC = FSC/(1.0+FSC);

[ omega ] @annotated
ECL  : 0.0391 : IIV on CL
EVC  : 0.0102 : IIV on VC
EVP1 : 0.0144 : IIV on VP1
EVP2 : 0.0333 : IIV on VP2
EQ1  : 0.0379 : IIV on Q1
EKA  : 0.0457 : IIV on KA
EFSC : 0.263  : IIV on FSC

[ omega ] @block @annotated
EKSYN : 0.281               : IIV on KSYN
EKDEG : 0.0867 0.0325       : IIV on KDEG
EIC50 : 0.0000 0.0000  1.18 : IIV on IC50
    
[ sigma ] @annotated 
ADDIV : 0.0193 : Additive error IV dose
ADDSC : 0.7330 : Additive error SC dose 
   
[ sigma ] @annotated
PDPROP : 0.0407 : Proportional error NTX
PDADD  :   20.7 : Additive error NTX

[ ode ]
double CLNL = VMAX/(CP+KM);
dxdt_SC     = -KA*SC;
dxdt_CENT   =  KA*SC - (CL+Q1+Q2+CLNL)*CENT/VC + Q1*P1/VP1 + Q2*P2/VP2;
dxdt_P1     =  CENT*Q1/VC - P1*Q1/VP1;
dxdt_P2     =  CENT*Q2/VC - P2*Q2/VP2;
dxdt_NTX    =  KSYN*(1.0 - CP/(IC50+CP)) - KDEG*NTX;

[ table ]
double IPRED = CP;
double PKEPS = IV==1 ? ADDIV : ADDSC;
double PKDV = exp(log(IPRED)+PKEPS);
double PDDV = NTX*(1+PDPROP) + PDADD;

[ capture ] @annotated
IPRED : Individual-predicted concenrtation (ng/mL)
PKDV  : Fc-OPG serum concentration (ng/mL)
PDDV  : NTX (nM BCE/mM creatinine per hr)
