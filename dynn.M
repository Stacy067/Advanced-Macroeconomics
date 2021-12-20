% Dimension of control vector (nc), state and co-state vectors (ns),
% and forcing (exogenous state) vector.

clc
clear
format compact


nc = 2;       % C, N 
ns = 1;        % K
ncs =1;       % lambda 
ne = 2;        % A, G                    

% ECONOMIC PARAMETER VALUES.

% First, the parameters of the production function.  These involve specifying
% factor shares (sk,sn) and the elasticity of substitution (zetakn).
% Recall that the Cobb-Douglas function implies zetakn=1, with sn and
% sk summing to unity
% Notice: zetakn is the inverse of the elasticity of substitution

alpha=0
sn = 0.6820
sk = 1-sn
zetakn = 1
sg=0.1452;
capomega=1;
xikk= -sn*zetakn
xikn=  sn*zetakn
xinn= -sk*zetakn
xink=  sk*zetakn

% and the growth rate of technical change and depreciation

gammax=1.0083
deltak=0.0175

% Second, the parameters of preferences.
% These are parameterized to match the three degrees of freedom implied
% by the steady state requirements.  Let sigma be the absolute value of the
% elasticity of marginal utility of consumption, it then follows that
% specifying sigma, i.e.,

r =0.0156
bstar=beta*(gammax^(1-sc/(sc+alpha*sg)));
si=(bstar*sk*capomega)/(gammax-bstar*(1-deltak))/(gammax/(gammax-(1-deltak)));
sc=(1-si-sg);
beta = (gammax^(sc/(sc+alpha*sg)))/(1+r)

bstar=beta*(gammax^(1+xicc)
sigma = sc/(sc+alpha*sg)


% implies a specification of two of the key elasticities of preferences.
xicc=-sigma
xilc=1-sigma

% The other two elasticities are free, so long as xicl=0 if sigma=1

xicl= 0 
xill= 0.4

% The steady state level of hours is 

nbar=0.3493

% We choose beta by requiring that there is a plasible value for the steady
% state interest rate



% Third, the parameters of government policy

omegak=0;
omegan=0;
omegaa=0;
omegag=0;

% STEADY-STATE CALCULATIONS
% The following steps calculate details of the steady-state of
% the transformed economy.

% First, the transformed economy's rate of time preference.

bstar=beta*(gammax^(1+xicc));

% Second, the ratio phi and the capital-output ratio are

phi=gammax/(gammax-(1-deltak));

kyratio=(bstar*sk*capomega)/(gammax-bstar*(1-deltak));

% Third, the steady-state consumption and investment ratios

si=kyratio*gammax/phi;
sc=(1-si-sg);


% Fourth, the steady-state elasticities of the aftertax marginal
% product of capital

etaa=(gammax-bstar*(1-deltak))*(omegaa+1)/gammax;
etak=(gammax-bstar*(1-deltak))*(omegak+xikk)/gammax;
etan=(gammax-bstar*(1-deltak))*(omegan+xikn)/gammax;
etag=(gammax-bstar*(1-deltak))*(omegag)/gammax;



% BASIC SYSTEM MATRICES:

% First, specify the matrices in control equations

Mcc=zeros(nc,nc)
%  C  N

Mcc= [ xicc  -xicl*nbar/(1-nbar)
       xilc  (-xinn-omegan)-xill*nbar/(1-nbar)];

Mcs = zeros(nc,ns+ncs);
% K lamda

Mcs = [ 0               1
       (xink+omegak)    1];

Mce = zeros(nc,ne);
% A G

Mce = [ 0         0
       (1+omegaa) omegag];

% Define matrices in state equations

Mss0 = zeros(ns+ncs,ns+ncs);
Mss1 = zeros(ns+ncs,ns+ncs);
% K lamda

Mss0 = [etak   1
        si*phi 0];
Mss1 = [   0           -1
        -sk-si*(phi-1)  0];

Msc0 = zeros(ns+ncs,nc);
Msc1 = zeros(ns+ncs,nc);
% C N

Msc0 = [0 -etan
        0     0];

Msc1 = [ 0    0
        -sc   sn];

Mse0 = zeros(ns+ncs,ne);
Mse1 = zeros(ns+ncs,ne);
% A G
Mse0 = [-etaa -etag
        0       0];

Mse1 = [ 0     0
         1   -sg];


% Some additional matrices links other observable flow variables (nf) to
% controls, states, costates and exogneous forcing variables.  
% In this application these are output, wages, and investment.  In particular,
% (A37-39) can be written as linear combinations of a stacked vector of states, controls,
% and forcing processes. (Shadow prices--or costates--do not enter in our case,
% but are included for completeness). The matrix formulation of these linear
% relations is as follows:

nf=3;

FVc = zeros(nf,nc);

FVc = [   0      sn   
          0      xinn+omegan 
       -(sc/si)  sn/si      ];

FVke= zeros(nf,ns+ne);

FVke = [sk              1              0
        xink+omegak     1+omegaa       0
        sk/si           1/si        -sg/si];

FVl = zeros(nf,ncs);

FVl = [0
       0
       0];


% This completes the specification of the model.

% At this point, we do some housekeeping

clear sn; clear sk; clear zetakn; 
clear xikk; clear xikn; clear xink; clear xinn;
clear gammax; clear deltak; clear sigma;
clear xicc; clear xilc; clear xicl; clear xill;
clear nbar; clear r; clear beta;
clear sg; clear capomega;
clear omegak; clear omegag; clear omegaa; clear omegan;
clear bstar; clear phi; clear kyratio;
clear si; clear sc;
clear etaa; clear etak; clear etan; clear etag;


% Some additional matrices links other observable flow variables (nf) to
% controls, states, costates and exogneous forcing variables.  
% In this application these are output, wages, and investment.  In particular,
% (A37-39) can be written as linear combinations of a stacked vector of states, controls,
% and forcing processes. (Shadow prices--or costates--do not enter in our case,
% but are included for completeness). The matrix formulation of these linear
% relations is as follows:
 



% FUNDAMENTAL STATE-COSTATE DIFFERENCE EQUATION


% Derive fundamental difference equation for states and co-states
% The memnonic MS is used to indicate the M* matrices in the text
% Notice that the symbol \ involves a particular concept of matrix
% division (see PCMatlab manual, p. 2-18) in which X=A\B is the solution
% to A*X=B.

MSss0 = Mss0 - Msc0*(Mcc\Mcs);

MSss1 = Mss1 - Msc1*(Mcc\Mcs);


MSse0 = Mse0 + Msc0*(Mcc\Mce);

MSse1 = Mse1 + Msc1*(Mcc\Mce);

%pause

% Put the fundamental difference equation in normal form

W = -(MSss0\MSss1);
R = MSss0\MSse0;
Q = MSss0\MSse1;

% Housekeeping

clear Mss0; clear Mss1; clear Msc0; clear Msc1; clear Mse0; clear Mse1;
clear MSss0; clear MSss1; clear MSse0; clear MSse1;


% EIGENVECTOR-EIGENVALUE DECOMPOSITION OF STATE TRANSITION MATRIX:

% Compute the eigenvalues and eigenvectors of the matrix W.

[P,MU] = eig(W);
MU

% Order the eigenvalues, in order of ascending absolute value.

AMU=abs(MU);

% Sort by the absolute value of the eigenvectors, then reorder
% the columns of the eigenvector matrix (P) by this sort. (See
% the discussion of sorting in the PCMATLAB manual, p. 3-113,
% where they use this example).

[MU,k] = sort(diag(AMU));
P=P(:,k);

% Next, we need to order the eigenvalue (diagonal) matrix in
% order of ascending absolute value.  As it was unclear how to
% accomplish this, the matrix is simply computed from the standard
% diagonalization formula given the sorted eigenvalue matrix and
% its inverse, denoted PS (matches the technical appendix notation,
% where the inverse of P is denoted with a star).

PS=inv(P);
MU=PS*W*P;


% PARTITIONING MATRICES

% The P, PS, MU, R and Q matrices must be partitioned, as discussed in
% the technical appendix section B.2, for the purpose of computing the
% different components of the solution.  Partitioning
% is undertaken with the colon notation in PCMATLAB, as discussed in
% the manual in the section starting p.2-31.  In view of potential
% multiple state variable problems, the partitioning uses ns as a
% parameter.

MU1 = MU(1:ns,1:ns);
MU2 = MU(ns+1:ns+ncs,ns+1:ns+ncs);

P11 = P(1:ns,1:ns);
P12 = P(1:ns,ns+1:ns+ncs);
P21 = P(ns+1:ns+ncs,1:ns);
P22 = P(ns+1:ns+ncs,ns+1:ns+ncs);


PS11 = PS(1:ns,1:ns);
PS12 = PS(1:ns,ns+1:ns+ncs);
PS21 = PS(ns+1:ns+ncs,1:ns);
PS22 = PS(ns+1:ns+ncs,ns+1:ns+ncs);


Rke = R(1:ns,1:ne);
Rle = R(ns+1:ns+ncs,1:ne);
Qke = Q(1:ns,1:ne);
Qle = Q(ns+1:ns+ncs,1:ne);

% COMPOSITE EXPRESSIONS

% The solution for the shadow price involves some composite expressions.
% These are:

SP1 = - MU2\(PS21*Rke+PS22*Rle);
SP2 = - MU2\(PS21*Qke+PS22*Qle);

% The solution for the capital stock involves some composite expressions.
% These are the following, where the memnonic is that the KLK gives the
% response of capital to lagged capital and the matrix KTL gives the
% response to the transformed shadow price.

KLK = P11*MU1/P11;
KTL = (P11*MU1*PS12+P12*MU2*PS22)/PS22;

% Note that these expressions involve the use of matrix right division
% in the multiple state variable case (PCMATLAB manual, p. 2-18,19).
% Thus, the formulas appear slightly different in these expressions
% than in the text.

% This terminates the program dyn.m:  It is necessary at this stage 
% to make decisions concerning the type of experiment that will be studied.
% The output of this program must be

% nc,ns,ne,nf
% Mcc,Mcs,Mce
% FVc,FVke,FVl
% W,R,Q
% MU2,PS21,PS22
% Rke,Qke,
% KTL,KEC,KLK
% SP1,SP2

% Thus, we eliminate all other material from the workspace

%clear P; clear PS; clear MU; clear AMU;
%clear MU1; clear PS11; clear PS12;
%clear P11; clear P12; clear P21; clear P22;
%clear Rle; clear Qle;
%clear k;

% and then save the worskpace contents on disk.  This enables multiple uses
% of the outputs of this program, since subsequent programs do housekeeping.

format

