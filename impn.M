% Matlab File: impn.m

% This is a prototype program for computing dynamic equiblibria
% in PC-Matlab. The economy analyzed is that discussed in appendix
% A of KPR "Production, Growth and Business Cycles," namely the
% neoclassical model with variable labor. The dimension of the
% state vector of controlled states is one in the present
% application, but we let it be ns to test more general procedures. 
% The dimension of the control vector is two, but we parameterize it
% by nc for similar reasons.  The dimension of the vector of
% uncontrolled states is two, but we denote it ne.  The base for
% this program and subsequent programs are the King/Rebelo program 
% #2 (Matlab file neovn.m). This program has been edited--all pauses
% were removed and only essential output is now printed.

% This program requires as input the output of Program #2, namely 
% the markov decision rules as summarized by Mke and H. 


% IMPULSE RESPONSES

% The calculation of impulse reponses involves two elements.  First, it is 
% necessary to specify initial conditions for the controlled state(k). Second,
% one needs to specify an values of the exogenous variables, i.e.,
% an (ns+ne x 1) vector of initial conditions. Consider just
% a total factor productivity shock from the steady-state position.

S = zeros(3,1);
answ1 = menu('Choose the impulse','K','A','G');
answ2 = input('Choose magnitude of impulse ');

S(answ1) = answ2;
clear answ1 answ2

% Specify the number of periods for which the impulse response is to calculated

nir = 30;

% Initialize the impulse response matrix

IR=zeros(nir,ns+ne+ncs+nc+nf+ncs+1);

% Calculate the impulse reponse for the nir periods.

for i=1:nir;
  %Index of  Impulse Response Period;
    IR(i,1)=i;
  %State and Exogneous Variables;
    IR(i,2:ns+ne+1)=S';
  %Shadow Price, Controls and Other Flows;
    IR(i,ns+ne+2:ns+ne+ncs+nc+nf+ncs+1)=(H*S)';
  % Update the state vector:;
    S = Mke*S;
end 
%'Impulse Response to Shock in Productivity'
%'Ordered as Follows'
%'Date K A G Lambda C N Y Y/N I R'

IR(1:nir,1:ns+ne+ncs+1);
IR(1:nir,ncs+ns+ne+2:ns+ne+ncs+nc+nf+ncs+1);



