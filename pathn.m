%  pathn.m
%  This program can be used after dyn.m to compute perfect foresight paths

%  To run this program ncs has to be defined and MU1 MU2 the partitions
%  of R, Q, P and PS cannot be erased from dyn


T = 20;   %  Number of periods during which the shocks vary
	  %  The Tth value has to be the new steady state value

TT = 30;  %  Number of periods for which paths are to be computed

e=zeros(ne,T);


% specify the shocks

for i=1:T-1
	e(1,i) = .9^(i-1);
end


k0 = zeros(ns,1);               % Initial Capital stock
 
f0 = 0;

% Fill in shocks

e = [e e(:,T)*ones(1,TT-T)];

c1 = PS21*Rke+PS22*Rle;
c2 = PS21*Qke+PS22*Qle;


for t = 1:T
	f0 = f0 + (MU2^(-t+1))*(c1*e(:,t+1) + c2*e(:,t));
end

f0 = -inv(MU2)*f0 - inv(eye(size(MU2))-inv(MU2))*(MU2^(-T-1))*(c1+c2)*e(:,T);
lamtil0 = f0;
lam0 = inv(PS22)*(lamtil0 - PS21*k0);

klam(1:ns+ncs,1) = [k0 
		   lam0];
f(:,1) = f0;

M1 = [P11*MU1*PS11   P11*MU1*PS12
     P21*MU1*PS11   P21*MU1*PS12];

M2 = [P11*(PS11*Rke+PS12*Rle)
     P21*(PS11*Rke+PS12*Rle)];


M3 = [P11*(PS11*Qke+PS12*Qle)
     P21*(PS11*Qke+PS12*Qle)];

M4 = [P12
      P22];


for t = 1:TT-1

	f(:,t+1) = MU2*f(:,t) + (c1*e(:,t+1) + c2*e(:,t));
	klam(:,t+1) = M1*klam(:,t) + M2*e(:,t+1) + M3*e(:,t) + M4*f(:,t+1);

end

%clear M1 M2 M3 M4 c1 c2 

controls = inv(Mcc)*(Mcs*klam + Mce*e);

ke = [klam(1:ns,:)
        e];

flows = FVc*controls + FVke*ke + FVl*klam(ns+1:ns+ncs,:);

%  Include variables in a single IR matrix so that the usual programs 
%  to plot impulse response functions can be used

nir = TT;
IR = [(1:TT)' klam(1:ns,:)' e' klam(ns+1:ns+ncs,:)' controls' flows' NaN*ones(TT,ncs)];
disp('paths computed; return to see the impulse responses')
pause
plotn
