r = 0.0156;g=0.0083;d=0.0175;a=0.6820;
lbda1_1=(1+r)/(1+g);
lbda1_2=a*(r+d)/((1+g)*(1-a));
lbda1_3=a*(r+d)/(1+r);

sig=[0;1;5;10000000]; phi=[0;0.5;0.95;1];
Q1 = zeros(4,1);
Q0 = zeros(4,1);

eta_ck_1 = zeros(4,4);
eta_kk_1 = zeros(4,4);
eta_ca_1 = zeros(4,4);
eta_ka_1 = zeros(4,4);

for i=1:4
for j=1:4
Q1(j,1)=lbda1_1-1+sig(j,1)*lbda1_3*(1-lbda1_1-lbda1_2);
Q2=1-lbda1_1-lbda1_2;
Q0(j,1)=sig(j,1)*lbda1_3*lbda1_1;

eta_ck_1(i,j) = (-Q1(j,1)-sqrt(Q1(j,1)^2-4*Q2*Q0(j,1)))/(2*Q2);
eta_kk_1(i,j) = lbda1_1-(Q1(j,1)+sqrt(Q1(j,1)^2-4*Q2*Q0(j,1)))/2;
eta_ca_1(i,j) = (-eta_ck_1(i,j)*lbda1_2+sig(j,1)*lbda1_3*(phi(i,1)-lbda1_2))/(phi(i,1)-1+Q2*(eta_ck_1(i,j)+sig(j,1)*lbda1_3));
eta_ka_1(i,j) = lbda1_2+(1-lbda1_1-lbda1_2)*eta_ca_1(i,j);
end
end

tab_1=zeros(2,2);
for i=1:4
for j=1:4
tab_1=[eta_ck_1(i,j) eta_ca_1(i,j) eta_kk_1(i,j) eta_ka_1(i,j)];
end
end

T_1 = zeros(4,1);
for i=1:4
T_1(i,1) = log(1/2)/log(eta_kk_1(1,i));
end



r = 0.0156;g=0.0083;d=0.0175;a=0.3640;
lbda2_1=(1+r)/(1+g);
lbda2_2=a*(r+d)/((1+g)*(1-a));
lbda2_3=a*(r+d)/(1+r);

sig=[0;1;5;10000000]; phi=[0;0.5;0.95;1];
Q1 = zeros(4,1);
Q0 = zeros(4,1);

eta_ck_2 = zeros(4,4);
eta_kk_2 = zeros(4,4);
eta_ca_2 = zeros(4,4);
eta_ka_2 = zeros(4,4);

for i=1:4
for j=1:4
Q1(j,1)=lbda2_1-1+sig(j,1)*lbda2_3*(1-lbda2_1-lbda2_2);
Q2=1-lbda2_1-lbda2_2;
Q0(j,1)=sig(j,1)*lbda2_3*lbda2_1;

eta_ck_2(i,j) = (-Q1(j,1)-sqrt(Q1(j,1)^2-4*Q2*Q0(j,1)))/(2*Q2);
eta_kk_2(i,j) = lbda2_1-(Q1(j,1)+sqrt(Q1(j,1)^2-4*Q2*Q0(j,1)))/2;
eta_ca_2(i,j) = (-eta_ck_2(i,j)*lbda2_2+sig(j,1)*lbda2_3*(phi(i,1)-lbda2_2))/(phi(i,1)-1+Q2*(eta_ck_2(i,j)+sig(j,1)*lbda2_3));
eta_ka_2(i,j) = lbda2_2+(1-lbda2_1-lbda2_2)*eta_ca_2(i,j);

end
end

tab_2=zeros(2,2);
for i=1:4
for j=1:4
tab_2=[eta_ck_2(i,j) eta_ca_2(i,j) eta_kk_2(i,j) eta_ka_2(i,j)];
end
end

T_2 = zeros(4,1);
for i=1:4
T_2(i,1) = log(1/2)/log(eta_kk_2(1,i));
end

