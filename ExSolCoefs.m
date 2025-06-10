function [C1,C2,C3,eetadel,ealdel,eta]=ExSolCoefs(om0,del,alpha) 
% del is the time step dt
% alpha is the half of the dampling: x''+2*alpha*x'...
eta2=om0^2-alpha^2; eta=sqrt(eta2);
eetadel=exp(1i*eta*del); a=1/eetadel; ealdel=exp(alpha*del);
I1=1i*(a-1)/eta; I2=(a*(1+1i*eta*del)-1)/eta2;
I3=(a*(del*eta*(2+1i*del*eta)-2*1i)+2*1i)/eta2/eta;
C1=(I3-I2*del)/2/del^2/ealdel;   C2=I1-I3/del^2;
C3=ealdel*(I2*del+I3)/2/del^2;
end