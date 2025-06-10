function [x,xd]=OneStep(x,xd,alpha,eta,eetaDel,ealDel,C1,C2,C3,sprev,s,snew)
A=x-1i*(xd+alpha*x)/eta;
A=A-1i*(C1*sprev+C2*s+C3*snew)/eta;
d=A*eetaDel;
y=real(d); yd=1i*eta*(d-conj(d))/2;
x=y/ealDel; xd=(yd-alpha*y)/ealDel;
end