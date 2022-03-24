function [D] = updateD(E,X,Y1,Y2,mu,Z,omega,gamma)
H = (X+Z-E-(Y1+Y2)/mu)/2;
temp_H = H;
[U,sigma,V] = svd(H,'econ');
D = U*(diag(max(diag(sigma)-1/(2*mu),0)))*V';   
temp_H(omega)=D(omega);
D = temp_H;
end
