function [E] = updateE(D,E,X,Y1,mu,alpha,omega)
G = X - D -Y1/mu;
temp_G=G;
E = max(abs(G)-alpha/mu,0).*sign(G);
temp_G(omega)=E(omega);
E=temp_G;
end
