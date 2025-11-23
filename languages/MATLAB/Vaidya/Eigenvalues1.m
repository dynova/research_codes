syms beta mu mu_H
format compact
J = ([mu+2*mu_H-beta-(mu+mu_H)^2/beta -(mu+mu_H)^2/beta; (beta-mu-mu_H)^2/beta (mu+mu_H-beta)*(mu+mu_H)/beta]);
eigs = eig(J);