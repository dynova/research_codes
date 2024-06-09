syms beta tau lambda mu mu_H mu_T
expr1 = 4*beta*lambda*mu*tau+4*beta*lambda*tau*mu_H+12*beta*mu*tau*mu_H-2*beta^3*lambda-2*beta^3*mu+2*beta^3*tau-4*beta^3*mu_H+4*beta^2*mu_H^2+4*mu^2*tau^2+4*tau^2*mu_H^2+4*beta^2*lambda*mu_H-6*beta^2*mu*tau+4*beta^2*mu*mu_H-8*beta^2*tau*mu_H+4*beta*mu^2*tau-4*beta*mu*tau^2-4*beta*tau^2*mu_H+8*beta*tau*mu_H^2+8*mu*tau^2*mu_H+4*beta*(beta+tau)*mu*tau+4*beta*(beta+tau)*tau*mu_H-4*mu*(beta+tau)*mu*tau-4*mu*(beta+tau)*tau*mu_H-4*mu_H*(beta+tau)*mu*tau-4*mu_H*(beta+tau)*tau*mu_H+beta^2*lambda^2+beta^2*mu^2+beta^2*tau^2+beta^4+2*beta^2*lambda*mu-2*beta^2*lambda*tau;
factor(expr1);
expr2 = beta^2*(mu+mu_T+lambda+beta+mu_T-tau)^2-expr1;
complicated = factor(expr2)/(-4*beta);
manuscript = beta*(tau-mu-mu_T)*(beta-mu-mu_H)+(beta+tau)*(lambda+mu_H)*(mu+mu_H-beta)...
            +beta*(mu+mu_T+lambda+mu_H)*(tau-mu-mu_T);
disp(factor(complicated-manuscript))