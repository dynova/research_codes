rvec = logspace(-2,4,1000);
[QGK1,QGK2] = Models_vs_Q(0,1e3,0,0,1e-4); 
[QGK3,QGK4] = Models_vs_Q(0,1e5,0,0,1e-4);
[QI1,QI2] = Models_vs_Q(2e-2,1e3,2e-5,2e-5,1e-4);
[QI3,QI4] = Models_vs_Q(2e-0,1e5,2e-5,2e-5,1e-4);
[QSD1,QSD2] = Models_vs_Q(2e-2,1e3,2e-5,2e-4,1e-4);
[QSD3,QSD4] = Models_vs_Q(2e-0,1e5,2e-5,2e-4,1e-4);
save 2a