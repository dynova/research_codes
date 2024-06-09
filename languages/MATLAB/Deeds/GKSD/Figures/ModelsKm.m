rvec = logspace(-2,4,1000);
[KMGK1,KMGK2] = Models_vs_Km(0,1e3,0,0,1e-4); 
[KMGK3,KMGK4] = Models_vs_Km(0,1e3,0,0,1e-2);
[KMI1,KMI2] = Models_vs_Km(2e-2,1e3,2e-5,2e-5,1e-4);
[KMI3,KMI4] = Models_vs_Km(2e-2,1e3,2e-5,2e-5,1e-2);
[KMSD1,KMSD2] = Models_vs_Km(2e-2,1e3,2e-5,2e-4,1e-4);
[KMSD3,KMSD4] = Models_vs_Km(2e-2,1e3,2e-5,2e-4,1e-2);
save 2b