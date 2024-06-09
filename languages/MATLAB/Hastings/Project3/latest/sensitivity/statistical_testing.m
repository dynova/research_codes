[v1,v2,v3,v4,v5,v6,v7,v8] = deal(cell(7,1));
a = 0.01;
for ind = 1:7
    n = ind + 2;
    NLine(n);
    v1{ind} = NLineC1m(n);
    v2{ind} = NLineC1e(n);
    v3{ind} = NLineC2m(n);
    v4{ind} = NLineC2e(n);
    v5{ind} = NLineC3(n);
    NTriangle(n);
    v6{ind} = NTriangleC1(n);
    v7{ind} = NTriangleC2(n);
    v8{ind} = NTriangleC3(n);
    disp(ind)
end

[x1,x2,x3,x4,x5,x6,x7,x8] = deal(zeros(6,5,3,3));
for ind = 1:6
    for i = 1:5
        for j = 1:3
            for k = 1:3
                t1 = nonzeros(squeeze(v1{ind}(i,j,k,:,:,:)));
                t2 = nonzeros(squeeze(v1{ind+1}(i,j,k,:,:,:)));
                x1(ind,i,j,k) = ttest2(t1,t2,'Vartype','unequal','Alpha',a);
                t3 = nonzeros(squeeze(v2{ind}(i,j,k,:,:,:)));
                t4 = nonzeros(squeeze(v2{ind+1}(i,j,k,:,:,:)));
                x2(ind,i,j,k) = ttest2(t3,t4,'Vartype','unequal','Alpha',a);
                t5 = nonzeros(squeeze(v3{ind}(i,j,k,:,:,:)));
                t6 = nonzeros(squeeze(v3{ind+1}(i,j,k,:,:,:)));
                x3(ind,i,j,k) = ttest2(t5,t6,'Vartype','unequal','Alpha',a);
                t7 = nonzeros(squeeze(v4{ind}(i,j,k,:,:,:)));
                t8 = nonzeros(squeeze(v4{ind+1}(i,j,k,:,:,:)));
                x4(ind,i,j,k) = ttest2(t7,t8,'Vartype','unequal','Alpha',a);
                t9 = nonzeros(squeeze(v5{ind}(i,j,k,:,:,:)));
                t10 = nonzeros(squeeze(v5{ind+1}(i,j,k,:,:,:)));
                x5(ind,i,j,k) = ttest2(t9,t10,'Vartype','unequal','Alpha',a);
                t11 = nonzeros(squeeze(v6{ind}(i,j,k,:,:,:)));
                t12 = nonzeros(squeeze(v6{ind+1}(i,j,k,:,:,:)));
                x6(ind,i,j,k) = ttest2(t11,t12,'Vartype','unequal','Alpha',a);
                t13 = nonzeros(squeeze(v7{ind}(i,j,k,:,:,:)));
                t14 = nonzeros(squeeze(v7{ind+1}(i,j,k,:,:,:)));
                x7(ind,i,j,k) = ttest2(t13,t14,'Vartype','unequal','Alpha',a);
                t15 = nonzeros(squeeze(v8{ind}(i,j,k,:,:,:)));
                t16 = nonzeros(squeeze(v8{ind+1}(i,j,k,:,:,:)));
                x8(ind,i,j,k) = ttest2(t15,t16,'Vartype','unequal','Alpha',a);                
            end
        end
    end
end

save statistical_testing.mat v1 v2 v3 v4 v5 v6 v7 v8 x1 x2 x3 x4 x5 x6 x7 x8;