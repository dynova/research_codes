folder = '/Users/abhishekmallela/Desktop/Programming/Mathematica/Hastings/files';
N = 25;
[array_b1, array_b2] = deal(linspace(0.01,0.99,2)*N);
array_d = linspace(0,1,2)*N;
len_b1 = length(array_b1);
len_b2 = length(array_b2);
len_d = length(array_d);
len_t = 4;
[surv, fptd] = deal(cell(len_d,len_b1,len_b2,len_t));
syms t;
for ind_d = 1:len_d
    for ind_b1 = 1:len_b1
        for ind_b2 = 1:len_b2
            for ind_t = 1:len_t
                surv{ind_d,ind_b1,ind_b2,ind_t} = eval(fileread(strcat(folder,...
                    ['/F' num2str(ind_d) num2str(ind_b1)...
                    num2str(ind_b2) num2str(ind_t) '.txt'])));
                fptd{ind_d,ind_b1,ind_b2,ind_t} = -diff(surv{ind_d,ind_b1,ind_b2,ind_t});
            end
        end
    end
end

save survival.mat surv fptd