
Results=[];

for i=0:20    
    mu=1+0.01/20*i;
    solver_file;
    solver_file_nonrobust;
    
    Results(i+1,1)=value_robust;
    Results(i+1,2)=value_nonrobust;    
end