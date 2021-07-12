% Bootstrap test

N_sample_bootstrap_outer = 1000;
N_sample_bootstrap_inner = 100;
counter = 0;
violated = 0;

for i=1:length(Error)

            Implemented_power = (x'*Diagrams + (Error(:,i).*x)'*Diagrams);

            Satisfies = prod(Implemented_power'<=[tau*ones(length(Indices_class_1),1) ; ones(length(Indices_class_2) + length(Indices_class_3),1) ])*prod( Implemented_power' >= [-tau*ones(length(Indices_class_1),1) ; -ones(length(Indices_class_2),1);  0.9*ones(length(Indices_class_3),1) ] )

            violated = violated + (1-Satisfies);
            
end