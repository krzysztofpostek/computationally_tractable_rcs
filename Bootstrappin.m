% Bootstrap test

counter = 0;
violated = 0;

for i = 1:length(Sample_indices)

            Implemented_power = (x'*Diagrams + (Error_matrix(:,Sample_indices(i)).*x)'*Diagrams);

            Satisfies = prod(Implemented_power'<=[tau*ones(length(Indices_class_1),1) ; ones(length(Indices_class_2) + length(Indices_class_3),1) ])*prod( Implemented_power' >= [-tau*ones(length(Indices_class_1),1) ; -ones(length(Indices_class_2),1);  0.9*ones(length(Indices_class_3),1) ] );

            violated = violated + (1-Satisfies);
            
end