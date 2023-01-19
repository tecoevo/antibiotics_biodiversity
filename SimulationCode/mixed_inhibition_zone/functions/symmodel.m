function F = symmodel(X,N,M,P,S,D,R,K_D,K_P,g,c_prod,c_deg,c_res)
% Symbolic expressions for the evolution equation of the mixed-inhibition zone model. 
% input and output are both rows.

r = zeros(1,N); % contains the growth rates corrected for metabolic costs 

% create symbolic variables for f and F
f = sym('f',[1 N], 'real'); % stores fitnesses
F = sym('F',[1 N], 'real'); % stores the functions on the RHS of the difference equation

for strain = 1:N 
    r(strain) = g(strain) - dot(P(strain,:),c_prod) - dot(D(strain,:),c_deg) - dot(R(strain,:),c_res);
end

for strain = 1:N
    Pkill_i_array = cell(1,M); % length M array of P_kill,ij for each j
    for antibiotic = 1:M 
        %effect of antibiotic j on strain i
        Pkill_ij = S(strain,antibiotic)*(exp(-K_D * dot(D(:,antibiotic) , X)))*(1 - exp(-K_P * dot(P(:,antibiotic) , X)));         
        Pkill_i_array{antibiotic} = Pkill_ij; 
    end
    Pkill_i = InclusionExclusion(Pkill_i_array); % combine effects of all antibiotics using inclusion-exclusion
    f(strain) = r(strain)*(1-Pkill_i);
end

F = (f .* X)/dot(f,X);

end
