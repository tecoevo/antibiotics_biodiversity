function F = model(X,N,M,P,S,D,R,K_D,K_P,g,c_prod,c_deg,c_res)
% this function does the same thing as symmodel, with the only difference being
% that it does not use symbolic variables. this is because symbolic computation 
% takes a significantly larger time. input and output are both rows. 

r = zeros(1,N); % contains the growth rates corrected for metabolic costs 

for strain = 1:N 
    r(strain) = g(strain) - dot(P(strain,:),c_prod) - dot(D(strain,:),c_deg) - dot(R(strain,:),c_res);
end

% The matrix P_kill contains the probabilities P_kill(i,j) of strain i being killed by the antibiotic j.
% The array f contains the fitnesses which are used to evolve the population abundances to the next time step.

P_kill = zeros(N,M);
f = zeros(1,N); % fitnesses

for strain = 1:N
    for antibiotic = 1:M 
        P_kill(strain,antibiotic) = S(strain,antibiotic)*exp(-K_D * dot(D(:,antibiotic) , X))*(1 - exp(-K_P * dot(P(:,antibiotic) , X)));         
    end
    Pkill_i = InclusionExclusion(num2cell(P_kill(strain,:)));
    f(strain) = r(strain)*(1-Pkill_i);
end

F = (f .* X)/dot(f,X);

end