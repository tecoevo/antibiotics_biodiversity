function trajectory = TimeEvolution(initialpt,N,M,P,S,D,R,K_D,K_P,g,c_prod,c_deg,c_res,timesteps)
% given initialpt, start the system at initialpt and iterate it for a
% duration equal to timesteps

trajectory = [initialpt zeros(N,timesteps-1)];
position = initialpt';

for time = 1:timesteps
    position = model(position,N,M,P,S,D,R,K_D,K_P,g,c_prod,c_deg,c_res);
    trajectory(:,time) = position';
end

end