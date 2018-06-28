% Calculates bottom boundary and side lithostatic pressures

% At first time step, calculates lithostatic pressure
if t==1
    p_lith = zeros([Ny,2]);
    k      = [2, Nx];
    for j=1:2
        for i=2:Ny
            p_lith(i,j) = p_lith(i-1,j) + rho(i,k(j))*PARAMS.gy*dy(i-1);
        end
        p_lith(2:end,j) = p_lith(1:end-1,j);
    end
    BC.p_lith_bot = mean(p_lith(end,:));
    BC.p_lith_top = 0;
else
    % Otherwise, calculate lithostatic pressure from the bottom up and top
    % down
    p_lith       = zeros([Ny,2]);
    p_lith(2,:)  = BC.p_lith_top;
    p_lith(Ny,:) = BC.p_lith_bot;
    k      = [2, Nx];
    
    for j=1:2
        for i=2:BC.sticky_bot_left+5
            p_lith(i,j) = p_lith(i-1,j) + rho(i,k(j))*PARAMS.gy*dy(i-1);
        end
        p_lith(2:BC.sticky_bot_left+5,j) = p_lith(1:BC.sticky_bot_left+4,j);
        l = fliplr(BC.lith_bot_left-5:Ny-1);
        for ll=1:length(l)
            i = l(ll);
            p_lith(i,j) = p_lith(i+1,j) - rho(i,k(j))*PARAMS.gy*dy(i-1);
        end
          p_lith(2:end,j) = p_lith(1:end-1,j);
    end
end

BC.p_lith = p_lith;
