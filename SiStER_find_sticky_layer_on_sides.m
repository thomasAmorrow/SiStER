% Calculates the depth of the phase1-phase2 transition on the
% boundaries. Gives the last element containing only phase1.

% Find depth of bottom of phase 1   
% Left side
phase1_bot_left     = topo_y(1);
BC.sticky_bot_left  = find((y-phase1_bot_left)<0,1,'last');


% Right side
phase1_bot_right    = topo_y(end);
BC.sticky_bot_right = find((y-phase1_bot_right)<0,1,'last');

if BC.sliding_lith(1)==1
   BC.sticky_bot_left = BC.sticky_bot_left-BC.sliding_lith(2);
   BC.sticky_bot_left(BC.sticky_bot_left<1)=1; 

   BC.sticky_bot_right = BC.sticky_bot_right-BC.sliding_lith(2);
   BC.sticky_bot_right(BC.sticky_bot_right<1)=1; 
end
