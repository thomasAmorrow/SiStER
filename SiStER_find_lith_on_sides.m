% function [BC] = SiStER_find_lith_on_sides(xm,ym,x,y,BC,etam,t)
if BC.mech_lith==0
    ym_left             = ym(xm>x(2)&xm<(x(3)-diff(x(2:3))/2));
    phase2_bot_left     = max(ym_left(im(xm>x(2)&xm<(x(3)-diff(x(2:3))/2))==2));
    BC.lith_bot_left  = find((y-phase2_bot_left)<0,1,'last');
    
    
    % Right side
    ym_right          = ym(xm>x(end-1)&xm<(x(end)-diff(x(end-1:end))/2));
    phase2_bot_right  = max(ym_right(im(xm>x(end-1)&xm<(x(end)-diff(x(end-1:end))/2))==2));
    BC.lith_bot_right = find((y-phase2_bot_right)<0,1,'last');
elseif BC.mech_lith==1
    % Method 1: Use viscosity to find bottom of mechanical lithosphere
    if ~isempty(etam(etam>=BC.min_lith_visc))
        if t==1
            % Find depth where viscosity drops below BC.min_lith_visc on sides
            % Left side
            lith_bot_left     = max(ym(etam==max(etam(:))));
            BC.lith_bot_left  = find((y-lith_bot_left)>0,1,'first');
            
            % Right side
            lith_bot_right    = max(ym(etam==max(etam(:))));
            BC.lith_bot_right = find((y-lith_bot_right)>0,1,'first')+1;
        else
            % Find depth where viscosity drops below BC.min_lith_visc on sides
            % Left side
            ym_left           = ym(xm>x(2)&xm<(x(3)-diff(x(2:3))/2));
            lith_bot_left     = max(ym_left(etam(xm>x(2)&xm<(x(3)-diff(x(2:3))/2))>=BC.min_lith_visc));
            BC.lith_bot_left  = find((y-lith_bot_left)>0,1,'first');
            
            % Right side
            ym_right          = ym(xm>x(end-1)&xm<(x(end)-diff(x(end-1:end))/2));
            lith_bot_right    = max(ym_right(etam(xm>x(end-1)&xm<(x(end)-diff(x(end-1:end))/2))>=BC.min_lith_visc));
            BC.lith_bot_right = find((y-lith_bot_right)>0,1,'first');
        end
    else % sometimes in the first timestep before visc. converge, this is needed
        BC.lith_bot_left  = length(y);
        BC.lith_bot_right = length(y);
    end
end

if BC.sliding_lith(1)==1
    BC.lith_bot_left = BC.lith_bot_left+BC.sliding_lith(2);
    BC.lith_bot_left(BC.lith_bot_left>length(y))=length(y);
    
    BC.lith_bot_right = BC.lith_bot_right+BC.sliding_lith(2);
    BC.lith_bot_right(BC.lith_bot_right>length(y))=length(y);
end
% end