function [im] = SiStER_initialize_marker_phases(Nphase,GEOM,xm,ym)

% assign material identity on markers
im=zeros(size(xm));

if GEOM(2).fzSW==0
    for kk = 1:Nphase
        
        if GEOM(kk).type==1 % layer
            
            im(ym>=GEOM(kk).top & ym<GEOM(kk).bot)=kk;
            
        elseif GEOM(kk).type==2 % circular inclusion
            
            rm=sqrt((xm-GEOM(kk).x0).^2 + (ym-GEOM(kk).y0).^2);
            im(rm<=GEOM(kk).rad)=kk;
            
        end
        
    end
end
 
if GEOM(2).fzSW==1
    for kk = 1:Nphase
        
        if GEOM(kk).type==1 % layer
            
            im(ym>=GEOM(kk).top & ym<GEOM(kk).bot)=kk;
            
            if kk==2 % FZ step added to layer kk=2
                if length(GEOM(kk).fz)>1
                    for bb=1:length(GEOM(kk).fz)
                        if bb==1
                            im(xm<=GEOM(kk).fz(bb)*(max(xm)) & ym>=GEOM(kk).top-GEOM(kk).step(bb) & ym<GEOM(kk).bot)=kk;
                        else
                            im(xm<=GEOM(kk).fz(bb)*(max(xm)) & xm>GEOM(kk).fz(bb-1)*(max(xm)) & ym>=GEOM(kk).top-GEOM(kk).step(bb) & ym<GEOM(kk).bot)=kk;
                        end
                    end
                                     
                else
                    im(xm<=GEOM(kk).fz(1)*(max(xm)) & ym>=GEOM(kk).top-GEOM(kk).step(1) & ym<GEOM(kk).bot)=kk;
                end
            end
            
            
            
        elseif GEOM(kk).type==2 % circular inclusion
            
            rm=sqrt((xm-GEOM(kk).x0).^2 + (ym-GEOM(kk).y0).^2);
            im(rm<=GEOM(kk).rad)=kk;
            
%         elseif GEOM(kk).type==3 % weak seed FZ
%                 for bb=1:length(GEOM(kk-1).fz)-1
%                             ep(xm<=GEOM(kk-1).fz(bb)*(max(xm))+0.5*GEOM(2).FZseed(1) & xm>GEOM(kk-1).fz(bb)*(max(xm))-0.5*GEOM(2).FZseed(1) ...
%                             & ym>=GEOM(kk-1).top-GEOM(kk-1).step(bb) & ym<GEOM(kk-1).top+GEOM(2).FZseed(2) ...
%                             & im>1)=MAT(im).ecrit;
%                 end
            end
            
        end
        
    end
end
