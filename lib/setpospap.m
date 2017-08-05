%%% 4-10-11: 
% function [] = setpospap(pos)
% Sets position to vector pos and paperpositionmode to auto

function [] = setpospap(pos)

if ~exist('pos','var')
    pos = get(gcf,'position');
end
set(gcf,'paperPositionMode','auto','position',pos);

end