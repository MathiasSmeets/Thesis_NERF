function jobj = findjobj(h, tag)
%FINDJOBJ short discription
%   long description
%   

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%
% Version History:
% - 2018/08/14   NV      Initial version
    
    fh = ancestor(h, 'Figure');
    if nargin == 1
        tag = get(h, 'Tag');
    end 
    state1  = warning('off', 'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
    state2  = warning('off', 'MATLAB:ui:javaframe:PropertyToBeRemoved');    
    fh = get(fh, 'JavaFrame');
    jobj = localfind(fh.getFigurePanelContainer.getComponent(0), tag);
    warning(state1);
    warning(state2);    
    
    function jobj = localfind(comp, tag)
        jobj = [];
        try 
            name = get(comp, 'Name');
            if strcmpi(char(name), tag)
                jobj = comp;
                return;
            end 
        end 
        for k = 0:comp.getComponentCount()-1
            jobj = localfind(comp.getComponent(k), tag);
            if ~isempty(jobj)
                return;
            end 
        end 
    end 
end
