function setminsize(fh)
%SETMINSIZE Set minimum size of figure.
%   SETMINSIZE(FH) sets the minimum size of the figure defined the handle FH to its current
%   size. 

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven.be)
%
% Version History:
% - 2018/08/13   NV      Initial version
    
    try 
        state1 = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
        state1 = warning('off','MATLAB:ui:javaframe:PropertyToBeRemoved');        
        jFrame = get(fh, 'JavaFrame');  
        warning(state1);        
        warning(state2);
        jClient = jFrame.fHG2Client;
        jWindow = jClient.getWindow;
        if isempty(jWindow)
            drawnow;
            pause(0.02);
            jWindow = jClient.getWindow;
        end
        jDim = jWindow.getSize;
        minsize = [jDim.getWidth, jDim.getHeight]; 
        jWindow.setMinimumSize(java.awt.Dimension(minsize(1), minsize(2)));
    end 
end
