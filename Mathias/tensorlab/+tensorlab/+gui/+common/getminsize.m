function minsize = getminsize(fh)
%SETMINSIZE Set minimum size of figure.
%   SETMINSIZE(FH, SZ) sets the minimum size of the figure defined the handle FH to the given
%   minimum size MINSIZE.
%   

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven.be)
%
% Version History:
% - 2018/08/13   NV      Initial version
    
    bakWarn1 = warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
    bakWarn2 = warning('off','MATLAB:ui:javaframe:PropertyToBeRemoved') ;   
    jFrame  = get(handle(fh), 'JavaFrame');  % Fails in Matlab < 7.0
    warning(bakWarn1);
    warning(bakWarn2);    
    jClient = jFrame.fHG2Client;
    jWindow = jClient.getWindow;
    if isempty(jWindow)
        drawnow;
        pause(0.02);
        jWindow = jClient.getWindow;
    end
    jDim = jWindow.getSize;
    minsize = [jDim.getWidth, jDim.getHeight];
end
