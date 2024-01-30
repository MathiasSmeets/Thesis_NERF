function varargout = gui_mlsv_inspect(T)
%GUI_MLSV_INSPECT Inspect multilinear singular values. 
%   GUI_MLSV_INSPECT(T) lets a user compute a MLSVD of a given tensor T and allows the
%   multilinear singular values to be inspected visually. 

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%
% Version History:
% - 2018/08/08   NV      Initial version

    if verLessThan('matlab', '9.0')
        warndlg('You are using an old matlab version, the GUI might not be fully functional')
    end

    import tensorlab.gui.mlsvd.*;
    name = inputname(1);
    model = ExplorationModel(T, name);
    ExplorationPresenter(model);

end 
