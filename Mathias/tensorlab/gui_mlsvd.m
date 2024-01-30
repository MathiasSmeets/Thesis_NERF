function gui_mlsvd()
% gui_mlsvd  Gui tool for the tensor cpd function
%
% gui_mlsvd() Starts the tensor cpd GUI tool
%
% Author(s):    Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%               Matthieu Vendeville (matthieu.vendeville@esat.kuleuven.be)
%
% Version History:
% - 2018/05/03   MV      Initial version

if verLessThan('matlab', '9.0')
    warndlg('You are using an old matlab version, the GUI might not be fully functional')
end

tensorlab.gui.mlsvd.MLSVDPresenter(tensorlab.gui.mlsvd.MLSVDModel());
