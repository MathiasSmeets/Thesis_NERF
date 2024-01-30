function gui_cpd()
% gui_cpd  Gui tool for the tensor cpd function
%
% gui_cpd() Starts the tensor cpd GUI tool
%
% Author(s):    Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%               Matthieu Vendeville (matthieu.vendeville@esat.kuleuven.be)
%
% Version History:
% - 2018/05/03   MV      Initial version
import tensorlab.gui.cpd.*

if verLessThan('matlab', '9.0')
    warndlg('You are using an old matlab version, the GUI might not be fully functional')
end

m = tensorlab.gui.cpd.CPDModel();
p = tensorlab.gui.cpd.CPDPresenter(m);
