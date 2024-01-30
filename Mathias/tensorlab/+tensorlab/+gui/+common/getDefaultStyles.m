function styles = getDefaultStyles()
%GETDEFAULTSTYLES Return default styles for objects.
%   STYLES = GETDEFAULTSTYLES returns a struct of default styles for Tensorlab GUIs.

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%
% Version History:
% - 2018/08/09   NV      Initial version
    
    styles = struct;
    
    label = struct;
    label.Style                 = 'text';
    label.HorizontalAlignment   = 'left';
    label.FontSize              = 10;
    label.FontName              = 'Arial';
    label.FontWeight            = 'bold';
    label.Units                 = 'pixels';
    styles.label = label;
    
    comment = label;
    comment.FontSize            = 9;
    comment.FontWeight          = 'normal';
    styles.comment = comment;
    
    value = struct;
    value.Style                 = 'text';
    value.HorizontalAlignment   = 'left';
    value.FontSize              = 10;
    styles.value = value;
    
    edit = struct;
    edit.Style                  = 'edit';
    edit.HorizontalAlignment    = 'center';
    edit.FontSize               = 10;            
    styles.edit = edit;

    panel = struct;
    panel.BorderType            = 'etchedin';
    panel.ForegroundColor       = [0.22 0.541 0.749];
    panel.FontSize              = 10;
    panel.FontName              = 'Arial';
    panel.FontWeight            = 'bold';
    panel.Units                 = 'Pixels';
    styles.panel = panel;
    
    buttongroup = struct;
    buttongroup.BorderType      = 'none';
    buttongroup.ForegroundColor = 'blue';
    buttongroup.FontSize        = 10;
    buttongroup.FontName        = 'Arial';
    buttongroup.Units           = 'Pixels';
    styles.buttongroup = buttongroup;
    
    button = struct;
    button.Style                = 'pushbutton';
    styles.button = button;

    accentbutton = button;
    accentbutton.BackgroundColor = [0.22 0.541 0.749];
    accentbutton.ForegroundColor = [1 1 1];
    accentbutton.FontWeight      = 'bold';
    styles.accentbutton = accentbutton;
    
    radio = struct;
    radio.Style                  = 'radiobutton';
    styles.radio = radio;
    
    checkbox = struct;
    checkbox.Style               = 'checkbox';
    styles.checkbox = checkbox;
    
    slider = struct;
    slider.Style                 = 'slider';
    styles.slider = slider;

    dropdown = struct;
    dropdown.Style               = 'popupmenu';
    styles.dropdown = dropdown;
    
    place = struct;
    place.BorderType             = 'none';
    place.Units                  = 'Pixels';
    styles.place = place;
    
    list = struct;
    list.Style                  = 'listbox';
    styles.list = list;
end
