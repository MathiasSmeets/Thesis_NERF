function helpel(a)




switch a
    
    case 1
    web https://www.tensorlab.net/doc 
           
    case 2
    d = dialog('Name','About CPD GUI Tensorlab','position', [100 200 600 550]);
    myicon = imread(fullfile(fileparts(mfilename('fullpath')), 'tensorlablogo.png'));
    
    txt = uicontrol('Parent',d,...
               'Style','text',...
               'Position',[25 60 550 315],...
               'fontweight', 'normal',...
               'String',{'CPD Tensorlab GUI BETA - 2018', '', '//////////','','','Uses Tensorslab V4.0:','','Tensorlab provides various tools for tensor computations, (coupled) tensor decompositions and complex optimization. In Tensorlab, datasets are stored as (possibly incomplete, sparse or structured) vectors, matrices and higher-order tensors, possibly obtained after tensorizing lower-order data. By mixing different types of tensor decompositions and factor transformations, a vast amount of factorizations can be computed. Users can choose from a library of preimplemented transformations and structures, including nonnegativity, orthogonality, Toeplitz and Vandermonde matrices to name a few, or even define their own factor transformations.','','','© Copyright 2018','','','Authors:', 'Nico Vervliet, Lieven De Lathauwer....etc'});
    
    btn = uicontrol('Parent',d,...
               'Position',[262 40 70 25],...
               'String','Close',...
               'Callback','delete(gcf)');
           
           hold on
         %  image(0,200,myicon);
           imshow(myicon);
           ylim([0 400]);
           axis off
          
          % Create and display the text label
url = 'https://www.tensorlab.net';
labelStr = ['<html>More info: <a href="">' url '</a></html>'];
jLabel = javaObjectEDT('javax.swing.JLabel', labelStr);
[hjLabel,hContainer] = javacomponent(jLabel, [10,10,250,20], gcf);
 
% Modify the mouse cursor when hovering on the label
hjLabel.setCursor(java.awt.Cursor.getPredefinedCursor(java.awt.Cursor.HAND_CURSOR));
 
% Set the label's tooltip
hjLabel.setToolTipText(['Visit the ' url ' website']);
 
% Set the mouse-click callback
set(hjLabel, 'MouseClickedCallback', @(h,e)web(['http://' url], '-browser'))
       
        
    case 3
    d = dialog('Position',[300 300 240 200],'Name','Referencing');

    txt = uicontrol('Parent',d,...
               'Style','text',...
               'Position',[10 20 225 180],...
               'String',{'', 'The toolbox can be cited as:', '', 'Vervliet N., Debals O., Sorber L., Van Barel M. and De Lathauwer L. Tensorlab 3.0, Available online, Mar. 2016. URL: https://www.tensorlab.net/','', 'Click text to copy to clipboard'});
            
           set(txt, 'Enable', 'Inactive');
           set(txt,'ButtonDownFcn',@copytxt);

% Create and display the text label
url = 'https://www.tensorlab.net';
labelStr = ['<html>More info: <a href="">' url '</a></html>'];
jLabel = javaObjectEDT('javax.swing.JLabel', labelStr);
[hjLabel,hContainer] = javacomponent(jLabel, [10,10,250,20], gcf);
 
% Modify the mouse cursor when hovering on the label
hjLabel.setCursor(java.awt.Cursor.getPredefinedCursor(java.awt.Cursor.HAND_CURSOR));
 
% Set the label's tooltip
hjLabel.setToolTipText(['Visit the ' url ' website']);
 
% Set the mouse-click callback
set(hjLabel, 'MouseClickedCallback', @(h,e)web(['http://' url], '-browser'))

    btn = uicontrol('Parent',d,...
               'Position',[90 40 70 25],...
               'String','Close',...
               'Callback','delete(gcf)');
           set(btn,'ButtonDownFcn',@copytxt);
        
    case 4
     web https://www.tensorlab.net -browser      
           
           
end

           function copytxt(h,~)
      clipboard('copy','Vervliet N., Debals O., Sorber L., Van Barel M. and De Lathauwer L. Tensorlab 3.0, Available online, Mar. 2016. URL: https://www.tensorlab.net/')
           end
       
end