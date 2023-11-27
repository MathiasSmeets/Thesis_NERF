function B = lfp_bandsum_comparison(varargin)

% sum power activity at different bands to compare between conditions
%
% Mattia D'Andola, Sep 2020

def_labs = {};
def_bands = 1:5;
def_norm = 1;
def_plt = 1;
%inputs control
parser = inputParser;
addRequired(parser, 'psd', @ismatrix)
addRequired(parser, 'f', @isvector)
addParameter(parser, 'bands', def_bands);
addParameter(parser, 'label_names', def_labs);
addParameter(parser, 'norm', def_norm);
addParameter(parser, 'plt',def_plt);

parse(parser,varargin{:})
psdmean = parser.Results.psd;
f = parser.Results.f;
bands = parser.Results.bands;
label_names = parser.Results.label_names;
norm = parser.Results.norm;
plt = parser.Results.plt;

if isempty(label_names)
    for i = 1:length(psdmean)
        label_names{i} = sprintf('Condition %d',i);
    end
end

%delta [0 4) Hz
B(:,1) = sum(psdmean(:,f>=0.5 & f<4),2);
% theta [4 8) Hz
B(:,2) = sum(psdmean(:,f>=4 & f<8),2);
%alpha [8 13) Hz
B(:,3) = sum(psdmean(:,f>=8 & f<13),2);
%beta [13 30) Hz
B(:,4) = sum(psdmean(:,f>=13 & f<30),2);
%sum activity in low gamma band [30 60)
B(:,5) = sum(psdmean(:,f>=30 & f<60),2);
%sum activity in low gamma band [60 100) Hz
B(:,6) = sum(psdmean(:,f>=60 & f<100),2);

B = abs(B');

if plt
    X = categorical({'0.5-4 Hz','4-8 Hz','8-13 Hz','13-30 Hz','30-60 Hz','60-100 Hz'});
    X = reordercats(X,{'0.5-4 Hz','4-8 Hz','8-13 Hz','13-30 Hz','30-60 Hz','60-100 Hz'});

    B = B(bands,:);
    X = X(bands);
    figure, hold on
    bar(X,B);
   
    legend(label_names)
    ylabel('$\sum$ (PSD) [au]','Interpreter','Latex')
end