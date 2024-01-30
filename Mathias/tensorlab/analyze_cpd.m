function varargout = analyze_cpd(T, U, varargin)
%ANALYZE_CPD short discription
%   long description
%   

% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%
% Version History:
% - 2018/08/23   NV      Initial version
    
    p = inputParser();
    p.addOptional('ShowPlots', true);
    p.addOptional('Display', true);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;
    
    R = size(U{1},2);
    
    output = struct;
    output.abserr = frobcpdres(T, U);
    output.relerr = output.abserr/frob(T);
    
    % Compute norms 
    nrmc = cellfun(@(u) sqrt(sum(abs(u).^2,1)), U, 'UniformOutput', false);
    nrm  = prod(vertcat(nrmc{:}),1);
    [output.norms, ind] = sort(nrm, 2, 'descend');
    U = cellfun(@(u) u(:,ind), U, 'UniformOutput', false);    

    % Compute congruence
    Utmp = cellfun(@(u,n) bsxfun(@rdivide, u, n), U, nrmc, 'UniformOutput', false);
    W    = cellfun(@(u) real(u'*u), Utmp, 'UniformOutput', false);
    output.congruence = prod(cat(3,W{:}),3);

    % Compute angles
    output.angles = acos(output.congruence)/pi*180;
    
    output.corcondia = corcondia(T, U);
    
    output.condition = cpd_cond(U);
    
    UHU = cellfun(@(u) u'*u, U, 'UniformOutput', false);
    [~,D] = eig(prod(cat(3,UHU{:}),3));
    D = sqrt(abs(diag(D)));
    output.ratiosv = max(D)/min(D);

    if options.Display
        fprintf('CPD analysis summary\r\n');
        fprintf('--------------------\r\n\r\n');
        fprintf('Absolute error:   %12.3e\r\n', output.abserr);
        fprintf('Relative error:   %12.3e\r\n', output.relerr);
        fprintf('Corcondia:        %12.1g %%\r\n', output.corcondia);
        fprintf('Condition number: %12.1g\r\n', output.condition);
        fprintf('Ratio min/max sv: %12.2g\r\n', output.ratiosv);
        fprintf('Congruence (max): %12.2f\r\n', max(max(output.congruence-eye(R))));
        minangle = acos(abs(output.congruence - eye(R)))/pi*180;
        fprintf('Angle (min):      %12.2f degrees\r\n', min(minangle(:)));
        fprintf('Norms (max):      %12.3e\r\n', output.norms(1));
        fprintf('Norms (min):      %12.3e\r\n', output.norms(end));
        fprintf('\r\n');
    end 
    
    if options.ShowPlots
        figure('Name', 'Angles between rank-1 terms');
        plot_angles(U);
        title('Angles between rank-1 terms');
        figure('Name', 'Error per slice');
        plot_errslice(T, U);
        title('Error per slice');
        figure('Name', 'Computed CPD versus original data');
        visualize(U, 'original', T);
        figure('Name', 'Norms of rank-1 terms'); 
        stem(output.norms, 'filled')
        title('Norms of rank-1 terms');
        figure('Name', 'Rank-1 terms'); 
        plot_rank1terms(U);
    end 
    
    if nargout == 1
        varargout = {output};
    else 
        varargout = {};
    end 
end
