function f = frob(T,squared,weights)
%FROB Frobenius norm
%   FROB(T) returns the Frobenius norm of the tensor T. In the case of
%   incomplete tensors, the weighted Frobenius norm error is used, i.e., the
%   norm of the known elements is computed. In the case T is an efficient
%   representation of a structured tensor, the norm is computed without
%   expanding the tensor.
%
%   FROB(T,'squared') returns the squared Frobenius norm of the tensor T.
%
%   See also norm.

% Authors: Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%          Otto Debals         (Otto.Debals@esat.kuleuven.be)
%          Laurent Sorber      (Laurent.Sorber@cs.kuleuven.be)
%          Marc Van Barel      (Marc.VanBarel@cs.kuleuven.be)
%          Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
% Version History:
% - 2018/04/05    NV      Added support weights
% - 2016/02/02    OD      Added support Hankel and Loewner
% - 2016/12/18    NV      Added structured tensor support 

    if nargin < 2, squared = false; end
    if nargin < 3, weights = {[]}; end
    widx = find(~cellfun(@isempty, weights));
    
    tomode = @(v,n) reshape(v, [ones(1,n-1), numel(v) 1]);

    switch getstructure(T)
      case 'full'
        for k = widx
            T = bsxfun(@times, T, tomode(weights{k},k));
        end 
        f = full(T(:)'*T(:));
      case {'sparse', 'incomplete'}
        for k = widx
            T.val = T.val .* weights{k}(T.sub{k});
        end 
        f = T.val(:)'*T.val(:);
      case 'cpd'
        if ~isempty(widx), 
            T(widx) = cellfun(@(w,z)bsxfun(@times,w,z), weights(widx), T(widx), 'UniformOutput', ...
                              false);
        end 
        W = cellfun(@(u) u'*u, T, 'UniformOutput', false);
        W = prod(cat(3, W{:}),3);
        f = abs(sum(W(:)));
      case 'lmlra'
        if ~isempty(widx), 
            T{1}(widx) = cellfun(@(w,z)bsxfun(@times,w,z), weights(widx), T{1}(widx), 'UniformOutput', ...
                                 false);
        end 
        W = cellfun(@(u) u'*u, T{1}, 'UniformOutput', false);
        S = tmprod(conj(T{2}), W, 1:length(W), 'T').*T{2};
        f = abs(sum(S(:)));
      case 'btd'
        % apply weights
        if ~isempty(widx)
            for r = 1:length(T)
                T{r}(widx) = cellfun(@(w,z)bsxfun(@times,w,z), weights(widx), T{r}(widx), ...
                                     'UniformOutput', false);
            end 
        end  
        % compute norm
        f = 0;
        for r = 1:length(T)
            W = cellfun(@(u) u'*u, T{r}(1:end-1), 'UniformOutput', false);
            S = tmprod(conj(T{r}{end}), W, 1:length(W), 'T').*T{r}{end};
            f = f + sum(S(:));
            for s = r+1:length(T)
                W = cellfun(@(u,v) u'*v, T{r}(1:end-1), T{s}(1:end-1), ...
                            'UniformOutput', false);
                S = tmprod(conj(T{r}{end}), W, 1:length(W), 'T').*T{s}{end};
                f = f + 2*real(sum(S(:)));
            end
        end
        f = abs(f);
      case 'tt'
        % apply weights
        if ~isempty(widx)
            if widx(1) == 1, T{1} = bsxfun(@times, weights(s), T{1}); end
            for k = widx(2:end)
                T{k} = bsxfun(@times, tomode(weights{k}, 2), T{k});
            end 
        end 
        % compute norm
        size_tens = cellfun('size', T, 2);
        ttr = [1 cellfun('size', T(2:end), 1) 1];
        tmp = T{1}'*T{1};
        f = tmp(:).';
        for n = 2:length(T)-1
            tmp = reshape(permute(T{n},[2 1 3]),size_tens(n),ttr(n)*ttr(n+1));
            tmp = tmp'*tmp;
            tmp = reshape(tmp, ttr([n n+1 n n+1]));
            tmp = permute(tmp, [1 3 2 4]);
            tmp = reshape(tmp, ttr(n)^2, ttr(n+1)^2);
            f = f * tmp;
        end
        tmp = conj(T{end}*T{end}');
        f = abs(f*tmp(:));
      case 'hankel'
        N = size(T.val,1);
        size_hankel = T.subsize.hankel;
        if ~isempty(widx)
            if T.ispermuted
                others = T.order+(1:numel(T.subsize.other));
            else 
                others = [1:T.dim-1 (T.dim+T.order):numel(T.size)];
            end 
            hdims = 1:numel(T.size);
            hdims(others) = [];
            for k = hdims
                if isempty(weights{k})
                    weights{k} = ones(size_hankel(k),1); 
                end
            end 
            w = fft(weights{hdims(1)}.^2,N,1);
            for k = hdims(2:end)
                w = w.*fft(weights{k}.^2,N,1);
            end
            w = ifft(w.');
            % weights in other modes
            weights = weights(others);
            widx = find(~cellfun(@isempty, weights));
            T.val = reshape(T.val, [size(T.val,1), T.subsize.other]);
            for k = widx
                T.val = bsxfun(@times, tomode(weights{k},k+1), T.val);
            end 
            T.val = reshape(T.val, size(T.val,1), []);
        else 
            if T.order==2,
                m = min(T.ind,N-T.ind+1);
                w = [1:m-1 m*ones(1,N-2*m+2) m-1:-1:1];
            else
                w = fft(ones(size_hankel(1),1),N);
                w = w(:);
                for i = 2:numel(size_hankel)
                    w = w.*fft(ones(size_hankel(i),1),N);
                end
                w = round(ifft(w).');
            end
        end 
        f = w*sum(T.val.*conj(T.val),2);
      case 'loewner'
        if T.order > 2
            % full case, as the structured case is not yet supported
            L = ful(T);
            f = frob(L,'squared',weights);
        else
            
            if ~isempty(widx)
                if T.ispermuted
                    others = T.order+(1:numel(T.subsize.other));
                else 
                    others = [1:T.dim-1 (T.dim+T.order):numel(T.size)];
                end 
                T.val = reshape(T.val, [size(T.val,1), T.subsize.other]);
                for k = 1:length(others)
                    if isempty(weights{others(k)}), continue; end
                    T.val = bsxfun(@times, tomode(weights{others(k)},k+1), T.val);
                end 
                T.val = reshape(T.val, size(T.val,1), []);
                weights = weights(T.dim:T.dim+1);
            end 
            
            if isempty(weights) || isempty(weights{1}), weights{1} = 1; end;
            if numel(weights) == 1 || isempty(weights{2}), weights{2} = 1; end;
           
            f = T.val(T.ind{1},:); 
            g = T.val(T.ind{2},:);
            
            if T.isequidistant
                f = bsxfun(@times, weights{1}, f); 
                g = bsxfun(@times, weights{2}, g); 
                                
                v = abs(T.structure.v).^2;
                tmp = ifft(fft(v,[],1).*fft(sum(f.*conj(f),2),numel(v),1));
                term1 = sum((weights{2}.^2) .' * tmp(end-numel(v)+size(f,1):end,:));
                
                tmp = ifft(fft(v(end:-1:1),[],1).*fft(sum(g.*conj(g),2),numel(v),1));
                term2 = sum((weights{1}.^2).'*tmp(end-numel(v)+size(g,1):end,:));

                f = bsxfun(@times, weights{1}, f); 
                g = bsxfun(@times, weights{2}, g); 
                tmp = ifft(bsxfun(@times,fft(v,[],1),fft(f,numel(v),1)),[],1);
                term3 = sum(sum(real(tmp(end-numel(v)+size(f,1):end,:).*conj(g))));
            else 
                XY = 1./bsxfun(@minus,T.t(T.ind{1}),T.t(T.ind{2}).');
                XY = bsxfun(@times, weights{1}, XY);
                XY = bsxfun(@times, weights{2}.', XY); 
                
                XYt2 = (XY.*conj(XY)).';
                term1 = sum(XYt2*sum(f.*conj(f),2));
                term2 = sum(sum(g.*conj(g),2).'*XYt2);
                term3 = sum(sum(real((XYt2*f).*conj(g))));
            end
            
            f = term1+term2-2*term3;
        end
      otherwise
        error('frob:notImplemented', 'Not yet implemented');
    end
    if ~ischar(squared), f = sqrt(f); end
end