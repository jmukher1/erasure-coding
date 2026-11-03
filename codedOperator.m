classdef codedOperator
% codedOperator - matrix-free operator: Bop * U evaluates fun(U).
% Used for the coded B' = P'BP so that code written for an explicit B
% (B*X) runs unchanged after a fault.
    properties
        fun
        n
    end
    methods
        function o = codedOperator(fun, n)
            o.fun = fun;  o.n = n;
        end
        function Y = mtimes(o, X)
            if isa(o, 'codedOperator')
                Y = o.fun(X);
            else
                error('codedOperator:mtimes', 'only codedOperator * matrix is supported');
            end
        end
        function varargout = size(o, d)
            if nargin < 2
                if nargout <= 1, varargout = {[o.n, o.n]}; else, varargout = {o.n, o.n}; end
            else
                varargout = {o.n};
            end
        end
    end
end
