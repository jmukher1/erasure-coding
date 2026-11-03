function faults = faultfunc_matlab2(i, err_iters, rowDeleteStart, rowDeleteEnd)
% rows erased at iteration i: rowDeleteStart:rowDeleteEnd 
    if (any(err_iters == i))
        faults = rowDeleteStart:rowDeleteEnd;
    else
        faults = [];
    end
end
