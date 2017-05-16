function out = padto(in, sizeout)
    [m, n] = size(in);
    if m > sizeout(1) || n > sizeout(2)
        fprintf('error')
    else
        out = [in zeros(m,sizeout(2)-n); zeros(sizeout(1)-m, n) zeros(sizeout(1) - m, sizeout(2) - n)];
    end
end
