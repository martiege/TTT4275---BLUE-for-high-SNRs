function est = BLUE(signal, H, C)

M1 = C \ H;
M2 = H' * M1;
M3 = M2 \ H';
M4 = M3 / C;

est = M4 * unwrap(angle(signal))';

end

