function aps = get_aps(ts)

times=[2:10];
aps=[0.1380    0.1553    0.1933    0.3018    0.7279    1.1917    1.5069    1.6594    1.7177];
aps=interp1(times,aps,ts);
% aps=aps.^2;
end

