function [R] = screw_exp(w)

nw2 = dot(w, w);
nw  = sqrt(nw2);

wx = so3_alg(w);

if nw < 1e-1
    R = eye(3) + wx + 0.5 * wx * wx;
else
    R = eye(3) + (sin(nw)/nw)*wx + (1-cos(nw))/nw2 * wx * wx;
end

end