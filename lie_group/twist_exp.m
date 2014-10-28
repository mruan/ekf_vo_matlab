function [T] = twist_exp(twist_coord)

assert(numel(twist_coord) == 6, 'twist coord must have 6 elements');

% Assume the first 3 elements are linear speed, then come the angular vel
nw = norm(twist_coord(4:6));
nw2 = nw*nw;

wx = so3_alg(twist_coord(4:6));

if nw < 1e-1
    R = eye(3) +     wx + 1/2 * wx * wx;
    V = eye(3) + 0.5*wx + 1/6 * wx * wx;
else
    R = eye(3) + (sin(nw)/nw)   *wx + (1 -cos(nw))/nw2 * wx * wx;
    V = eye(3) + (1-cos(nw))/nw2*wx + (nw-sin(nw))/nw2/nw * wx * wx;
end

T = [R V*twist_coord(1:3); 0 0 0 1];

end