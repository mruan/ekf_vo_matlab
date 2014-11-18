function [ w ] = screw_log( R )
% Work based on Edward Rosten's C++ implemention of SO3 class
% http://www.edwardrosten.com/cvd/toon/html-user/so3_8h_source.html

cos_t = 0.5*(trace(R)-1);

if cos_t > 1.0
    cos_t = 1.0;
elseif cos_t < -1.0
    cos_t =-1.0;
end

w = 0.5*[R(3,2) - R(2,3); R(1,3) - R(3,1); R(2,1) - R(1,2)];

sin_t_abs = norm(w);

if cos_t > 1/sqrt(2) 
    % a ~ [0 ~ Pi/4] => use asin
    if sin_t_abs >0
        w = w * asin(sin_t_abs)/sin_t_abs;
        % else w is returned directly -> e.g. w = [0 0 0] for R==I
    end
elseif cos_t > -1/sqrt(2) 
    % a ~ [Pi/4 ~ 3Pi/4] => use acos, w/o antisymmetric part
    t = acos(cos_t);
    w = w * t / sin_t_abs;
else % use symmetric part
   % antisymmetric part vanishes, but still large rotation, need
   % information from symmetric part
    t = pi - asin(sin_t_abs);
    d = [R(1,1) - cos_t; R(2,2) - cos_t; R(3,3)-cos_t]; % diag(R)-cos_t
    D = d.^2;
    if D(1) > D(2) && D(1) > D(3) % 1st is largest, fill 1st column
        r2 = [d(1); ...
              0.5*(R(1,2)+R(2,1));...
              0.5*(R(1,3)+R(3,1))];
    elseif D(1) > D(2)            % 2nd is largest, fill 2nd column
        r2 = [0.5*(R(2,1)+R(1,2));...
              d(2);...
              0.5*(R(2,3)+R(3,2))];
    else
        r2 = [0.5*(R(3,1)+R(1,3));...
              0.5*(R(3,2)+R(2,3));...
              d(3)];
    end
    
    if dot(r2, w) < 0 % Flip if point in the wrong direction
       r2 = -r2; 
    end
    r2 = r2/norm(r2);
    w = t*r2;
end

end

