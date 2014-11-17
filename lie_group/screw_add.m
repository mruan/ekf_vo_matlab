function [result] = screw_add(addent1, addent2)
% Also works as screw minus

result = screw_log(screw_exp(addent1)*screw_exp(addent2));

end

