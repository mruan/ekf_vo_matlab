delta = 0.0001;
x = obj.x;

z0 = proj_P3(x, K, R, t);

z1 = proj_P3(x+[delta 0 0 0]', K, R, t);

z2 = proj_P3(x+[0 delta 0 0]', K, R, t);

z3 = proj_P3(x+[0 0 delta 0]', K, R, t);

z4 = proj_P3(x+[0 0 0 delta]', K, R, t);

H = [z1-z0 z2-z0 z3-z0 z4-z0]/delta