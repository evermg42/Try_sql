function [v11,v12,v21,v22]= gradient_x(u)
    v11 = diff(u.M{1},1,1);v12 = diff(u.M{1},1,2);
    v21 = diff(u.M{2},1,1);v22 = diff(u.M{2},1,2);