function p_i = func_pi(theta)

    alpha = theta.alpha;
    d1 = theta.d1;
    d2 = theta.d2;
    d3 = theta.d3;
    s12 = theta.s12;
    s23 = theta.s23;
    s13 = theta.s13;
    t1 = theta.t1;
    t2 = theta.t2;
    t3 = theta.t3;
    
    e = exp(1);

    p_i = zeros(10, 1);

    p_i(1) = t1^2 + t2^2 + t3^2 + 1;
    p_i(2) = t1*(e^d1);
    p_i(3) = t1*s12 + t2*(e^d2);
    p_i(4) = t1*s13 + t2*s23 + t3*(e^d3);
    p_i(5) = s12^2 + s13^2 + s23^2 + e^(2*d2) + e^(2*d3);
    p_i(6) = s13^2 + s23^2 + e^(2*d1) + e^(2*d3);
    p_i(7) = s12^2 + e^(2*d1) + e^(2*d2);
    p_i(8) = -s12*e^(d1);
    p_i(9) = -s12*s13 - s23*e^(d2);
    p_i(10) = -s13*e^(d1);

    p_i = e^(2*alpha)*p_i;

end