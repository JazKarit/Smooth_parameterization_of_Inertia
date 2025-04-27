function p_i = func_pi(theta)

    alpha = theta(1);
    d1 = theta(2);
    d2 = theta(3);
    d3 = theta(4);
    s12 = theta(5);
    s23 = theta(6);
    s13 = theta(7);
    t1 = theta(8);
    t2 = theta(9);
    t3 = theta(10);
    
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