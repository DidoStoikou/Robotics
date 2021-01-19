function solution = Create_Trajectory (initial, final, t0, T, D)

% equation variables
syms a0 a1 a2 a3 a4 a5 
syms b0 b1 
syms c0 c1 c2 c3 c4 c5 

% time intervals
t1 = t0 + D;
t2 = T - D;

% equations
eqn1 = a5*t0^5 + a4*t0^4 + a3*t0^3 + a2*t0^2 + a1*t0 + a0 == initial;
eqn2 = 5*a5*t0^4 + 4*a4*t0^3 + 3*a3*t0^2 + 2*a2*t0 + a1 == 0;
eqn3 = 20*a5*t0^3 + 12*a4*t0^2 + 6*a3*t0 + 2*a2 == 0;

eqn4 = a5*t1^5 + a4*t1^4 + a3*t1^3 + a2*t1^2 + a1*t1 + a0 == b1*t1 + b0;
eqn5 = 5*a5*t1^4 + 4*a4*t1^3 + 3*a3*t1^2 + 2*a2*t1 + a1 == b1;
eqn6 = 20*a5*t1^3 + 12*a4*t1^2 + 6*a3*t1 + 2*a2 == 0;

eqn7 = c5*t2^5 + c4*t2^4 + c3*t2^3 + c2*t2^2 + c1*t2 + c0 == b1*t2 + b0;
eqn8 = 5*c5*t2^4 + 4*c4*t2^3 + 3*c3*t2^2 + 2*c2*t2 + c1 == b1;
eqn9 = 20*c5*t2^3 + 12*c4*t2^2 + 6*c3*t2 + 2*c2 == 0;

eqn10 = c5*T^5 + c4*T^4 + c3*T^3 + c2*T^2 + c1*T + c0 == final;
eqn11 = 5*c5*T^4 + 4*c4*T^3 + 3*c3*T^2 + 2*c2*T + c1 == 0;
eqn12 = 20*c5*T^3 + 12*c4*T^2 + 6*c3*T + 2*c2 == 0;

eqn13 = (b1*t1 + b0 - (a5*t0^5 + a4*t0^4 + a3*t0^3 + a2*t0^2 + a1*t0 + a0))/D == 0.5*b1;
eqn14 = (c5*T^5 + c4*T^4 + c3*T^3 + c2*T^2 + c1*T + c0 - b1*t2 - b0)/D == 0.5*b1;

% solving system
[A,B] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7, eqn8, eqn9, eqn10, eqn11, eqn12, eqn13, eqn14], [a5, a4, a3, a2, a1, a0, b1, b0, c5, c4, c3, c2, c1, c0]);

solution = linsolve(A,B);

end