% Rocky_closed_loop_poles.m
%
% 1) Symbolically calculates closed loop transfer function of PI disturbannce
% rejection control system for Rocky. 
% Currently no motor model (M =1).Placeholder for motor model (1st order TF)
%
% 2) Specify location of (target)poles based on desired reponse. The number of
% poles = denominator polynomial of closed loop TF
%
% 3) Extract the closed loop denomiator poly and set = polynomial of target
% poles
%
% 4) Solve for Ki and Kp to match coefficients of polynomials. In general,
% this will be underdefined and will not be able to place poles in exact
% locations. 
%
% 5) Plot impulse response to see closed-loop behavior. 
%
% based on code by SG. last modified 3/12/21 CL

clear all; 
close all;

syms s a b l g Kp Ki Jp JiCp Ci Cp   % define symbolic variables

Hvtheta = -s/l/(s^2-g/l);       % TF from velocity to angle of pendulum

% M = a*b/(s+a)                 % TF of motor
K = Kp + Ki/s
M = (a*b)/(s+a);
Nested_loops = M/(1+((Jp+JiCp/s+(Ci/s^2)))*M);
%M = 1;                          % TF without motor  
%  
%closed loop transfer function from disturbance d(t)totheta(t)
Hcloop = 1/(1-Hvtheta*Nested_loops*K);  

pretty(simplify(Hcloop))       % to display the total transfer function

% Substitute parameters and solve
% system parameters
g = 9.81;
l = 0.4053;  %effective length 
a = 13.89178;           %nomical motor parameters
b = 0.003304;        %nomical motor parameters

Hcloop_sub = subs(Hcloop); % sub parameter values into Hcloop

% specify locations of the target poles,
% choose # based on order of Htot denominator
% e.g., want some oscillations, want fast decay, etc. 

% Old, tested poles
% p1 = -1 + 1*1i;
% p2 = -1 - 1*1i;
% p3 = -4 + 1*1i;
% p4 = -4 - 1*1i;
% p5 = -5;

% Cara's calculated poles V1
% p1 = -4.42557 + 2.1434*1i;
% p2 = -4.42557 - 2.1434*1i;
% p3 = -4.1797 + 2.5903*1i;
% p4 = -4.1797 - 2.5903*1i;
% p5 = -4.9173;

% Example code from past student
  p1 = -4.521;
  p2 = -12.413;
  p3 = -3.4212;
  p4 = -12.7681;
  p5 = -9.3469;

% Cara's calculated poles V2
%  p1 = -4.9173;
%  p2 = -2.2821;
%  p3 = -6.5689;
%  p4 = -1.5893;
%  p5 = -6.77;

% Cara's calculated poles V3
%  p1 = -4.9173;
%  p2 = 4.9173*(-cosd(25) + sind(25));
%  p3 = 4.9173*(-cosd(25) - sind(25));
%  p4 = 4.9173*(-cosd(30) + sind(30));
%  p5 = 4.9173*(-cosd(30) - sind(30));


% target characteristic polynomial
% if motor model (TF) is added, order of polynomial will increases
tgt_char_poly = (s-p1)*(s-p2)*(s-p3)*(s-p4)*(s-p5);

% get the denominator from Hcloop_sub
[n d] = numden(Hcloop_sub);

% find the coefficients of the denominator polynomial TF
coeffs_denom = coeffs(d, s);

% divide though the coefficient of the highest power term
coeffs_denom = coeffs(d, s)/(coeffs_denom(end));

% find coefficients of the target charecteristic polynomial
coeffs_tgt = coeffs(tgt_char_poly, s);

% solve the system of equations setting the coefficients of the
% polynomial in the target to the actual polynomials
solutions = solve(coeffs_denom == coeffs_tgt,  Kp, Ki, Ci,JiCp, Jp);

% display the solutions as double precision numbers
Kp = double(solutions.Kp)
Ki = double(solutions.Ki)
Ci = double(solutions.Ci)
JiCp =  double(solutions.JiCp)
Jp = double(solutions.Jp)





% Location of the poles of the closed-loop TF.
% NOTE there are only 2 unknowns but 3 polynomial coefficients so 
% the problem is underdetermined and the closed loop poles don't exact
% match the target poles. 
% use trial-and-error to tune response
closed_loop_poles = vpa (roots(subs(coeffs_denom)), 4)

% Plot impulse response of closed-loop system
    TFstring = char(subs(Hcloop));
    % Define 's' as transfer function variable
    s = tf('s');
    % Evaluate the expression
    eval(['TFH = ',TFstring]);
    figure (1)
    impulse(TFH);   %plot the impulse reponse
    figure(2)
    step(TFH)
    
    







