function p = penalty(ksmall,klarge,g1,g2)
% input
% ksmall: parameter that is smaller
% klarge: parameter that is larger
% g1: penalty weight one
% g2: penalty weight two

p = g1*heaviside(ksmall-klarge)*((ksmall-klarge)/klarge)^2 + g2*heaviside(ksmall-klarge);

end