function [ w ] = Lambert_W_Fritsch( x )
% This function evaluates the upper branch of the Lambert-W function for
% real input x. Function is not optimized for vectorial processing.

% Values
% 
% x - Input Range: [0,inf)
% w - W(x)

% Function written by F. Esqueda 2/10/17 based on implementation presented 
% by Darko Veberic - "Lambert W Function for Applications in Physics"
% Available at https://arxiv.org/pdf/1209.0735.pdf

%% Compute Lambert-W

% Error threshold
thresh = 10e-12;

if(x)==0
    w=0;
    return
end
% Initial guess
if x<0.14546954290661823
    
    num = 1 + 5.931375839364438*x + 11.39220550532913*x^2 +  7.33888339911111*x^3 + 0.653449016991959*x^4;
    den = 1 + 6.931373689597704*x + 16.82349461388016*x^2 + 16.43072324143226*x^3 + 5.115235195211697*x^4;
    w = x * num / den;

elseif x<8.706658967856612
    
    num = 1 + 2.4450530707265568*x + 1.3436642259582265*x^2 + 0.14844005539759195*x^3 + 0.0008047501729130*x^4;
    den = 1 + 3.4447089864860025*x + 3.2924898573719523*x^2 + 0.9164600188031222*x^3  + 0.05306864044833221*x^4;
    w = x * num / den;
        
else

    a = log(x);
    b = log(a);    
    ia = 1/a;
    w = a - b + (b*ia) + 0.5*b*(b - 2)*(ia*ia) + (1/6)*(2*b*b - 9*b + 6)*(ia*ia*ia);

end

% Comput Lambert-W (limited to 20 iterations)
for m=1:20
    
    w1 = w + 1;
    z = log(x) - log(w) - w;
    q = 2*w1*(w1 + (2/3)*z);
    e = (z/w1)*((q - z)/(q - 2*z));
    w = w * (1 + e);
    
    if abs(e)<thresh
        break; 
    end

end

end

