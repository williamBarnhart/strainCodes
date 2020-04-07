function  b=ZeroTwoPi(a)

b=a;
twopi=2*pi;
if b<0
    b=b+twopi;
elseif b>=twopi
    b=b-twopi;
end
