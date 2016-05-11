function [c] = minmod(a,b) % calculate minmod for tvd scheme
if (a*b>0) 
    if (abs(a)<abs(b)) 
     c = a;
    elseif (abs(a)>abs(b))
        c=b;
    end;
else c=0;

end;
if (a==b)
    c=a;
end;
    