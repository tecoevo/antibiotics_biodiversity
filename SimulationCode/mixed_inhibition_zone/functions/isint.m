function truthvalue = isint(num)
%Returns 1 if num is an integer, 0 otherwise. 
    if mod(num,1) == 0
        truthvalue=1;
    else
        truthvalue=0;
    end    
end