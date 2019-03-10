function value = comp2dec(bi,LSB)

if bi(end) == 0
    value = bin2dec_GPSSDR(bi(1:end-1)) * 2^LSB;
    
else
    temp = bi(1:end-1);
    for i=1:length(temp)
        if temp(i) == 1;
            temp2(i) = 0;
        else
            temp2(i) = 1;
        end
    end
    value = -1 * (1+bin2dec_GPSSDR(temp2)) * 2^LSB;
    
end