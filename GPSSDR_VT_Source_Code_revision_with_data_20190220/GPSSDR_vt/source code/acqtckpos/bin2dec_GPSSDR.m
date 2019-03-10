%Written by Jordan Krcmaric

function d = bin2dec_GPSSDR(b)

    for i=1:size(b,1) 
        d(i,:) = polyval(fliplr(b(i,:)),2);
    end
    
end
