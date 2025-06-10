function [array] = shiftArray(array,newValue,maxLength)
    if length(array) >= maxLength
        array = [array(2:end),newValue];
    else
        array = [array, newValue];
    end
end

