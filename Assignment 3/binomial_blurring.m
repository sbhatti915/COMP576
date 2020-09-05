function [out] = binomial_blurring(img)


    [h, w] = size(img); 
    out = im2double(img);
    for i = 1:h
        for j = 1:w
            if j == w
                out(i,j) = out(i,j)*2 / 2;
            else
                out(i,j) = (out(i,j) + out(i, j+1)) / 2;
            end
        end
        for j = w:-1:1
            if j == 1
                out(i,j) = out(i,j)*2 / 2;
            else
                out(i,j) = (out(i,j) + out(i, j-1)) / 2;
            end
        end
    end

        
    for j = 1:w
        for i = 1:h
            if i == h
                out(i,j) = out(i,j)*2 / 2;
            else
                out(i,j) = (out(i,j) + out(i+1, j)) / 2;
            end
        end

        for i = h:-1:1
            if i == 1
                out(i,j) = out(i,j)*2 / 2;
            else
                out(i,j) = (out(i,j) + out(i-1, j)) / 2;
            end
        end
    end
        
end 
