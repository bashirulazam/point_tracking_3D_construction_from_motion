% Relative 2D coordinates
function [c, r] = relative2D(c_original, r_original)
    c_centroid = mean(c_original);
    r_centroid = mean(r_original);
    c = c_original - c_centroid;
    r = r_original - r_centroid;
end 
