function [p] = min_max_wf(Pt,vector)
    %vector is reciprocal ,so vector is in decreasing order
    len = length(vector);
    allo_set = 1:len;
    %initialize 
    %assuming that level>all bottom ,so s.t. level > max(vector.^(1/2))
    level = (Pt+sum(vector))/sum(vector.^(1/2));
    while (level<vector(allo_set(1))^(1/2))
        %delet the first element
        allo_set(1) = []; 
        level = (Pt+sum(vector(allo_set)))/sum(vector(allo_set).^(1/2));
    end
    p = max(level*vector.^(1/2)-vector,0);
end
