function rgb = GenerateRandomColorHue(color, sensitivity)

% sensitivity = 0.75;

if(strcmpi(color,'red'))
    eps = sensitivity*rand(1,3) - 0.5*sensitivity;
    rgb = clip([1 0 0] + eps);
end

if(strcmpi(color,'green'))
    eps = sensitivity*rand(1,3)- 0.5*sensitivity;
    rgb = clip([0 1 0] + eps);
end

if(strcmpi(color,'blue'))
    eps = sensitivity*rand(1,3)- 0.5*sensitivity;
    rgb = clip([0 0 1] + eps);
end

if(strcmpi(color,'black'))
    eps = sensitivity*rand(1,3)- 0.5*sensitivity;
    rgb = clip([0 0 0] + eps);
end

if(strcmpi(color,'yellow'))
    eps = sensitivity*rand(1,3)- 0.5*sensitivity;
    rgb = clip([1 1 0] + eps);
end

if(strcmpi(color,'magenta'))
    eps = sensitivity*rand(1,3)- 0.5*sensitivity;
    rgb = clip([1 0 1] + eps);
end



end

% function to clip values to within 0 and 1
function v = clip(v)
    for i =1:length(v)
        if(v(i) > 1)
            v(i) = 1;
        elseif(v(i) < 0)
            v(i) = 0;
        end
    end
end