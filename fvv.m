function ovec = fvv(fun,ivec)
len = length(ivec); siz = size(ivec);
ovec = zeros(siz);
for i = 1:len
    ovec(i) = fun(ivec(i));
end

end
