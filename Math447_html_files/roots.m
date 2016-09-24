function root=roots(f);

cl=0;
for i=-pi:.01:pi;
    if f(i)+f(i+1)==1;
        cl=cl+1;
    end
end
root=cl