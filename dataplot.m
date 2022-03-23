margin = importdata("newfile.txt");
force = importdata("newfile1.txt");
trajectory = importdata("newfile2.txt");
[a,b] = size(margin);
for i = 1:a
    for j = 1:b
        if margin(i,j) >= 100000
            margin(i,j) = 0;
        end
    end
end

plot(margin)
legend