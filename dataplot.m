margin = importdata("newfile.txt");
force = importdata("newfile1.txt");
trajectory = importdata("newfile2.txt");
[a,b] = size(margin);
index = 1;
index2 = 3;
check_violation = false;
for i = 1:a
    for j = 1:b-1
        if margin(i,j) > 0
            violation_ticks(index,index2) = j;
            index2 = index2 + 1;
            violation_ticks(index,1) = i;
            violation_ticks(index,2) = margin(i,11);
            check_violation = true;
        end
        if margin(i,j) >= 100000
            margin(i,j) = 0;
        end
    end
    if check_violation
        index = index +1;
        check_violation = false;
    end
    index2 = 3;
end

plot(margin)
legend
