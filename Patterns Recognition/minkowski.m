clear all
clc

% Calculate the Minkowski distance for 2 spots
n_List = [1, 2, 3, 10];
r_list = 1 : 10;
for k = 1 : length(n_List)

    n = n_List(k);
    a = randn(n, 1);
    b = randn(n, 1);
    distances_r = zeros(1, length(r_list));
    for i = 1 : length(r_list)
        r = r_list(i);
        SUM = 0;
        for dim = 1 : n
            SUM = SUM + abs(a(dim) - b(dim)) ^ r;
        end
        distance_r = (SUM)^(1/r);
        distances_r(i) = distance_r;
    end
    plot(distances_r);
    hold on
    
end

xlabel("r - Minkowski");
ylabel("Distance(r)");
legend("1D", "2D", "3D", "10D");
