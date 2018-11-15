function G = create_img_graph(prob)

    n = size(prob.imgs, 1);
    A = zeros(n, n);
    for i  = 1 : n
        pts_i = prob.img_pts(prob.img_pts(:, 4) == prob.imgs(i, 1), :);

        for j  = i+1 : n
            pts_j = prob.img_pts(prob.img_pts(:, 4) == prob.imgs(j, 1), :);
            [~, idx_i, idx_j] = intersect(pts_i(:,5), pts_j(:,5));
            %if length(idx_i) > 0
                %fprintf('%i %i n = %i\n', i, j, length(idx_i));
                A(i, j) = length(idx_i);
                A(j, i) = length(idx_i);
            %end

        end
    end
    G = graph(A);