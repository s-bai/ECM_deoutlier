function m_swap = row_column_swap(m, index_target)
    %%  This function swap the targeted row and column to the 1st and 2nd places
    % % The rest rows and columns are sorted descendingly afterwards.

    size_m = size(m, 1);

    index = (1:size_m).';

    for ii = 1:2

        switch ii
            case 1
                m_aug = [index, m];
            case 2
                m_aug = m_aug.';

                m_aug_temp = m_aug(2:end, :);

                m_aug = [index, m_aug_temp];

                clear m_aug_temp;
        end

        m_aug(index_target, 1) = [1; 2];

        m_aug(setdiff(index, index_target), 1) = setdiff(index, [1; 2]);

        m_aug = sortrows(m_aug, 1);
    end

    m_swap = m_aug(:, 2:end).';

end
