
function filter_string = filter_lookup(filter_type, filter_tag)

switch filter_type
    case 'cellfilter'
        switch filter_tag
            case 'non_mu_clusters'
                filter_string = '($numspikes > 100) && (all(cellfun(''isempty'',(arrayfun(@(x) strfind(x,''mua''), $tags, ''un'', 0)))))';
        end
end
end