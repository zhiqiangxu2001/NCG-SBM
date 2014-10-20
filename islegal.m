function legal = islegal(v)

    legal = ~any(imag(v(:))) & ~isnan(v(:)) & ~isinf(v(:));
    
end