function[ vectors ] = plotProperHeight(vectors,...
                                       otherVectors,...
                                       ref)

offset = 0;
for kk = 1:size(vectors,1)
   vectors(kk,:) = vectors(kk,:) - min(vectors(kk,:)) + offset;
   if nargin > 1
       if  any(ref == kk)
        idx    = find(ref == kk);
        otherVectors(idx,:) = otherVectors(idx,:) -...
                              min(otherVectors(idx,:)) +...
                              offset;
        offset = 1.05*max(max(vectors(kk,:)),...
                          max(otherVectors(idx,:)));
       end
   else
        offset = max( vectors(kk,:) ) * 1.05;   
   end
end

end