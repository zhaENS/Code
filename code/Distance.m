function D = Distance(FirstBeed,LastBeed,MetBeed)


D = zeros(3,1);

D(1) = sqrt(sum((FirstBeed-LastBeed).^2));
D(2) = sqrt(sum((FirstBeed-MetBeed).^2));
D(3) = sqrt(sum((MetBeed-LastBeed).^2));





end