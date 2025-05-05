function v = get_v(t,D,TD,TI)

v=zeros(size(t));

for i=1:length(D)
    j=find( (t>=(i-1)*TD) & (t<(i-1)*TD+TI));
    v(j)=D(i)/TI;
end

end