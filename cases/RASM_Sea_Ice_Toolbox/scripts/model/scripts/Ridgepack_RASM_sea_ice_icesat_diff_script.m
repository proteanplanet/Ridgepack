


rasmcases={'R1009RBRceap01a','R2100aRBRcaaa01a'};
quicknames={'RASM 1.1','RASM 2.1'};

pubdir='/Users/aroberts/science/publications/2018_RASM/Figures'


%springfall={'spring','fall'};
springfall={'fall'};
columns=5;

pub=false;

for j=1:length(springfall)

 Ridgepack_RASM_sea_ice_icesat_diff(rasmcases,char(springfall{j}),quicknames,columns,pub)

end


