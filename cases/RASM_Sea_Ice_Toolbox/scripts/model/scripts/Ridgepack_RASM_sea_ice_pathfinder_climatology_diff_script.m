


%rasmcases={'R1009RBRceap01a','R1009Gceap01d','R1009RBRcevp01a','R1009Gcevp01d'};
%quicknames={'EAP','EAP-CORE2','EVP','EVP-CORE2'};

%pubdir='/Users/aroberts/science/publications/2015_Anisotropic/Figures' 

%rasmcases={'R1009RBRceap01a','R2100aRBRcaaa01a','R2100bRBRcaaa01a'};
%quicknames={'RASM 1.1','RASM 2.0.a','RASM 2.0.b'};

%rasmcases={'R1009RBRceap01a','R2100aRBRcaaa01a','R2100aRBRcpon01a','R2100aRBRctop01a'};
%quicknames={'RASM 1.1','RASM 2.1.a','RASM 2.1.pon','RASM 2.1.top'};

rasmcases={'R1009RBRceap01a','R2100aRBRcaaa01a'};
quicknames={'RASM 1.1','RASM 2.1'};

pubdir='/Users/aroberts/science/publications/2018_RASM/Figures'

startyear=1990
endyear=2009

format=2;

pub=false;
%pub=true;

Ridgepack_RASM_sea_ice_pathfinder_stream_diff(rasmcases,quicknames,startyear,endyear,format,pub,pubdir)

