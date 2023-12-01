To use the resp dataset for the macros, load resp.sas7bdat into SAS and run the following SAS code:

data resp;
set resp;
SexNum = (Sex = "M");
UniqID = Center*100 + ID;
Trt = (Treatment = "A");
run;

This code produces
- a numeric variable for sex
- a unique patient identifier for each patient in the dataset
- a numeric variable for treatment

To use the skin dataset for the macros, load skin.sas7bdat into SAS and run the following SAS code:

data skin;
set skin;
Stage4 = (STAGE = 4);
Stage5 = (STAGE = 5);
ID = _N_;
if INV = 5 then center = 1;
else if INV = 6 then center = 2;
else if INV = 8 then center = 3;
else if INV = 9 then center = 4;
else if INV = 10 then center = 5;
else center = 6;
center2 = center;
if center = 4 then center2 = 3;
run;

This code produces
- two numeric variables for stage
- a unique patient identifier for each patient in the dataset
- a center variable
- a center2 variable which pools two of the original centers
