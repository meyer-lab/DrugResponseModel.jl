""" Organize the combination experimental data into dataframes. """

### 1. control, Palbociclib 50, Palbociclib 50 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]:
### 2. control, Palbociclib 50, Palbociclib 50 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
### 3. control, Gemcitabine 10, Gemcitabine 10 nM + palbociclibs [25 nM, 50 nM, 100 nM, 250 nM]
### 4. control, Gemcitabine 10, Gemcitabine 10 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]
### 5. control, Lapatinib  100, Lapatinib  100 nM + palbociclibs [25 nM, 50 nM, 100 nM, 250 nM]
### 6. control, Lapatinib  100, Lapatinib  100 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
### 7. control, Doxorubicin 20, Doxorubicin 20 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
### 8. control, Paclitaxel 2, paclitaxel 2 nM + [palbo50, lap50, lap100, gem10 nM]
### 9. control, Doxorubicin 20, Dox 20 nM + [pax2, , gem10, palbo50, lap100]


# data import
gt1, gt2 = DrugResponseModel.import_combination("AU01001");
gt1_2, gt2_2 = DrugResponseModel.import_combination("AU01101");
gt1_3, gt2_3 = DrugResponseModel.import_combination("AU00901");

GS1 = cat(gt1, gt1_2, gt1_3, dims = 4);
GS2 = cat(gt2, gt2_2, gt2_3, dims = 4);

meanGS1 = mean(GS1, dims = 4);
meanGS2 = mean(GS2, dims = 4);
meanGS2[:, :, 19] .= mean(cat(gt2[:, :, 19], gt2_2[:, :, 19], dims = 3), dims = 3)[:, :, 1]

g = zeros(3, 193, 6, 9)
g[:, :, 1, :] .= meanGS2[:, :, 1, 1]
########### 1. Palbociclib 50 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]
g[:, :, 2, 1] .= meanGS1[:, :, 8, 1]
g[:, :, 3:6, 1] .= meanGS2[:, :, 3:6, 1]
# well 2: 3,4,5,6
# palbo50 alone: well 1: 8

########### Palbociclib 50 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
g[:, :, 2, 2] .= meanGS1[:, :, 8, 1]
g[:, :, 3:6, 2] .= meanGS2[:, :, 21:24, 1]
# well 2: 21, 22, 23, 24

########### Gemcitabine 10 nM + palbociclibs [25 nM, 100 nM, 250 nM]
g[:, :, 2, 3] .= meanGS1[:, :, 20, 1]
g[:, :, 3, 3] .= meanGS1[:, :, 18, 1]
g[:, :, 4, 3] .= meanGS1[:, :, 22, 1]
g[:, :, 5:6, 3] .= meanGS2[:, :, 23:24, 1]
# well 1: 18, 23, 24
# gem 10 alone: well 1: 20

########### Gemcitabine 10 nM + lapatinibs [25 nM, 50 nM, 100 nM, 250 nM]
g[:, :, 2, 4] .= meanGS1[:, :, 20, 1]
g[:, :, 3:6, 4] .= meanGS2[:, :, 9:12, 1]
# well 2: 9, 10, 11, 12
# gem 10 alone: well 1: 20


########### Lap 100 nM + gemcitabines [5, 10, 17 nM, 30 nM]
g[:, :, 2, 5] .= meanGS1[:, :, 4, 1]
g[:, :, 3, 5] .= meanGS2[:, :, 7, 1]
g[:, :, 4, 5] .= meanGS2[:, :, 11, 1]
g[:, :, 5, 5] .= meanGS2[:, :, 13, 1]
g[:, :, 6, 5] .= meanGS2[:, :, 19, 1]
# well 2: 7, 11, 13, 19
# lap 100 alone well 1 : 4

########### Lap 100 nM + palbociclibs [25 nM, 50, 100 nM, 250 nM]
g[:, :, 2, 6] .= meanGS1[:, :, 4, 1]
g[:, :, 3, 6] .= meanGS2[:, :, 8, 1]
g[:, :, 4, 6] .= meanGS2[:, :, 5, 1]
g[:, :, 5, 6] .= meanGS2[:, :, 14, 1]
g[:, :, 6, 6] .= meanGS2[:, :, 20, 1]
# well 2: 8, 5, 14, 20

########## Dox 20 nM + gemcitabines [5 nM, 10 nM, 17 nM, 30 nM]
g[:, :, 2, 7] .= meanGS1[:, :, 6, 1]
g[:, :, 3:6, 7] .= meanGS2[:, :, 15:18, 1]
# well 2: 15, 16, 17, 18
# dox 20 alone well 1 : 6 

########### Pax 2 nM + [palbo50, lap50 nM, lap100 nM, gem10]
g[:, :, 2, 8] .= meanGS1[:, :, 13, 1]
g[:, :, 3, 8] .= meanGS1[:, :, 14, 1]
g[:, :, 4, 8] .= meanGS2[:, :, 2, 1]
g[:, :, 5:6, 8] .= meanGS1[:, :, 16:17, 1]
# well 1: 14, well 2: 2, well 1: 16, 17 
# pax 2 alone well 1: 13

########### Dox 20 nM + [pax2, gem10, palbo50, lap100]
g[:, :, 2, 9] .= meanGS1[:, :, 6, 1]
g[:, :, 3, 9] .= meanGS1[:, :, 15, 1]
g[:, :, 4, 9] .= meanGS2[:, :, 16, 1]
g[:, :, 5, 9] .= meanGS1[:, :, 12, 1]
g[:, :, 6, 9] .= meanGS1[:, :, 11, 1]
# well 1: 15, well 2: 16, well 1: 12, 11 
