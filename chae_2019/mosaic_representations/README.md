DATA DESCRIPTION


Property strength vectors are located in the ”dragon_properties” folder:

“dragon_property_descriptions.mat” contains the list of the physical-chemical properties names (abbreviations).

“dragon_property_matrix.mat” contains the values taken by these properties across the set of 49 monomolecular odorants used in this study. Rows denote properties and columns denote odorants used.




All throughout for neural data, rows denote odors and columns ROIs (glomeruli, mitral or tufted cells). 

Glomerular data has been split into two folders:

 "glomerular_responses__main" 
contains the glomerular responses from awake head fixed mice as described in the text.

"glomerular_responses__awake_vs_anesthetized_control"
contains the data used to compare tuning of glomerular activity to physical-chemical properties in awake vs. anaesthetized mice (Supplementary Figure 5a-d). 

In both cases, each folder contains one folder per animal using the following convention: "animal_#(number)" or "animal_#(number)_(condition)", where (number) refers to the corresponding animal identifier number, and (condition) ---when applicable--- describes if the animal was awake or anesthetized. 

Each folder, in turn, contains two subfolders, one per each hemibulb:
1) The file "(side)_hemibulb_ROI_responses.mat", where (side) corresponds to "left" or "right", contains the dR/R responses for each odor and ROI (i.e each putative glomerulus). The values of each cell correspond to the dR/R values; those are expected to be negative (<0) since intrinsic imaging measures the decrease in reflectance to IR light correlated with neural activity. For example, if for the matrix "right_hemibulb_ROI_responses" the 8th row, 3rd column value is -0.0014, that means that the dR/R response of ROI #3 to odor #8 was -0.0014. The identities of each odor can be found in the accompanying file "glomerular_odor_panel.xlsx" and respectively “MT_cell.xlsx”. 

    The dR/R responses here are provided as obtained. For data analysis, non-significant responses were set to zero using a threshold value of t = -0.00045 (i.e. all the odor responses whose amplitude was ABOVE this threshold were set to 0). For more details, see "Odor responses" in the Methods section of the manuscript. 

 	2) The file "(side)_hemibulb_ROI_descriptors.mat", where (side) corresponds to "left" or "right", contains the physical location descriptors for each ROI. Each column corresponds to the ROI identifier numbers, which are also shared with the "(side)_hemibulb_ROI_responses.mat" file. 
    Numbers in the "(side)_hemibulb_ROI_descriptors" matrix represent measurements in pixels.
 
    *) the 1st row represents the average of X-coordinate values of all pixels from a given ROI; 
    *) the 2nd row represents the minimum value of X-coordinate values of all pixels from a given ROI (i.e. is the X-coordinate of the upper-leftmost pixel of that ROI);  
    *) the 3rd row represents the maximum value of X-coordinate values of all pixels from a given ROI (i.e. is the X-coordinate of the lower-rightmost pixel of that ROI);  
    *) the 4th row represents the minimum value of Y-coordinate values of all pixels from a given ROI (i.e. is the Y-coordinate of the upper-leftmost pixel of that ROI);  
    *) the 5th row represents the maximum value of Y-coordinate values of all pixels from a given ROI (i.e. is the Y-coordinate of the lower-rightmost pixel of that ROI); 
    *) the 6th row represents the average of the Y-coordinate values of all the pixels from a given ROI; 
    
    As values were obtained from digitized versions of the images of the intrinsic optical imaging responses to odorants, the X-axis points from left to right, and the Y-axis points from top to bottom. In all cases, the rostral side of the animal was on the right side of the image and the caudal side on the left. For the X-axis, the pixel's width is 9.03 +/- 0.05 µm; for the Y-axis, the pixel's height is 11.28 +/- 0.06 µm.   


ADDITIONAL COMMENTS: 
*) The data in "glomerular_responses__main_text" corresponds to 5 mice, while the data in "glomerular_responses__awake_vs_anesthetized_control" corresponds to 4 additional mice. 

ROIs identification numbers, locations and sizes were assigned independently for all the datasets, and thus ROI identifier numbers across datasets are unrelated.



Mitral and tufted cell data:

M/T cell data has been split into two folders: “mitral_cells” and “tufted_cells”
“mitral_cells” contains two data sets sampling 55 odorants across 2 different nominal oil dilutions, 1:100 and 1:3,000, and 33 odorants with 1:100 oil dilution. The data is split into two folders, “55_odors” and “33_odors”. 
“tufted_cells” contains tufted cell responses sampling 55 odorants across 2 different nominal oil dilution, at 1:100 and 1:3,000. 

 “mitral_cells/55_odors“ and ”tufted_cells” contain .mat files for each individual field of view. 
The file name “(cell identity-m or t_concentration identity-h or l)_field_of_view_(field number).mat” , where  “(cell_concentration identity)” corresponds to mitral cells (m), tufted cells (t), 1:100 oil dilution (h) and 1:3,000 oil dilution (l).  The Odorants list, animal numbers, and which fields of view were sampled across 2 different concentrations can be found in the file MT_cell.xlsx. 

“mitral_cell/33_odors“ contains mitral cell responses to 33 odorants (Supplementary Figure 2.a).  The Odorants list, animal numbers can be found in the file MT_cell.xlsx. 


Data description (mat file)
Each folder contains mat files for each individual field of view.  

R : mean odor responses (dF/F) across repeats. 
P: p-value for a balanced two-way ANOVA for comparing between the mean responses of air period across repeats and the mean responses of odor period across repeats.   We used p<0.1 for threshold to determine responsive odor responses in further analysis. Odor responses (R) are replaced by 0 if p≥0.1 (see Methods).
Pixed_size_x: Pixel size (x) in µm
Pixel_size_y: Pixel size (y) in µm
x: centroid of ROI the x coordinate
y: centroid of ROI the y coordinate
xy  
1) 1st row: ROI x left border of FOV 
2) 2nd row: ROI x right border of FOV
3) 3rd row: ROI y top border of FOV
4) 4th row: ROI y bottom border of FOV
