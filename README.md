# Double-K-S-potential-using-MATLAB-MTEX
MATLAB/MTEX codes for calculating the double Kurdjumov–Sachs (K-S) potential along the main BCC rolled texture fibers (α- and γ-fibers) to evaluate texture influence on austenite nucleation propensity (selective nucleation).

The codes were developed to evaluate the double K-S potential along the main cold-rolled and partially recrystallized BCC texture fibers (α- and γ-fibers) and to examine the relationship between ferrite crystallographic texture and the propensity for austenite nucleation.

The methodology and mathematical definition of the double K-S potential are described in detail in the associated paper. This repository is mainly provided to make the computational procedure used in the study available and easier to apply to other research.

## Codes

Two MATLAB scripts are provided:

double_KS_potential_alpha_fiber.m — Calculates the double K-S potential along the α-fiber (<110> || RD).
double_KS_potential_gamma_fiber.m — Calculates the double K-S potential along the γ-fiber ({111} || ND).

The codes calculate the ferrite (BCC) ODF from the imported EBSD data and then evaluate the double K-S potential along the selected texture fiber.
The output includes the double K-S potential for each reference α₁ orientation along the fiber and CSV files containing the corresponding candidate α₂ orientations and their corresponding calculations.

### Before Running the Codes

The EBSD data import is not included directly in these scripts because the crystal symmetry, phase order, and specimen coordinate system can be different for different EBSD datasets.
The EBSD data should first be imported using the MTEX Import Wizard. The Import Wizard automatically generates an import script that can be used to load the EBSD dataset.
After importing the data, make sure that the ferrite phase is named "BCC" in MATLAB. The code uses ebsd('BCC') to access the ferrite orientations.

## The general workflow is:

Import the EBSD dataset using the MTEX Import Wizard.
Check the crystal phases and assign the ferrite phase the name BCC.
Check the specimen reference frame and apply the required rotation, if required.
Run the appropriate double K-S potential code.

## Output

the double K-S potential along the selected texture fiber;
CSV files containing the candidate α₂ orientations;
ODF-based statistics for the candidate α₂ orientations;
Plot showing the variation of the double K-S potential along the selected fiber.

The main output files are saved in a results folder created by the script.


Note on the Code

The original calculation code was developed for the analysis presented in the associated study. The code was subsequently revised with an AI assistant to improve its organization, comments, and readability for sharing with other researchers.

## Data Availability

The experimental EBSD datasets used in the study are not included in this repository.

## Citation

## Contact
M. Milad Jafarzad-Shayan (miladjshayan@postech.ac.kr)

For questions regarding the double K-S potential methodology or the associated study, please refer to the corresponding author of the paper.

For questions regarding the double K-S potential methodology or the
associated study, please refer to the corresponding author of the
paper.
