# Load necessary libraries
library(mrgsolve)    # For PK/PD modeling and simulation
library(dplyr)       # For data manipulation
library(ggplot2)     # For plotting
library(tidyr)       # For data tidying (e.g., crossing)
library(cowplot)     # For arranging plots
library(scales)      # For plot scaling functions (e.g., label_number)
# library(RColorBrewer) # Optional, if using scale_color_brewer

# --- 1. Define the Model ---
# Model now includes calculation and capture of CP
model_string <- '
$PARAM
KA = 0,       // Set KA=0 for IV infusion
WT = 70,      // Default Weight (kg) - will be overridden by idata
CLcr = 30,    // Default Creatinine Clearance (mL/min) - will be overridden by idata
CL_pop = 3.09,  // Population typical value for CL (L/hr)
V1_pop = 32.1,  // Population typical value for V1 (L)
V2_pop = 113,   // Population typical value for V2 (L)
Q_pop = 39.7,   // Population typical value for Q (L/hr)
MIC = 2,      // Default MIC value (mg/L) - will be overridden by idata or param()
ALB = 25.8,   // Default Albumin (g/L) - will be overridden by idata
TBIL = 10,    // Default Total Bilirubin (umol/L ?) - will be overridden by idata
GGT = 30,     // Default Gamma-Glutamyl Transferase (U/L ?) - will be overridden by idata
// --- ELF: Add parameter for individual ELF/Plasma ratio (will be set in idata) ---
ELF_RATIO_indiv = 2.41, // Default value, will be overridden

$OMEGA @name OM @labels eta_CL eta_V1 eta_Q // Renamed eta_V2 in original model to eta_Q for consistency
// Variances for inter-individual variability (IIV)
0.073   // IIV on CL (Variance) -> Corresponds to CV% ~ sqrt(exp(0.073)-1) * 100 = 27.6%
0.525   // IIV on V1 (Variance) -> Corresponds to CV% ~ sqrt(exp(0.525)-1) * 100 = 83.1%
0.035   // IIV on Q (Variance) -> Corresponds to CV% ~ sqrt(exp(0.035)-1) * 100 = 18.9%

$MAIN
// Individual parameter calculations incorporating IIV and covariates
double CLcr_Lhr = CLcr * 0.06;
double CL = (CL_pop + CLcr/77.0 * 3.28) * exp(eta_CL); // CL in L/hr

// V1 calculation with WT covariate
double V1 = V1_pop * pow(WT/61.0, 1.95) * exp(eta_V1) ; // V1 in L

// V2 calculation with WT and ALB covariates
double V2 = V2_pop * pow(WT/61.0, 1.61) * pow(log(ALB)/log(25.8), 4.52) * exp(0); // V2 in L

// Q calculation with GGT and TBIL covariates
double Q = Q_pop * pow(log(GGT+1e-6)/log(30+1e-6), 0.956) * pow(log(TBIL+1e-6)/log(10+1e-6), -0.912) * exp(eta_Q); // Q in L/hr

$SIGMA
1 // Placeholder, residual error not used for PTA or these PK plots

$CMT CENT PERIPH // Define compartments (Central, Peripheral)
$ODE
dxdt_CENT = - (CL/V1) * CENT - (Q/V1) * CENT + (Q/V2) * PERIPH;
dxdt_PERIPH = (Q/V1) * CENT - (Q/V2) * PERIPH;

$TABLE
// Calculate concentration in the central compartment
double CP = CENT / V1;
// Capture individual CL for PTA calculation
double CL_indiv = CL;

// --- ELF: Calculate estimated ELF concentration using the INDIVIDUAL ratio ---
double CP_ELF_est = CP * ELF_RATIO_indiv; // Estimated ELF concentration

// Capture CP, CP_ELF_est, CL, individual ELF ratio, and other parameters
$CAPTURE CP CP_ELF_est CL_indiv ELF_RATIO_indiv WT CLcr ALB TBIL GGT CL V1 V2 Q MIC
'

# --- 2. Compile the Model ---
mod <- mcode("antibiotic_model", model_string)

# --- 3. Simulation Setup ---
set.seed(1234) # Set seed for reproducibility
n_sim <- 5000   # Number of individuals to simulate

# Define MIC values to test for PTA
mic_values_pta <- c(0.25, 0.5, 1, 2, 4, 8, 16)

# Define Dosing Regimens (for PTA and PK plots)
# List for PTA calculation (Daily Dose)
dosing_regimens_pta <- list(
    "TGC 100mg LD, 50mg q12h" = 50 * 2,  # 100 mg/day
    "TGC 150mg LD, 75mg q12h" = 75 * 2,  # 150 mg/day
    "TGC 200mg LD, 100mg q12h" = 100 * 2, # 200 mg/day
    "TGC 2.4mg/kg LD, 1.2mg/kg q12h" = NA # Weight-based dosing
)
regimen_names <- names(dosing_regimens_pta) # For consistent naming

# Define the target AUC/MIC ratio for PTA (Assuming same target applies to AUCelf/MIC)
auc_mic_target <- 6.93

# --- ELF: Define the ELF/Plasma ratio range for sampling ---
elf_ratio_min <- 2.41
elf_ratio_max <- 2.41

# Define the specific MIC to use for this plot
eff_tox_mic <- 2.0 # Example: Use MIC=1 for the plot

# Define the toxicity threshold (AUC > threshold)
eff_tox_auc_threshold <- 20.76 # Example: Toxicity if AUC > 100 mg*h/L

# --- 4. Create Individual Data (idata) ---
# Add individual ELF ratio sampled from uniform distribution
sim_data <- data.frame(
    ID = 1:n_sim,
    WT = 70,                
    CLcr = 120,               # Fixed Creatinine Clearance (mL/min)
    ALB = runif(n_sim, min = 25.8, max = 30.3),                # Fixed Albumin (g/L)
    TBIL = runif(n_sim, min = 10, max = 30.4), # Simulated Total Bilirubin
    GGT = runif(n_sim, min = 30, max = 103),  # Simulated Gamma-Glutamyl Transferase
    # --- ELF: Assign individual ELF ratio from uniform distribution ---
    ELF_RATIO_indiv = runif(n_sim, min = elf_ratio_min, max = elf_ratio_max)
)
