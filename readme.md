# Simulation Workflow Instructions

This document provides instructions for running simulations in a SLURM-enabled environment. SLURM is required to execute the following commands as they use job scheduling and parallelization features. For users without SLURM, an alternative approach using standard R scripts or shell loops is outlined at the end of this document.

---

## Prerequisites/Initial Setup

1. **SLURM Environment**: Ensure SLURM is installed and properly configured on your system.
2. **R Environment**: Install R version 4.2.2 and all required libraries.
3. **Shell Access**: Have access to a command-line terminal.  
4. **Directory Structure**: Ensure that all folders mentioned in this document are placed in the same parent directory to maintain consistent file paths. Change file paths in `run.R` files and `compile_result` files accordingly.

---

## Simulations Overview

### Simulation 1: Standard Compositional Covariate Settings (Dirichlet Distribution)
1. Navigate to the directory `run103`:
   ```bash
   cd run103
   ```
2. Submit jobs using SLURM:
   ```bash
   sbatch --array=1-3000 run.sh
   ```
3. Compile results:
   ```bash
   Rscript compile_result_presentation.R
   ```
4. Navigate to `run106`:
   ```bash
   cd run106
   ```
5. Submit jobs:
   ```bash
   sbatch --array=1-500 run.sh
   ```
6. Compile results:
   ```bash
   Rscript compile_result_presentation.R
   ```

### Simulation 2: Standard Compositional Covariate Settings (Logistic-Normal Distribution)
1. Navigate to `run105`:
   ```bash
   cd run105
   ```
2. Submit jobs:
   ```bash
   sbatch --array=1-6500 run.sh
   ```
3. Compile results:
   ```bash
   Rscript compile_result_presentation.R
   ```
4. Navigate to `run108`:
   ```bash
   cd run108
   ```
5. Submit jobs:
   ```bash
   sbatch --array=1-500 run.sh
   ```
6. Compile results:
   ```bash
   Rscript compile_result_presentation.R
   ```

### Simulation 3: Non-Compositional Covariate Settings (Normal Distribution)
1. Navigate to `run103`:
   ```bash
   cd run103
   ```
2. Submit jobs:
   ```bash
   sbatch --array=3000-6000 run.sh
   ```
3. Compile results:
   ```bash
   Rscript compile_result_presentation.R
   ```
4. Navigate to `run106_1`:
   ```bash
   cd run106_1
   ```
5. Submit jobs:
   ```bash
   sbatch --array=1-500 run.sh
   ```
6. Compile results:
   ```bash
   Rscript compile_result_presentation.R
   ```

### Simulation 4: Sparse Compositional Covariate Settings
1. Navigate to `run109`:
   ```bash
   cd run109
   ```
2. Submit jobs:
   ```bash
   sbatch --array=1-500 run.sh
   ```
3. Compile results:
   ```bash
   Rscript compile_result_presentation.R
   ```
4. Navigate to `run107`:
   ```bash
   cd run107
   ```
5. Submit jobs:
   ```bash
   sbatch --array=1-50 run.sh
   ```
6. Compile results:
   ```bash
   Rscript compile_result.R
   ```

### Simulation 5: Robustness Study
1. Navigate to `run94`:
   ```bash
   cd run94
   ```
2. Submit jobs:
   ```bash
   sbatch --array=1-7500 run.sh
   ```
3. Compile results:
   ```bash
   Rscript combined_plot.R
   ```
4. Navigate to `run107`:
   ```bash
   cd run107
   ```
5. Submit jobs:
   ```bash
   sbatch --array=1-500 run.sh
   ```
6. Compile results:
   ```bash
   Rscript compile_result_presentation.R
   ```

### Simulation 6: Bonferroni p-values
#### Simulation 6.1: Standard Compositional Covariate Settings (Dirichlet Distribution)
1. Navigate to `run110`:
   ```bash
   cd run110
   ```
2. Submit jobs:
   ```bash
   sbatch --array=1-7500 run.sh
   ```
3. Compile results:
   ```bash
   Rscript compile_result_bonferroni.R
   ```
4. Navigate to `run111`:
   ```bash
   cd run111
   ```
5. Submit jobs:
   ```bash
   sbatch --array=1-1000 run.sh
   ```
6. Compile results:
   ```bash
   Rscript compile_result_bonferroni.R
   ```

#### Simulation 6.2: Standard Compositional Covariate Settings (Logistic-Normal Distribution)
1. Navigate to `run112`:
   ```bash
   cd run112
   ```
2. Run R script:
   ```bash
   Rscript Run_from_pval.R
   ```
3. Compile results:
   ```bash
   Rscript compile_result_bonferroni.R
   ```
4. Navigate to `run108`:
   ```bash
   cd run108
   ```
5. Submit jobs:
   ```bash
   sbatch --array=1-1000 run.sh
   ```
6. Compile results:
   ```bash
   Rscript compile_result_bonferroni.R
   ```

#### Simulation 6.3: Non-Compositional Covariate Settings (Normal Distribution)
1. Navigate to `run113`:
   ```bash
   cd run113
   ```
2. Run R script:
   ```bash
   Rscript run_from_pval.R
   ```
3. Compile single condition results:
   ```bash
   Rscript compile_result_normal_bonferroni_single.R
   ```
4. Compile multiple condition results:
   ```bash
   Rscript compile_result_normal_bonferroni_multiple.R
   ```

#### Simulation 6.4: Sparse Compositional Covariate Settings
1. Navigate to `run114`:
   ```bash
   cd run114
   ```
2. Run R script:
   ```bash
   Rscript run_from_pval.R
   ```
3. Compile results:
   ```bash
   Rscript compile_result_bonferroni_conditioning_multiple.R
   ```
4. Navigate to `run107`:
   ```bash
   cd run107
   ```
5. Submit jobs:
   ```bash
   sbatch --array=1-50 run.sh
   ```
6. Compile results:
   ```bash
   Rscript compile_result_bonferroni.R
   ```

### Simulation 7: Performance of adaFilter
1. Navigate to `run103`:
   ```bash
   cd run103
   ```
2. Submit jobs:
   ```bash
   sbatch --array=1-3000 run.sh
   ```
3. Compile results:
   ```bash
   Rscript compile_result_presentation_adafilter.R
   ```

### Simulation 8: Computational Speedups
1. Navigate to `run104_1`:
   ```bash
   cd run104_1
   ```
2. Submit jobs:
   ```bash
   sbatch --array=1-1500 run.sh
   ```
3. Compile results:
   ```bash
   Rscript compile_result_presentation.R
   ```

### Real data analysis
**Data sources:**
- Paper: [Biostatistics Supplementary Data](https://academic.oup.com/biostatistics/article/22/4/687/5689688#supplementary-data)  
- Data: [HMP dataset on GitHub](https://github.com/diegotomassi/sdr4comp/tree/master/data/HMP)

1. Navigate to `run127`:
   ```bash
   cd run127
   ```
2. Run:
   ```bash
   Rscript Run.R
   ```

---

## Running Without SLURM

To run simulations without SLURM, ensure you input the number of cores available on your device in the script and replace the `run.sh` file with the following script:

```bash
#!/bin/bash

# Parse custom range argument
if [[ $1 =~ --array=([0-9]+)-([0-9]+) ]]; then
  start=${BASH_REMATCH[1]}
  end=${BASH_REMATCH[2]}
else
  echo "Error: Please provide a valid range in the format --array=start-end"
  exit 1
fi

# Set the number of cores manually
n_cores=8
max_jobs=$((n_cores - 1))

echo "Using $max_jobs cores for parallel execution"
echo "Running jobs for range $start to $end"

# Function to run jobs with concurrency control
run_jobs() {
  for i in $(seq $start $end); do
    Rscript Run.R --array_id=$i &
    # Limit background jobs to max_jobs at a time
    if (( $(jobs -r | wc -l) >= max_jobs )); then
      wait
    fi
  done
  wait  # Wait for all background jobs to finish
}

# Execute the job function
run_jobs
```

### Steps to Execute

1. Replace the `run.sh` file in the relevant folder with the above script.
2. Navigate to the current folder in the terminal. For example:
   ```bash
   cd run105
   ```
3. Modify the `Run.R` file to handle command-line arguments instead of SLURM environment variables. Replace this line:
   ```r
   taskid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
   ```
   with this line:
   ```r
   args <- commandArgs(trailingOnly = TRUE)
   taskid <- as.numeric(sub("--array_id=", "", args[1]))
   ```
4. Make the script executable:
   ```bash
   chmod +x run.sh
   ```
5. Run the command with your desired array range:
   ```bash
   ./run.sh --array=1-500
   ```

This will execute the jobs specified in the given range without SLURM.

---

### Important Considerations

1. This non-SLURM workaround is just one of many approaches. There may exist more efficient or advanced alternatives.

2. Using this script may lead to slower execution compared to SLURM. In SLURM, tasks are scheduled dynamically, and the next task starts as soon as one finishes. However, with this script, the entire set of tasks in the current batch must complete before scheduling the next batch of jobs.

3. A common issue that might arise is a conflict with the R version in a Conda environment if present. Running `conda deactivate` in the terminal before executing the script might resolve this issue.

---
