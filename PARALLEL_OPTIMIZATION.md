# Parallel Processing Optimization for RunIGPModel

## Overview

The `RunIGPModel` function has been optimized to use multi-core processing by default to reduce execution time while maintaining same model outputs. (Tested on default dataset with 2 chains, original: ~18.5mins; optimized: ~12mins)

## Key Improvements

### 1. New Parameters

- `use.parallel = TRUE` (default): Enable/disable parallel processing
- `n.cores = NULL` (default): Number of cores to use (auto-detected if NULL)

### 2. Automatic Setup

- Automatically detects available CPU cores
- Sets up parallel backend using `doParallel` and `foreach`
- Properly manages cluster cleanup
- Installs missing packages automatically

### 3. Core Allocation

- Uses `min(number_of_chains, available_cores - 1)` cores
- Leaves one core free for system processes

## Usage Examples

### Basic Usage (Parallel by Default)

```r
# This will automatically use parallel processing
RunIGPModel(data.raw = data.raw, interval = 50)
```

### Custom Core Configuration

```r
# Use specific number of cores
RunIGPModel(data.raw = data.raw,
            interval = 50,
            n.cores = 4)
```

### Disable Parallel Processing

```r
# Force sequential processing
RunIGPModel(data.raw = data.raw,
            interval = 50,
            use.parallel = FALSE)
```


### Backward Compatibility

- All existing function calls continue to work unchanged
- Default behavior now uses parallel processing
- Original sequential mode available via `use.parallel = FALSE`


### Model Output Verification

The optimized function produces **identical results** to the original:

- Same MCMC chains with identical seeds
- Same parameter estimates and credible intervals
- Same diagnostic outputs
- Same file structure and naming

Only the execution time is improved - all scientific outputs remain unchanged.

### Memory Considerations

- Ensure adequate RAM for multiple simultaneous JAGS processes
