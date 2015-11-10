# CAMB Test/Benchmark/Profile (TBP) suite build system
#
# Developed by:
#   Marco Raveri (mraveri@sissa.it)
#

# directory structure of the test suite:
TEST_DIR=$(abspath .)
CAMB_DIR=$(TEST_DIR)/..
SCRIPT_DIR=$(TEST_DIR)/script
RESULTS_DIR=$(TEST_DIR)/results
LOG_DIR=$(TEST_DIR)/logs
LEGACY_DIR=$(TEST_DIR)/legacy_results

# standard targets:

default: all

all: dir create_spectra run_comparison run_benchmarks run_profile

clean: clean_spectra clean_comparison clean_benchmarks clean_profile

deep_clean:
	@rm -rf $(LEGACY_DIR) $(LOG_DIR) $(RESULTS_DIR)

# directory structure:

dir: $(LEGACY_DIR) $(LOG_DIR) $(RESULTS_DIR)

$(LEGACY_DIR):
	@mkdir -p $(LEGACY_DIR)
$(LOG_DIR):
	@mkdir -p $(LOG_DIR)
$(RESULTS_DIR):
	@mkdir -p $(RESULTS_DIR)

# target to create the legacy results:

create_legacy: $(LEGACY_DIR) create_spectra
	@echo "********************************************"
	@echo "CAMB TBP: created new legacy results"
	@echo "********************************************"
	@mv -f $(RESULTS_DIR)/* $(LEGACY_DIR)/
	@rm -rf $(RESULTS_DIR)

clean_legacy:
	@rm -rf $(LEGACY_DIR)/*
	@rm -rf $(LOG_DIR)/create_legacy.log

# target to run camb to get the spectra:

create_spectra: $(RESULTS_DIR) $(LOG_DIR)
	@./script/create_spectra.sh

clean_spectra:
	@rm -rf $(RESULTS_DIR)/*
	@rm -rf $(LOG_DIR)/create_spectra.log

# target to run the actual comparison between spectra and legacy spectra:

run_comparison: dir
	@./script/compare_legacy.sh

clean_comparison:
	@rm -rf $(LOG_DIR)/compare_legacy.log

# target to run the benchmarks:

run_benchmarks: $(LOG_DIR)
	@./script/run_benchmarks.sh

clean_benchmarks:
	@rm -rf $(LOG_DIR)/run_benchmarks.log

# target to run the profiler:

run_profile: $(LOG_DIR)
	@./script/run_profile.sh

clean_profile:
	@rm -rf $(LOG_DIR)/run_profile.log
