
# Entries to modify as needed

INSTALLDIR := /usr/local/bin

# For ATLAS BLAS
LIBS := -lf77blas -llapack -latlas -lm -lgsl -lgslcblas
# For Intel MKL BLAS
# LIBS := -lmkl_gf_lp64 -lmkl_sequential -lmkl_lapack -lmkl_core \
# -lm -lgsl -lgslcblas

INCLUDES := -I/usr/include/gsl -I/usr/include

CFLAGS := -O3 -Wall


# Do not modify below this line

BUILDDIR := build
SRCDIR := src
TESTDIR := test

RM := rm -rf

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS := $(wildcard $(SRCDIR)/*.c)

OBJS := $(C_SRCS:$(SRCDIR)/%.c=$(BUILDDIR)/%.o)

C_DEPS := $(C_SRCS:$(SRCDIR)/%.c=$(BUILDDIR)/%.d)

BINARY := rowavedt


$(BUILDDIR)/%.o: $(SRCDIR)/%.c
	mkdir -p $(BUILDDIR)
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc $(INCLUDES) -std=gnu99 \
$(CFLAGS) -c -fmessage-length=0 -MMD -MP \
	-MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


USER_OBJS :=

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif


# All Target
all: rowavedt

# Tool invocations
rowavedt: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C Linker'
	gcc  -o"$(BINARY)" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

install: rowavedt
	@echo "Installing to $(INSTALLDIR)"
	cp $(BINARY) $(INSTALLDIR)

# Basic test
.PHONY : test
test: rowavedt
	mkdir -p test/
	Rscript scripts/mk_wavelet_basis.R -c 128 > test/basis.dat
	./rowavedt test/basis.dat 2048 128 \
		data/y.dat `wc -l data/y.dat` \
		> test/output.txt
	./rowavedt test/basis.dat 2048 128 \
		data/yMissing.dat `wc -l data/yMissing.dat` \
		>> test/output.txt
	Rscript scripts/screen_time_series.R --alpha=0.0001 test/output.txt \
		test/detections.txt test/detection_stats.txt
	Rscript scripts/compute_features.R --detections=test/detections.txt \
		test/output.txt test/basis.dat \
		> test/features.txt

# Other Targets
.PHONY : clean
clean:
	-$(RM) $(OBJS)$(C_DEPS)$(EXECUTABLES) rowavedt
	-@echo ' '

