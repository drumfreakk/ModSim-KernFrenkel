level := 3

filename := kernfrenkel

CFLAGS := 

# Variables to run through
OUT_BASE ?= out

# Those satisfy the SBPP condition
COS_LIST ?= 0.97 0.98 0.99
LAMBDA_LIST ?= 1.02 1.12 1.22 1.32 1.42 1.52 1.62 1.72 1.82 1.92 2.02

# # Random place holder currently
# TEMPERATURES ?= 0.02 0.04
# DENSITIES ?= 0.3 0.5 0.7


# F.Romano tetrahedral structure, delta 0.07
TEMPERATURES ?= 0.01
DENSITIES ?= 0.1 0.5 0.72

PRESSURES ?= 10

run: out
	./build/out

out: $(filename).o
	gcc -Wall -g3 -O$(level) $(CFLAGS) -lm -I../inc/   ./build/$(filename).o -o ./build/out

$(filename).o: $(filename).c
	mkdir -p build
	mkdir -p out
	gcc -Wall -g3 -O$(level) $(CFLAGS) -lm -I../inc/ -c $(filename).c -o ./build/$(filename).o

dev: clean debug_compile

prod: clean build_compile

time: time_compile

debug_compile: CFLAGS += -DVERBOSE
debug_compile: level := g
debug_compile: run

build_compile: level := 3
build_compile: run

time_compile: filename:=timing
time_compile: level := 0
time_compile: run


clean:
	rm -f build/*


# Run through: patchdistance + coshalfangle
lambda_cos: out
	@for lambda in $(LAMBDA_LIST); do \
		for coshalfangle in $(COS_LIST); do \
			echo "lambda=$$lambda cos=$$coshalfangle"; \
			out_dir=$(OUT_BASE)/lambda_$${lambda}_coshalfangle_$${coshalfangle}; \
			mkdir -p $$out_dir; \
			./build/out --patchdistance $$lambda --coshalfangle $$coshalfangle --outdir $$out_dir; \
		done; \
	done

# Run through: temperature + density
temp_density: out
	@for temp in $(TEMPERATURES); do \
        for rho in $(DENSITIES); do \
            out_dir=$(OUT_BASE)/T_$${temp}_rho_$${rho}; \
            mkdir -p $$out_dir; \
            ./build/out --temperature $$temp --density $$rho --outdir $$out_dir; \
        done; \
    done


# Run through: temperature + pressure
temp_pressure: CFLAGS += -DNPT
temp_pressure: out
	@for temp in $(TEMPERATURES); do \
        for p in $(PRESSURES); do \
            out_dir=$(OUT_BASE)/T_$${temp}_P_$${p}; \
            mkdir -p $$out_dir; \
            ./build/out --temperature $$temp --pressure $$p --outdir $$out_dir; \
        done; \
    done
