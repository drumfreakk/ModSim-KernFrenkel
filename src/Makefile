level := 3

filename := kernfrenkel

CFLAGS := 

# Variables to run through
OUT_BASE ?= out

# Those satisfy the SBPP condition
# COS_LIST ?= 0.97 0.98 0.99
# LAMBDA_LIST ?= 1.02 1.12 1.22 1.32 1.42 1.52 1.62 1.72 1.82 1.92 2.02

# # F.Romano tetrahedral structure, delta 0.07
# TEMPERATURES ?= 0.01
# DENSITIES ?= 0.1 0.5 0.72

# PRESSURES ?= 10


# # current setup for finding phase diagrams
# TEMPERATURES ?= 0.12
# DENSITIES ?= 0.75



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

generate_coords: filename:=generate_coordinates
generate_coords: level := 3
generate_coords: run 

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
lambda_coshalfangle: out
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





STARTING_CRYSTALS = FCC BCC HCP DIAMOND
TEMPERATURES ?= 0.08 0.16
PRESSURES ?= 4 12
# Run through: crystals, temperature, pressure
# I added a possibility of adding simple cubic but we need to make those files to run it, its an interesting thing to look at tho
temp_pressure_measurement: CFLAGS += -DNPT
temp_pressure_measurement: out
	@for crystal in $(STARTING_CRYSTALS); do \
        case "$$crystal" in \
            FCC) initfile="../fcc/500_rot.snap" ;; \
            BCC) initfile="../bcc/432_rot.snap" ;; \
            HCP) initfile="../hcp/252_rot.snap" ;; \
            SC) initfile="../sc/500_rot.snap" ;; \
            DIAMOND) initfile="../dc/512_rot.snap" ;; \
            *) echo "Configuration for $$crystal not found!"; exit 1 ;; \
        esac; \
        echo "Starting simulation for $$crystal (INIT: $$initfile)"; \
        for temp in $(TEMPERATURES); do \
            for p in $(PRESSURES); do \
                out_dir=$(OUT_BASE)/$$crystal/T_$${temp}_P_$${p}; \
                mkdir -p $$out_dir; \
                ./build/out --initfile $$initfile --temperature $$temp --pressure $$p --outdir $$out_dir; \
            done; \
        done; \
    done



STARTING_CRYSTALS = FCC BCC HCP DIAMOND
COS_LIST ?= 0.97 0.98 0.99
LAMBDA_LIST ?= 1.02 1.12

# Run through: crystals, lambda interaction lenght, coshlafangle
lambda_coshalfangle_measurement: CFLAGS += -DNPT
lambda_coshalfangle_measurement: out
	@for crystal in $(STARTING_CRYSTALS); do \
        case "$$crystal" in \
            FCC) initfile="../fcc/500_rot.snap" ;; \
            BCC) initfile="../bcc/432_rot.snap" ;; \
            HCP) initfile="../hcp/252_rot.snap" ;; \
            SC) initfile="../sc/500_rot.snap" ;; \
            DIAMOND) initfile="../dc/512_rot.snap" ;; \
            *) echo "Configuration for $$crystal not found!"; exit 1 ;; \
        esac; \
		echo "Starting simulation for $$crystal (INIT: $$initfile)"; \
		for lambda in $(LAMBDA_LIST); do \
			for coshalfangle in $(COS_LIST); do \
				echo "lambda=$$lambda cos=$$coshalfangle"; \
				out_dir=$(OUT_BASE)/$$crystal/lambda_$${lambda}_coshalfangle_$${coshalfangle}; \
				mkdir -p $$out_dir; \
				./build/out --initfile $$initfile --patchdistance $$lambda --coshalfangle $$coshalfangle --outdir $$out_dir; \
			done; \
		done; \
    done


STARTING_CRYSTALS = FCC BCC HCP DIAMOND
TEMPERATURES ?= 0.12
DENSITIES ?= 0.75
# Run through: crystals, temperature, density
temp_density_measurement: out
	@for crystal in $(STARTING_CRYSTALS); do \
        case "$$crystal" in \
            FCC) initfile="../fcc/500_rot.snap" ;; \
            BCC) initfile="../bcc/432_rot.snap" ;; \
            HCP) initfile="../hcp/252_rot.snap" ;; \
            SC) initfile="../sc/500_rot.snap" ;; \
            DIAMOND) initfile="../dc/512_rot.snap" ;; \
            *) echo "Configuration for $$crystal not found!"; exit 1 ;; \
        esac; \
        echo "Starting simulation for $$crystal (INIT: $$initfile)"; \
        for temp in $(TEMPERATURES); do \
			for rho in $(DENSITIES); do \
				out_dir=$(OUT_BASE)/$$crystal/T_$${temp}_rho_$${rho}; \
				mkdir -p $$out_dir; \
				./build/out --initfile $$initfile --temperature $$temp --density $$rho --outdir $$out_dir; \
			done; \
    	done; \
    done