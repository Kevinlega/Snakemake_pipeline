################################################################################
#                                                                              #
# Copyright 2017 Megaprobe-Lab                                                 #
#                                                                              #
# This is software created by the megaprobe lab under the GPL3 license.        #
#                                                                              #
# This program reproduces the pipeline for de novo rna sequencing research. To #
# run the program, just move into the folder that contains the makefile and    #
# make sure the sources folder contains the source files for the compilation   #
# of Trinity, khmer, Trimmomatic and ________.                                 #
#                                                                              #
# Then utilize the command: "make" to run the desired operation.               #
################################################################################
# Flags
SHELL:=/bin/bash
# Working variables :)
CPU=$(grep -c processor /proc/cpuinfo)
MEM=$(free available -h | grep "Mem:" | awk '{ print $3 }')


# Make all targets
all: QualityControl Trinity Comparison clean specs
	@echo "Currently making the $@ recipee."
	@echo "Hello $(USER), if unsure how to run this file try make help for specific options."

# Testing area:
# Currently testing:
#		Variable expansion in makefile to figure out cpu arquitecture and operating system.
PHONY:
	@echo "Currently making the $@ recipee."
	@echo "Cpu: " $(value CPU)
	@echo "Mem: " $(value MEM)

# Non-working computer specs recipee
specs:
	@echo "Currently making the $@ recipee."
	@echo "You have $(value CPU) available cores and $(value MEM) of RAM.\n"

# Choose preffered style help vs help2
help:
	@echo "Currently making the $@ recipee."
	@echo "Usage: 'make target' \n"
	@echo "Possible make targets are:"
	@echo "	all"
	@echo " 	Makes all the following targets."
	@echo	"	QualityControl"
	@echo " 	Proper Description"
	@echo "	Trinity"
	@echo " 	Proper Description"
	@echo "	Comparison"
	@echo " 	Proper Description"
	@echo "	clean"
	@echo " 	Proper Description\n"

# Choose preffered style help vs help2
help2:
	@echo "Currently making the $@ recipee."
	@echo "Usage: 'make target' \n"
	@echo "Possible make targets are:"
	@echo "	all QualityControl Trinity Comparison clean\n"
	@echo "\nFor more help with a specific target try help_target"

# Specific help for QualityControl recipee
help_QualityControl:
	@echo "Currently making the $@ recipee."
	@echo "To perform Quality Control by normalizing and trimming the rna seq file data."
	@echo "\nThe required inputs are: \n"

QualityControl:
	@echo "Currently making the $@ recipee."
# 	Trimmomatic
# 	Normalize
# 	the other script

# 	script to calculate # of lines and decide how many CPU is needed
# 	run it
# LIST = $(shell lscpu)
# for i in $(LIST);
#		do echo "\t $$i \n";
# done
CALCCPU=5
CALCMEM=3
Trinity: specs
	@echo "Currently making the $@ recipee."
	ifeq ($(CPU), $(CALCCPU))
		CPUARGS:= CALCCPU
	else
		CPUARGS:= CPU
	endif
	@echo CPUARGS
ifeq ($CPU==$CALCMEM )
	MEMARGS:= CALCRAM
else
	MEMARGS:= MEM
endif
	@echo CPUARGS

# All fasta and fastq files
	# LIST = one two three
	# for i in $(LIST); do \
	# 	echo $$i; \
	# done

Comparison:
	@echo "Currently making the $@ recipee."
	@echo "Looking for fata files:"
	@find -O3 . -name "*.fasta"
	@find -O3 . -name "*.fa"
	@find -O3 . -name "*.seq"

clean:
	@echo "Currently making the $@ recipee."
	@echo "Cleaning Options"
