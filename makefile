#compiler
FC = gfortran

UNAME := $(shell uname)

# debugger flag
CFLAGS += -g
# netCDF flags
ifeq ($(UNAME), Linux) #netCDF libraries for WSL/Linux
ncdfCFLAGS += -I/usr/include
ncdfCFLAGS += -L/opt/local/lib
endif
ifeq ($(UNAME), Darwin) #netCDF libraries for MacOS
ncdfCFLAGS += -I/opt/homebrew/Cellar/netcdf-fortran/4.6.1/include
ncdfCFLAGS += -L//opt/homebrew/Cellar/netcdf-fortran/4.6.1/lib
endif
#ncdfCFLAGS += -lnetcdf
ncdfCFLAGS += -lnetcdff

# source files
SRCS = m_csv_io base_utils m_zrand m_random main
OBJS = $(SRCS:=.o)

# executable
MAIN = main.out

# compile the project
all : $(MAIN)
	@echo Model compiled succesfully

# build executable from object files
$(MAIN) : $(OBJS)
	$(FC) $(CFLAGS) -o $(MAIN) $(OBJS) $(ncdfCFLAGS)

.SUFFIXES : .o .f90

# compile fortran files to object files
.f90.o :
	$(FC) $(CFLAGS) -c $< $(ncdfCFLAGS)

# commands
clean :
	$(RM) *.o *.mod -r output/*.csv output/*.dat output/*.nc
ofiles :
	$(RM) *.o
run : 
	./$(MAIN)
rmdata :
	$(RM) output/*.nc output/*.dat output/*.csv
rmplots :
	$(RM) R/RPlots/*.pdf
warn :
	$(FC) $(CFLAGS) -c -Wall main.f90 $(ncdfCFLAGS)