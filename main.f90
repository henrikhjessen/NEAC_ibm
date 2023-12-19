!This is the code for the individual based model, being adapted for NEA cod.
!Currently being written as part of my PhD thesis.
!Version control here: https://tegsvn.uib.no/svn/tegsvn/branches/henrik/ecoevo_ibm
!SVN ignore: "svn propset svn:ignore -RF svnignore.txt ."
!To compile use: gfortran -g m_csv_io.f90 base_utils.f90 m_zrand.f90 m_random.f90 main.f90
!-I/usr/include -L/opt/local/lib -L/opt/local/lib -lnetcdff    <- netCDF flags, place AFTER above line ^
!To check data after running: ncdump -h output/DATANAME.nc

PROGRAM evoeco_ibm
USE netcdf
USE csv_io
USE base_utils
USE m_zrand
USE m_random

IMPLICIT NONE
!----------------------------------------------------------------------------------------------------------------------------------
!SET INPUT PARAMETERS & FLAGS
!----------------------------------------------------------------------------------------------------------------------------------

!FLAGS
integer, parameter  :: WriteRes = 1                  !Write results? 1=yes, 0=no
integer, parameter  :: WriteCSV = 0                  !Do we also want the final population in csv format? 1=yes, 0=no
integer, parameter  :: WriteInterval = 20             !Write results every XXX years.
integer, parameter  :: NoRuns = 20                   !Number of runs to do. Each runs results are saved in separate files.

integer, parameter  :: ReadInitialPopulation = 2     !Read from file? 0=no, 1=continue last (finalpop), 2=use "startpop"
integer, parameter  :: SmallPop = 0                  !Smaller population for debugging, testing etc? 1=yes, 0=no
integer, parameter  :: pop_evolve = 1                !Use inherited values or base values? 1=evolve, 0=no evolve
integer, parameter  :: GearType = 4                  !1 = Gillnet, 2 = Trawl, 3 = FlatMort, 4 = Mixed fishing, 0 = No fishing
integer, parameter  :: dens_dep = 1                  !Include density dependance in FoodEnv. 1=yes, 0=no

integer, parameter  :: ClimScen = 0                  !Climate scenario. Base these on IPCC scenarios. 

!GENERAL PARAMETERS
integer, parameter  :: srp = SELECTED_REAL_KIND(8)   !Precision of reals
integer, parameter  :: Horizon = 2000                 !Amount of years we run the model for.
integer             :: last_run_time                 !The amount of years the mode ran previously, used when reading earlier Pop()

!FISH POPULATION PARAMETERS
integer, parameter  :: MaxAge       = 20             !Maximum age of any given individual
integer, parameter  :: TargetInd    = 150000         !Target nr of individuals, adjusts the population size
integer, parameter  :: MaxInd       = 800000         !Maximum number of individuals, gives the size of the population array.

!STOCHASTICITY PARAMETERS
real(srp), parameter    ::  phenotypic_deviance = 0.05 !Spread around the genotypic value for expressed traits (0.05 PLACEHOLDER)
real(srp), parameter    ::  inheritance_deviance = 0.14 !Spread around midparental values when inheriting traits (0.05 PLACEHOLDER)
real(srp), parameter    ::  FoodEnv_deviance = 0.03  !Spread around average food environment, represented by "1" (0.05 PLACEHOLDER)
real(srp), parameter    ::  temp_deviance = 0.05     !Spread around average temperature (0.05 PLACEHOLDER)
real(srp), parameter    ::  recruitment_deviance = 0.05 !Spread around number of recruits (0.05 PLACEHOLDER)

!MATURATION PARAMETERS
real(srp), parameter    ::  PMRN_width = 20          !Width of reaction norm            PLACEHOLDER
real(srp), parameter    ::  PMRN_slope = -2.5           !Slope of reaction norm            PLACEHOLDER

!ENERGETICS PARAMETERS
real(srp), parameter    ::  k            = 0.01      !Condition factor of somatic weight, W(g) = K * L(cm)^3
real(srp), parameter    ::  EnDensGon    = 6.93E6    !Energy density of gonads, Holdway and Beamish 1984. J kg^-1
real(srp), parameter    ::  EnDensSom    = 4.62E6    !Energy density of somatic tissue, Holdway and Beamish 1984. J kg^-1
real(srp), parameter    ::  GSIThreshold = 0.05      !Minimum GSI available after migration in order to NOT skip spawning
real(srp), parameter    ::  beta = 0.70              !Allometric scaling exponent from Quince 2006 and Boukal 2014
real(srp), parameter    ::  c_R = 0.5                !Conversion efficiency energy -> tissue. Holt and Jørgensen 2015
real(srp), parameter    ::  c_phi = 0.15             !Energetic cost of foraging coefficient. Holt and Jørgensen 2015
real(srp), parameter    ::  c_SMR = 4.67E6           !Metabolic rate coefficient. Unit is J kg^-1. Clarke and Johnston 1999.
real(srp), parameter    ::  c_SDA = 0.17             !Cost of digestion. From Holt and Jørgensen 2016, citing Hansson et al. 1996.

real(srp), parameter    ::  D_M = 780                !Distance to migrate, single way, kilometers. Jørgensen and Fiksen 2006.
real(srp), parameter    ::  c_COT = 4.18E1           !Cost of transport coefficient, J*km^-1, Ware 1978
real(srp), parameter    ::  c_u = 0.138              !Optimal swimming speed constant for pelagic fish, s^-1, Ware 1978
real(srp), parameter    ::  b_2 = 0.43               !Scaling factor for optimal cruising speed for a pelagic fish, Ware 1978
real(srp), parameter    ::  b_3 = 1.02               !Length scaling factor, Ware 1978
real(srp), parameter    ::  b_4 = 2.42               !Swimming speed scaling factor, Ware 1978

real(srp), parameter    ::  v1 = 4.11E6              !Maximal oxygen uptake parameter, adapted from Claireaux et al 2000
real(srp), parameter    ::  v2 = 0.015               !Maximal oxygen uptake parameter, adapted from Claireaux et al 2000
real(srp), parameter    ::  v3 = 1.062               !Maximal oxygen uptake parameter, adapted from Claireaux et al 2000
real(srp), parameter    ::  v4 = 7.13E6              !Maximal oxygen uptake parameter, adapted from Claireaux et al 2000

real(srp), parameter    ::  InvestDecay = 0.8        !Rate of decay of allocation to somatic growth, Quince et al. 2006, Eq. 7

real(srp), parameter    ::  cF1 = 2.4                !Foraging -> energy scaling factor, see OneNote 2022-2. cF1
real(srp), parameter    ::  cF2 = 0.29               !Foraging -> energy scaling factor, see OneNote 2022-2  cF2

!MORTALITY PARAMETERS
real(srp), parameter    ::  c_predation = 0.66       !Predation mortality coefficient, Holt and Jørgensen 2014
real(srp), parameter    ::  c_foraging = 0.03        !Foraging mortality coefficient, Holt and Jørgensen 2014
real(srp), parameter    ::  c_respiration = 11       !Respiration mortality coefficient, Holt and Jørgensen 2014
real(srp), parameter    ::  M_fixed = 0.07           !Size-independant mortality, eg. sickness, Holt and Jørgensen 2014
real(srp), parameter    ::  pred_exp = 0.75          !Predation exponent, McGurk 1986
real(srp), parameter    ::  for_risk = 1.8           !Foraging exponent, Holt and Jørgensen 2014
real(srp), parameter    ::  repro_exp = 2.5          !Reproduction exponent, Holt and Jørgensen 2014
real(srp), parameter    ::  resp_exp = 3.            !Respiration exponent, Holt and Jørgensen 2014

!RECRUITMENT PARAMETERS
real(srp), parameter    ::  small_tau = 1.4877E-6    !Bevertonholt parameter, typically alpha
real(srp), parameter    ::  small_eta = 7.046E-11    !Bevertonholt parameter, typically beta
real(srp), parameter    ::  large_tau = 1.877E-6     !Bevertonholt parameter, typically alpha, Enberg et al 2009 EQ 13
real(srp), parameter    ::  large_eta = 2.346E-11    !Bevertonholt parameter, typically beta, Enberg et al 2009 EQ 13
real(srp), parameter    ::  EggWeight = 4E-7         !The weight of a single egg, kg, Enberg et al. 2009, Kjesbu et al. 1998

!FISHING PARAMETERS
real(srp),  parameter   :: PropTrawl = 1.0           !Proportion that experiences trawling, rest experience gillnet (Only mixed)
real(srp),  parameter   :: Fmax = 0.2                !Maximum annual harvest rate
real(srp),  parameter   :: Lmax = 130.                !Size of maximum selectivity in Jørgensen 2009 selectivity equation
real(srp),  parameter   :: FlatMort = 0.2            !Used for unselective fisheries
real(srp),  parameter   :: FishSelWidth = 0.28*Lmax  !Sigma in Jørgensen 2009 selectivity equation
integer, parameter      :: HarvestStart = 0          !Year when harvesting starts
integer, parameter      :: HarvestDuration = Horizon-HarvestStart           !Duration of harvesting period
integer, parameter      :: HarvestStop = HarvestStart + HarvestDuration     !Year when harvesting stops


!----------------------------------------------------------------------------------------------------------------------------------
!DEFINE ARGUMENTS TO BE USED
!----------------------------------------------------------------------------------------------------------------------------------

!POPULATION VARIABLES AND TRAITS
real(srp), allocatable, dimension(:,:)   :: Pop, PopCopy, PopAgeSorted !Population array
real(srp), allocatable, dimension(:)     :: AnnualPopSummary           !The sum of all columns in Pop(), calculated yearly
real(srp), allocatable, dimension(:,:)   :: IntermediateResults, IMResultsAgeSorted !Storage for results not found in Pop()
real(srp), allocatable, dimension(:,:)   :: Heritability_holder        !Holds the heritabilities for saving purposes
integer     :: t    !Using for looping through time
integer     :: Ind  !Individual number. Used for looping through individuals
integer     :: NoAlive !Amount of living fish in the population
integer     :: NoRecr  !Number of recruits added during recruitment
integer     :: run_looper !Used for looping entire model for several runs


real(srp)   ::  Length, SomWeight, GonWeight    !Length, somatic weight and gonad weight
real(srp)   ::  CumGonWeight                    !Cumulative gonad weight, for whole life
real(srp)   ::  Pop_biomass                     !Total population biomass
real(srp)   ::  Maturity                        !Maturity. 0 is immature, any other number is age at maturation
real(srp)   ::  Age                             !Age
real(srp)   ::  Status                          !1=Alive, -1=Dead, -2=Fished

!MATURATION/PMRN
real(srp)   ::  Len50, LenMat                           !Length at 50% maturity, length at maturation
real(srp)   ::  PhenIntercept, PhenSlope, pmrnWidth     !PMRN intercept, slope and width
real(srp)   ::  ProbMat                                 !Probability of maturing
real(srp)   ::  pmrnDelta                               !Delta in ProbMat equation

!ENERGETICS
real(srp)   ::  B_SMR, N_y, N_y_max  !Standard metabolic rate, targeted energy availability and max energy availability
real(srp)   ::  V, V_max !Energy burned through metabolism, and maximum value of V
real(srp)   ::  FoodEnv !The amount of food in the environment.
real(srp)   ::  Foraging, Foraging_max, Foraging_min !Foraging behavior/Intensity and maximum values based on oxygen limitation.
real(srp)   ::  B_M !Energetic cost of migrating                                                  
real(srp)   ::  SomGro, GonGro !Somatic growth and gonad growth
real(srp)   ::  p_t, InitInvest !Allocation to somatic growth with accesories, Quince 2006a.
real(srp)   ::  Appetite !Desired energy intake, scales Foraging Intensity based on FoodEnv.
real(srp)   ::  Arrhenius !The Arrhenius temperature function, keeping it seperate from the SMR calculation
real(srp)   ::  SkipSpawnThreshold !Total energy needed after somatic growth in order to migrate and spawn (instead of skipping)
real(srp)   ::  SpawnSkip !Whether spawning was skipped or not. 1=yes, 0=no.
real(srp)   ::  Gross_N_in !Total energy in, before subtracting metabolism etc., compared to N_y
real(srp)   ::  B_SDA, B_phi !Cost of specific dynamic action and cost of foraging
real(srp)   ::  weight_watcher !Range ~0-1, relative energy intake decreaser with weight, see logbook 02-06-2023

!MORTALITY
real(srp)   ::  M_predation, M_foraging, M_reproduction, M_respiration  !Types of natural mortality contribution
real(srp)   ::  TrawlMort, GillnetMort !Instantaneous mortality rates in different fishing scenarios
real(srp)   ::  FishMort !Fisheries mortality, the one that's implemented
real(srp)   ::  Z, Surv !Total mortality and annual survival

!RECRUITMENT
real(srp)   ::  TEP !Total egg production, summarised every year for all individuals
real(srp)   ::  tau, eta !Beverton-Holt parameters, typically alpha and beta. Refer to parameters further up for values.
real(srp), allocatable, dimension(:,:) :: PropFec !Proportionate fecundity array, to calculate chance of being a parent
real(srp)   ::  ran_num1, ran_num2 !Random numbers used to select parents for recruits
integer     :: Parent1, Parent2 !The individual ID of the two weighted-randomly selected parents for each recruit.
integer     :: NoMature !Number of mature individuals during recruitment

!CLIMATE VARIABLES
real(srp)   ::  Temp  !Temperature


!FISHERIES VARIABLES
real(srp)   ::  FishSelCurve    !Fisheries selectivity curve, length based function dependant on gear type

!RESULTS VARIABLES
real(srp), allocatable, dimension(:,:,:) :: ncdf_PopResults !Storage array for results to be saved.
integer   :: ResultsAge !Used to write for every age in ncdf_PopResults()
character(26), dimension(:) :: results_list(NoRuns) !List of filenames that are saved, printed at the end of all runs.

!netCDF SAVING VARIABLES
integer     :: ncid, varid_val
integer     :: array_length, array_width, a_length, a_width, a_time, run_pars_dim, run_pars_val, run_heritability_val, ranseed_val
integer     :: ncdf_status, heritability_dim, heritability_size, ran_seed_dim
!ACCESSORY VARIABLES
character(*), parameter :: output_folder = 'output' !Folder where all output is stored
character(10) :: char_lastruntime, char_horizon !Character versions, used to save in meaningful filenames
character(26)  :: ncdf_datafile !Name of results-file untimately saved
character(8)  :: DATE !The date of model initiation, based on system datetime
character(10) :: TIME !The time of model initiation, based on system datetime
logical :: is_error_log !Logical indicator flag for errors
logical :: WriteResThisYear !Used to determine which years results are saved
integer :: ran_seed(33) !To store the random seed for the run
integer, parameter :: ran_seed_length = SIZE(ran_seed)
integer, dimension(MaxInd)  ::  index_holder !Temporary vector used when sorting an array
real(srp)   :: time_spent, time_remain !Processor time spent on the model in total
real(srp)   :: pct_progress !How far is the model in its current run
real(srp), dimension(:), allocatable :: parameter_holder !Array used to save all parameters
integer, parameter :: SaveSteps = (Horizon / WriteInterval) +1 !Amount of timesteps actually saved
integer     :: parameter_holder_length !Length of parameter_holder()
integer :: u, o !Used for saving the final population as binary
real(srp) :: count_popagesorted, count_intermediateresults !Used when saving results
integer :: loop_counter


CALL FS_MKDIR(output_folder) !Creates a folder named output. From m_csv_io.
CALL FS_CHDIR(output_folder) !Changes the directory to 'output' From m_csv_io.

IF (SmallPop .EQ. 1) THEN !If we are aiming for a smaller population for playing around
    tau = small_tau
    eta = small_eta
ELSEIF (SmallPop .EQ. 0) THEN !Standard run
    tau = large_tau
    eta = large_eta
ENDIF
!----------------------------------------------------------------------------------------------------------------------------------
!ALLOCATION OF NEEDED MEMORY
!----------------------------------------------------------------------------------------------------------------------------------

ALLOCATE(Pop(1:MaxInd,20),PopCopy(1:MaxInd,20))         !Whole population, and copy for temporarily holding it
ALLOCATE(AnnualPopSummary(1:20))                        !Vector for summin columns in Pop()
ALLOCATE(ncdf_PopResults(1:SaveSteps,1:MaxAge,48))      !Summarised results for all years divided into age groups.
ALLOCATE(PopAgeSorted(1:MaxInd,20))                     !Pop(), but with only 1 age group remaining
ALLOCATE(IntermediateResults(1:MaxInd,10))              !Similar to Pop(), but for keeping track of e.g. mortalities
ALLOCATE(IMResultsAgeSorted(1:MaxInd,10))               !Age sorted version of IntermediateResults()
ALLOCATE(parameter_holder(1:51))                        !For saving all parameter values
ALLOCATE(PropFec(1:Maxind,3))                           !Proportion of fecundity, for use in calculating inheritance
ALLOCATE(Heritability_holder(1:SaveSteps,3))            !Holds heritabilities on select traits to save them


!----------------------------------------------------------------------------------------------------------------------------------
!STORE PARAMETERS FOR CURRENT RUNS
!This is for saving the parameters at the end of the code
!----------------------------------------------------------------------------------------------------------------------------------
    parameter_holder( 1) = Horizon
    parameter_holder( 2) = MaxAge
    parameter_holder( 3) = TargetInd
    parameter_holder( 4) = MaxInd
    parameter_holder( 5) = k 
    parameter_holder( 6) = EnDensGon
    parameter_holder( 7) = EnDensSom
    parameter_holder( 8) = HarvestDuration
    parameter_holder( 9) = HarvestStart
    parameter_holder(10) = beta
    parameter_holder(11) = c_R
    parameter_holder(12) = c_phi
    parameter_holder(13) = c_SMR
    parameter_holder(14) = c_SDA
    parameter_holder(15) = D_M
    parameter_holder(16) = c_COT
    parameter_holder(17) = c_u
    parameter_holder(18) = b_2
    parameter_holder(19) = b_3
    parameter_holder(20) = b_4
    parameter_holder(21) = v1
    parameter_holder(22) = v2
    parameter_holder(23) = v3
    parameter_holder(24) = v4
    parameter_holder(25) = InvestDecay
    parameter_holder(26) = c_predation
    parameter_holder(27) = c_foraging
    parameter_holder(28) = c_respiration
    parameter_holder(29) = M_fixed
    parameter_holder(30) = pred_exp
    parameter_holder(31) = for_risk
    parameter_holder(32) = repro_exp
    parameter_holder(33) = resp_exp
    parameter_holder(34) = tau
    parameter_holder(35) = eta
    parameter_holder(36) = EggWeight
    parameter_holder(37) = ClimScen
    parameter_holder(38) = GearType
    parameter_holder(39) = PropTrawl
    parameter_holder(40) = Fmax
    parameter_holder(41) = Lmax
    parameter_holder(42) = FlatMort
    parameter_holder(43) = FishSelWidth
    parameter_holder(44) = dens_dep
    parameter_holder(45) = cF1
    parameter_holder(46) = cF2
    parameter_holder(47) = phenotypic_deviance
    parameter_holder(48) = inheritance_deviance
    parameter_holder(49) = FoodEnv_deviance
    parameter_holder(50) = temp_deviance
    parameter_holder(51) = recruitment_deviance
    parameter_holder_length = SIZE(parameter_holder)

!----------------------------------------------------------------------------------------------------------------------------------
!INITIATE THE RUN
!----------------------------------------------------------------------------------------------------------------------------------
DO run_looper=1, NoRuns
PRINT '(5x,A,i2,A2,i3)', "Now starting run number:",run_looper,"/",NoRuns
!----------------------------------------------------------------------------------------------------------------------------------
!CREATE NETCDF FILE FOR STORING STUFF
!Documentation at: https://www.unidata.ucar.edu/software/netcdf/docs/netcdf_documentation.html
!----------------------------------------------------------------------------------------------------------------------------------
    CALL date_and_time(DATE=DATE, TIME=TIME)
    CALL RANDOM_SEED(get=ran_seed)
    WRITE(ncdf_datafile, fmt='(A,A,A,A,A)') DATE,"_", TIME(1:6),"_","popdata.nc"
    array_length        = size(ncdf_PopResults, DIM=2)
    array_width         = size(ncdf_PopResults, DIM=3)
    heritability_size   = size(Heritability_holder, DIM=2)
    ncid = 1
    ncdf_status = nf90_create(path = ncdf_datafile, cmode = nf90_NETCDF4, ncid = ncid)

    ncdf_status = nf90_def_dim(ncid, "Age_groups", array_length, a_length)
    ncdf_status = nf90_def_dim(ncid, "Result_cols", array_width, a_width)
    ncdf_status = nf90_def_dim(ncid, "Years", SaveSteps, a_time)
    ncdf_status = nf90_def_dim(ncid, "Heritability_dim", heritability_size, heritability_dim)
    ncdf_status = nf90_def_dim(ncid, "Parameter_dim", parameter_holder_length, run_pars_dim)
    ncdf_status = nf90_def_dim(ncid, "Ran_seed_dim", ran_seed_length, ran_seed_dim)

    ncdf_status = nf90_def_var(ncid, "PopResults", nf90_float, [a_time, a_length, a_width], varid_val)
    ncdf_status = nf90_def_var(ncid, "Parameters", nf90_float, [run_pars_dim], run_pars_val)
    ncdf_status = nf90_def_var(ncid, "Heritabilities", nf90_float, [a_time, heritability_dim], run_heritability_val)
    ncdf_status = nf90_def_var(ncid, "RanSeed", nf90_int, [ran_seed_dim], ranseed_val)

    ncdf_status = nf90_put_att(ncid, nf90_global, "Note", "Contains the data from the runs of my IBM")
    ncdf_status = nf90_put_att(ncid, varid_val, "Note", "Contains results averaged for every age group, every year")
    ncdf_status = nf90_put_att(ncid, run_pars_val, "Note", "Contains parameters for the run")
    ncdf_status = nf90_put_att(ncid, run_heritability_val, "Note", "Contains heritabilities of select traits for year")
    ncdf_status = nf90_put_att(ncid, run_heritability_val, "Legend", "Age-at-mat")
    ncdf_status = nf90_put_att(ncid, ranseed_val, "Note", "Contains the random seed string for this run")
    ncdf_status = nf90_put_att(ncid, run_pars_val, "Legend", "Horizon, MaxAge, TargetInd, MaxInd, k, EnDensGon, EnDensSom, &
    & HarvestDuration, HarvestStart, beta, c_R, c_phi, c_SMR, c_SDA, D_M, c_COT, c_u, b_2, b_3, b_4, v1, v2, v3, v4, InvestDecay, &
    & c_predation, c_foraging- c_respiration, M_fixed, pred_exp, for_risk, repro_exp, resp_exp, tau, eta, EggWeight, ClimScen &
    & GearType, PropTrawl, Fmax, Lmax, FlatMort, FishSelWidth, density_dependance, cF1, cF2, phenotypic_deviance, &
    & inheritance_deviance, FoodEnv_deviance, temp_deviance, recruitment_deviance")
    ncdf_status = nf90_put_att(ncid, nf90_global, "Traits_in_order", "year,age,temp,FoodEnv,n_recr,n_survived,n_died,SomWeight,& 
    &SomWeight_SD, GonWeight,GonWeight_SD,length,length_SD,M_predation,Mpred_SD,M_foraging, Mfor_SD,M_reproduction,Mrepr_SD,&
    &M_respiration,Mresp_SD,FishMort,FishMort_SD, f_int,fint_SD,n_mature,n_spawnskip,pop_biomass,som_growth, somgro_SD,&
    &gon_growth, gongro_SD, appetite, appetite_SD, somatic allocation, allocation_SD, PMRN intercept, intercept_SD, age-at-mat, &
    &age-at-mat_SD, length-at-mat, length-at-mat_SD, cumulative-gonad-weight, CumGonWeight_SD")

    ncdf_status = nf90_put_var(ncid, run_pars_val, parameter_holder)
    ncdf_status = nf90_put_var(ncid, ranseed_val, ran_seed)

    ncdf_status = nf90_close(ncid = ncid)

!----------------------------------------------------------------------------------------------------------------------------------
!SET ARGUMENTS AND ARRAYS TO 0
!----------------------------------------------------------------------------------------------------------------------------------
Pop = 0.
PopCopy = 0.
PopAgeSorted = 0.
IntermediateResults = 0.
IMResultsAgeSorted = 0.
ncdf_PopResults = 0.
NoAlive = 0
NoRecr = 0
FoodEnv = 0.
TEP = 0.
Temp = 0.
last_run_time = 0

!----------------------------------------------------------------------------------------------------------------------------------
!INITIALIZE POPULATION
!----------------------------------------------------------------------------------------------------------------------------------
 !The population takes the form of an array, with each line corresponding to an individual fish, and each column corresponding
 !to a trait of interest, such as age, length and maturity status. 
 !While currently not functional, there will be an option to read in an initial population rather than initializing from scratch
 !in every run.
IF (ReadInitialPopulation .EQ. 1) THEN
    INQUIRE(FILE='finalpop.dat', exist=is_error_log ) ! Check file exists
    IF ( is_error_log .EQV. .FALSE. ) then ! Note: inquire exist is TRUE when no error, i.e. the reverse of is error
        PRINT *, "ERROR: data file for initial population cannot be found"
        STOP
    ENDIF
        OPEN(16, file='last_run_time.dat', form='unformatted',status='old') !When loading old population, also keep track of how
            read (16) last_run_time                                         !many years the population has been run
        CLOSE(16)
        PRINT '(5x,A,i4,x,A)', "Last run finished after", last_run_time, "years"   
        PRINT '(5x,A,i4)', "This program will now run an additional:", Horizon    
        PRINT '(5x,A)', "Reading population from finalpop.dat"

        OPEN(15, file='finalpop.dat', form='unformatted',status='old')  !Read binary data, store in Pop()
            read (15) Pop
        CLOSE(15)

ELSEIF (ReadInitialPopulation .EQ. 2) THEN
    INQUIRE(FILE='startpop.dat', exist=is_error_log ) ! Check file exists
    IF ( is_error_log .EQV. .FALSE. ) then ! Note: inquire exist is TRUE when no error, i.e. the reverse of is error
        PRINT *, "ERROR: data file for initial population cannot be found"
        STOP
    ENDIF
    OPEN(15, file='startpop.dat', form='unformatted',status='old')  !Read binary data, store in Pop()
            read (15) Pop
    CLOSE(15)
ELSE

    DO Ind = 1, TargetInd !Pop() initializing loop
        Pop(Ind, 0) = dble(ind)  !Invidual number, needs to be REAL
        Pop(Ind, 1) = 1.    !Status; alive=1, dead=-1
        Pop(Ind, 2) = 3.    !Age
        Pop(Ind, 3) = 20.   !Length
        Pop(Ind, 4) = 0.25  !Somatic Weight
        Pop(Ind, 5) = 0.    !Stored energy (may be energy or weight, revisit) Liver weight?
        Pop(Ind, 6) = 0.    !Maturity, 0=immature, number=age at maturation
        Pop(Ind, 7) = 0.    !Gonad weight
        Pop(Ind, 8) = 0.    !CumGonadweight
        Pop(Ind, 9) = 112.   !PMRN Intercept                                     - Evolving
        Pop(Ind,10) = PMRN_slope   !PMRN Slope
        Pop(Ind,11) = PMRN_width   !40. + n()*10. !PMRN Width
        Pop(Ind,12) = 0.68   !Gen. Initial somatic investment                    - Evolving
        Pop(Ind,13) = 1.60E7    !Gen. appetite, from Quince                      - Evolving
        Pop(Ind,14) = 0.    !Midparental age at maturation for calculating the heritability
        Pop(Ind,15) = Pop(ind, 9) * (random_normal()*phenotypic_deviance+1) !Phenotypic Intercept
        Pop(Ind,16) = Pop(ind,10) * (random_normal()*phenotypic_deviance+1) !Phenotypic Slope
        Pop(Ind,17) = Pop(ind,12) * (random_normal()*phenotypic_deviance+1) !Phenotypic initial investment
        Pop(Ind,18) = Pop(ind,13) * (random_normal()*phenotypic_deviance+1) !Phenotypic appetite
        Pop(Ind,19) = 0.    !Length at maturation
        Pop(Ind,20) = 0.    !Midparental Length at maturation
    ENDDO !Pop() initializing loop
ENDIF


!----------------------------------------------------------------------------------------------------------------------------------
!TIME & INDIVIDUAL LOOPS
!----------------------------------------------------------------------------------------------------------------------------------
 !The model runs over 'Horizon' years, by performing a loop over all individuals for every timestep (year)
 !The individual loop grabs the traits from the Pop() array, gives them meaningful names, performes operations using them,
 !such as calculating maturation probability, mortality rates etc., and then stores them back in the array. 
 !This is performed at every time step, after which upkeep is performed: recruitment occurs, results are saved, dead individuals
 !are removed etc., and the model then progresses to the next year and does this again.

!Time loop, annual time steps from year=1 to Horizon
DO t = 1, Horizon


    Pop_biomass = SUM(Pop(:,4),DIM=1) !Total population biomass

    FoodEnv = 1. * (random_normal()*FoodEnv_deviance+1)
    IF (dens_dep .EQ. 1) FoodEnv = FoodEnv + 0.15 - 0.00012 * (Pop_biomass/1000)
    IF (dens_dep .EQ. 2) FoodEnv = FoodEnv + 0.15 - 0.00040 * (Pop_biomass/1000)
    IF (dens_dep .EQ. 3) FoodEnv = FoodEnv + 0.15 - 0.00045 * (Pop_biomass/1000)
 

    Temp    = 4. * (random_normal()*temp_deviance+1)   !Temperature in environment, degrees celcius.
    IF (ClimScen .EQ. 1) Temp = Temp + (0.25 * t) / (10 + 0.3 * t)
    !IF (ClimScen .EQ. 2) Temp = Temp + (0.18 * t) / (8 + 0.023 * t) !Higher scenario
    IF (ClimScen .EQ. 2) Temp = Temp + (0.16 * t) / (5.6 + 0.05 * t) !Lower scenario
    IF (ClimScen .EQ. 3) Temp = Temp + (9 * exp( -2.7 * exp( -0.011 * t ) )) - 0.6
    IF (ClimScen .EQ. 4) Temp = Temp + (11.5 * exp( -2.68 * exp( -0.011 * t ) )) - 0.7

    NoAlive = COUNT(Pop(:,1)/=0.,DIM=1)  !counts the number of individuals alive in the population

    IF (mod(t,5)==0) THEN !Print update to terminal to see how things are going
        pct_progress = (real(t)/real(Horizon))*100.
        CALL cpu_time(time_spent)
        time_remain = ((time_spent / pct_progress) * (100.-pct_progress))/60.
        PRINT '(a8,a13,5a11,2a14,a18)', "Run", "Year", "Temp", "FoodEnv", "#Alive", "#Mature", "#Recruits","Biomass (t)", &
        & "%Completion", "Mins. remaining"
        PRINT '(i3,a2,i3,i6,a2,i5,2f11.3,3i11,2f14.2,f18.2)',run_looper,"/",NoRuns, last_run_time,'+', t, Temp, FoodEnv, NoAlive,&
        & NoMature, NoRecr, Pop_biomass/1000, pct_progress, time_remain
        PRINT *, "===============================================================================================================&
        &============="
    ENDIF


    IF (NoAlive .EQ. 0) THEN !Stop the run if there are no individuals left alive
        PRINT '(5x,A,i4)', "No individuals left alive, ending run now, made it to year ", t
        STOP 
    ENDIF

    !Individual loop, each individual treated seperately
    DO Ind=1, NoAlive
        !Set arguments to 0 to avoid errors.
        Length = 0. 
        SomWeight = 0.
        GonWeight = 0.
        Maturity = 0.
        Age = 0.
        Status = 0.
        FishMort = 0.
        M_predation = 0.
        M_foraging = 0.
        M_reproduction = 0. 
        M_respiration = 0.
        Len50 = 0.
        LenMat = 0.                      
        PhenIntercept = 0. 
        PhenSlope = 0.
        pmrnWidth = 0.
        ProbMat = 0.
        pmrnDelta = 0.
        B_SMR = 0.
        B_SDA = 0.
        B_phi = 0. 
        N_y = 0.
        V = 0.
        V_max = 0.
        Foraging = 0. 
        Foraging_max = 0.
        B_M = 0.                                                
        SomGro = 0. 
        GonGro = 0.
        p_t = 0.
        InitInvest = 0.
        Appetite = 0.
        Arrhenius = 0.
        Z = 0. 
        Surv = 0.
        FishSelCurve = 0.
        TrawlMort = 0.
        GillnetMort = 0.
        SpawnSkip = 0.
        Gross_N_in = 0.
        CumGonWeight = 0.

        !Store traits from Pop() in more meaningful names here
        Status          = Pop(Ind, 1)
        Age             = Pop(Ind, 2)
        Length          = Pop(Ind, 3)
        SomWeight       = Pop(Ind, 4)
        Maturity        = Pop(Ind, 6)
        GonWeight       = Pop(Ind, 7)
        CumGonWeight    = Pop(Ind, 8)
        pmrnWidth       = Pop(Ind, 11)
        PhenIntercept   = Pop(Ind, 15)
        PhenSlope       = Pop(Ind, 16)
        InitInvest      = Pop(Ind, 17)
        LenMat          = Pop(Ind, 19)
        Appetite        = Pop(Ind, 18)
        !Etc, etc
        
        
        !**************************************************************************************************************************
        !MATURATION
        !Maturation is based on probabilistic maturation reaction norms (PMRN), where each individual has a chance of maturing
        !based on their age and size. Maturation age is saved in 'Maturity' for matured individuals.
        !https://www.int-res.com/abstracts/meps/v335/p253-269/ <- More info on PMRN's
        !**************************************************************************************************************************
        IF (Maturity .eq. 0.) THEN
            Len50 = PhenIntercept + Age-1 * PMRN_slope
            pmrnDelta = PMRN_width/(log(0.75_SRP/0.25_SRP) - log(0.25_SRP/0.75_SRP))

            ProbMat = 1. / ( 1. + exp(-(Length - Len50) / pmrnDelta) )
            
            IF (zrand() .lt. ProbMat) THEN
                Maturity = Age      !Maturation age is saved
                LenMat   = Length      !Maturation length    
            ENDIF
    
        ENDIF !Maturation loop

        !**************************************************************************************************************************
        !ENERGETICS
        !Energetics in the model combines principles from different previous models. Foragin behavior is scaled by 'Appetite', a
        !evolving trait set for each individual (Quince 2006) setting energy availability. Energy is divided between processes
        !based on the Wisconsin Bioenergetics Framework as parametised by Holt and Jørgensen (2014), limited by an oxygen budget
        !adapted from Claireaux (2000).
        !**************************************************************************************************************************
        Arrhenius = (exp( 15.7 - ( 5020 / (Temp + 273.15) )) * 434 * 24 * 365) / (c_SMR * 0.05**beta) !Temperature function
        !Temperature function was originally standardised for a 50g fish by Clarke and Johnston, hence the division.
        !Further notes on these values are found in OneNote 2020 week 43-45
        !Standard metabolic rate
        B_SMR   = c_SMR * Arrhenius * (SomWeight+GonWeight)**beta !Holt and Jørgensen 2014 Eq. 5

        Foraging = (Appetite * SomWeight**beta) / &
        & ( FoodEnv * B_SMR * cF1 - SomWeight**beta * Appetite * cF2 )

        ! Legacy code, failed attempt at regulating larger body sizes
        ! Foraging_min = (-cF1 * c_SDA * FoodEnv * weight_watcher + cF1 * FoodEnv * weight_watcher - cF2 - c_phi - &
        ! & SQRT(cF1**2 * c_SDA**2 * FoodEnv**2 * weight_watcher**2 - 2 * cF1**2 * c_SDA * FoodEnv**2 * weight_watcher**2 + &
        ! & cF1**2 * FoodEnv**2 * weight_watcher**2 + 2 * cF1 * cF2 * c_SDA * FoodEnv * weight_watcher + &
        ! & 2* cF1 * c_SDA * FoodEnv * c_phi * weight_watcher - 2 * cF1 * cF2 * FoodEnv * weight_watcher - &
        ! & 2* cF1 * FoodEnv * c_phi * weight_watcher + cF2 - 2* cF2 * c_phi + c_phi**2)) &
        ! & / (2 * cF2 * c_phi)
        ! Foraging_max = (-cF1 * c_SDA * FoodEnv * weight_watcher + cF1 * FoodEnv * weight_watcher - cF2 - c_phi + &
        ! & SQRT(cF1**2 * c_SDA**2 * FoodEnv**2 * weight_watcher**2 - 2 * cF1**2 * c_SDA * FoodEnv**2 * weight_watcher**2 + &
        ! & cF1**2 * FoodEnv**2 * weight_watcher**2 + 2 * cF1 * cF2 * c_SDA * FoodEnv * weight_watcher + &
        ! & 2* cF1 * c_SDA * FoodEnv * c_phi * weight_watcher - 2 * cF1 * cF2 * FoodEnv * weight_watcher - &
        ! & 2* cF1 * FoodEnv * c_phi * weight_watcher + cF2 - 2* cF2 * c_phi + c_phi**2)) &
        ! & / (2 * cF2 * c_phi)

        
        !Minimum foraging - Below this energy cannot cover basic metabolism
        Foraging_min = -(c_SDA *FoodEnv* cF1 -FoodEnv *cF1 + SQRT(cF1**2 * c_SDA**2 * FoodEnv**2- 2* cF1**2 * c_SDA * FoodEnv**2& 
        & + cF1**2 * FoodEnv**2 +2*cF1*cF2*c_SDA*FoodEnv +2*cF1*c_SDA*FoodEnv*c_phi-2*cF1*cF2*FoodEnv-2*cF1*FoodEnv*c_phi &
        & + cF2**2 - 2*cF2*c_phi + c_phi**2)+cF2 + c_phi) / (2*cF2*c_phi)
        !Maximum foraging - Above this energy no longer covers the cost of foraging
        Foraging_max = (-c_SDA *FoodEnv* cF1 +FoodEnv *cF1 + SQRT(cF1**2 * c_SDA**2 * FoodEnv**2- 2* cF1**2 * c_SDA * FoodEnv**2& 
        & + cF1**2 * FoodEnv**2 +2*cF1*cF2*c_SDA*FoodEnv +2*cF1*c_SDA*FoodEnv*c_phi-2*cF1*cF2*FoodEnv-2*cF1*FoodEnv*c_phi &
        & + cF2**2 - 2*cF2*c_phi + c_phi**2)-cF2 - c_phi) / (2*cF2*c_phi)

        IF (Foraging .GE. Foraging_max) Foraging = Foraging_max - 0.05
        IF (Foraging .LE. Foraging_min) Foraging = Foraging_min + 0.05

        !weight_watcher = ww1 + (ww2/ (1 + ww3**(ww4*(SomWeight-ww5)) ) ) !See logbook 02-06-2023
        Gross_N_in = FoodEnv * B_SMR * ((cF1*Foraging)/(1+cF2*Foraging)) 
        B_SDA = c_SDA * Gross_N_in
        B_phi = c_phi * Foraging * B_SMR
        N_y = (Gross_N_in - B_SDA - B_SMR - B_phi) * c_R

        !Respiration, total oxygen consumption and maximum oxygen consumption
        V = c_SDA * Gross_N_in + B_SMR + c_phi * Foraging * B_SMR + (1-c_R) * &
        & (Gross_N_in - c_SDA * Gross_N_in - B_SMR - c_phi * Foraging * B_SMR)
        V_max = ( v1 * Temp**(-v2 * Temp + v3) + v4 ) * SomWeight**beta

        IF (Maturity .EQ. 0.) THEN !Immature individuals only dedicate energy to somatic growth
            SomGro = (N_y / EnDensSom)                                              !Somatic growth
            GonGro = 0.                                                             !Gonad growth (=0)
         
         ELSE IF (Maturity .GT. 0.1) THEN !Mature individuals may dedicate energy to gonads as well

            p_t = InitInvest * InvestDecay **(Age - Maturity)                       !Calculate somatic investment
            B_M = (c_COT * Length**b_3 *(c_u * Length**b_2)) *2 * D_M               !Cost of migration

            SkipSpawnThreshold = B_M + ((GSIThreshold*SomWeight*EnDensGon)/(1-GSIThreshold))!Energy threshold, based on required 
            ! GSI post-spawning migrations, GSIThreshold, set earlier. 

            IF ((N_y * (1-p_t)) .GT. SkipSpawnThreshold) THEN 
                SpawnSkip = 0.
                SomGro = ((N_y * p_t) / EnDensSom)                                  !Somatic growth
                GonGro = ((N_y * (1-p_t) - B_M ) / EnDensGon)                       !Gonad growth
            ELSE
                SpawnSkip = 1.
                SomGro = (N_y / EnDensSom)                                          !Somatic growth
                GonGro = 0.                                                         !Gonad growth (=0)
            ENDIF
            
        ENDIF !Growth loop

        !Actual growth
        SomWeight = SomWeight + SomGro                                          !Somatic weight in kg
        GonWeight = GonGro                                                      !Gonad weight in kg
        Length = ( ( (SomWeight + GonWeight) * 1000 ) / k ) ** ( 1. / 3. )      !Weight needed to be in grams

        CumGonWeight = CumGonWeight + GonWeight

        !**************************************************************************************************************************
        !MORTALITY
        !Mortality calculations taken from Holt and Jørgensen 2014. Fisheries mortality can come from either trawl fishing or
        !gillnet fishing. How to assign these has not yet been decided.
        !**************************************************************************************************************************
        !Natural mortality
        M_predation      = c_predation * Length **(-pred_exp)                               !Size-dependent predation mortality
        M_foraging       = c_foraging * Foraging**(for_risk) * M_predation                   !Foraging mortality
        IF (Maturity .EQ. 0.) THEN
            M_reproduction  = 0.                                                            !Reproductive mortality
        ELSEIF (Maturity .GT. 0.) THEN
            M_reproduction  = ( ( GonWeight / SomWeight ) / 0.10 )**repro_exp * M_predation
        ENDIF
        M_respiration    = c_respiration * (V / V_max)**(resp_exp) * M_predation            !Respiration mortality

        !Fisheries mortality, calculate for both methods
        IF (Length .LT. Lmax) THEN
            TrawlMort = (exp((-( Length - Lmax )**2) / (2 * FishSelWidth**2))) * Fmax 
         ELSE
            TrawlMort = 1. * Fmax
        ENDIF
        GillnetMort = (exp((-( Length - Lmax )**2) / (2 * FishSelWidth**2))) * Fmax

        !Which fisheries mortality gets implemented
        IF (( t .GE. HarvestStart) .AND. (t .LE. HarvestStop)) THEN
            IF (GearType .EQ. 1) THEN !Gillnet
                    FishMort = GillnetMort   
                ELSEIF (GearType .EQ. 2) THEN !Trawling
                    FishMort = TrawlMort
                ELSEIF (GearType .EQ. 3) THEN !Flat mortality
                    FishMort = FlatMort
                ELSEIF (GearType .EQ. 4) THEN !Mixed fishing
                    IF (zrand() .LT. PropTrawl) THEN
                        FishMort = TrawlMort
                    ELSE
                        FishMort = GillnetMort
                    ENDIF
                ELSEIF (GearType .EQ. 0) THEN
            FishMort = 0.  
            ENDIF
        ENDIF

        !Total mortality is the sum of component mortalities
        Z = M_predation + M_foraging + M_reproduction + M_respiration + M_fixed + FishMort
        Surv = exp(-Z)

        !Flag individuals for death, which will happen after reproduction, during "upkeep"
        IF (Age .GT. 0.0) THEN
            IF (zrand() .GT. Surv) THEN                         !If chance of survival is low, probably die
                Status = -1.                                    !Flag for death
                IF (zrand() .LT. (FishMort/Z)) Status = -2.
            ENDIF
         ENDIF
        

        !Put arguments back into Pop() here
        Pop(Ind, 1)  = Status
        Pop(Ind, 2)  = Age
        Pop(Ind, 3)  = Length
        Pop(Ind, 4)  = SomWeight
        Pop(Ind, 6)  = Maturity
        Pop(Ind, 7)  = GonWeight
        Pop(Ind, 8)  = CumGonWeight
        Pop(Ind, 15) = PhenIntercept 
        Pop(Ind, 16) = PhenSlope
        Pop(Ind, 17) = InitInvest
        Pop(Ind, 19) = LenMat


        !Save intermediate values into IntermediateResults()
        IntermediateResults(Ind, 0)  = dble(Ind)
        IntermediateResults(Ind, 1)  = Status
        IntermediateResults(Ind, 2)  = M_predation
        IntermediateResults(Ind, 3)  = M_foraging
        IntermediateResults(Ind, 4)  = M_reproduction
        IntermediateResults(Ind, 5)  = M_respiration
        IntermediateResults(Ind, 6)  = FishMort
        IntermediateResults(Ind, 7)  = Foraging
        IntermediateResults(Ind, 8)  = SpawnSkip
        IntermediateResults(Ind, 9)  = SomGro
        IntermediateResults(Ind, 10) = GonGro

    ENDDO !Individual loop


    NoAlive = COUNT(Pop(:,1)/=0.,DIM=1)  !Counts the number of individuals alive in the population
    AnnualPopSummary = SUM(Pop,DIM=1)
    !******************************************************************************************************************************
    !RECRUITMENT
    !Recruitment is based on Beverton-Holt, using the total fecundity to calculate Total Egg Production (TEP). The 'tau' and 'eta'
    !variables are typicalled called 'alpha' and 'beta' in B-H functions. In this model, they are scaled based on TargetInd, to
    !alter recruitment so the population-size equilibrium can be adjusted to allow the model to run faster when tinkering.
    !******************************************************************************************************************************
    ! DO Ind = 1, MaxInd                  !Remove dead individuals (set all columns to 0.)
    !     IF (Pop(Ind,1) .EQ. -1.) THEN   !If individuals are dead
    !         Pop(Ind, :) = 0.
    !     ENDIF
    ! ENDDO                               !Remove dead individuals
    
    NoMature = COUNT(Pop(:,6).NE.0., DIM=1)
    PropFec = 0.
    DO Ind = 1, MaxInd !Calculate how much each individual in Pop() contribute to the total fecundity of the population
        PropFec(Ind, 1) = dble(ind)
        PropFec(Ind, 2) = (Pop(Ind,7))/AnnualPopSummary(7)
        !IF(Pop(Ind,1) .LT. 0.) PropFec(Ind,:) = 0.
    ENDDO

    CALL ARRAY_INDEX( PropFec( : , 2), index_holder )   
    index_holder(1:MaxInd) = index_holder(MaxInd:1:-1)
    PropFec(:,:) = PropFec( index_holder , :)           !Sort PropFec according to relative contribution

    DO Ind = 1, MaxInd !Calculate cumulative contribution to fecundity
        IF (Ind .EQ. 1) THEN
            PropFec(Ind, 3) = PropFec(Ind, 2)
        ELSE
            PropFec(Ind, 3) = PropFec(Ind-1, 3) + PropFec(Ind, 2)
        ENDIF
    ENDDO
    TEP = ( AnnualPopSummary(7) / (EggWeight) ) !Total egg production, based on total gonad weight in the population
    NoRecr = nint(( ( tau * TEP ) / ( 1 + eta * TEP ) ) * (zrand()*recruitment_deviance+1) ) !Number of recruits to generate, B-H 

    IF (NoRecr .LT. 1) THEN
        PRINT *, "It happened again"
    ENDIF


    IF (NoAlive + NoRecr .LT. MaxInd) THEN !If population tabe is full, no recruitment
        IF (pop_evolve .EQ. 1) THEN !Evolving recruitment
        DO Ind = NoAlive+1, NoAlive+NoRecr !ADD RECRUITS TO THE ARRAY, this is done for every recruit
            !Select two parents
            loop_counter = 0
            100 CONTINUE
            ran_num1 = zrand()
            ran_num2 = zrand()
            Parent1 = random_pull(PropFec, 3, ran_num1)
            Parent2 = random_pull(PropFec, 3, ran_num2)
            IF ( Parent1 .EQ. Parent2 ) THEN
                loop_counter = loop_counter + 1
                IF (loop_counter .GE. 20) THEN
                    NoRecr = 0
                    PRINT *, "Not enough mature fish to select parents, skipping recruitment in year", t
                    GO TO 200
                ENDIF
                GO TO 100
            ENDIF

            Pop(Ind, 0) = dble(ind)  !Invidual number, needs to be REAL
            Pop(Ind, 1) = 1.    !Status; alive=1, dead=-1
            Pop(Ind, 2) = 0.    !Age
            Pop(Ind, 3) = 12.   !Length
            Pop(Ind, 4) = 0.01   !Somatic weight
            Pop(Ind, 5) = 0.    !Stored energy (may be energy or weight, revisit) Liver weight?
            Pop(Ind, 6) = 0.    !Maturity, 0=immature, number=age at maturation
            Pop(Ind, 7) = 0.    !Gonad Weight
            Pop(Ind, 8) = 0.    !CumGonadweight
            Pop(Ind, 9) = ((Pop(Parent1, 9)+Pop(Parent2, 9))/2)*(random_normal()*inheritance_deviance+1) !PMRN Intercept - Evolving 
            Pop(Ind,10) = PMRN_slope    !PMRN Slope
            Pop(Ind,11) = PMRN_width    !PMRN Width
            Pop(Ind,12) = 0.50   !Gen. Initial investment
            ! IF (Pop(Ind,12) .LE. 0.0) Pop(Ind,12) = 0.01
            ! IF (Pop(Ind,12) .GT. 1.0) Pop(Ind,12) = 1.0
            Pop(Ind,13) = ((Pop(Parent1,13)+Pop(Parent2,13))/2)*(random_normal()*inheritance_deviance+1)!Gen. appetite, from Quince
            Pop(Ind,14) = ((Pop(Parent1,6)+Pop(Parent2,6))/2)    !Midparental age at maturation
            Pop(Ind,15) = Pop(ind, 9) * (random_normal()*phenotypic_deviance+1) !Phenotypic Intercept
            Pop(Ind,16) = Pop(ind,10) * (random_normal()*phenotypic_deviance+1) !Phenotypic Slope
            Pop(Ind,17) = Pop(ind,12) * (random_normal()*phenotypic_deviance+1) !Phenotypic initial somatic investment
            IF (Pop(Ind,17) .LE. 0.0) Pop(Ind,12) = 0.01
            IF (Pop(Ind,17) .GT. 1.0) Pop(Ind,12) = 1.0
            Pop(Ind,18) = Pop(ind,13) * (random_normal()*phenotypic_deviance+1) !Phenotypic appetite
            Pop(Ind,19) = 0.           !Length at maturation
            Pop(Ind,20) = ((Pop(Parent1,19)+Pop(Parent2,19))/2)           !Midparental Length at maturation
        ENDDO !Recruitment loop
        ELSE  !Non-evolving recruitment
            DO Ind = NoAlive+1, NoAlive+NoRecr !ADD RECRUITS TO THE ARRAY, non-evolving
                Pop(Ind, 0) = dble(ind)  !Invidual number, needs to be REAL
                Pop(Ind, 1) = 1.    !Status; alive=1, dead=-1
                Pop(Ind, 2) = 0.    !Age
                Pop(Ind, 3) = 10.   !Length
                Pop(Ind, 4) = 0.1   !Somatic weight
                Pop(Ind, 5) = 0.    !Stored energy (may be energy or weight, revisit) Liver weight?
                Pop(Ind, 6) = 0.    !Maturity, 0=immature, number=age at maturation
                Pop(Ind, 7) = 0.    !Gonad Weight
                Pop(Ind, 8) = 0.    !CumGonadweight
                Pop(Ind, 9) = 80.   !PMRN Intercept
                Pop(Ind,10) = PMRN_slope    !PMRN Slope
                Pop(Ind,11) = PMRN_width    !PMRN Width
                Pop(Ind,12) = 0.5   !Gen. Initial somatic investment
                Pop(Ind,13) = 1.8E7    !Gen. appetite, from Quince
                Pop(Ind,14) = 0.           !Midparental age at maturation for calculating the heritability
                Pop(Ind,15) = Pop(ind, 9)  !Phenotypic Intercept
                Pop(Ind,16) = Pop(ind,10)  !Phenotypic Slope
                Pop(Ind,17) = Pop(ind,12)  !Phenotypic allocation propensity
                Pop(Ind,18) = Pop(ind,13)  !Phenotypic appetite
                Pop(Ind,19) = 0.           !Length at maturation
                Pop(Ind,20) = 0.           !Midparental Length at maturation
            ENDDO
        ENDIF
    ELSE
        PRINT *, "Population array full, no recruitment in year", t
    ENDIF
    200 CONTINUE
    !******************************************************************************************************************************
    !UPKEEP
    !Results are written, individuals are aged, and dead individuals are removed from the array, after which it gets resorted.
    !******************************************************************************************************************************
    NoAlive = COUNT(Pop(:,1)/=0.,DIM=1)

    !Here we write results for all years of the run, divided into age groups. These are the main results, and will be saved
    !alongside heritabilities and parameters within the same netCDF file.
    IF (WriteRes .EQ. 1)THEN
        WriteResThisYear = .FALSE.
        !IF (ReadInitialPopulation .EQ. 2 .AND. t .EQ. 1) WriteResThisYear = .TRUE.
        IF (t .EQ. 1) WriteResThisYear = .TRUE.
        IF (MOD(t, WriteInterval) .EQ. 0) WriteResThisYear = .TRUE.

        IF (WriteResThisYear .EQV. .TRUE.) THEN
        DO ResultsAge=1, MaxAge
            
            !Sort out individuals of the right age, essentially Subset() from R
            PopAgeSorted = 0.
            IMResultsAgeSorted = 0.

            !Age sort Pop()
            DO, Ind = 1, MaxInd
                IF (Pop(Ind, 2) .EQ. ResultsAge) THEN
                    PopAgeSorted(Ind,:) = Pop(Ind,:)
                ENDIF
            ENDDO !Subsetting Pop()

            !Age sort IntermediateResults()
            DO, Ind = 1, MaxInd
                IF (Pop(Ind, 2) .EQ. ResultsAge) THEN
                    IMResultsAgeSorted(Ind,:) = IntermediateResults(Ind,:)
                ENDIF
            ENDDO !Subsetting IntermediateResults()
            
            IF (t .EQ. 1) THEN
                !Put together interesting results
                count_popagesorted = COUNT(PopAgeSorted(:,1)/=0.,DIM=1)
                count_intermediateresults = COUNT(IMResultsAgeSorted(:,1)/=0.,DIM=1)

                ncdf_PopResults(1, ResultsAge, 1)  = dble(t) + dble(last_run_time)                                  !Year
                ncdf_PopResults(1, ResultsAge, 2)  = dble(ResultsAge)                                                !Age
                ncdf_PopResults(1, ResultsAge, 3)  = Temp                                                    !Temperature
                ncdf_PopResults(1, ResultsAge, 4)  = FoodEnv                                            !Food environment
                ncdf_PopResults(1, ResultsAge, 5)  = dble(NoRecr)                                     !Number of recruits 
                ncdf_PopResults(1, ResultsAge, 6)  = COUNT(PopAgeSorted(:,1)==1.,DIM=1)              !n survived age/year
                ncdf_PopResults(1, ResultsAge, 7)  = COUNT(PopAgeSorted(:,1) .LT. 0.,DIM=1)              !n died in age/year
                ncdf_PopResults(1, ResultsAge, 8)  = SUM(PopAgeSorted(:,4),DIM=1)&                             !SomWeight
                                                                &/ count_popagesorted
                ncdf_PopResults(1, ResultsAge, 9)  = std_dev(PopAgeSorted(:,4),PopAgeSorted(:,1))           !Somweight_SD
                ncdf_PopResults(1, ResultsAge, 10) = SUM(PopAgeSorted(:,7),DIM=1)&                             !GonWeight
                                                                &/ count_popagesorted
                ncdf_PopResults(1, ResultsAge, 11) = std_dev(PopAgeSorted(:,7),PopAgeSorted(:,1))           !Gonweight_SD
                ncdf_PopResults(1, ResultsAge, 12) = SUM(PopAgeSorted(:,3),DIM=1)&                                !Length
                                                                &/ count_popagesorted
                ncdf_PopResults(1, ResultsAge, 13) = std_dev(PopAgeSorted(:,3),PopAgeSorted(:,1))              !Length_SD
                ncdf_PopResults(1, ResultsAge, 14) = SUM(IMResultsAgeSorted(:,2),DIM=1)&                     !M_predation
                                                                &/ count_intermediateresults
                ncdf_PopResults(1, ResultsAge, 15) = std_dev(IMResultsAgeSorted(:,2),IMResultsAgeSorted(:,1))   !Mpred_SD
                ncdf_PopResults(1, ResultsAge, 16) = SUM(IMResultsAgeSorted(:,3),DIM=1)&                      !M_foraging
                                                                &/ count_intermediateresults
                ncdf_PopResults(1, ResultsAge, 17) = std_dev(IMResultsAgeSorted(:,3),IMResultsAgeSorted(:,1))    !Mfor_SD
                ncdf_PopResults(1, ResultsAge, 18) = SUM(IMResultsAgeSorted(:,4),DIM=1)&                  !M_reproduction
                                                                &/ count_intermediateresults
                ncdf_PopResults(1, ResultsAge, 19) = std_dev(IMResultsAgeSorted(:,4),IMResultsAgeSorted(:,1))   !Mrepr_SD
                ncdf_PopResults(1, ResultsAge, 20) = SUM(IMResultsAgeSorted(:,5),DIM=1)&                   !M_respiration
                                                                &/ count_intermediateresults
                ncdf_PopResults(1, ResultsAge, 21) = std_dev(IMResultsAgeSorted(:,5),IMResultsAgeSorted(:,1))   !Mresp_SD
                ncdf_PopResults(1, ResultsAge, 22) = SUM(IMResultsAgeSorted(:,6),DIM=1)&                        !FishMort
                                                                &/ count_intermediateresults
                ncdf_PopResults(1, ResultsAge, 23) = std_dev(IMResultsAgeSorted(:,6),IMResultsAgeSorted(:,1))!FishMort_SD
                ncdf_PopResults(1, ResultsAge, 24) = SUM(IMResultsAgeSorted(:,7),DIM=1)&              !Foraging intensity
                                                                &/ count_intermediateresults
                ncdf_PopResults(1, ResultsAge, 25) = std_dev(IMResultsAgeSorted(:,7),IMResultsAgeSorted(:,1))!Foraging_SD
                ncdf_PopResults(1, ResultsAge, 26) = COUNT(PopAgeSorted(:,6)>0.,DIM=1)                          !n mature
                ncdf_PopResults(1, ResultsAge, 27) = SUM(IMResultsAgeSorted(:, 8),DIM=1)              !n skipped_spawning
                ncdf_PopResults(1, ResultsAge, 28) = SUM(PopAgeSorted(:,4),DIM=1)                      !Age-group biomass
                ncdf_PopResults(1, ResultsAge, 29) = SUM(IMResultsAgeSorted(:, 9),DIM=1)&                 !Somatic Growth
                                                                &/ count_intermediateresults
                ncdf_PopResults(1, ResultsAge, 30) = std_dev(IMResultsAgeSorted(:,9),IMResultsAgeSorted(:,1))  !SomGro_SD
                ncdf_PopResults(1, ResultsAge, 31) = SUM(IMResultsAgeSorted(:,10),DIM=1)&                   !Gonad Growth
                                                                &/ count_intermediateresults
                ncdf_PopResults(1, ResultsAge, 32) = std_dev(IMResultsAgeSorted(:,10),IMResultsAgeSorted(:,1)) !GonGro_SD
                ncdf_PopResults(1, ResultsAge, 33) = SUM(PopAgeSorted(:,13),DIM=1)&                     !Genetic appetite
                                                                &/ count_popagesorted
                ncdf_PopResults(1, ResultsAge, 34) = std_dev(PopAgeSorted(:,13),PopAgeSorted(:,1))           !Appetite_SD
                ncdf_PopResults(1, ResultsAge, 35) = SUM(PopAgeSorted(:,12),DIM=1)&                   !Somatic allocation
                                                                &/ count_popagesorted
                ncdf_PopResults(1, ResultsAge, 36) = std_dev(PopAgeSorted(:,12),PopAgeSorted(:,1))         !Allocation_SD
                ncdf_PopResults(1, ResultsAge, 37) = SUM(PopAgeSorted(:,9),DIM=1)&                        !PMRN intercept
                                                                &/ count_popagesorted
                ncdf_PopResults(1, ResultsAge, 38) = std_dev(PopAgeSorted(:,9),PopAgeSorted(:,1))           !Intercept_SD
                ncdf_PopResults(1, ResultsAge, 39) = SUM(PopAgeSorted(:,6),DIM=1,MASK=PopAgeSorted(:,6).NE.0.)& !Ageatmat
                                                                &/ count(PopAgeSorted(:,6).NE.0)
                ncdf_PopResults(1, ResultsAge, 40) = std_dev(PopAgeSorted(:,6),PopAgeSorted(:,6))          !Age-at-mat_SD
                ncdf_PopResults(1, ResultsAge, 41) = SUM(PopAgeSorted(:,19),DIM=1,MASK=PopAgeSorted(:,6).NE.0.)&!Lenatmat
                                                                &/ count(PopAgeSorted(:,6).NE.0)
                ncdf_PopResults(1, ResultsAge, 42) = std_dev(PopAgeSorted(:,19),PopAgeSorted(:,6))         !Len-at-mat_SD
                ncdf_PopResults(1, ResultsAge, 43) = SUM(PopAgeSorted(:,8),DIM=1,&                          !CumGonW
                                                                & MASK=PopAgeSorted(:,1).EQ.-1. .AND. PopAgeSorted(:,6) .GT. 0.)& 
                                                                &/ count(PopAgeSorted(:,1).EQ.-1 .AND. PopAgeSorted(:,6) .GT. 0.)
                !ncdf_PopResults(1, ResultsAge, 43) = SUM(PopAgeSorted(:,8),DIM=1,MASK=PopAgeSorted(:,1).EQ.-1.)&!CumGonW
                                                                !&/ count(PopAgeSorted(:,1).EQ. -1)
                ncdf_PopResults(1, ResultsAge, 44) = std_dev(PopAgeSorted(:,8),PopAgeSorted(:,1))        !CumGonWeight_SD
                ncdf_PopResults(1, ResultsAge, 45) = COUNT(PopAgeSorted(:,1) .EQ. -2.,DIM=1)            !Number of fish caught
                ncdf_PopResults(1, ResultsAge, 46) = SUM(PopAgeSorted(:,4),DIM=1,MASK=PopAgeSorted(:,1).EQ.-2.)!Fisheries yield
                ncdf_PopResults(1, ResultsAge, 47) = SUM(PopAgeSorted(:,4),DIM=1, MASK=PopAgeSorted(:,6).GT.0.)!SSB
                ncdf_PopResults(1, ResultsAge, 48) = SUM(PopAgeSorted(:,4),DIM=1)!SSB
                !Checklist: If adding any additional traits, also remember to update the allocation of ncdf_PopResults(),
                !IntermediateResults() and IMResultsAgeSorted().
            ELSE
            !Put together interesting results
            count_popagesorted = COUNT(PopAgeSorted(:,1)/=0.,DIM=1)
            count_intermediateresults = COUNT(IMResultsAgeSorted(:,1)/=0.,DIM=1)

            IF (count_popagesorted .GE. 1.) THEN !If we actually have living individuals in this age group

            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 1)  = dble(t) + dble(last_run_time)                              !Year
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 2)  = dble(ResultsAge)                                            !Age
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 3)  = Temp                                                !Temperature
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 4)  = FoodEnv                                        !Food environment
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 5)  = dble(NoRecr)                                 !Number of recruits 
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 6)  = COUNT(PopAgeSorted(:,1)==1.,DIM=1)          !n survived age/year
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 7)  = COUNT(PopAgeSorted(:,1) .LT. 0.,DIM=1)       !n died in age/year
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 8)  = SUM(PopAgeSorted(:,4),DIM=1)&                         !SomWeight
                                                            &/ count_popagesorted
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 9)  = std_dev(PopAgeSorted(:,4),PopAgeSorted(:,1))       !Somweight_SD
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 10) = SUM(PopAgeSorted(:,7),DIM=1)&                         !GonWeight
                                                            &/ count_popagesorted
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 11) = std_dev(PopAgeSorted(:,7),PopAgeSorted(:,1))       !Gonweight_SD
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 12) = SUM(PopAgeSorted(:,3),DIM=1)&                            !Length
                                                            &/ count_popagesorted
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 13) = std_dev(PopAgeSorted(:,3),PopAgeSorted(:,1))          !Length_SD
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 14) = SUM(IMResultsAgeSorted(:,2),DIM=1)&                 !M_predation
                                                            &/ count_intermediateresults
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 15) =std_dev(IMResultsAgeSorted(:,2),IMResultsAgeSorted(:,1))!Mpred_SD
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 16) = SUM(IMResultsAgeSorted(:,3),DIM=1)&                  !M_foraging
                                                            &/ count_intermediateresults
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 17) = std_dev(IMResultsAgeSorted(:,3),IMResultsAgeSorted(:,1))!Mfor_SD
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 18) = SUM(IMResultsAgeSorted(:,4),DIM=1)&              !M_reproduction
                                                            &/ count_intermediateresults
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 19) =std_dev(IMResultsAgeSorted(:,4),IMResultsAgeSorted(:,1))!Mrepr_SD
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 20) = SUM(IMResultsAgeSorted(:,5),DIM=1)&               !M_respiration
                                                            &/ count_intermediateresults
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 21) =std_dev(IMResultsAgeSorted(:,5),IMResultsAgeSorted(:,1))!Mresp_SD
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 22) = SUM(IMResultsAgeSorted(:,6),DIM=1)&                    !FishMort
                                                            &/ count_intermediateresults
            ncdf_PopResults((t/WriteInterval)+1,ResultsAge,23)=std_dev(IMResultsAgeSorted(:,6),IMResultsAgeSorted(:,1))!FishMort_SD
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 24) = SUM(IMResultsAgeSorted(:,7),DIM=1)&          !Foraging intensity
                                                            &/ count_intermediateresults
            ncdf_PopResults((t/WriteInterval)+1,ResultsAge,25)=std_dev(IMResultsAgeSorted(:,7),IMResultsAgeSorted(:,1))!Foraging_SD
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 26) = COUNT(PopAgeSorted(:,6)>0.,DIM=1)                      !n mature
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 27) = SUM(IMResultsAgeSorted(:, 8),DIM=1)          !n skipped_spawning
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 28) = SUM(PopAgeSorted(:,4),DIM=1)                  !Age-group biomass
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 29) = SUM(IMResultsAgeSorted(:, 9),DIM=1)&             !Somatic Growth
                                                            &/ count_intermediateresults
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 30)=std_dev(IMResultsAgeSorted(:,9),IMResultsAgeSorted(:,1))!SomGro_SD
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 31) = SUM(IMResultsAgeSorted(:,10),DIM=1)&               !Gonad Growth
                                                            &/ count_intermediateresults
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge,32)=std_dev(IMResultsAgeSorted(:,10),IMResultsAgeSorted(:,1))!GonGro_SD
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 33) = SUM(PopAgeSorted(:,13),DIM=1)&                 !Genetic appetite
                                                            &/ count_popagesorted
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 34) = std_dev(PopAgeSorted(:,13),PopAgeSorted(:,1))       !Appetite_SD
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 35) = SUM(PopAgeSorted(:,12),DIM=1)&               !Somatic allocation
                                                            &/ count_popagesorted
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 36) = std_dev(PopAgeSorted(:,12),PopAgeSorted(:,1))     !Allocation_SD
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 37) = SUM(PopAgeSorted(:,9),DIM=1)&                    !PMRN intercept
                                                            &/ count_popagesorted
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 38) = std_dev(PopAgeSorted(:,9),PopAgeSorted(:,1))       !Intercept_SD
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge,39)=SUM(PopAgeSorted(:,6),DIM=1,MASK=PopAgeSorted(:,6).NE.0.)&!Ageatmat
                                                            &/ count(PopAgeSorted(:,6).NE.0)
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 40) = std_dev(PopAgeSorted(:,6),PopAgeSorted(:,6))      !Age-at-mat_SD
            ncdf_PopResults((t/WriteInterval)+1,ResultsAge,41)=SUM(PopAgeSorted(:,19),DIM=1,MASK=PopAgeSorted(:,6).NE.0.)&!Lenatmat
                                                            &/ count(PopAgeSorted(:,6).NE.0)
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 42) = std_dev(PopAgeSorted(:,19),PopAgeSorted(:,6))     !Len-at-mat_SD
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 43) = SUM(PopAgeSorted(:,8),DIM=1,&                           !CumGonW
                                                            & MASK=PopAgeSorted(:,1).EQ.-1. .AND. PopAgeSorted(:,6) .GT. 0.)& 
                                                            &/ count(PopAgeSorted(:,1).EQ.-1 .AND. PopAgeSorted(:,6) .GT. 0.)
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 44) = std_dev(PopAgeSorted(:,8),PopAgeSorted(:,1))    !CumGonWeight_SD
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 45) = COUNT(PopAgeSorted(:,1) .EQ. -2.,DIM=1)   !Number of fish caught
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 46) = SUM(PopAgeSorted(:,4),DIM=1,MASK=PopAgeSorted(:,1).EQ.-2.)!Yield
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 47) = SUM(PopAgeSorted(:,4),DIM=1, MASK=PopAgeSorted(:,6).GT.0.)  !SSB
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 48) = SUM(PopAgeSorted(:,4),DIM=1)                      !Total biomass
            !Checklist: If adding any additional traits, also remember to update the allocation of ncdf_PopResults(),
            !IntermediateResults() and IMResultsAgeSorted().
            ELSE !if we don't have any living individuals in this age group
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 1)  = dble(t) + dble(last_run_time)                              !Year
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 2)  = dble(ResultsAge)                                            !Age
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 3)  = Temp                                                !Temperature
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 4)  = FoodEnv                                        !Food environment
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 5)  = dble(NoRecr)                                 !Number of recruits 
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 6)  = COUNT(PopAgeSorted(:,1)==1.,DIM=1)          !n survived age/year
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 7)  = COUNT(PopAgeSorted(:,1) .LT. 0.,DIM=1)       !n died in age/year
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 8)  = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 9)  = 0.  
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 10) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 11) = 0.   
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 12) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 13) = 0. 
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 14) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 15) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 16) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 17) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 18) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 19) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 20) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 21) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 22) = 0.
            ncdf_PopResults((t/WriteInterval)+1,ResultsAge,23)= 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 24) = 0.
            ncdf_PopResults((t/WriteInterval)+1,ResultsAge,25)= 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 26) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 27) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 28) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 29) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 30)= 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 31) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge,32)= 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 33) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 34) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 35) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 36) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 37) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 38) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge,39)= 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 40) = 0.
            ncdf_PopResults((t/WriteInterval)+1,ResultsAge,41)= 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 42) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 43) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 44) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 45) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 46) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 47) = 0.
            ncdf_PopResults((t/WriteInterval)+1, ResultsAge, 48) = 0.
            !Checklist: If adding any additional traits, also remember to update the allocation of ncdf_PopResults(),
            !IntermediateResults() and IMResultsAgeSorted().
            ENDIF
            ENDIF
            !WRITE PopResults() TABLE TO NETCDF FILE 1
            ncdf_status = nf90_open(ncdf_datafile, nf90_WRITE, ncid)
            ncdf_status = nf90_put_var(ncid, varid_val, ncdf_PopResults)
            ncdf_status = nf90_close(ncid = ncid)
        ENDDO
        ENDIF
    ENDIF   



    IF (t .NE. Horizon) THEN !Don't age in final year, for FinalPop() saving puprposes
    !Ageing and new status
        DO Ind = 1, NoAlive             !Aging loop
            Status      = Pop(Ind, 1)
            Age         = Pop(Ind, 2)
            
            Age = Age + 1. !Add 1 to age
            
            IF (Age .GT. MaxAge) THEN   !If above max age, set status to -1 (dead)
                Status = -1.
            END IF

            Pop(Ind, 1) = Status
            Pop(Ind, 2) = Age 
        ENDDO                           !Aging loop
    ENDIF !Don't age in final year for saving puprposes

    !KILL
    DO Ind = 1, MaxInd                  !Remove dead individuals (set all columns to 0.)
        IF (Pop(Ind,1) .LT. 0.) THEN   !If individuals are dead
            Pop(Ind, :) = 0.
        ENDIF
    ENDDO                               !Remove dead individuals

    IF (MOD(t, WriteInterval)==0) THEN
        Heritability_holder (t/WriteInterval, 1)    = dble(t)
        Heritability_holder (t/WriteInterval, 2)    = p_corr(Pop(:,14),Pop(:, 6),Pop(:, 6)) !Age at mat
        Heritability_holder (t/WriteInterval, 3)    = p_corr(Pop(:,20),Pop(:,19),Pop(:, 6)) !Length at mat

        ncdf_status = nf90_open(ncdf_datafile, nf90_WRITE, ncid)
        ncdf_status = nf90_put_var(ncid, run_heritability_val, Heritability_holder)
        ncdf_status = nf90_close(ncid = ncid)
    ENDIF

    !RE-INDEX goes here                                     !Place empty rows below non-empty rows. This is done so "1:NoAlive" 
    CALL ARRAY_INDEX( Pop( : , 1), index_holder )           !will only affect non-empty rows, i.e. living individuals.
    index_holder(1:MaxInd) = index_holder(MaxInd:1:-1)
    Pop(:,:) = Pop( index_holder , :)

ENDDO !Time loop
!----------------------------------------------------------------------------------------------------------------------------------
!WRITE TO FILES
!----------------------------------------------------------------------------------------------------------------------------------
WRITE(char_lastruntime, fmt='(i10)') last_run_time
WRITE(char_horizon, fmt='(i10)') Horizon

CALL cpu_time(time_spent)
PRINT '(/,5x,A, f8.3, a4)',  "Run exited with no errors after:", time_spent/60, "min"
!IF (WriteRes .EQ. 1) 
IF (WriteRes .EQ. 1) PRINT '(5x,A,A,/)', "Results for this run are saved in: ", ncdf_datafile

!Write binary files, only used when we wish to continue the population where last run left off.
open(newunit=u,file='finalpop.dat',form='unformatted')
    write(u) Pop(:,:) !Writes entire population table
close(u)
open(newunit=o,file='last_run_time.dat',form='unformatted')
    write(o) Horizon + last_run_time !Keeps track of the end year of current run
close(o)
IF (WriteCSV .EQ. 1) CALL CSV_MATRIX_WRITE(Pop(:,:),"finalpop.csv")
results_list(run_looper) = ncdf_datafile

!----------------------------------------------------------------------------------------------------------------------------------
!END THE CURRENT RUN
!----------------------------------------------------------------------------------------------------------------------------------
ENDDO
PRINT '(/,5x,A)', "Successfully executed all runs."
PRINT '(5x, A,/,5x,A)', "Results for runs, in order, are: ", results_list
PRINT '(5x,A,/)', "    <・ )))><<         ♥   >^)))<～～"
!----------------------------------------------------------------------------------------------------------------------------------
!EXTRA FUNCTIONS USED WITHIN THE CODE
!Non-instrinsic functions that I needed to build myself are included below
!----------------------------------------------------------------------------------------------------------------------------------
CONTAINS
recursive integer function random_pull (sorted_array, col_number, target_num) result(ran_pull_result)
    !This function is intended to find the first occurence that is equal to or greater than a targeted number (target_num) in a
    !specified column (col_number) in an array where the specified column in sorted ascending (sorted_array).
    !The function then returns the value of the first column, assumed to be an ID-number, corresponding to the desired value in
    !col_number as an integer.
    IMPLICIT NONE 
    real(srp), dimension(:,:), intent(in)   :: sorted_array
    real(srp), intent(in)                   :: target_num
    integer, intent(in)                     :: col_number
    real(srp)                               :: testing_val
    integer :: array_size, array_midpoint
    ran_pull_result = 0
    array_size = SIZE(sorted_array, DIM=1)
    IF (array_size .LE. 2) THEN !When the two lines surrounding the targeted number is reached
        ran_pull_result = nint(sorted_array(2,1))
        RETURN
    ENDIF

    array_midpoint = nint(real(array_size)/2)
    testing_val = sorted_array(array_midpoint,col_number)

    !If the midpoint value is identical to the desired number (unlikely) 
    IF (testing_val .EQ. target_num) THEN
        ran_pull_result = nint(sorted_array(array_midpoint,1))
        RETURN
    ENDIF
    !Otherwise, proceed with narrowing down
    IF (target_num .LT. testing_val) THEN
        ran_pull_result = random_pull(sorted_array(1:array_midpoint,:), col_number, target_num)
        RETURN
    ELSEIF (target_num .GT. testing_val) THEN
        ran_pull_result = random_pull(sorted_array(array_midpoint:array_size,:), col_number, target_num)
        RETURN
    ENDIF
end function random_pull

real(srp) function std_dev (array_in, mask_array)
    !Returns the (population) standard deviation of a one-dimensional array (array_in) containing real values.
    !Optionally, a masking array (mask_array) can be included, which must be of the same size as array_in. 
    !Any value of "0." within the masking array tells the function to ignore the corresponding entry in array_in.
    real(srp), dimension(:), intent(in)             :: array_in
    real(srp), dimension(:), intent(in), optional   :: mask_array
    real(srp)                                       :: sumsq_holder(size(array_in))
    real(srp)                                       :: sumsq, array_mean, array_length
    integer                                         :: stdev_looper, nonzero_count

    IF( present( mask_array )) THEN !Masked standard deviation
        IF( SIZE(array_in) .NE. SIZE(mask_array) ) PRINT *, "Warning: Masking array was not the same size as stdev array."
        array_length = SIZE( array_in )
        nonzero_count = COUNT( mask_array .NE. 0)
        array_mean = SUM( array_in, MASK=mask_array .NE. 0) / nonzero_count

        DO stdev_looper=1, nint(array_length)
            sumsq_holder(stdev_looper) = (array_in(stdev_looper) - array_mean)**2
        ENDDO
        sumsq = SUM(sumsq_holder, MASK=mask_array .NE. 0)

        std_dev = sqrt(sumsq / nonzero_count)
    ELSE !Un-masked standard deviation
        array_length = SIZE(array_in)
        array_mean = SUM(array_in) / array_length

        DO stdev_looper=1, nint(array_length)
            sumsq_holder(stdev_looper) = (array_in(stdev_looper) - array_mean)**2
        ENDDO
        sumsq = SUM(sumsq_holder)

        std_dev = sqrt(sumsq / array_length)
    ENDIF
end function std_dev

real(srp) function p_corr (array_in1, array_in2, mask_array)
    !Returns the Pearson correlation of two 1D arrays, which must be of the same size. Includes an optional masking array, with 
    !all entires where mask_array=0 being ignored.
    real(srp), dimension(:), intent(in)             :: array_in1, array_in2
    real(srp), dimension(:), intent(in), optional   :: mask_array
    real(srp)   :: x_sqr(size(array_in1)), y_sqr(size(array_in1)), xy(size(array_in1))
    integer     :: n_count, pcorr_looper

    IF( SIZE(array_in1) .NE. SIZE(array_in2) ) PRINT *, "Warning: arrays not same size"
    IF( present( mask_array )) THEN
        n_count = count(mask_array .NE. 0.)
        DO pcorr_looper=1, size(array_in1)
            x_sqr(pcorr_looper) = array_in1(pcorr_looper)**2
            y_sqr(pcorr_looper) = array_in2(pcorr_looper)**2
            xy(pcorr_looper)    = array_in1(pcorr_looper)*array_in2(pcorr_looper)
        ENDDO
        
        p_corr = ( n_count*SUM(xy,MASK=mask_array.NE.0.) - SUM(array_in1,MASK=mask_array.NE.0.)* &
        & SUM(array_in2,MASK=mask_array.NE.0.) ) / (SQRT( ( n_count*SUM(x_sqr,MASK=mask_array.NE.0.) - &
        & SUM(array_in1,MASK=mask_array.NE.0.)**2 ) * (n_count*SUM(y_sqr,MASK=mask_array.NE.0.) - &
        & SUM(array_in2,MASK=mask_array.NE.0.)**2) ) )
     ELSE
        n_count = size(array_in1)
        DO pcorr_looper=1, n_count
            x_sqr(pcorr_looper) = array_in1(pcorr_looper)**2
            y_sqr(pcorr_looper) = array_in2(pcorr_looper)**2
            xy(pcorr_looper)    = array_in1(pcorr_looper)*array_in2(pcorr_looper)
        ENDDO
        
        p_corr = ( n_count*SUM(xy) - SUM(array_in1)*SUM(array_in2) ) / &
        & (SQRT( ( n_count*SUM(x_sqr) - SUM(array_in1)**2 ) * (n_count*SUM(y_sqr) - SUM(array_in2)**2) ) )
    ENDIF
end function p_corr

END PROGRAM evoeco_ibm