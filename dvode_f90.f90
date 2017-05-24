MODULE DVODE_F90_M
! This version is the December 2005 release.
! Last change: 12/06/05
! _____________________________________________________________________
! Working Precision
IMPLICIT NONE
! Define the working precision for DVODE_F90. Change D0 to E0 in the
! next statement to convert to single precision.
INTEGER, PARAMETER, PRIVATE :: WP = KIND(1.0D0)
! ______________________________________________________________________
! Overview

! The f77 ordinary differential equation solver VODE.f is applicable to
! nonstiff systems of odes and to stiff systems having dense or banded
! Jacobians. DVODE_F90 is a Fortran 90 extension of VODE.f. While
! retaining all of the features available in VODE.f, we have
! incorporated several new options in DVODE_F90 including:
!   1. the ability to solve stiff systems with sparse Jacobians
!   2. internal management of storage and work arrays
!   3. specification of options via optional keywords
!   4. the ability to perform root finding or "event detection"
!   5. various new diagnostic and warning messages
!   6. the ability to impose solution bounds
!   7. several specialized options for dealing with sparsity
! ______________________________________________________________________
! Version Information

! This is DVODE_F90, the double precision FORTRAN 90 extension of the
! f77 DVODE.f ordinary differential equation solver. This version uses
! MA28 for sparse Jacobians. This file and related information can be
! obtained at the following support page:
!
!     http://www.radford.edu/~thompson/vodef90web/
!
! We are indebted to Richard Cox (ORNL) for providing us with his
! implementation of MA28 in LSOD28.f (a variant of Alan Hindmarsh's
! lsodes.f). We are indebted to Alan Hindmarsh for numerous contributions.
! In particular, we borrowed liberally from the f77 solvers VODE.f,
! LSODAR.f, and LSODES.f while developing DVODE_F90. We are indebted
! to Doug Salane for providing us with his JACSP Jacobian routines.
!
! If you find a bug or encounter a problem with DVODE_F90, please
! contact one of us:
!    G.D. Byrne (gbyrne@wi.rr.com)
!    S. Thompson (thompson@radford.edu)
! A set of quick start instructions is provided below.
! ______________________________________________________________________
! Note to LAPACK and MA48 Users

! We can provide a version of DVODE_F90 that allows the use of LAPACK
! for dense and banded Jacobians and the use of MA48 for sparse Jacobians.
! (Alternatively, this file contains the changes necessary. The changes
! may be located by searching on the string 'build change'.
! Using LAPACK with DVODE_F90 requires the file 'lapackd_f90_m.f90' which
! we can provide. Alternatively, you can use another LAPACK F90 module.
! If the module has a name other than lapackd_f90_m.f90, you will need
! to modify a USE statement below. You may wish to link to LAPACK routines
! that are not contained in an F90 module by deleting the USE statement,
! provided that your F90 compiler allows you do this. Should you do this,
! then depending on the version of LAPACK to which you are linking, it may
! be necessary to change the array declarations in subroutines DVJAC and
! DVSOL below. We are very interested in hearing from anyone who uses
! LAPACK from a source other than the file lapackd_f90_m.f90.
! Using MA48 DVODE_F90 requires the file 'jacobian_for_ma48.f90' which
! we can provide. In addition, you must have licensed access to Harwell's
! HSL library which contains MA48. Please note that MA48 is not a part
! of DVODE_F90 and is not distributed with DVODE_F90. We are grateful
! to Harwell for allowing us the opportunity to provide an MA48 based
! version of DVODE_F90 for individuals who have access to HSL. We are
! very interested in learning of your experience and suggestions for
! improvement if you decide to try the MA48 based solution option.
! ______________________________________________________________________
! Note on F90/F95 Compilers

! To date we have used DVODE_F90 successfully with all F90/F95 compilers
! to which we have access. In particular, we have used it with the Lahey
! F90 and Lahey-Fujitsu F95 compilers, the Compaq Visual F90 compiler,
! the g90 compiler, and the INTEL, Portland, Salford, and SUN compilers.
! It should be noted that compilers such as Salford's FTN95 complain
! about uninitialized arrays passed as subroutine arguments and the use of
! slices of two dimensional arrays as one dimensional vectors, and will
! not run using the strictest compiler options. It is perfectly safe to
! use the /-CHECK compiler option to avoid these FTN95 runtime checks.
! DVODE_F90 does not use any variable for numerical purposes until it
! has been assigned an appropriate value.
! ______________________________________________________________________
! Quick Start Instructions

! (1) Compile this file. Then compile, link, and execute the program
!     example1.f90. The output is written to the file example1.dat.
!     Verify that the last line of the output is the string
!     'No errors occurred.'
! (2) Repeat this process for the program example2.f90.
!
! Other test programs you may wish to run to verify your installation
! of DVODE_F90 are:
!
! (3) Run the test programs nonstiffoptions.f90 and stiffoptions.f90
!     and verify that the last line in the output files produced is
!     'No errors occurred.' They solve the problems in the Toronto
!     test suites using several different error tolerances and various
!     solution options. Note that stiffoptions.f90 takes several
!     minutes to run because it performs several thousand separate
!     integrations.
! (4) Locate the file robertson.f90 in the demo programs and look at
!     how options are set using SET_OPTS, how DVODE_F90 is called to
!     obtain the solution at desired output times, and how the
!     derivative and Jacobian routines are supplied. Note too the
!     manner in which the solution is constrained to be nonnegative.
! (5) Locate demoharmonic.f90 and look at how root finding options
!     are set and how the event residual routine is supplied to
!     DVODE_F90.
! (6) The other demo programs available from the DVODE_F90 support
!     page illustrate various other solution options available in
!     DVODE_F90. The demo programs may be obtained from
!
!        http://www.radford.edu/~thompson/vodef90web/index.html
! ______________________________________________________________________
! DVODE_F90 Full Documentation Prologue

! Section 1.  Setting Options in DVODE_F90
! Section 2.  Calling DVODE_F90
! Section 3.  Choosing Error Tolerances
! Section 4.  Choosing the Method of Integration
! Section 5.  Interpolation of the Solution and Derivatives
! Section 6.  Handling Events (Root Finding)
! Section 7.  Gathering Integration Statistics
! Section 8.  Determining Jacobian Sparsity Structure Arrays
! Section 9.  Original DVODE Documentation Prologue
! Section 10. Example Usage

! Note: Search on the string 'Section' to locate these sections. You
! may wish to refer to the support page which has the sections broken
! into smaller pieces.
! ______________________________________________________________________
! Section 1.  Setting Options in DVODE_F90
!
! You can use any of three options routines:
!
! SET_NORMAL_OPTS
! SET_INTERMEDIATE_OPTS
! SET_OPTS

! OPTIONS = SET_NORMAL_OPTS(DENSE_J, BANDED_J, SPARSE_J,               &
!   USER_SUPPLIED_JACOBIAN, LOWER_BANDWIDTH, UPPER_BANDWIDTH,          &
!   RELERR, ABSERR, ABSERR_VECTOR, NEVENTS)

! OPTIONS = SET_INTERMEDIATE_OPTS(DENSE_J, BANDED_J, SPARSE_J,         &
!   USER_SUPPLIED_JACOBIAN,LOWER_BANDWIDTH, UPPER_BANDWIDTH,           &
!   RELERR, ABSERR, ABSERR_VECTOR,TCRIT, H0, HMAX, HMIN, MAXORD,       &
!   MXSTEP, MXHNIL, NZSWAG, USER_SUPPLIED_SPARSITY, MA28_RPS,          &
!   NEVENTS, CONSTRAINED, CLOWER, CUPPER, CHANGE_ONLY_f77_OPTIONS)     &

! OPTIONS = SET_OPTS(METHOD_FLAG, DENSE_J, BANDED_J, SPARSE_J,         &
!   USER_SUPPLIED_JACOBIAN, SAVE_JACOBIAN, CONSTANT_JACOBIAN,          &
!   LOWER_BANDWIDTH, UPPER_BANDWIDTH, SUB_DIAGONALS, SUP_DIAGONALS,    &
!   RELERR, RELERR_VECTOR, ABSERR, ABSERR_VECTOR, TCRIT, H0, HMAX,     &
!   HMIN, MAXORD, MXSTEP, MXHNIL, YMAGWARN, SETH, UPIVOT, NZSWAG,      &
!   USER_SUPPLIED_SPARSITY, NEVENTS, CONSTRAINED, CLOWER, CUPPER,      &
!   MA28_ELBOW_ROOM, MC19_SCALING, MA28_MESSAGES, MA28_EPS,            &
!   MA28_RPS, CHANGE_ONLY_f77_OPTIONS, JACOBIAN_BY_JACSP)

! Please refer to the documentation prologue for each of these functions
! to see what options may be used with each. Note that input to each is
! via keyword and all variables except error tolerances are optional.
! Defaults are used for unspecified options. If an option is available
! in SET_NORMAL OPTS, it is available and has the same meaning in
! SET_INTERMEDIATE_OPTS and SET_OPTS. Similarly, if an option is available
! in SET_INTERMEDIATE_OPTS, it is available and has the same meaning in
! SET_OPTS.

! The first two functions are provided merely for convenience.
! SET_NORMAL_OPTS is available simply to relieve you of reading the
! documentation for SET_OPTS and to use default values for all but
! the most common options. SET_INTERMEDIATE_OPTS is available to allow
! you more control of the integration while still using default values
! for less commonly used options. SET_OPTS allows you to specify any
! of the options available in DVODE_F90.

! Roughly, SET_NORMAL_OPTS is intended to provide for dense, banded,
! and numerical sparse Jacobians without the need to specify other
! specialized options. SET_INTERMEDIATE_OPTIONS is intended to allow
! more general sparse Jacobian options. SET_OPTS is intended to provide
! access to all options in DVODE_F90.

! Please note that SET_INTERMEDIATE_OPTS can be invoked using the same
! arguments as SET_NORMAL_OPTS; and SET_OPTS can be invoked using the
! same arguments as either SET_NORMAL_OPTS or SET_INTERMEDIATE_OPTS.
! If you wish you can simply delete SET_NORMAL_OPTS as well as
! SET_INTERMEDIATE_OPTS and use only SET_OPTS for all problems. If you
! do so, you need only include the options that you wish to use when
! you invoke SET_OPTIONS.

! In the following description any reference to SET_OPTS applies equally
! to SET_NORMAL_OPTS and SET_INTERMEDIATE OPTS.

! Before calling DVODE_F90 for the first time, SET_OPTS must be invoked.
! Typically, SET_OPTS is called once to set the desired integration
! options and parameters. DVODE_F90 is then called in an output loop to
! obtain the solution for the desired times. A detailed description of
! the DVODE_F90 arguments is given in a section below. Detailed descriptions
! of the options available via SET_OPTS are given in the documentation prologue.
! Although each option available in the f77 version of DVODE as well as
! several additional ones are available in DVODE_F90 via SET_OPTS,
! several of the available options are not relevant for most problems
! and need not be specified. Refer to the accompanying demonstration
! programs for specific examples of each usage. Note that after any call
! to DVODE_F90, you may call GET_STATS to gather relevant integration
! statistics. After your problem has completed, you may call
! RELEASE_ARRAYS to deallocate any internal arrays allocated by
! DVODE_F90 and to determine how much storage was used by DVODE_F90.
!
! To communicate with DVODE_F90 you will need to include the following
! statement in your calling program:
!    USE DVODE_F90_M
! and include the following statement in your type declarations section:
!    TYPE(VODE_OPTS) :: OPTIONS
! Below are brief summaries of typical uses of SET_OPTS.
! Nonstiff Problems:
! OPTIONS = SET_OPTS(RELERR=RTOL,ABSERR=ATOL)
!    The above use of SET_OPTS will integrate your system of odes
!    using the nonstiff Adams methods while using a relative error
!    tolerance of RTOL and an absolute error tolerance of ATOL.
!    Your subsequent call to DVODE_F90 might look like:
!    CALL DVODE_F90(F,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS)
! OPTIONS = SET_OPTS(RELERR=RTOL,ABSERR=ATOL,NEVENTS=NG)
!    If you wish to do root finding, SET_OPTS can be used as above.
!    Here, NEVENTS is the desired number of root finding functions.
!    Your subsequent call to DVODE_F90 might look like:
!    CALL DVODE_F90(F,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,G_FCN=G)
!    Here F is the name of your derivative subroutine and G is the
!    name of your subroutine to evaluate residuals of the root
!    finding functions.
! OPTIONS = SET_OPTS(RELERR=RTOL,ABSERR_VECTOR=ATOL)
!    This use of SET_OPTS indicates that a scalar relative error
!    tolerance and a vector of absolute error tolerances will be
!    used.
! Stiff Problems, internally generated dense Jacobian:
! OPTIONS = SET_OPTS(DENSE_J=.TRUE.,RELERR=RTOL,ABSERR=ATOL)
!    This use of DENSE_J=.TRUE. indicates that DVODE_F90 will
!    use the stiffly stable BDF methods and will approximate
!    the Jacobian, considered to be a dense matrix, using
!    finite differences. Your subsequent call to DVODE_F90
!    might look like:
!    CALL DVODE_F90(F,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS)
! OPTIONS = SET_OPTS(DENSE_J=.TRUE.,ABSERR=ATOL,RELERR=RTOL, &
!                    USER_SUPPLIED_JACOBIAN=.TRUE.)
!    If you know the Jacobian and wish to supply subroutine JAC
!    as described in the documentation for DVODE_F90, the options
!    call could look like the above.
!    Your subsequent call to DVODE_F90 might look like:
!    CALL DVODE_F90(F1,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,J_FCN=JAC)
!    Here, JAC is the name of the subroutine that you provide to
!    evaluate the known Jacobian.
! Stiff Problems, internally generated banded Jacobian:
! OPTIONS = SET_OPTS(BANDED_J=.TRUE.,RELERR=RTOL,ABSERR=ATOL, &
!                       LOWER_BANDWIDTH=ML,UPPER_BANDWIDTH=MU)
!    This use of BANDED_J=.TRUE. indicates that DVODE_F90 will
!    use the stiffly stable BDF methods and will approximate the
!    Jacobian, considered to be a banded matrix, using finite
!    differences. Here ML is the lower bandwidth of the Jacobian
!    and ML is the upper bandwidth of the Jacobian.
! Stiff Problems, internally generated sparse Jacobian:
! OPTIONS = SET_OPTS(SPARSE_J=.TRUE.,ABSERR=ATOL,RELERR=RTOL)
!    This use of SET_OPTS indicates that the Jacobian is a sparse
!    matrix. Its structure will be approximated internally by
!    making calls to your derivative routine. If you know the
!    structure before hand, you may provide it directly in a
!    variety of ways as described in the documentation prologue
!    for SET_OPTS. In addition, several other options related
!    to sparsity are available.
! More complicated common usage:
!    Suppose you need to solve a stiff problem with a sparse Jacobian.
!    After some time, the structure of the Jacobian changes and you
!    wish to have DVODE_F90 recalculate the structure before continuing
!    the integration. Suppose that initially you want to use an absolute
!    error tolerance of 1.0D-5 and that when the Jacobian structure is
!    changed you wish to reduce the error tolerance 1.0D-7. Your calls
!    might look like this.
!    RTOL = ...
!    ATOL = 1.0D-5
!    OPTIONS = SET_OPTS(SPARSE_J=.TRUE.,ABSERR=ATOL,RELERR=RTOL)
!    Output loop:
!       CALL DVODE_F90(FCN,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS)
!       At desired time:
!       ISTATE = 3
!       ATOL = 1.0D-7
!       OPTIONS = SET_OPTS(SPARSE_J=.TRUE.,ABSERR=ATOL,RELERR=RTOL)
!    End of output loop

! In the following we have summarized and described how some of the demonstration
! programs set options and call DVODE_F90. In each case the necessary parameters
! are defined before invoking SET_OPTS. The call to DVODE_F90 is in a loop in
! which the output time is successively updated. The actual programs are available
! from the web support page
!
!    http://www.radford. edu/~thompson/vodef90web/index.html/
!
!                              Problem Summary
!
! Problem                  NEQ      Jacobian            Features Illustrated
!
! Prologue Example 1        3        Dense               Basic
!
! Prologue Example 2        3        Dense               Root finding
!
! Robertson                 3        Dense               Solution bounds
!
! Harmonic Oscillator       4        Nonstiff            Root finding
!
! Flow Equations         5-1800      Sparse              Automatic determination
!                                                        of sparsity arrays
!
! Diurnal Kinetics      50-5000      Sparse or banded    Sparsity options
!
!                       Options Used and DVODE_F90 Call
!
! Prologue Example 1
!
!    OPTIONS = SET_NORMAL_OPTS(DENSE_J=.TRUE., ABSERR_VECTOR=ATOL, RELERR=RTOL,     &
!              USER_SUPPLIED_JACOBIAN=.TRUE.)
!    CALL DVODE_F90(FEX,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,J_FCN=JEX)
!
!    The problem consists of a stiff system of NEQ=3 equations. The dense
!    Jacobian option (DENSE_J) is used. A vector ATOL(*) of error tolerances
!    is used. A scalar relative error tolerance RTOL is used. Subroutine JEX
!    is provided to evaluate the analytical Jacobian. If the last argument
!    J_FCN=JEX is omitted (as in Example 2), a numerical Jacobian will
!    be used.
!
! Prologue Example 2
!
!    OPTIONS = SET_NORMAL_OPTS(DENSE_J=.TRUE., RELERR=RTOL, ABSERR_VECTOR=ATOL,     &
!              NEVENTS=NG)
!    CALL DVODE_F90(FEX,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,G_FCN=GEX)
!
!    The system in Example 1 is used to illustrate root finding. It is
!    desired to locate the times at which two of the solution components
!    attain prescribed values. NEVENTS=2 informs the solver that two such
!    functions are used. Subroutine GEX is used to calculate the residuals
!    for these two functions. A dense numerical Jacobian is used.
!
! Robertson Problem
!
!    OPTIONS = SET_INTERMEDIATE_OPTS(DENSE_J=.TRUE., RELERR_VECTOR=RTOL,            &
!              ABSERR_VECTOR=ABSERR_TOLERANCES, CONSTRAINED=BOUNDED_COMPONENTS,     &
!              CLOWER=LOWER_BOUNDS, CUPPER=UPPER_BOUNDS)
!
!    CALL DVODE_F90(DERIVS,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,J_FCN=JACD)
!    The system used in Examples 1 and 2 is solved over a much larger
!    time interval. The solution is constrained to be nonegative. This
!    is done by prescribing the components to be constrained (BOUNDED_COMPONENTS).
!    Artificially large values are used to impose upper bounds (UPPER_BOUNDS)
!    and lower bounds of zero are used to force a nonnegative solution.
!
! Harmonic Oscillator Problem
!
!    OPTIONS = SET_NORMAL_OPTS(RELERR=RTOL, ABSERR=ATOL, NEVENTS=NG)
!    CALL DVODE_F90(DERIVS,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,G_FCN=GEVENTS)
!
!    A nonstiff system of NEQ=4 equations is solved. The nonstiff option is
!    used because neither DENSE_ nor BANDED_J nor SPARSE_J is present. It is
!    desired to find the times at which Y(2) or Y(3) is equal to 0. Residuals
!    for the two corresponding event functions are calculated in subroutine
!    GEVENTS.
!
! Flow Equations Problem
!
!    OPTIONS = SET_OPTS(SPARSE_J=SPARSE, ABSERR=ATOL(1), RELERR=RTOL(1),            &
!              MXSTEP=100000, NZSWAG=20000)
!    CALL DVODE_F90(FCN,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS)
!
!    This is a stiff system of equations resulting a method of lines
!    discretization. The Jacobian is sparse. Scalar absolute and relative
!    error tolerances are used. The Jacobian structure and a numerical
!    Jacobian are used. The solver is limited to a maximum of MXSTEP steps.
!    NZSWAG is the amount by which allocated array sizes will be increased.
!    The accompanying test program may be used to illutrate several other
!    solution options.
!
! Diurnal Kinetics Problem
!
!    OPTIONS = SET_OPTS(SPARSE_J=SPARSE, BANDED_J=BANDED, DENSE_J=DENSE,            &
!              ABSERR_VECTOR=ATOL(1:NEQ), RELERR=RTOL(1), MXSTEP=100000,            &
!              NZSWAG=50000, HMAX=MAXH, LOWER_BANDWIDTH=ML, UPPER_BANDWIDTH=MU,     &
!              MA28_ELBOW_ROOM=10, MC19_SCALING=.TRUE., MA28_MESSAGES=.FALSE.,      &
!              MA28_EPS=1.0D-4, MA28_RPS=.TRUE.,                                    &
!              USER_SUPPLIED_SPARSITY=SUPPLY_STRUCTURE)
!   CALL USERSETS_IAJA(IA, IADIM, JA, JADIM)
!   CALL DVODE_F90(FCN, NEQ, Y, T, TOUT, ITASK, ISTATE, OPTIONS)
!
!   This problem can be used to illustrate most solution options. Here, dense,
!   banded, or sparse Jacobians are used depending on the values of the first
!   three parameters. A vector error tolerance is used and a scalar relative
!   error tolerance is used. If a banded solution is desired, it is necessary
!   to supply the bandwidths ML and MU. If a sparse solution is desired,
!   several special options are used. The most important one is MA28_RPS to
!   force the solver to update the partial pivoting sequence when singular
!   iteration matrices are encountered. The sparsity pattern is determined
!   numerically if SUPPLY_STRUCTURE is FALSE. Otherwise the user will supply
!   the pattern by calling subroutine USERSETS_IAJA.
! ______________________________________________________________________
! Section 2.  Calling DVODE_F90
!
! DVODE_F90 solves the initial value problem for stiff or nonstiff
! systems of first order ODEs,
!     dy/dt = f(t,y), or, in component form,
!     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
! DVODE_F90 is a package based on the EPISODE and EPISODEB packages,
! and on the ODEPACK user interface standard. It was developed from
! the f77 solver DVODE developed by Brown, Byrne, and Hindmarsh.
! DVODE_F90 also provides for the solution of sparse systems in a
! fashion similar to LSODES and LSOD28. Currently, MA28 is used
! to perform the necessary sparse linear algebra. DVODE_F90 also
! contains the provision to do root finding in a fashion similar
! to the LSODAR solver from ODEPACK.

! Communication between the user and the DVODE_F90 package, for normal
! situations, is summarized here. This summary describes only a subset
! of the full set of options available. See the full description for
! details, including optional communication, nonstandard options, and
! instructions for special situations.
!    CALL DVODE_F90(F,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,J_FCN=JAC,G_FCN=GEX)
!    The arguments in the call list to DVODE_F90 have the following
!    meanings.
! F       = The name of the user-supplied subroutine defining the
!           ODE system. The system must be put in the first-order
!           form dy/dt = f(t,y), where f is a vector-valued function
!           of the scalar t and the vector y. Subroutine F is to
!           compute the function f. It is to have the form
!                SUBROUTINE F(NEQ,T,Y,YDOT)
!                DOUBLE PRECISION T,Y(NEQ),YDOT(NEQ)
!           where NEQ, T, and Y are input, and the array YDOT = f(t,y)
!           is output. Y and YDOT are arrays of length NEQ.
!           Subroutine F should not alter Y(1),...,Y(NEQ).
!           If F (and JAC) are not contained in a module available to
!           your calling program, you must declare F to be EXTERNAL
!           in the calling program.
! NEQ     = The size of the ODE system (number of first order
!           ordinary differential equations).
! Y       = A double precision array for the vector of dependent variables,
!           of length NEQ or more. Used for both input and output on the
!           first call (ISTATE = 1), and only for output on other calls.
!           On the first call, Y must contain the vector of initial
!           values. In the output, Y contains the computed solution
!           evaluated at T.
! T       = The independent variable. In the input, T is used only on
!           the first call, as the initial point of the integration.
!           In the output, after each call, T is the value at which a
!           computed solution Y is evaluated (usually the same as TOUT).
!           On an error return, T is the farthest point reached.
! TOUT    = The next value of t at which a computed solution is desired.
!           TOUT is Used only for input. When starting the problem
!           (ISTATE = 1), TOUT may be equal to T for one call, then
!           should not equal T for the next call. For the initial T,
!           an input value of TOUT unequal to T is used in order to
!           determine the direction of the integration (i.e. the
!           algebraic sign of the step sizes) and the rough scale
!           of the problem. Integration in either direction (forward
!           or backward in t) is permitted. If ITASK = 2 or 5 (one-step
!           modes), TOUT is ignored after the first call (i.e. the
!           first call with TOUT \= T). Otherwise, TOUT is required
!           on every call. If ITASK = 1, 3, or 4, the values of TOUT
!           need not be monotone, but a value of TOUT which backs up
!           is limited to the current internal t interval, whose
!           endpoints are TCUR - HU and TCUR. (Refer to the description
!           of GET_STATS for a description of TCUR and HU.)
! ITASK   = An index specifying the task to be performed.
!           Input only. ITASK has the following values and meanings.
!           1  means normal computation of output values of y(t) at
!              t = TOUT (by overshooting and interpolating).
!           2  means take one step only and return.
!           3  means stop at the first internal mesh point at or
!              beyond t = TOUT and return.
!           4  means normal computation of output values of y(t) at
!              t = TOUT but without overshooting t = TCRIT.
!              TCRIT must be specified in your SET_OPTS call. TCRIT
!              may be equal to or beyond TOUT, but not behind it in
!              the direction of integration. This option is useful
!              if the problem has a singularity at or beyond t = TCRIT.
!           5  means take one step, without passing TCRIT, and return.
!              TCRIT must be specified in your SET_OPTS call.
!           If ITASK = 4 or 5 and the solver reaches TCRIT (within
!           roundoff), it will return T = TCRIT(exactly) to indicate
!           this (unless ITASK = 4 and TOUT comes before TCRIT, in
!           which case answers at T = TOUT are returned first).
! ISTATE  = an index used for input and output to specify the
!           the state of the calculation.
!           In the input, the values of ISTATE are as follows.
!           1  means this is the first call for the problem
!              (initializations will be done). See note below.
!           2  means this is not the first call, and the calculation
!              is to continue normally, with no change in any input
!              parameters except possibly TOUT and ITASK.
!           3  means this is not the first call, and the
!              calculation is to continue normally, but with
!              a change in input parameters other than
!              TOUT and ITASK. Desired changes require SET_OPTS
!              be called prior to calling DVODE_F90 again.
!           A preliminary call with TOUT = T is not counted as a
!           first call here, as no initialization or checking of
!           input is done. (Such a call is sometimes useful to
!           include the initial conditions in the output.)
!           Thus the first call for which TOUT is unequal to T
!           requires ISTATE = 1 in the input.
!           In the output, ISTATE has the following values and meanings.
!            1  means nothing was done, as TOUT was equal to T with
!               ISTATE = 1 in the input.
!            2  means the integration was performed successfully.
!            3  means a root of one of your root finding functions
!               has been located.
!           A negative value of ISTATE indicates that DVODE_F90
!           encountered an error as described in the printed error
!           message. Since the normal output value of ISTATE is 2,
!           it does not need to be reset for normal continuation.
!           Also, since a negative input value of ISTATE will be
!           regarded as illegal, a negative output value requires
!           the user to change it, and possibly other input, before
!           calling the solver again.
! OPTIONS = The options structure produced by your call to SET_OPTS.
! JAC     = The name of the user-supplied routine (MITER = 1 or 4 or 6)
!           If you do not specify that a stiff method is to be used
!           in your call to SET_OPTS, you need not include JAC in
!           your call to DVODE_F90. If you specify a stiff method and
!           that a user supplied Jacobian will be supplied, JAC must
!           compute the Jacobian matrix, df/dy, as a function of the
!           scalar t and the vector y. It is to have the form:
!              SUBROUTINE JAC(NEQ, T, Y, ML, MU, PD, NROWPD)
!              DOUBLE PRECISION T, Y(NEQ), PD(NROWPD,NEQ)
!           where NEQ, T, Y, ML, MU, and NROWPD are input and the array
!           PD is to be loaded with partial derivatives (elements of the
!           Jacobian matrix) in the output. PD must be given a first
!           dimension of NROWPD. T and Y have the same meaning as in
!           Subroutine F.
!           In the full matrix case (MITER = 1), ML and MU are
!           ignored, and the Jacobian is to be loaded into PD in
!           columnwise manner, with df(i)/dy(j) loaded into PD(i,j).
!           In the band matrix case (MITER = 4), the elements
!           within the band are to be loaded into PD in columnwise
!           manner, with diagonal lines of df/dy loaded into the rows
!           of PD. Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).
!           ML and MU are the half-bandwidth parameters. (See IUSER).
!           The locations in PD in the two triangular areas which
!           correspond to nonexistent matrix elements can be ignored
!           or loaded arbitrarily, as they are overwritten by DVODE_F90.
!           In the sparse matrix case the elements of the matrix
!           are determined by the sparsity structure given by the
!           IA and JA pointer arrays. Refer to the documentation
!           prologue for SET_OPTS for a description of the arguments
!           for JAC since they differ from the dense and banded cases.
!           JAC need not provide df/dy exactly. A crude
!           approximation (possibly with a smaller bandwidth) will do.
!           In either case, PD is preset to zero by the solver,
!           so that only the nonzero elements need be loaded by JAC.
!           In the sparse matrix case, JAC has a different form:
!                SUBROUTINE JAC (N, T, Y, IA, JA, NZ, PD)
!           Given the number of odes N, the current time T, and the
!           current solution vector Y, JAC must do the following:
!              If NZ = 0 on input:
!              Replace NZ by the number of nonzero elements in the
!              Jacobian. The diagonal of the Jacobian must be included.
!              Do NOT define the arrays IA, JA, PD at this time.
!              Once JAC has been called with NZ = 0 and you have
!              defined the value of NZ, future calls to JAC will use
!              this value of NZ.
!              When a call is made with NZ unequal to 0, you must
!              define the sparsity structure arrays IA and JA, and
!              the sparse Jacobian PD.
!                 IA defines the number of nonzeros including the
!                 diagonal in each column of the Jacobian. Define
!                 IA(1) = 1 and for J = 1,..., N,
!                 IA(J+1) = IA(J) + number of nonzeros in column J.
!                 Diagonal elements must be included even if they are
!                 zero. You should check to ensure that IA(N+1)-1 = NZ.
!                 JA defines the rows in which the nonzeros occur.
!                 For I = 1,...,NZ, JA(I) is the row in which the Ith
!                 element of the Jacobian occurs. JA must also include
!                 the diagonal elements of the Jacobian.
!                 PD defines the numerical value of the Jacobian
!                 elements. For I = 1,...,NZ, PD(I) is the numerical
!                 value of the Ith element in the Jacobian. PD must
!                 also include the diagonal elements of the Jacobian.
! GFUN    = the name of the subroutine to evaluate the residuals for
!           event functions. If you do not specify that events are
!           present (by specifying NEVENTS > 0 in SET_OPTS), you
!           need not include GFUN in your call list for DVODE_F90.
!           If GFUN is not contained in a module available to your
!           calling program, you must declare GFUN to be EXTERNAL
!           in your calling program.
! To continue the integration after a successful return, simply
! reset TOUT and call DVODE_F90 again. No other parameters need
! be reset unless ISTATE=3 in which case, reset it to 2 before
! calling DVODE_F90 again.
! ______________________________________________________________________
! Section 3.  Choosing Error Tolerances
!
! This is the most important aspect of solving odes numerically.
! You may supply any of four keywords and values. If you wish to
! use scalar error tolerances, you should supply ABSERR and RELERR.
! For a good many problems, it is advisable to supply a vector of
! absolute error tolerances ABSERR_VECTOR = desired vector. This
! allows you to use different tolerances for each component of
! the solution. If ABSERR_VECTOR is supplied, it must be a vector
! of length NEQ where NEQ is the number of odes in your system.
! Similarly, you may supply a vector of relative error tolerances,
! RELERR_VECTOR. If no tolerances are specified, DVODE_F90 will use
! default error tolerances ABSERR=1D-6 and RELERR=1D-4; but it is
! strongly recommended that you supply values that are appropriate
! for your problem. In the event you do not supply error tolerances,
! DVODE_F90 will print a reminder that the default error tolerances
! are not appropriate for all problems.
!
! RELERR can be set to a scalar value as follows.
! Determine the number of significant digits of accuracy desired,
! which will be a positive integer, say, N.
! Then set RELERR = 10**-(N+1).
! The authors recommend that RELERR be no larger than 10**-4.
! The authors recommend a vector valued absolute error tolerance,
! which can be set as follows.
! For the I-th component of the solution vector, Y(I), determine
! the positive number FLOOR(I) at which ABS(Y(I)) becomes
! negligible for the problem at hand. FLOOR(I) is sometimes called
! the problem zero or the floor value for the I-th component and is
! problem dependent. For a given problem that is not scaled, these
! floor values may well vary by up to 9 orders of magnitude.
! Set ABSERR(I) = FLOOR(I) or to be conservative
! ABSERR_VECTOR(I) = 0.1*FLOOR(I). There is no variable FLOOR in
! DVODE_F90. If it is difficult to divine the components of ABSERR,
! (or FLOOR) make a reasonable guess, run the problem, then set that
! ABSERR_VECTOR so for I = 1, 2,...NEQ,
! ABSERR_VECTOR(I) = 1D-6*RELERR*MAX{ABS(Y(I,TOUT): for all TOUT}.
! The correct choices for RELERR and ABSERR can and do have
! significant impact on both the quality of the solution and run
! time. Counter intuitively, error tolerances that are too loose
! can and do increase run time significantly and the quality of
! the solution may be seriously compromised.
! Examples:
! 1. OPTIONS = SET_OPTS(DENSE_J=.TRUE., ABSERR=1D-8,RELERR=1D-8)
!    This will yield MF = 22. Both the relative error tolerance
!    and the absolute error tolerance will equal 1D-8.
! 2. OPTIONS = SET_OPTS(DENSE_J=.TRUE.,RELERR=1D-5, &
!                       ABSERR_VECTOR=(/1D-6,1D-8/))
!    For a system with NEQ=2 odes, this will yield MF = 22. A scalar
!    relative error tolerance equal to 1D-5 will be used. Component 1
!    of the solution will use an absolute error tolerance of 1D-6
!    while component 2 will use an absolute error tolerance of 1D-8.
! ______________________________________________________________________
! Section 4.  Choosing the Method of Integration
!
! If you wish to specify a numerical value for METHOD_FLAG, it can equal
! any of the legal values of MF for DVODE or LSODES. (See below.) If
! you do not wish to specify a numerical value for METHOD_FLAG, you
! may specify any combination of the five logical keywords DENSE_J,
! BANDED_J, SPARSE_J, USER_SUPPLIED_JACOBIAN, SAVE_JACOBIAN that you
! wish. Appropriate values will be used in DVODE_F90 for any variables
! that are not present. The first three flags indicate the type of
! Jacobian, dense, banded, or sparse. If USER_SUPPLIED_JACOBIAN=.TRUE.,
! the Jacobian will be determined using the subroutine JAC you supply
! in your call to DVODE_F90. Otherwise, an internal Jacobian will be
! generated using finite differences.
! Examples:
! 1. OPTIONS = SET_OPTS(METHOD_FLAG=22,...)
!    DVODE will use the value MF=22 as described in the documentation
!    prologue. In this case, the stiff BDF methods will be used with
!    a dense, internally generated Jacobian.
! 2. OPTIONS = SET_OPTS(METHOD_FLAG=21,...)
!    This is the same an Example 1 except that a user supplied dense
!    Jacobian will be used and DVODE will use MF=21.
! 3. OPTIONS = SET_OPTS(DENSE_J=.TRUE.,...)
!    This will yield MF = 22 as in Example 1, provided
!    USER_SUPPLIED_JACOBIAN and SAVE_JACOBIAN are not present, or if
!    present are set to their default
!     values of .FALSE. and .TRUE., respectively.
! 4. OPTIONS = SET_OPTS(DENSE_J=.TRUE.,&
!                       USER_SUPPLIED_JACOBIAN=.TRUE.,...)
!    This will yield MF = 21 as in Example 2, provided SAVE_JACOBIAN
!    is not present, or if present is set to its default value .FALSE.
! Notes:
! 1. If you specify more than one of DENSE_J, BANDED_J, and SPARSE_J,
!    DENSE_J takes precedence over BANDED_J which in turn takes
!    precedence over SPARSE_J.
! 2. By default, DVODE_F90 saves a copy of the Jacobian and reuses the
!    copy when necessary. For problems in which storage requirements
!    are acute, you may wish to override this default and have
!    DVODE_F90 recalculate the Jacobian rather than use a saved copy.
!    You can do this by specifying SAVE_JACOBIAN=.FALSE. It is
!    recommended that you not do this unless necessary since it can
!    have a significant impact on the efficiency of DVODE_F90. (For
!    example, when solving a linear problem only one evaluation of
!    the Jacobian is required with the default option.)
! 3. If you choose BANDED_J = .TRUE. or if you supply a value of MF
!    that corresponds to a banded Jacobian, you must also supply the
!    lower  bandwidth ML and the upper bandwidth of the Jacobian MU
!    by including the keywords
!    LOWER_BANDWIDTH = value of ML and UPPER_BANDWIDTH = value of M
!                   More on Method Selection
! The keyword options available in SET_OPTS are intended to replace
! the original method indicator flag MF. However, if you wish to
! retain the flexibility of the original solver, you may specify MF
! directly in your call to SET_OPTS. This is done by using the
! keyword METHOD_FLAG=MF in your SET_OPTS where MF has any of the
! values in the following description. Refer to the demonstration
! program demosp.f90 for an example in which this is done.
! MF     = The method flag. Used only for input. The legal values of
!          MF are:
!          10, 11, 12, 13, 14, 15, 16, 17, 20, 21, 22, 23, 24, 25, 26,
!          27, -11, -12, -14, -15, -21, -22, -24, -25, -26, -27.
!          MF is a signed two-digit integer, MF = JSV*(10*METH+MITER).
!          JSV = SIGN(MF) indicates the Jacobian-saving strategy:
!            JSV =  1 means a copy of the Jacobian is saved for reuse
!                     in the corrector iteration algorithm.
!            JSV = -1 means a copy of the Jacobian is not saved
!                     (valid only for MITER = 1, 2, 4, or 5).
!          METH indicates the basic linear multistep method:
!            METH = 1 means the implicit Adams method.
!            METH = 2 means the method based on backward
!                     differentiation formulas (BDF-s).
!          MITER indicates the corrector iteration method:
!            MITER = 0 means functional iteration (no Jacobian matrix
!                      is involved).
!            MITER = 1 means chord iteration with a user-supplied
!                      full (NEQ by NEQ) Jacobian.
!            MITER = 2 means chord iteration with an internally
!                      generated (difference quotient) full Jacobian
!                      (using NEQ extra calls to F per df/dy value).
!            MITER = 4 means chord iteration with a user-supplied
!                      banded Jacobian.
!            MITER = 5 means chord iteration with an internally
!                      generated banded Jacobian (using ML+MU+1 extra
!                      calls to F per df/dy evaluation).
!            MITER = 6 means chord iteration with a user-supplied
!                      sparse Jacobian.
!            MITER = 7 means chord iteration with an internally
!                      generated sparse Jacobian
!          If MITER = 1, 4, or 6 the user must supply a subroutine
!          JAC(the name is arbitrary) as described above under JAC.
!          For other values of MITER, JAC need not be provided.
! ______________________________________________________________________
! Section 5.  Interpolation of the Solution and Derivative
!
! Following a successful return from DVODE_F90, you may call
! subroutine DVINDY to interpolate the solution or derivative.
! SUBROUTINE DVINDY(T, K, DKY, IFLAG)
! DVINDY computes interpolated values of the K-th derivative of the
! dependent variable vector y, and stores it in DKY. This routine
! is called with K = 0 or K = 1 and T = TOUT. In either case, the
! results are returned in the array DKY of length at least NEQ which
! must be declared and dimensioned in your calling program. The
! computed values in DKY are obtained by interpolation using the
! Nordsieck history array.
! ______________________________________________________________________
! Section 6.  Handling Events (Root Finding)
!
!    DVODE_F90 contains root finding provisions. It attempts to
!    locates the roots of a set of functions
!         g(i) = g(i,t,y(1),...,y(NEQ))  (i = 1,...,ng).
!    To use root finding include NEVENTS=NG in your call to SET_OPTS
!    where NG is the number of root finding functions. You must then
!    supply subroutine GFUN in your call to DVODE_F90 using
!    G_FCN=GFUN as the last argument. GFUN must has the form
!               SUBROUTINE GFUN(NEQ, T, Y, NG, GOUT)
!    where NEQ, T, Y, and NG are input, and the array GOUT is output.
!    NEQ, T, and Y have the same meaning as in the F routine, and
!    GOUT is an array of length NG. For i = 1,...,NG, this routine is
!    to load into GOUT(i) the value at (T,Y) of the i-th constraint
!    function g(i). DVODE_F90 will find roots of the g(i) of odd
!    multiplicity (i.e. sign changes) as they occur during the
!    integration. GFUN must be declared EXTERNAL in the calling
!    program. Note that because of numerical errors in the functions
!    g(i) due to roundoff and integration error, DVODE_F90 may return
!    false roots, or return the same root at two or more nearly equal
!    values of t. This is particularly true for problems in which the
!    integration is restarted (ISTATE = 1) at a root. If such false
!    roots are suspected, you should consider smaller error tolerances
!    and/or higher precision in the evaluation of the g(i). Note
!    further that if a root of some g(i) defines the end of the
!    problem, the input to DVODE_F90 should nevertheless allow
!    integration to a point slightly past that root, so that DVODE_F90
!    can locate the root by interpolation. Each time DVODE_F90 locates
!    a root of one of your event functions it makes a return to the
!    calling program with ISTATE = 3. When such a return is made and
!    you have processed the results, simply change ISTATE = 2 and call
!    DVODE_F90 again without making other changes.
! ______________________________________________________________________
! Section 7.  Gathering Integration Statistics
!
! SUBROUTINE GET_STATS(RSTATS, ISTATS, NUMEVENTS, JROOTS)
! Caution:
! RSTATS and ISTATS must be declared and dimensioned in your
! main program. The minimum dimensions are:
! DIMENSION RSTATS(22), ISTATS(31)
! This subroutine returns the user portions of the original DVODE
! RUSER and IUSER arrays, and if root finding is being done, it
! returns the original LSODAR JROOT vector. NUMEVENTS and JROOTS
! are optional parameters. NUMEVENTS is the number of root functions
! and JROOTS is an integer array of length NUMEVENTS.
! Available Integration Statistics:
! HU      RUSER(11) The step size in t last used (successfully).
! HCUR    RUSER(12) The step size to be attempted on the next step.
! TCUR    RUSER(13) The current value of the independent variable
!                   which the solver has actually reached, i.e. the
!                   current internal mesh point in t. In the output,
!                   TCUR will always be at least as far from the
!                   initial value of t as the current argument T,
!                   but may be farther (if interpolation was done).
! TOLSF   RUSER(14) A tolerance scale factor, greater than 1.0,
!                   computed when a request for too much accuracy was
!                   detected (ISTATE = -3 if detected at the start of
!                   the problem, ISTATE = -2 otherwise). If ITOL is
!                   left unaltered but RTOL and ATOL are uniformly
!                   scaled up by a factor of TOLSF for the next call,
!                   then the solver is deemed likely to succeed.
!                   (The user may also ignore TOLSF and alter the
!                   tolerance parameters in any other way appropriate.)
! NST     IUSER(11) The number of steps taken for the problem so far.
! NFE     IUSER(12) The number of f evaluations for the problem so far.
! NJE     IUSER(13) The number of Jacobian evaluations so far.
! NQU     IUSER(14) The method order last used (successfully).
! NQCUR   IUSER(15) The order to be attempted on the next step.
! IMXER   IUSER(16) The index of the component of largest magnitude in
!                   the weighted local error vector (E(i)/EWT(i)),
!                   on an error return with ISTATE = -4 or -5.
! LENRW   IUSER(17) The length of RUSER actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
! LENIW   IUSER(18) The length of IUSER actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
! NLU     IUSER(19) The number of matrix LU decompositions so far.
! NNI     IUSER(20) The number of nonlinear (Newton) iterations so far.
! NCFN    IUSER(21) The number of convergence failures of the nonlinear
!                   solver so far.
! NETF    IUSER(22) The number of error test failures of the integrator
!                   so far.
! MA28AD_CALLS      IUSER(23) The number of calls made to MA28AD
! MA28BD_CALLS      IUSER(24) The number of calls made to MA28BD
! MA28CD_CALLS      IUSER(25) The number of calls made to MA28CD
! MC19AD_CALLS      IUSER(26) The number of calls made to MC19AD
! IRNCP             IUSER(27) The number of compressions done on array JAN
! ICNCP             IUSER(28) The number of compressions done on array ICN
! MINIRN            IUSER(29) Minimum size for JAN array
! MINICN            IUSER(30) Minimum size for ICN array
! MINNZ             IUSER(31) Number of nonzeros in sparse Jacobian
! JROOTS  JROOTS    Optional array of component indices for components
!                   having a zero at the current time
! ______________________________________________________________________
! Section 8.  Determining Jacobian Sparsity Structure Arrays
!
! If you are solving a problem with a sparse Jacobian, the arrays
! that define the sparsity structure are needed. The arrays may
! be determined in any of several ways.
! 1. If you choose the default mode by indicating SPARSE=.TRUE.,
!    the sparsity arrays will be determined internally by DVODE_F90
!    by making calls to your derivative subroutine. This mode is
!    equivalent to using the integration method flag MF = 227.
! 2. The DVODE_F90 method flag MF is defined to be
!    MF = 100*MOSS + 10*METH + MITER. If you supply MF = 227 (or 217),
!    the sparse Jacobian will be determined using finite differences;
!    and the sparsity arrays will be determined by calling your
!    derivative subroutine.
! 3. If you supply MF = 126 (or 116), you must supply the Jacobian
!    subroutine JAC to define the exact Jacobian. JAC must have the
!    following form:
!           SUBROUTINE JAC (N, T, Y, IA, JA, NZ, PD)
!    Given the number of odes N, the current time T, and the current
!    solution vector Y, JAC must do the following:
!    -  If NZ = 0 on input:
!       Replace NZ by the number of nonzero elements in the Jacobian.
!       The diagonal of the Jacobian must be included.
!       Do NOT define the arrays IA, JA, PD at this time.
!       Once JAC has been called with NZ = 0 and you have defined the
!       value of NZ, future calls to JAC will use this value of NZ.
!    -  When a call is made with NZ unequal to 0, you must define the
!       sparsity structure arrays IA and JA, and the sparse Jacobian
!       PD.
!         - IA defines the number of nonzeros including the diagonal
!           in each column of the Jacobian. Define IA(1) = 1 and for
!           J = 1,..., N,
!           IA(J+1) = IA(J) + number of nonzeros in column J.
!           Diagonal elements must be include even if they are zero.
!           You should check to ensure that IA(N+1)-1 = NZ.
!         - JA defines the rows in which the nonzeros occur. For
!           I = 1,...,NZ, JA(I) is the row in which the Ith element
!           of the Jacobian occurs. JA must also include the diagonal
!           elements of the Jacobian.
!         - PD defines the numerical value of the Jacobian elements.
!           For I = 1,...,NZ, PD(I) is the numerical value of the
!           Ith element in the Jacobian. PD must also include the
!           diagonal elements of the Jacobian.
! 4. If you wish to supply the IA and JA arrays directly, use
!    MF = 27. In this case, after calling SET_OPTS, you must call
!    SET_IAJA supplying the arrays IAUSER and JAUSER described in
!    the documentation prologue for SET_IAJA. These arrays will be
!    used when approximate Jacobians are determined using finite
!    differences.
! There are two user callable sparsity structure subroutines:
! USERSETS_IAJA may be used if you wish to supply the sparsity
! structure directly.
! SUBROUTINE USERSETS_IAJA(IAUSER,NIAUSER,JAUSER,NJAUSER)
!     Caution:
!     If it is called, USERSETS_IAJA must be called after the
!     call to SET_OPTS.
!     Usage:
!     CALL SET_IAJA(IAUSER,NIAUSER,JAUSER,NJAUSER)
!       In this case, IAUSER of length NIAUSER will be used for
!       IA; and JAUSER of length NJAUSER will be used for JA.
!     Arguments:
!     IAUSER  = user supplied IA array
!     NIAUSER = length of IAUSER array
!     JAUSER  = user supplied JA vector
!     NJAUSER = length of JAUSER array
! The second subroutine allows you to approximate the sparsity
! structure using derivative differences. It allows more flexibility
! in the determination of perturbation increments used.
! SUBROUTINE SET_IAJA(DFN,NEQ,T,Y,FMIN,NTURB,DTURB,IAUSER, &
!   NIAUSER, JAUSER, NJAUSER)
!     Caution:
!     If it is called, SET_IAJA must be called after the call to
!     SET_OPTS.
!     Usage:
!     SET_IAJA may be called in one of two ways:

!     CALL SET_IAJA(DFN,NEQ,T,Y,FMIN,NTURB,DTURB)
!       In this case IA and JA will be determined using calls
!       to your derivative routine DFN.
!     CALL SET_IAJA(DFN,NEQ,T,Y,FMIN,NTURB,DTURB,IAUSER,NIAUSER, &
!       JAUSER, NJAUSER)
!       In this case, IAUSER of length NIAUSER will be used for
!       IA; and JAUSER of length NJAUSER will be used for JA.
!       T, Y, FMIN, NTURB, and DTURB will be ignored (though
!       they must be present in the argument list).
!     Arguments:
!     DFN     = DVODE derivative subroutine
!     NEQ     = Number of odes
!     T       = independent variable t
!     Y       = solution y(t)
!     FMIN    = Jacobian threshold value. Elements of the Jacobian
!               with magnitude smaller than FMIN will be ignored.
!               FMIN will be ignored if it is less than or equal
!               to ZERO.
!     NTURB   = Perturbation flag. If NTURB=1, component I of Y
!               will be perturbed by 1.01D0.
!               If NTURB=NEQ, component I of Y will be perturbed
!               by ONE + DTURB(I).
!     DTURB   = perturbation vector of length 1 or NEQ.
!     If these four optional parameters are present, IAUSER and JAUSER
!     will be copied to IA and JA rather than making derivative calls
!     to approximate IA and JA:
!        IAUSER  = user supplied IA array
!        NIAUSER = length of IAUSER array
!        JAUSER  = user supplied JA vector
!        NJAUSER = length of JAUSER array
! ______________________________________________________________________
! Section 9.  Original DVODE.F Documentation Prologue
!
! SUBROUTINE DVODE(F, NEQ, Y, T, TOUT, ITASK, ISTATE, OPTS, JAC, GFUN)
! DVODE: Variable-coefficient Ordinary Differential Equation solver,
! with fixed-leading-coefficient implementation.
! Note:
! Numerous changes have been made in the documentation and the code
! from the original Fortran 77 DVODE solver. With regard to the new
! F90 version, if you choose options that correspond to options
! available in the original f77 solver, you should obtain the same
! results. In all testing, identical results have been obtained
! between this version and a simple F90 translation of the original
! solver.
! DVODE solves the initial value problem for stiff or nonstiff
! systems of first order ODEs,
!     dy/dt = f(t,y), or, in component form,
!     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
! DVODE is a package based on the EPISODE and EPISODEB packages, and
! on the ODEPACK user interface standard, with minor modifications.
! This version is based also on elements of LSODES and LSODAR.
! Authors:
!               Peter N. Brown and Alan C. Hindmarsh
!               Center for Applied Scientific Computing, L-561
!               Lawrence Livermore National Laboratory
!               Livermore, CA 94551
!               George D. Byrne
!               Illinois Institute of Technology
!               Chicago, IL 60616
! References:
! 1. P. N. Brown, G. D. Byrne, and A. C. Hindmarsh, "VODE: A Variable
!    Coefficient ODE Solver," SIAM J. Sci. Stat. Comput., 10 (1989),
!    pp. 1038-1051. Also, LLNL Report UCRL-98412, June 1988.
! 2. G. D. Byrne and A. C. Hindmarsh, "A Polyalgorithm for the
!    Numerical Solution of Ordinary Differential Equations,"
!    ACM Trans. Math. Software, 1 (1975), pp. 71-96.
! 3. A. C. Hindmarsh and G. D. Byrne, "EPISODE: An Effective Package
!    for the Integration of Systems of Ordinary Differential
!    Equations," LLNL Report UCID/30112, Rev. 1, April 1977.
! 4. G. D. Byrne and A. C. Hindmarsh, "EPISODEB: An Experimental
!    Package for the Integration of Systems of Ordinary Differential
!    Equations with Banded Jacobians," LLNL Report UCID/30132, April
!    1976.
! 5. A. C. Hindmarsh, "ODEPACK, a Systematized Collection of ODE
!    Solvers," in Scientific Computing, R. S. Stepleman et al., eds.,
!    North-Holland, Amsterdam, 1983, pp. 55-64.
! 6. K. R. Jackson and R. Sacks-Davis, "An Alternative Implementation
!    of Variable Step-Size Multistep Formulas for Stiff ODEs," ACM
!    Trans. Math. Software, 6 (1980), pp. 295-318.
!                     Summary of Usage
! Communication between the user and the DVODE package, for normal
! situations, is summarized here. This summary describes only a subset
! of the full set of options available. See the full description for
! details, including optional communication, nonstandard options,
! and instructions for special situations. See also the example
! problem (with program and output) following this summary.
! A. First provide a subroutine of the form:
!           SUBROUTINE F(NEQ, T, Y, YDOT)
!           REAL(KIND=WP) T, Y(NEQ), YDOT(NEQ)
! which supplies the vector function f by loading YDOT(i) with f(i).
! B. Next determine (or guess) whether or not the problem is stiff.
! Stiffness occurs when the Jacobian matrix df/dy has an eigenvalue
! whose real part is negative and large in magnitude, compared to the
! reciprocal of the t span of interest. If the problem is nonstiff,
! use a method flag MF = 10. If it is stiff, there are four standard
! choices for MF(21, 22, 24, 25), and DVODE requires the Jacobian
! matrix in some form. In these cases (MF > 0), DVODE will use a
! saved copy of the Jacobian matrix. If this is undesirable because of
! storage limitations, set MF to the corresponding negative value
! (-21, -22, -24, -25). (See full description of MF below.)
! The Jacobian matrix is regarded either as full (MF = 21 or 22),
! or banded (MF = 24 or 25). In the banded case, DVODE requires two
! half-bandwidth parameters ML and MU. These are, respectively, the
! widths of the lower and upper parts of the band, excluding the main
! diagonal. Thus the band consists of the locations (i,j) with
! i-ML <= j <= i+MU, and the full bandwidth is ML+MU+1.
! C. If the problem is stiff, you are encouraged to supply the Jacobian
! directly (MF = 21 or 24), but if this is not feasible, DVODE will
! compute it internally by difference quotients (MF = 22 or 25).
! If you are supplying the Jacobian, provide a subroutine of the form:
!           SUBROUTINE JAC(NEQ, T, Y, ML, MU, PD, NROWPD)
!           REAL(KIND=WP) T, Y(NEQ), PD(NROWPD,NEQ)
! which supplies df/dy by loading PD as follows:
!     For a full Jacobian (MF = 21), load PD(i,j) with df(i)/dy(j),
! the partial derivative of f(i) with respect to y(j). (Ignore the
! ML and MU arguments in this case.)
!     For a banded Jacobian (MF = 24), load PD(i-j+MU+1,j) with
! df(i)/dy(j), i.e. load the diagonal lines of df/dy into the rows of
! PD from the top down.
!     In either case, only nonzero elements need be loaded.
! D. Write a main program which calls subroutine DVODE once for
! each point at which answers are desired. This should also provide
! for possible use of logical unit 6 for output of error messages
! by DVODE. On the first call to DVODE, supply arguments as follows:
! F      = Name of subroutine for right-hand side vector f.
!          This name must be declared external in calling program.
! NEQ    = Number of first order ODEs.
! Y      = Array of initial values, of length NEQ.
! T      = The initial value of the independent variable.
! TOUT   = First point where output is desired (/= T).
! ITOL   = 1 or 2 according as ATOL(below) is a scalar or array.
! RTOL   = Relative tolerance parameter (scalar).
! ATOL   = Absolute tolerance parameter (scalar or array).
!          The estimated local error in Y(i) will be controlled so as
!          to be roughly less (in magnitude) than
!             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
!             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
!          Thus the local error test passes if, in each component,
!          either the absolute error is less than ATOL(or ATOL(i)),
!          or the relative error is less than RTOL.
!          Use RTOL = 0.0 for pure absolute error control, and
!          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
!          control. Caution: Actual (global) errors may exceed these
!          local tolerances, so choose them conservatively.
! ITASK  = 1 for normal computation of output values of Y at t = TOUT.
! ISTATE = Integer flag (input and output). Set ISTATE = 1.
! IOPT   = 0 to indicate no optional input used.
! JAC    = Name of subroutine for Jacobian matrix (MF = 21 or 24).
!          If used, this name must be declared external in calling
!          program. If not used, pass a dummy name.
! MF     = Method flag. Standard values are:
!          10 for nonstiff (Adams) method, no Jacobian used.
!          21 for stiff (BDF) method, user-supplied full Jacobian.
!          22 for stiff method, internally generated full Jacobian.
!          24 for stiff method, user-supplied banded Jacobian.
!          25 for stiff method, internally generated banded Jacobian.
! E. The output from the first call (or any call) is:
!      Y = Array of computed values of y(t) vector.
!      T = Corresponding value of independent variable (normally TOUT).
! ISTATE = 2  if DVODE was successful, negative otherwise.
!          -1 means excess work done on this call. (Perhaps wrong MF.)
!          -2 means excess accuracy requested. (Tolerances too small.)
!          -3 means illegal input detected. (See printed message.)
!          -4 means repeated error test failures. (Check all input.)
!          -5 means repeated convergence failures. (Perhaps bad
!             Jacobian supplied or wrong choice of MF or tolerances.)
!          -6 means error weight became zero during problem. (Solution
!             component I vanished, and ATOL or ATOL(I) = 0.)
! F. To continue the integration after a successful return, simply
! reset TOUT and call DVODE again. No other parameters need be reset.
!         Full Description of User Interface to DVODE
! The user interface to DVODE consists of the following parts.
! i.  The call sequence to subroutine DVODE, which is a driver
!      routine for the solver. This includes descriptions of both
!      the call sequence arguments and of user-supplied routines.
!      Following these descriptions is
!        * a description of optional input available through the
!          call sequence,
!        * a description of optional output (in the work arrays), and
!        * instructions for interrupting and restarting a solution.
! ii. Descriptions of other routines in the DVODE package that may be
!      (optionally) called by the user. These provide the ability to
!      alter error message handling, save and restore the internal
!      PRIVATE variables, and obtain specified derivatives of the
!      solution y(t).
! iii. Descriptions of PRIVATE variables to be declared in overlay
!      or similar environments.
! iv. Description of two routines in the DVODE package, either of
!      which the user may replace with his own version, if desired.
!      these relate to the measurement of errors.
! Part i. Call Sequence.
! The call sequence parameters used for input only are
!     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, JAC, MF,
! and those used for both input and output are
!     Y, T, ISTATE.
! The work arrays RUSER and IUSER are used for conditional and
! optional input and optional output. (The term output here refers
! to the return from subroutine DVODE to the user's calling program.)
! The legality of input parameters will be thoroughly checked on the
! initial call for the problem, but not checked thereafter unless a
! change in input parameters is flagged by ISTATE = 3 in the input.
! The descriptions of the call arguments are as follows.
! F      = The name of the user-supplied subroutine defining the
!          ODE system. The system must be put in the first-order
!          form dy/dt = f(t,y), where f is a vector-valued function
!          of the scalar t and the vector y. Subroutine F is to
!          compute the function f. It is to have the form
!               SUBROUTINE F(NEQ, T, Y, YDOT)
!               REAL(KIND=WP) T, Y(NEQ), YDOT(NEQ)
!          where NEQ, T, and Y are input, and the array YDOT = f(t,y)
!          is output. Y and YDOT are arrays of length NEQ.
!          Subroutine F should not alter Y(1),...,Y(NEQ).
!          F must be declared EXTERNAL in the calling program.

!          If quantities computed in the F routine are needed
!          externally to DVODE, an extra call to F should be made
!          for this purpose, for consistent and accurate results.
!          If only the derivative dy/dt is needed, use DVINDY instead.
! NEQ    = The size of the ODE system (number of first order
!          ordinary differential equations). Used only for input.
!          NEQ may not be increased during the problem, but
!          can be decreased (with ISTATE = 3 in the input).
! Y      = A real array for the vector of dependent variables, of
!          length NEQ or more. Used for both input and output on the
!          first call (ISTATE = 1), and only for output on other calls.
!          On the first call, Y must contain the vector of initial
!          values. In the output, Y contains the computed solution
!          evaluated at T. If desired, the Y array may be used
!          for other purposes between calls to the solver.
!          This array is passed as the Y argument in all calls to
!          F and JAC.
! T      = The independent variable. In the input, T is used only on
!          the first call, as the initial point of the integration.
!          In the output, after each call, T is the value at which a
!          computed solution Y is evaluated (usually the same as TOUT).
!          On an error return, T is the farthest point reached.
! TOUT   = The next value of t at which a computed solution is desired.
!          Used only for input.
!          When starting the problem (ISTATE = 1), TOUT may be equal
!          to T for one call, then should /= T for the next call.
!          For the initial T, an input value of TOUT /= T is used
!          in order to determine the direction of the integration
!          (i.e. the algebraic sign of the step sizes) and the rough
!          scale of the problem. Integration in either direction
!          (forward or backward in t) is permitted.
!          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
!          the first call (i.e. the first call with TOUT /= T).
!          Otherwise, TOUT is required on every call.
!          If ITASK = 1, 3, or 4, the values of TOUT need not be
!          monotone, but a value of TOUT which backs up is limited
!          to the current internal t interval, whose endpoints are
!          TCUR - HU and TCUR. (See optional output, below, for
!          TCUR and HU.)
! ITOL   = An indicator for the type of error control. See
!          description below under ATOL. Used only for input.
! RTOL   = A relative error tolerance parameter, either a scalar or
!          an array of length NEQ. See description below under ATOL.
!          Input only.
! ATOL   = An absolute error tolerance parameter, either a scalar or
!          an array of length NEQ. Input only.
!          The input parameters ITOL, RTOL, and ATOL determine
!          the error control performed by the solver. The solver will
!          control the vector e = (e(i)) of estimated local errors
!          in Y, according to an inequality of the form
!                      rms-norm of (E(i)/EWT(i)) <= 1,
!          where       EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
!          and the rms-norm (root-mean-square norm) here is
!          rms-norm(v) = sqrt(sum v(i)**2 / NEQ). Here EWT
!          is a vector of weights which must always be positive, and
!          the values of RTOL and ATOL should all be nonnegative.
!          The following table gives the types (scalar/array) of
!          RTOL and ATOL, and the corresponding form of EWT(i).
!           ITOL    RTOL       ATOL          EWT(i)
!            1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
!            2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
!            3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
!            4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
!          When either of these parameters is a scalar, it need not
!          be dimensioned in the user's calling program.
!          If none of the above choices (with ITOL, RTOL, and ATOL
!          fixed throughout the problem) is suitable, more general
!          error controls can be obtained by substituting
!          user-supplied routines for the setting of EWT and/or for
!          the norm calculation. See Part iv below.
!          If global errors are to be estimated by making a repeated
!          run on the same problem with smaller tolerances, then all
!          components of RTOL and ATOL(i.e. of EWT) should be scaled
!          down uniformly.
! ITASK  = An index specifying the task to be performed.
!          Input only. ITASK has the following values and meanings.
!          1  means normal computation of output values of y(t) at
!             t = TOUT(by overshooting and interpolating).
!          2  means take one step only and return.
!          3  means stop at the first internal mesh point at or
!             beyond t = TOUT and return.
!          4  means normal computation of output values of y(t) at
!             t = TOUT but without overshooting t = TCRIT.
!             TCRIT must be input as RUSER(1). TCRIT may be equal to
!             or beyond TOUT, but not behind it in the direction of
!             integration. This option is useful if the problem
!             has a singularity at or beyond t = TCRIT.
!          5  means take one step, without passing TCRIT, and return.
!             TCRIT must be input as RUSER(1).
!          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
!          (within roundoff), it will return T = TCRIT(exactly) to
!          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
!          in which case answers at T = TOUT are returned first).
! ISTATE = an index used for input and output to specify the
!          the state of the calculation.
!          In the input, the values of ISTATE are as follows.
!          1  means this is the first call for the problem
!             (initializations will be done). See note below.
!          2  means this is not the first call, and the calculation
!             is to continue normally, with no change in any input
!             parameters except possibly TOUT and ITASK.
!             (If ITOL, RTOL, and/or ATOL are changed between calls
!             with ISTATE = 2, the new values will be used but not
!             tested for legality.)
!          3  means this is not the first call, and the
!             calculation is to continue normally, but with
!             a change in input parameters other than
!             TOUT and ITASK. Changes are allowed in
!             NEQ, ITOL, RTOL, ATOL, IOPT, MF, ML, MU,
!             and any of the optional input except H0.
!             (See IUSER description for ML and MU.)
!          Caution:
!          If you make a call to DVODE_F90 with ISTATE=3, you will
!          first need to call SET_OPTS again, supplying the new
!          necessary option values.
!          Note:  A preliminary call with TOUT = T is not counted
!          as a first call here, as no initialization or checking of
!          input is done. (Such a call is sometimes useful to include
!          the initial conditions in the output.)
!          Thus the first call for which TOUT /= T requires
!          ISTATE = 1 in the input.
!          In the output, ISTATE has the following values and meanings.
!           1  means nothing was done, as TOUT was equal to T with
!              ISTATE = 1 in the input.
!           2  means the integration was performed successfully.
!          -1  means an excessive amount of work (more than MXSTEP
!              steps) was done on this call, before completing the
!              requested task, but the integration was otherwise
!              successful as far as T. (MXSTEP is an optional input
!              and is normally 5000.)  To continue, the user may
!              simply reset ISTATE to a value > 1 and call again.
!              (The excess work step counter will be reset to 0.)
!              In addition, the user may increase MXSTEP to avoid
!              this error return. (See optional input below.)
!          -2  means too much accuracy was requested for the precision
!              of the machine being used. This was detected before
!              completing the requested task, but the integration
!              was successful as far as T. To continue, the tolerance
!              parameters must be reset, and ISTATE must be set
!              to 3. The optional output TOLSF may be used for this
!              purpose. (Note: If this condition is detected before
!              taking any steps, then an illegal input return
!              (ISTATE = -3) occurs instead.)
!          -3  means illegal input was detected, before taking any
!              integration steps. See written message for details.
!              Note:  If the solver detects an infinite loop of calls
!              to the solver with illegal input, it will cause
!              the run to stop.
!          -4  means there were repeated error test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as T.
!              The problem may have a singularity, or the input
!              may be inappropriate.
!          -5  means there were repeated convergence test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as T.
!              This may be caused by an inaccurate Jacobian matrix,
!              if one is being used.
!          -6  means EWT(i) became zero for some i during the
!              integration. Pure relative error control (ATOL(i)=0.0)
!              was requested on a variable which has now vanished.
!              The integration was successful as far as T.
!          Note:  Since the normal output value of ISTATE is 2,
!          it does not need to be reset for normal continuation.
!          Also, since a negative input value of ISTATE will be
!          regarded as illegal, a negative output value requires the
!          user to change it, and possibly other input, before
!          calling the solver again.
! IOPT   = An integer flag to specify whether or not any optional
!          input is being used on this call. Input only.
!          The optional input is listed separately below.
!          IOPT = 0 means no optional input is being used.
!                   Default values will be used in all cases.
!          IOPT = 1 means optional input is being used.
! RUSER  = A real working array (real(wp)).
!          The length of RUSER must be at least 22
!             20 + NYH * (MAXORD + 1) where
!          NYH    = the initial value of NEQ,
!          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
!                   smaller value is given as an optional input),
!          The first 22 words of RUSER are reserved for conditional
!          and optional input and optional output.

!          The following word in RUSER is a conditional input:
!            RUSER(1) = TCRIT = critical value of t which the solver
!                       is not to overshoot. Required if ITASK is
!                       4 or 5, and ignored otherwise. (See ITASK.)
! IUSER  = An integer work array. The length of IUSER must be at least 30.
!             30        if MITER = 0 or 3 (MF = 10, 13, 20, 23), or
!             30 + NEQ  otherwise (ABS(MF) = 11,12,14,15,16,17,21,22,
!             24,25,26,27).
!          The first 30 words of IUSER are reserved for conditional and
!          optional input and optional output.

!          The following 2 words in IUSER are conditional input:
!            IUSER(1) = ML  These are the lower and upper
!            IUSER(2) = MU  half-bandwidths, respectively, of the
!                       banded Jacobian, excluding the main diagonal.
!                       The band is defined by the matrix locations
!                       (i,j) with i-ML <= j <= i+MU. ML and MU
!                       must satisfy  0 <= ML,MU  <= NEQ-1.
!                       These are required if MITER is 4 or 5, and
!                       ignored otherwise. ML and MU may in fact be
!                       the band parameters for a matrix to which
!                       df/dy is only approximately equal.
! Note:  The work arrays must not be altered between calls to DVODE
! for the same problem, except possibly for the conditional and
! optional input, and except for the last 3*NEQ words of RUSER.
! The latter space is used for internal scratch space, and so is
! available for use by the user outside DVODE between calls, if
! desired (but not for use by F or JAC).
! JAC    = The name of the user-supplied routine (MITER = 1 or 4 or 6)
!          to compute the Jacobian matrix, df/dy, as a function of
!          the scalar t and the vector y. It is to have the form
!               SUBROUTINE JAC(NEQ, T, Y, ML, MU, PD, NROWPD)
!               REAL(KIND=WP) T, Y(NEQ), PD(NROWPD,NEQ)
!          where NEQ, T, Y, ML, MU, and NROWPD are input and the array
!          PD is to be loaded with partial derivatives (elements of the
!          Jacobian matrix) in the output. PD must be given a first
!          dimension of NROWPD. T and Y have the same meaning as in
!          Subroutine F.
!               In the full matrix case (MITER = 1), ML and MU are
!          ignored, and the Jacobian is to be loaded into PD in
!          columnwise manner, with df(i)/dy(j) loaded into PD(i,j).
!               In the band matrix case (MITER = 4), the elements
!          within the band are to be loaded into PD in columnwise
!          manner, with diagonal lines of df/dy loaded into the rows
!          of PD. Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).
!          ML and MU are the half-bandwidth parameters. (See IUSER).
!          The locations in PD in the two triangular areas which
!          correspond to nonexistent matrix elements can be ignored
!          or loaded arbitrarily, as they are overwritten by DVODE.
!               In the sparse matrix case the elements of the matrix
!          are determined by the sparsity structure given by the
!          IA and JA pointer arrays. Refer to the documentation
!          prologue for SET_OPTS for a description of the arguments
!          for JAC since they differ from the dense and banded cases.
!               JAC need not provide df/dy exactly. A crude
!          approximation (possibly with a smaller bandwidth) will do.
!               In either case, PD is preset to zero by the solver,
!          so that only the nonzero elements need be loaded by JAC.
!          Each call to JAC is preceded by a call to F with the same
!          arguments NEQ, T, and Y. Thus to gain some efficiency,
!          intermediate quantities shared by both calculations may be
!          saved in a user common block by F and not recomputed by JAC,
!          if desired. Also, JAC may alter the Y array, if desired.
!          JAC must be declared external in the calling program.
! MF     = The method flag. Used only for input. The legal values of
!          MF are 10, 11, 12, 13, 14, 15, 16, 17, 20, 21, 22, 23, 24,
!          25, 26, 27, -11, -12, -14, -15, -21, -22, -24, -25, -26,
!          -27.
!          MF is a signed two-digit integer, MF = JSV*(10*METH+MITER).
!          JSV = SIGN(MF) indicates the Jacobian-saving strategy:
!            JSV =  1 means a copy of the Jacobian is saved for reuse
!                     in the corrector iteration algorithm.
!            JSV = -1 means a copy of the Jacobian is not saved
!                     (valid only for MITER = 1, 2, 4, or 5).
!          METH indicates the basic linear multistep method:
!            METH = 1 means the implicit Adams method.
!            METH = 2 means the method based on backward
!                     differentiation formulas (BDF-s).
!          MITER indicates the corrector iteration method:
!            MITER = 0 means functional iteration (no Jacobian matrix
!                      is involved).
!            MITER = 1 means chord iteration with a user-supplied
!                      full (NEQ by NEQ) Jacobian.
!            MITER = 2 means chord iteration with an internally
!                      generated (difference quotient) full Jacobian
!                      (using NEQ extra calls to F per df/dy value).
!            MITER = 3 means chord iteration with an internally
!                      generated diagonal Jacobian approximation
!                      (using 1 extra call to F per df/dy evaluation).
!            MITER = 4 means chord iteration with a user-supplied
!                      banded Jacobian.
!            MITER = 5 means chord iteration with an internally
!                      generated banded Jacobian (using ML+MU+1 extra
!                      calls to F per df/dy evaluation).
!            MITER = 6 means chord iteration with a user-supplied
!                      sparse Jacobian.
!            MITER = 7 means chord iteration with an internally
!                      generated sparse Jacobian
!          If MITER = 1, 4, or 6 the user must supply a subroutine
!          JAC(the name is arbitrary) as described above under JAC.
!          For other values of MITER, a dummy argument can be used.
!                         Optional Input
! The following is a list of the optional input provided for in the
! call sequence. (See also Part ii.)  For each such input variable,
! this table lists its name as used in this documentation, its
! location in the call sequence, its meaning, and the default value.
! The use of any of this input requires IOPT = 1, and in that
! case all of this input is examined. A value of zero for any of
! these optional input variables will cause the default value to be
! used. Thus to use a subset of the optional input, simply preload
! locations 5 to 10 in RUSER and IUSER to 0.0 and 0, respectively,
! and then set those of interest to nonzero values.
! NAME    LOCATION      MEANING AND DEFAULT VALUE
! H0      RUSER(5)  The step size to be attempted on the first step.
!                   The default value is determined by the solver.
! HMAX    RUSER(6)  The maximum absolute step size allowed.
!                   The default value is infinite.
! HMIN    RUSER(7)  The minimum absolute step size allowed.
!                   The default value is 0. (This lower bound is not
!                   enforced on the final step before reaching TCRIT
!                   when ITASK = 4 or 5.)
! MAXORD  IUSER(5)  The maximum order to be allowed. The default
!                   value is 12 if METH = 1, and 5 if METH = 2.
!                   If MAXORD exceeds the default value, it will
!                   be reduced to the default value.
!                   If MAXORD is changed during the problem, it may
!                   cause the current order to be reduced.
! MXSTEP  IUSER(6)  Maximum number of (internally defined) steps
!                   allowed during one call to the solver.
!                   The default value is 5000.
! MXHNIL  IUSER(7)  Maximum number of messages printed (per problem)
!                   warning that T + H = T on a step (H = step size).
!                   This must be positive to result in a non-default
!                   value. The default value is 10.
!                          Optional Output
! As optional additional output from DVODE, the variables listed
! below are quantities related to the performance of DVODE
! which are available to the user. These are communicated by way of
! the work arrays, but also have internal mnemonic names as shown.
! Except where stated otherwise, all of this output is defined
! on any successful return from DVODE, and on any return with
! ISTATE = -1, -2, -4, -5, or -6. On an illegal input return
! (ISTATE = -3), they will be unchanged from their existing values
! (if any), except possibly for TOLSF, LENRW, and LENIW.
! On any error return, output relevant to the error will be defined,
! as noted below.
! NAME    LOCATION      MEANING
! HU      RUSER(11) The step size in t last used (successfully).
! HCUR    RUSER(12) The step size to be attempted on the next step.
! TCUR    RUSER(13) The current value of the independent variable
!                   which the solver has actually reached, i.e. the
!                   current internal mesh point in t. In the output,
!                   TCUR will always be at least as far from the
!                   initial value of t as the current argument T,
!                   but may be farther (if interpolation was done).
! TOLSF   RUSER(14) A tolerance scale factor, greater than 1.0,
!                   computed when a request for too much accuracy was
!                   detected (ISTATE = -3 if detected at the start of
!                   the problem, ISTATE = -2 otherwise). If ITOL is
!                   left unaltered but RTOL and ATOL are uniformly
!                   scaled up by a factor of TOLSF for the next call,
!                   then the solver is deemed likely to succeed.
!                   (The user may also ignore TOLSF and alter the
!                   tolerance parameters in any other way appropriate.)
! NST     IUSER(11) The number of steps taken for the problem so far.
! NFE     IUSER(12) The number of f evaluations for the problem so far.
! NJE     IUSER(13) The number of Jacobian evaluations so far.
! NQU     IUSER(14) The method order last used (successfully).
! NQCUR   IUSER(15) The order to be attempted on the next step.
! IMXER   IUSER(16) The index of the component of largest magnitude in
!                   the weighted local error vector (e(i)/EWT(i)),
!                   on an error return with ISTATE = -4 or -5.
! LENRW   IUSER(17) The length of RUSER actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
! LENIW   IUSER(18) The length of IUSER actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
! NLU     IUSER(19) The number of matrix LU decompositions so far.
! NNI     IUSER(20) The number of nonlinear (Newton) iterations so far.
! NCFN    IUSER(21) The number of convergence failures of the nonlinear
!                   solver so far.
! NETF    IUSER(22) The number of error test failures of the integrator
!                   so far.
! The following two arrays are segments of the RUSER array which
! may also be of interest to the user as optional output.
! For each array, the table below gives its internal name,
! its base address in RUSER, and its description.
!                  Interrupting and Restarting
! If the integration of a given problem by DVODE is to be interrupted
! and then later continued, such as when restarting an interrupted run
! or alternating between two or more ODE problems, the user should save,
! following the return from the last DVODE call prior to the
! interruption, the contents of the call sequence variables and
! internal PRIVATE variables, and later restore these values before the
! next DVODE call for that problem. To save and restore the PRIVATE
! variables, use subroutine DVSRCO, as described below in part ii.
! In addition, if non-default values for either LUN or MFLAG are
! desired, an extra call to XSETUN and/or XSETF should be made just
! before continuing the integration. See Part ii below for details.
! Part ii. Other Routines Callable.
! The following are optional calls which the user may make to
! gain additional capabilities in conjunction with DVODE.
! (The routines XSETUN and XSETF are designed to conform to the
! SLATEC error handling package.)
!     FORM OF CALL                  FUNCTION
!  CALL XSETUN(LUN)           Set the logical unit number, LUN, for
!                             output of messages from DVODE, if
!                             the default is not desired.
!                             The default value of LUN is 6.
!  CALL XSETF(MFLAG)          Set a flag to control the printing of
!                             messages by DVODE.
!                             MFLAG = 0 means do not print. (Danger:
!                             This risks losing valuable information.)
!                             Either of the above calls may be made at
!                             any time and will take effect immediately.
!  CALL DVINDY(...)           Provide derivatives of y, of various
!                             orders, at a specified point T, if
!                             desired. It may be called only after
!                             a successful return from DVODE.
! The detailed instructions for using DVINDY are as follows.
! The form of the call is:
!      CALL DVINDY(T,K,DKY,IFLAG)
! The input parameters are:
! T         = Value of independent variable where answers are desired
!             (normally the same as the T last returned by DVODE).
!             For valid results, T must lie between TCUR - HU and TCUR.
!             (See optional output for TCUR and HU.)
! K         = Integer order of the derivative desired. K must satisfy
!             0 <= K <= NQCUR, where NQCUR is the current order
!             (see optional output). The capability corresponding
!             to K = 0, i.e. computing y(T), is already provided
!             by DVODE directly. Since NQCUR >= 1, the first
!             derivative dy/dt is always available with DVINDY.
! The output parameters are:
! DKY       = A real array of length NEQ containing the computed value
!             of the K-th derivative of y(t).
! IFLAG     = Integer flag, returned as 0 if K and T were legal,
!             -1 if K was illegal, and -2 if T was illegal.
!             On an error return, a message is also written.
! Part iii. Optionally Replaceable Solver Routines.
! Below are descriptions of two routines in the DVODE package which
! relate to the measurement of errors. Either routine can be
! replaced by a user-supplied version, if desired. However, since such
! a replacement may have a major impact on performance, it should be
! done only when absolutely necessary, and only with great caution.
! (Note: The means by which the package version of a routine is
! superseded by the user's version may be system-dependent.)
! (a) DEWSET.
! The following subroutine is called just before each internal
! integration step, and sets the array of error weights, EWT, as
! described under ITOL/RTOL/ATOL above:
!     SUBROUTINE DEWSET(NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
! where NEQ, ITOL, RTOL, and ATOL are as in the DVODE call sequence,
! YCUR contains the current dependent variable vector, and
! EWT is the array of weights set by DEWSET.
! If the user supplies this subroutine, it must return in EWT(i)
! (i = 1,...,NEQ) a positive quantity suitable for comparison with
! errors in Y(i). The EWT array returned by DEWSET is passed to the
! DVNORM function (See below.), and also used by DVODE in the
! computation of the optional output IMXER, the diagonal Jacobian
! approximation, and the increments for difference quotient Jacobians.
! In the user-supplied version of DEWSET, it may be desirable to use
! the current values of derivatives of y. Derivatives up to order NQ
! are available from the history array YH, described above under
! Optional Output. In DEWSET, YH is identical to the YCUR array,
! extended to NQ + 1 columns with a column length of NYH and scale
! factors of h**j/factorial(j). On the first call for the problem,
! given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
! NYH is the initial value of NEQ. Thus, for example, the current
! value of dy/dt can be obtained as YCUR(NYH+i)/H  (i=1,...,NEQ)
! (and the division by H is unnecessary when NST = 0).
! (b) DVNORM.
! The following is a function which computes the weighted
! root-mean-square norm of a vector v:
!     D = DVNORM(N, V, W)
! where:
!   N = the length of the vector,
!   V = real array of length N containing the vector,
!   W = real array of length N containing weights,
!   D = sqrt((1/N) * sum(V(i)*W(i))**2).
! DVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
! EWT is as set by subroutine DEWSET.
! If the user supplies this routine, it should return a nonnegative
! value of DVNORM suitable for use in the error control in DVODE.
! None of the arguments should be altered by DVNORM.
! For example, a user-supplied DVNORM function might:
!   -substitute a max-norm of (V(i)*W(i)) for the rms-norm, or
!   -ignore some components of V in the norm, with the effect of
!    suppressing the error control on those components of Y.
!_______________________________________________________________________
! Other Routines in the DVODE Package

! In addition to subroutine DVODE, the DVODE package includes the
! following subroutines and function routines (not user callable):
!  DVHIN       computes an approximate step size for the initial step.
!  DVINDY_CORE computes an interpolated value of the y vector at t=TOUT.
!  DVINDY      computes an interpolated value of the y vector at t=TOUT.
!              (user callable)
!  DVSTEP      is the core integrator, which does one step of the
!              integration and the associated error control.
!  DVSET       sets all method coefficients and test constants.
!  DVNLSD,     solves the underlying nonlinear system -- the corrector.
!  DVNLSS28
!  DVJAC,      computes and preprocesses the Jacobian matrix J = df/dy
!  DVJACS28    and the Newton iteration matrix P = I - (h/l1)*J.
!  DVSOL,      manages solution of linear system in chord iteration.
!  DVSOLS28
!  DVJUST      adjusts the history array on a change of order.
!  DEWSET      sets the error weight vector EWT before each step.
!  DVNORM      computes the weighted r.m.s. norm of a vector.
!  DACOPY      is a routine to copy a two-dimensional array to another.
!  DGEFA_F90 and DGESL_F90 are routines from LINPACK for solving full
!              systems of linear algebraic equations.
!  DGBFA_F90 and DGBSL_F90 are routines from LINPACK for solving banded
!              linear systems.
!  DAXPY_F90, DSCAL_F90, and DCOPY_F90 are basic linear algebra modules
!              (BLAS).
!  DVCHECK     does preliminary checking for roots, and serves as an
!              interface between subroutine DVODE_F90 and subroutine
!              DVROOTS.
!  DVROOTS     finds the leftmost root of a set of functions.
! ______________________________________________________________________
! Section 10.  Example Usage
!
! MODULE example1
! The following is a simple example problem, with the coding
! needed for its solution by DVODE_F90. The problem is from
! chemical kinetics, and consists of the following three rate
! equations:
!     dy1/dt = -.04d0*y1 + 1.d4*y2*y3
!     dy2/dt = .04d0*y1 - 1.d4*y2*y3 - 3.d7*y2**2
!     dy3/dt = 3.d7*y2**2
! on the interval from t = 0.0d0 to t = 4.d10, with initial
! conditions y1 = 1.0d0, y2 = y3 = 0.0d0. The problem is stiff.
! The following coding solves this problem with DVODE_F90,
! using a user supplied Jacobian and printing results at
! t = .4, 4.,...,4.d10. It uses ITOL = 2 and ATOL much smaller
! for y2 than y1 or y3 because y2 has much smaller values. At
! the end of the run, statistical quantities of interest are
! printed. (See optional output in the full DVODE description
! below.) Output is written to the file example1.dat.
! CONTAINS
!     SUBROUTINE FEX(NEQ, T, Y, YDOT)
!     IMPLICIT NONE
!     INTEGER NEQ
!     DOUBLE PRECISION T, Y, YDOT
!     DIMENSION Y(NEQ), YDOT(NEQ)
!     YDOT(1) = -.04D0*Y(1) + 1.D4*Y(2)*Y(3)
!     YDOT(3) = 3.E7*Y(2)*Y(2)
!     YDOT(2) = -YDOT(1) - YDOT(3)
!     RETURN
!     END SUBROUTINE FEX
!     SUBROUTINE JEX(NEQ, T, Y, ML, MU, PD, NRPD)
!     IMPLICIT NONE
!     INTEGER NEQ,ML,MU,NRPD
!     DOUBLE PRECISION PD, T, Y
!     DIMENSION Y(NEQ), PD(NRPD,NEQ)
!     PD(1,1) = -.04D0
!     PD(1,2) = 1.D4*Y(3)
!     PD(1,3) = 1.D4*Y(2)
!     PD(2,1) = .04D0
!     PD(2,3) = -PD(1,3)
!     PD(3,2) = 6.E7*Y(2)
!     PD(2,2) = -PD(1,2) - PD(3,2)
!     RETURN
!     END SUBROUTINE JEX
! END MODULE example1
!******************************************************************

!     PROGRAM runexample1
!     USE DVODE_F90_M
!     USE example1
!     IMPLICIT NONE
!     DOUBLE PRECISION ATOL, RTOL, T, TOUT, Y, RSTATS
!     INTEGER NEQ, ITASK, ISTATE, ISTATS, IOUT, IERROR, I
!     DIMENSION Y(3), ATOL(3), RSTATS(22), ISTATS(31)
!     TYPE(VODE_OPTS) :: OPTIONS
!     OPEN(UNIT=6, FILE = 'example1.dat')
!     IERROR = 0
!     NEQ = 3
!     Y(1) = 1.0D0
!     Y(2) = 0.0D0
!     Y(3) = 0.0D0
!     T = 0.0D0
!     TOUT = 0.4D0
!     RTOL = 1.D-4
!     ATOL(1) = 1.D-8
!     ATOL(2) = 1.D-14
!     ATOL(3) = 1.D-6
!     ITASK = 1
!     ISTATE = 1
!     OPTIONS = SET_OPTS(DENSE_J=.TRUE.,ABSERR_VECTOR=ATOL, &
!       RELERR=RTOL, USER_SUPPLIED_JACOBIAN=.TRUE.)
!     DO IOUT = 1,12
!       CALL DVODE_F90(FEX,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,J_FCN=JEX)
!       CALL GET_STATS(RSTATS,ISTATS)
!       WRITE(6,63)T,Y(1),Y(2),Y(3)
!       DO I = 1, NEQ
!          IF (Y(I) < 0.0D0) IERROR = 1
!       END DO
!       IF (ISTATE < 0) THEN
!          WRITE(6,64)ISTATE
!          STOP
!       END IF
!       TOUT = TOUT*10.0D0
!     END DO
!     WRITE(6,60) ISTATS(11),ISTATS(12),ISTATS(13),ISTATS(19), &
!                 ISTATS(20),ISTATS(21),ISTATS(22)
!     IF (IERROR == 1) THEN
!        WRITE(6,61)
!     ELSE
!        WRITE(6,62)
!     END IF
! 60  FORMAT(/'  No. steps =',I4,'   No. f-s =',I4,        &
!             '  No. J-s =',I4,'   No. LU-s =',I4/         &
!             '  No. nonlinear iterations =',I4/           &
!             '  No. nonlinear convergence failures =',I4/ &
!             '  No. error test failures =',I4/)
! 61  FORMAT(/' An error occurred.')
! 62  FORMAT(/' No errors occurred.')
! 63  FORMAT(' At t =',D12.4,'   y =',3D14.6)
! 64  FORMAT(///' Error halt: ISTATE =',I3)
!     STOP
!     END PROGRAM runexample1
!
! MODULE example2
! The following is a modification of the previous example
! program to illustrate root finding. The problem is from
! chemical kinetics, and consists of the following three
! rate equations:
!     dy1/dt = -.04d0*y1 + 1.d4*y2*y3
!     dy2/dt = .04d0*y1 - 1.d4*y2*y3 - 3.d7*y2**2
!     dy3/dt = 3.d7*y2**2
! on the interval from t = 0.0d0 to t = 4.d10, with initial
! conditions y1 = 1.0d0, y2 = y3 = 0.0d0. The problem is stiff.
! In addition, we want to find the values of t, y1, y2,
! and y3 at which:
!   (1) y1 reaches the value 1.d-4, and
!   (2) y3 reaches the value 1.d-2.
! The following coding solves this problem with DVODE_F90
! using an internally generated dense Jacobian and
! printing results at t = .4, 4., ..., 4.d10, and at the
! computed roots. It uses ITOL = 2 and ATOL much smaller
! for y2 than y1 or y3 because y2 has much smaller values.
! At the end of the run, statistical quantities of interest
! are printed (see optional outputs in the full description
! below). Output is written to the file example2.dat.
! CONTAINS
!     SUBROUTINE FEX (NEQ, T, Y, YDOT)
!     IMPLICIT NONE
!     INTEGER NEQ
!     DOUBLE PRECISION T, Y, YDOT
!     DIMENSION Y(3), YDOT(3)
!     YDOT(1) = -0.04D0*Y(1) + 1.0D4*Y(2)*Y(3)
!     YDOT(3) = 3.0D7*Y(2)*Y(2)
!     YDOT(2) = -YDOT(1) - YDOT(3)
!     RETURN
!     END SUBROUTINE FEX
!     SUBROUTINE GEX (NEQ, T, Y, NG, GOUT)
!     IMPLICIT NONE
!     INTEGER NEQ, NG
!     DOUBLE PRECISION T, Y, GOUT
!     DIMENSION Y(3), GOUT(2)
!     GOUT(1) = Y(1) - 1.0D-4
!     GOUT(2) = Y(3) - 1.0D-2
!     RETURN
!     END SUBROUTINE GEX
! END MODULE example2
!******************************************************************
!     PROGRAM runexample2
!     USE DVODE_F90_M
!     USE example2
!     IMPLICIT NONE
!     INTEGER ITASK, ISTATE, NG, NEQ, IOUT, JROOT, ISTATS, &
!     IERROR, I
!     DOUBLE PRECISION ATOL, RTOL, RSTATS, T, TOUT, Y
!     DIMENSION Y(3), ATOL(3), RSTATS(22), ISTATS(31), JROOT(2)
!     TYPE(VODE_OPTS) :: OPTIONS
!     OPEN (UNIT=6, FILE='example2.dat')
!     IERROR = 0
!     NEQ = 3
!     Y(1) = 1.0D0
!     Y(2) = 0.0D0
!     Y(3) = 0.0D0
!     T = 0.0D0
!     TOUT = 0.4D0
!     RTOL = 1.0D-4
!     ATOL(1) = 1.0D-8
!     ATOL(2) = 1.0D-12
!     ATOL(3) = 1.0D-8
!     ITASK = 1
!     ISTATE = 1
!     NG = 2
!     OPTIONS = SET_OPTS(DENSE_J=.TRUE.,RELERR=RTOL, &
!       ABSERR_VECTOR=ATOL,NEVENTS=NG)
!     DO 40 IOUT = 1,12
! 10    CONTINUE
!       CALL DVODE_F90(FEX,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,G_FCN=GEX)
!       CALL GET_STATS(RSTATS, ISTATS, NG, JROOT)
!       WRITE(6,20) T, Y(1), Y(2), Y(3)
!       DO I = 1, NEQ
!          IF (Y(I) < 0.0D0) IERROR = 1
!       END DO
! 20    FORMAT(' At t =',D12.4,'   Y =',3D14.6)
!       IF (ISTATE < 0) GOTO 60
!       IF (ISTATE == 2) GOTO 40
!       WRITE(6,30) JROOT(1),JROOT(2)
! 30    FORMAT(5X,' The above line is a root, JROOT =',2I5)
!       ISTATE = 2
!       GOTO 10
! 40  TOUT = TOUT*10.0D0
!     WRITE(6,50) ISTATS(11), ISTATS(12), ISTATS(13), ISTATS(10)
!     IF (IERROR == 1) THEN
!        WRITE(6,61)
!     ELSE
!        WRITE(6,62)
!     END IF
! 50  FORMAT(/' No. steps =',I4,'  No. f-s =',I4,'  No. J-s =',I4, &
!     '  No. g-s =',I4/)
!     STOP
! 60  WRITE(6,70) ISTATE
! 61  FORMAT(/' An error occurred.')
! 62  FORMAT(/' No errors occurred.')
! 70  FORMAT(///' Error halt.. ISTATE =',I3)
!     STOP
!     END PROGRAM runexample2
!_______________________________________________________________________
! BEGINNING OF DVODE_F90_M PRIVATE SECTION.
! Note: This global information is used throughout DVODE_F90.
!_______________________________________________________________________
! JACSPDB arrays and parameters.
LOGICAL, PRIVATE :: USE_JACSP, LIKE_ORIGINAL_VODE
INTEGER, PRIVATE :: INFODS, LIWADS, MAXGRPDS, MINGRPDS, NRFJACDS,    &
  NCFJACDS, LWKDS, LIWKDS
INTEGER, ALLOCATABLE, PRIVATE :: INDROWDS(:), INDCOLDS(:),           &
  NGRPDS(:), IPNTRDS(:), JPNTRDS(:), IWADS(:), IWKDS(:), IOPTDS(:)
REAL (WP), ALLOCATABLE, PRIVATE :: YSCALEDS(:), WKDS(:), FACDS(:)
REAL (WP), PRIVATE :: U125, U325
!_______________________________________________________________________
LOGICAL, PARAMETER, PRIVATE :: USE_MA48_FOR_SPARSE=.FALSE.
!_______________________________________________________________________
! *****MA48 build change point. Replace the above statement.
! LOGICAL, PARAMETER, PRIVATE :: USE_MA48_FOR_SPARSE=.TRUE.
!_______________________________________________________________________
! *****MA48 build change point. Insert these statements.
! MA48 type declarations:
! TYPE(ZD01_TYPE) MATRIX
! TYPE(MA48_CONTROL) CONTROL
! TYPE(MA48_FACTORS) FACTORS
! TYPE(MA48_AINFO) AINFO
! TYPE(MA48_FINFO) FINFO
! TYPE(MA48_SINFO) SINFO
!_______________________________________________________________________
! .. Parameters ..
!     IPCUTH_MAX - maximum number of times the solver will halve the
!                  stepsize to prevent an infeasible prediction if
!                  solution bounds are used
!     KFC        - maximum number of consecutive convergence failures
!                  before crashing the order
!     KFH        - maximum number of consecutive error test failures
!                  before giving up (changed from 7 to 15)
!     MAXCOR     - maximum number of corrections
!     MSBP       - maximum number of steps before forming new P matrix
!     MXNCF      - maximum number of consecutive convergence failures
!                  before giving up
!     MXHNLO     - maximum number of T+H=T messages
!     MXSTP0     - maximum number of integration steps
!     L*****     - lengths and pointers for some internal arrays
!     INTEGER, PARAMETER, PRIVATE :: KFC = -3, KFH = -7, LENIV1 = 33,    &
    INTEGER, PARAMETER, PRIVATE :: IPCUTH_MAX = 100, KFC = -3,         &
      KFH = -15, LENIV1 = 33,                                          &
      LENIV2 = 8, LENRV1 = 48, LENRV2 = 1, LIWUSER = 30, LRWUSER = 22, &
      MAXCOR = 3, MAX_ARRAY_SIZE = 900000000, MSBP = 20, MXHNL0 = 10,  &
      MXNCF = 10, MXSTP0 = 5000
!_______________________________________________________________________
! *****LAPACK build change point. Use .TRUE. for LAPACK.
!     LOGICAL, PARAMETER, PRIVATE :: USE_LAPACK = .TRUE.
!_______________________________________________________________________
    REAL (WP), PARAMETER, PRIVATE :: ADDON = 1.0E-6_WP
    REAL (WP), PARAMETER, PRIVATE :: BIAS1 = 6.0_WP
    REAL (WP), PARAMETER, PRIVATE :: BIAS2 = 6.0_WP
    REAL (WP), PARAMETER, PRIVATE :: BIAS3 = 10.0_WP
    REAL (WP), PARAMETER, PRIVATE :: CCMAX = 0.3_WP
    REAL (WP), PARAMETER, PRIVATE :: CORTES = 0.1_WP
    REAL (WP), PARAMETER, PRIVATE :: CRDOWN = 0.3_WP
    REAL (WP), PARAMETER, PRIVATE :: ETACF = 0.25_WP
    REAL (WP), PARAMETER, PRIVATE :: ETAMIN = 0.1_WP
    REAL (WP), PARAMETER, PRIVATE :: ETAMX1 = 1.0E4_WP
    REAL (WP), PARAMETER, PRIVATE :: ETAMX2 = 10.0_WP
    REAL (WP), PARAMETER, PRIVATE :: ETAMX3 = 10.0_WP
    REAL (WP), PARAMETER, PRIVATE :: ETAMXF = 0.2_WP
    REAL (WP), PARAMETER, PRIVATE :: FIVE = 5.0_WP
    REAL (WP), PARAMETER, PRIVATE :: FOUR = 4.0_WP
    REAL (WP), PARAMETER, PRIVATE :: HALF = 0.5_WP
    REAL (WP), PARAMETER, PRIVATE :: HUN = 100.0_WP
    REAL (WP), PARAMETER, PRIVATE :: HUNDRETH = 0.01_WP
    REAL (WP), PARAMETER, PRIVATE :: ONE = 1.0_WP
    REAL (WP), PARAMETER, PRIVATE :: ONEPSM = 1.00001_WP
    REAL (WP), PARAMETER, PRIVATE :: PT1 = 0.1_WP
    REAL (WP), PARAMETER, PRIVATE :: PT2 = 0.2_WP
    REAL (WP), PARAMETER, PRIVATE :: RDIV = 2.0_WP
    REAL (WP), PARAMETER, PRIVATE :: SIX = 6.0_WP
    REAL (WP), PARAMETER, PRIVATE :: TEN = 10.0_WP
    REAL (WP), PARAMETER, PRIVATE :: TENTH = 0.1_WP
    REAL (WP), PARAMETER, PRIVATE :: THOU = 1000.0_WP
    REAL (WP), PARAMETER, PRIVATE :: THRESH = 1.5_WP
    REAL (WP), PARAMETER, PRIVATE :: TWO = 2.0_WP
    REAL (WP), PARAMETER, PRIVATE :: ZERO = 0.0_WP

! Beginning of DVODE_F90 interface.
! ..
! .. Generic Interface Blocks ..
    INTERFACE DVODE_F90

!         VODE_F90 is the interface subroutine that is actually invoked
!         when the user calls DVODE_F90. It in turn calls subroutine
!         DVODE which is the driver that directs all the work.
        MODULE PROCEDURE VODE_F90

!         GET_STATS can be called to gather integration statistics.
        MODULE PROCEDURE GET_STATS

!         DVINDY can be called to interpolate the solution and derivative.
        MODULE PROCEDURE DVINDY

!         RELEASE_ARRAYS can be called to release/deallocate the work arrays.
        MODULE PROCEDURE RELEASE_ARRAYS

!         SET_IAJA can be called to set sparse matrix descriptor arrays.
        MODULE PROCEDURE SET_IAJA

!         USERSETS_IAJA can be called to set sparse matrix descriptor arrays.
        MODULE PROCEDURE USERSETS_IAJA

!         CHECK_STAT can be called to stop if a storage allocation or
!         deallocation error occurs.
        MODULE PROCEDURE CHECK_STAT

!         JACSP can be called to calculate a Jacobian using Doug Salane's
!         algoritm
        MODULE PROCEDURE JACSP

!         DVDSM can be called to calculate sparse pointer arrays needed
!         by JACSP
        MODULE PROCEDURE DVDSM

    END INTERFACE
! ..
! .. Derived Type Declarations ..
    TYPE, PUBLIC :: VODE_OPTS
      REAL (WP), DIMENSION (:), POINTER :: ATOL, RTOL
      INTEGER :: MF, METH, MITER, MOSS, ITOL, IOPT, NG
      LOGICAL :: DENSE, BANDED, SPARSE
    END TYPE VODE_OPTS
! ..
! .. Local Scalars ..
!_______________________________________________________________________
! *****MA48 build change point. Insert these statements.
!     For communication with subroutine ma48_control_array:
!     REAL (WP), PUBLIC :: COPY_OF_U_PIVOT
!_______________________________________________________________________
    REAL (WP), PRIVATE :: ACNRM, ALPHA, BIG, BIG1, CCMXJ, CGCE, CONP, CRATE,  &
      DRC, DRES, DXMAX, EPS, ERRMAX, ETA, ETAMAX, FRACINT, FRACSUB, H, HMIN,  &
      HMXI, HNEW, HSCAL, HU, MEPS, MRESID, MRMIN, PRL1, RC, RESID, RL1,       &
      RMIN, SETH, T0ST, THEMAX, TLAST, TN, TOL, TOL1, TOUTC, UMAX, UROUND,    &
      U_PIVOT, X2, WM1, WM2
    INTEGER, PRIVATE :: ADDTOJA, ADDTONNZ, CONSECUTIVE_CFAILS,                &
      CONSECUTIVE_EFAILS, ELBOW_ROOM, IADIM, IANPIV, IAVPIV,                  &
      ICF, ICNCP, IFAIL, IMAX, IMIN, INEWJ, INIT, IPUP, IRANK, IRFND, IRNCP,  &
      ISTART, ISTATC, ITASKC, JADIM, JCUR, JMIN, JSTART, JSV, KFLAG, KOUNTL,  &
      KUTH, L, LARGE, LAST, LENIGP, LICN_ALL, LIRN_ALL, LIW, LIWM, LMAX,      &
      LOCJS, LP, LRW, LWM, LWMDIM, LWMTEMP, LYH, LYHTEMP, MANPIV, MAPIV, MAXG,&
      MAXIT, MAXORD, MB28, MB48, METH, MICN, MICNCP, MINICN, MINIRN, MIRANK,  &
      MIRN, MIRNCP, MITER, MLP, MOSS, MP, MSBG, MSBJ, MXHNIL, MXSTEP, N, NZB, &
      NCFN, NDROP, NDROP1, NDX, NETF, NEWH, NEWQ, NFE, NGC, NGE, NGP, NHNIL,  &
      NJE, NLP, NLU, NNI, NNZ, NOITER, NQ, NQNYH, NQU, NQWAIT, NSLG, NSLJ,    &
      NSLP, NSRCH, NSRCH1, NST, NSUBS, NSUPS, NUM, NUMNZ, NYH, NZ_ALL,        &
      NZ_SWAG, PREVIOUS_MAXORD, WPD, WPS, MA28AD_CALLS, MA28BD_CALLS,         &
      MA28CD_CALLS, MC19AD_CALLS, MAX_MINIRN, MAX_MINICN, MAX_NNZ, BNGRP
!       MA48AD_CALLS, MA48BD_CALLS, MA48CD_CALLS
! *****MA48 build change point. Insert the above line.
    LOGICAL, PRIVATE :: ABORT, ABORT1, ABORT2, ABORT3, ABORTA, ABORTB,        &
      ALLOW_DEFAULT_TOLS, BUILD_IAJA, BOUNDS, CHANGED_ACOR, GROW, IAJA_CALLED,&
      J_HAS_BEEN_COMPUTED, J_IS_CONSTANT, LBIG, LBIG1, LBLOCK, MA48_WAS_USED, &
      OK_TO_CALL_MA28, SUBS, SUPS, OPTS_CALLED, REDO_PIVOT_SEQUENCE,          &
      SCALE_MATRIX, SPARSE, USE_FAST_FACTOR, YMAXWARN
! ..
! .. Local Arrays ..
    REAL (WP), ALLOCATABLE, PRIVATE :: ACOR(:), CSCALEX(:), EWT(:),           &
      FPTEMP(:), FTEMP(:), FTEMP1(:), G0(:), G1(:), GX(:), JMAT(:),           &
      LB(:), PMAT(:), RSCALEX(:), RWORK(:), SAVF(:), UB(:), WM(:),            &
      WMTEMP(:), WSCALEX(:,:), YHNQP2(:), YHTEMP(:), YMAX(:), YNNEG(:),       &
      YTEMP(:), DTEMP(:)
    REAL (WP), PRIVATE :: EL(13), RUSER(22), TAU(13), TQ(5)
    INTEGER, ALLOCATABLE, PRIVATE :: BIGP(:), BJGP(:), IA(:), IAB(:), IAN(:), &
      ICN(:), IDX(:), IGP(:), IKEEP28(:,:), IW28(:,:), IWORK(:), JA(:),       &
      JAB(:), JAN(:), JATEMP(:), JGP(:), JROOT(:), JVECT(:), SUBDS(:), SUPDS(:)
    INTEGER, PRIVATE :: IDISP(2), IUSER(30), LNPIV(10), LPIV(10)
    INTEGER, PRIVATE :: MORD(2) = (/ 12, 5 /)
! ..
! .. Public Subroutines and Functions ..
PUBLIC ::                                                        &
DAXPY_F90, DCOPY_F90, DDOT_F90, DGBFA_F90, DGBSL_F90, DGEFA_F90, &
DGESL_F90, DSCAL_F90, IDAMAX_F90
! ..
! .. Private Subroutines and Functions ..
PRIVATE ::                                                       &
CHECK_DIAG    , DACOPY        , DEWSET        , DGROUP        ,  &
DGROUPDS      , DVCHECK       , DVHIN         , DVINDY_BNDS   ,  &
DVINDY_CORE   , DVJAC         , DVJACS28      , DVJUST        ,  &
DVNLSD        , DVNLSS28      , DVNORM        , DVNRDN        ,  &
DVNRDP        , DVNRDS        , DVODE         , DVPREPS       ,  &
DVROOTS       , DVSET         , DVSOL         , DVSOLS28      ,  &
DVSRCO        , DVSTEP        , GDUMMY        , IUMACH        ,  &
IXSAV         , JACSPDB       , JDUMMY        , MA28AD        ,  &
MA28BD        , MA28CD        , MA28DD        , MA28ID        ,  &
MA30AD        , MA30BD        , MA30CD        , MA30DD        ,  &
MC13E         , MC19AD        , MC20AD        , MC20BD        ,  &
MC21A         , MC21B         , MC22AD        , MC23AD        ,  &
MC24AD        , SET_ICN       , XERRDV        , XSETF         ,  &
XSETUN        , DEGR          , IDO           , NUMSRT        ,  &
SEQ           , SETR          , SLO           , SRTDAT        ,  &
FDJS
! DVJACS48      , DVNLSS48      , DVPREPS48     , DVSOLS48
!_______________________________________________________________________
! *****MA48 build change point. Insert the above line.
!_______________________________________________________________________
! ..
! .. Intrinsic Functions ..
    INTRINSIC KIND
! ..
! .. Data Statements ..
    DATA OPTS_CALLED/ .FALSE./
    DATA MP/6/, NLP/6/, MLP/6/, NSRCH/32768/, ISTART/0/, MAXIT/16/, &
      LBIG/ .FALSE./, LBLOCK/ .TRUE./, GROW/ .TRUE./,               &
      TOL/0.0_WP/, CGCE/0.5_WP/, BIG/0.0_WP/, ABORT1/ .TRUE./,      &
      ABORT2/ .TRUE./, ABORT3/ .FALSE./, ABORT/ .FALSE./, MIRN/0/,  &
      MICN/0/, MIRNCP/0/, MICNCP/0/, MIRANK/0/, NDROP1/0/,          &
      MRMIN/0.0D0/, MRESID/0/, OK_TO_CALL_MA28/.FALSE./
! ..
! END OF DVODE_F90 PRIVATE SECTION.
!_______________________________________________________________________

  CONTAINS

    SUBROUTINE VODE_F90(F,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTS,J_FCN,G_FCN)
! ..
! This is an interface for DVODE to allow JAC and GFUN to be
! OPTIONAL arguments.
! ..
   IMPLICIT NONE
! ..
! .. Structure Arguments ..
      TYPE (VODE_OPTS) :: OPTS
! ..
! .. Scalar Arguments ..
      REAL (WP), INTENT (INOUT) :: T, TOUT
      INTEGER, INTENT (INOUT) :: ISTATE
      INTEGER, INTENT (IN) :: ITASK, NEQ
! ..
! .. Array Arguments ..
      REAL (WP), INTENT (INOUT) :: Y(*)
! ..
! .. Subroutine Arguments ..
      OPTIONAL :: G_FCN, J_FCN
      EXTERNAL J_FCN
! ..
! .. Subroutine Interfaces ..
   INTERFACE
     SUBROUTINE F(NEQ,T,Y,YDOT)
       INTEGER, PARAMETER :: WP = KIND(1.0D0)
       INTEGER NEQ
       REAL(WP) T
       REAL(WP), DIMENSION(NEQ) :: Y, YDOT
       INTENT(IN)  :: NEQ, T, Y
       INTENT(OUT) :: YDOT
     END SUBROUTINE F
   END INTERFACE

   INTERFACE
     SUBROUTINE G_FCN(NEQ,T,Y,NG,GROOT)
       INTEGER, PARAMETER :: WP = KIND(1.0D0)
       INTEGER NEQ, NG
       REAL(WP) T
       REAL(WP), DIMENSION(NEQ) :: Y
       REAL(WP), DIMENSION(NG) :: GROOT(NG)
       INTENT(IN)  :: NEQ, T, Y, NG
       INTENT(OUT) :: GROOT
     END SUBROUTINE G_FCN
   END INTERFACE

!    Note:
!    The best we can do here is to declare J_FCN to be
!    external. The interface for a sparse problem differs
!    from that for a banded or dense problem. The following
!    would suffuce for a banded or dense problem.
!    INTERFACE
!      SUBROUTINE J_FCN(NEQ,T,Y,ML,MU,PD,NROWPD)
!        INTEGER, PARAMETER :: WP = KIND(1.0D0)
!        INTEGER NEQ, ML, MU, NROWPD
!        REAL(WP) T
!        REAL(WP), DIMENSION(NEQ) :: Y
!        REAL(WP), DIMENSION(NEQ) :: PD(NROWPD,NEQ)
!        INTENT(IN)  :: NEQ, T, Y, ML, MU, NROWPD
!        INTENT(INOUT) :: PD
!      END SUBROUTINE J_FCN
!    END INTERFACE
!    The following would suffice for a sparse problem.
!    INTERFACE
!      SUBROUTINE J_FCN(NEQ,T,Y,IA,JA,NZ,P)
!        INTEGER, PARAMETER :: WP = KIND(1.0D0)
!        INTEGER NEQ, NZ
!        REAL(WP) T
!        REAL(WP), DIMENSION Y(*), P(*)
!        INTEGER, DIMENSION IA(*), JA(*)
!        INTENT(IN) :: NEQ, T, Y
!        INTENT(INOUT) IA, JA, NZ, P
!      END SUBROUTINE J_FCN
!    END INTERFACE
! ..
! .. Local Scalars ..
      INTEGER :: HOWCALL, METH, MFA, MITER, MOSS, NG
      CHARACTER (80) :: MSG
! ..
! .. Intrinsic Functions ..
      INTRINSIC ABS, PRESENT
! ..
! .. FIRST EXECUTABLE STATEMENT VODE_F90
! ..
!       Check that SET_OPTS has been called.
      IF (.NOT.OPTS_CALLED) THEN
        MSG = 'You have not called SET_OPTS before'
        CALL XERRDV(MSG,10,1,0,0,0,0,ZERO,ZERO)
        MSG = 'calling DVODE_F90 the first time.'
        CALL XERRDV(MSG,10,2,0,0,0,0,ZERO,ZERO)
      END IF

!       Check that JAC is present if it is needed.
      IF (PRESENT(J_FCN)) THEN
      ELSE
!         Note:
!         MOSS is irrelevant. OPTS%MF is two digits after the
!         call to SET_OPTS.
        MFA = ABS(OPTS%MF)
        MOSS = MFA/100
        METH = (MFA-100*MOSS)/10
        MITER = MFA - 100*MOSS - 10*METH
        IF (MITER==1 .OR. MITER==4 .OR. MITER==6) THEN
          MSG = 'You have specified a value of the integration'
          CALL XERRDV(MSG,20,1,0,0,0,0,ZERO,ZERO)
          MSG = 'method flag MF which requires that you supply'
          CALL XERRDV(MSG,20,1,0,0,0,0,ZERO,ZERO)
          MSG = 'a Jacobian subroutine JAC; but FAC is not'
          CALL XERRDV(MSG,20,1,0,0,0,0,ZERO,ZERO)
          MSG = 'present in the argument list.'
          CALL XERRDV(MSG,20,2,0,0,0,0,ZERO,ZERO)
        END IF
      END IF

!       Check that GFUN is present if it is needed.

      IF (PRESENT(G_FCN)) THEN
      ELSE
        NG = OPTS%NG
        IF (NG>0) THEN
          MSG = 'You have indicated that events are present but'
          CALL XERRDV(MSG,30,1,0,0,0,0,ZERO,ZERO)
          MSG = 'you have not supplied a GFUN subroutine.'
          CALL XERRDV(MSG,30,2,0,0,0,0,ZERO,ZERO)
        END IF
      END IF

!     Determine how DVODE will be called.

!     HOWCALL = 1: JDUMMY, GDUMMY
!               2: JAC, GFUN
!               3: JAC, GDUMMY
!               4: JDUMMY, GFUN
      HOWCALL = 1
      IF (PRESENT(J_FCN)) THEN
        IF (PRESENT(G_FCN)) THEN
          HOWCALL = 2
        ELSE
          HOWCALL = 3
        END IF
      ELSE
        IF (PRESENT(G_FCN)) THEN
          HOWCALL = 4
        ELSE
          HOWCALL = 1
        END IF
      END IF

!       Call DVODE to do the integration.

      IF (HOWCALL==1) THEN
        CALL DVODE(F,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTS,JDUMMY,GDUMMY)
      ELSE IF (HOWCALL==2) THEN
        CALL DVODE(F,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTS,J_FCN,G_FCN)
      ELSE IF (HOWCALL==3) THEN
        CALL DVODE(F,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTS,J_FCN,GDUMMY)
      ELSE IF (HOWCALL==4) THEN
        CALL DVODE(F,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTS,JDUMMY,G_FCN)
      END IF
      RETURN

    END SUBROUTINE VODE_F90
!_______________________________________________________________________

    SUBROUTINE JDUMMY(NEQ,T,Y,ML,MU,PD,NROWPD)
! ..
! This is a dummy Jacobian subroutine for VODE_F90 (never called).
! ..
   IMPLICIT NONE
! ..
! .. Scalar Arguments ..
      REAL (WP) :: T
      INTEGER :: ML, MU, NEQ, NROWPD, I
      LOGICAL DUMMY
! ..
! .. Array Arguments ..
      REAL (WP) :: PD(NROWPD,*), Y(*)
! ..
    INTENT(IN) T, Y, ML, MU, NROWPD
    INTENT(INOUT) PD
! ..
! .. FIRST EXECUTABLE STATEMENT JDUMMY
! ..
!       Get rid of some needless compiler warning messages.
      DUMMY = .FALSE.
      IF (DUMMY) THEN
        I = NEQ
        I = ML
        I = MU
        I = NROWPD
        PD(1,1) = T
        PD(1,1) = Y(1)
        PD(1,1) = DBLE(REAL(I))
      END IF
      RETURN

    END SUBROUTINE JDUMMY
!_______________________________________________________________________

    SUBROUTINE GDUMMY(NEQ,T,Y,NG,GOUT)
! ..
! This is a dummy event subroutine for VODE_F90 (never called).
! ..
   IMPLICIT NONE
! ..
! .. Scalar Arguments ..
      REAL (WP) :: T
      INTEGER :: NEQ, NG, I
      LOGICAL DUMMY
! ..
! .. Array Arguments ..
      REAL (WP) :: GOUT(*), Y(*)
! ..
    INTENT(IN) NEQ, T, Y, NG
    INTENT(OUT) GOUT
! ..
! .. FIRST EXECUTABLE STATEMENT JDUMMY
! ..
!       Get rid of some needless compiler warning messages.
      DUMMY = .FALSE.
      IF (DUMMY) THEN
        I = NEQ
        I = NG
        GOUT(1) = T
        GOUT(1) = Y(1)
        GOUT(1) = DBLE(REAL(I))
      END IF
      RETURN

    END SUBROUTINE GDUMMY
!_______________________________________________________________________

SUBROUTINE SET_OPTS_2(HMAX,HMIN,MXSTEP)
! ..
! Allow the maximum step size, the minimum step size, and the maximum
! number of steps to be changed without restarting the integration.
! ..
!                     Quick Summary of Options
! HMAX                   - Maximum step size in DVODE
! HMIN                   - Minimum step size in DVODE
! MXSTEP                 - Maximum number of integration steps in DVODE
! ..
   IMPLICIT NONE
! ..
! .. Scalar Arguments ..
      REAL (WP), OPTIONAL, INTENT (IN) :: HMAX, HMIN
      INTEGER, OPTIONAL, INTENT (IN) :: MXSTEP
! ..
! .. Local Scalars ..
      CHARACTER (80) :: MSG
! ..
! .. Intrinsic Functions ..
      INTRINSIC ALLOCATED, PRESENT
! ..
! .. FIRST EXECUTABLE STATEMENT SET_OPTS_2
! ..
!       Check that SET_OPTS has been called:
      IF (.NOT.OPTS_CALLED) THEN
        MSG = 'You have not called SET_OPTS before'
        CALL XERRDV(MSG,40,1,0,0,0,0,ZERO,ZERO)
        MSG = 'calling subroutine SET_OPTS_2.'
        CALL XERRDV(MSG,40,2,0,0,0,0,ZERO,ZERO)
      END IF
      IF (PRESENT(HMAX)) THEN
        RUSER(6) = HMAX
        MSG = 'HMAX changed in SET_OPTS_2.'
        CALL XERRDV(MSG,50,1,0,0,0,1,HMAX,ZERO)
      END IF
      IF (PRESENT(HMIN)) THEN
        RUSER(7) = HMIN
        MSG = 'HMIN changed in SET_OPTS_2.'
        CALL XERRDV(MSG,60,1,0,0,0,1,HMIN,ZERO)
      END IF
      IF (PRESENT(MXSTEP)) THEN
        IUSER(6) = MXSTEP
        MSG = 'MXSTEP changed in SET_OPTS_2.'
        CALL XERRDV(MSG,70,1,1,MXSTEP,0,0,ZERO,ZERO)
      END IF

END SUBROUTINE SET_OPTS_2
!_______________________________________________________________________

FUNCTION SET_NORMAL_OPTS(DENSE_J, BANDED_J, SPARSE_J,                &
  USER_SUPPLIED_JACOBIAN, LOWER_BANDWIDTH, UPPER_BANDWIDTH,          &
  RELERR, ABSERR, ABSERR_VECTOR, NEVENTS) RESULT(OPTS)

! FUNCTION SET_NORMAL_OPTS:
!    Jacobian type:
!       DENSE_J, BANDED_J, SPARSE_J
!    Analytic Jacobian:
!       USER_SUPPLIED_JACOBIAN
!    If banded Jacobian:
!       LOWER_BANDWIDTH,UPPER_BANDWIDTH
!    Error tolerances:
!       RELERR, ABSERR, ABSERR_VECTOR
!    Rootfinding:
!       NEVENTS
! RESULT(OPTS)

! Note:
! Invoking SET_NORMAL_OPTS causes an integration restart. A common
! situation is one in which all you wish to change is one of the
! vode.f77 optional parameters HMAX, HMIN, or MXSTEP. Once the
! integration is started and this is all you wish to do, you can
! change any of these parameters without restarting the integration
! simply by calling subroutine SET_OPTS_2:
!      CALL SET_OPTS_2(HMAX,HMIN,MXSTEP)
! Each of the three arguments is optional and only the ones actually
! supplied will be used. Changes will take effect in the same manner
! as in the VODE.f77 solver.
!
! NORMAL_OPTIONS sets user parameters for DVODE via keywords.
! Values that are defined herein will be used internally by
! DVODE. All option keywords are OPTIONAL and order is not
! important. These options should be adequate for most problems.
! If you wish to use more specialized options, you must use
! SET_INTERMEDIATE_OPTS or SET_OPTS rather than NORMAL_OPTS.
! If you wish to use SET_INTERMEDIATE_OPTS or SET_OPTS, you
! may use any of the SET_NORMAL_OPTS keywords or any of the
! keywords available for these two functions. Of course, you
! may opt to simply use SET_OPTS for all problems.

! Note that DVODE_F90 requires that one of SET_NORMAL_OPTS or
! SET_INTERMEDIATE_OPTS or SET_OPTS is called before the first
! time DVODE_F90 is called.
!
! Important Note:
! If feasible, you should use the dense or banded option; but
! SET_NORMAL_OPTS allows you to use a sparse internal Jacobian
! (i.e., one that is determined using finite differences) and
! structure pointere array that are determined internally
! using finite differences. If any of the following are true
!    (1) DVODE_F90 doesn't perform satisfactorily for
!        your problem,
!    (2) you are solving a very large problem,
!    (3) you wish to supply the sparse pointer arrays
!        directly,
!    (4) you wish to supply an analytical sparse Jacobian,
! or
!    (5) you wish to use one of the specialized sparse
!        Jacobian options,
! you are encouraged to use SET_OPTS which contains several
! provisions for solving sparse problems more efficiently.

! Option Types:
! DENSE_J                - logical
! BANDED_J               - logical
! SPARSE_J               - logical
! USER_SUPPLIED_JACOBIAN - logical
! LOWER_BANDWIDTH        - integer
! UPPER_BANDWIDTH        - integer
! RELERR                 - real(wp) scalar
! ABSERR                 - real(wp) scalar
! ABSERR_VECTOR          - real(wp) vector
! NEVENTS                - integer
! Options:
! ABSERR                 = Absolute error tolerance
! ABSERR_VECTOR          = Vector of absolute error tolerances
! RELERR                 = Scalar relative error tolerance
! NEVENTS                = Number of event functions (requires
!                          user-supplied GFUN)
! DENSE_J                = Use dense linear algebra if .TRUE.
! BANDED_J               = Use banded linear algebra if .TRUE.
!   LOWER_BANDWIDTH      = Lower bandwidth of the Jacobian
!                          (required if BANDED_J = .TRUE.)
!   UPPER_BANDWIDTH      = Upper bandwidth of the Jacobian
!                          (required if BANDED_J = .TRUE.)
! SPARSE_J               = Use sparse linear algebra if .TRUE.
! USER_SUPPLIED_JACOBIAN = Exact Jacobian option
!                          (requires user-supplied JAC;
!                          ignored for SPARSE_J=.TRUE.)
!
! Note: DENSE_J takes precedence over BANDED_J which in turn
! takes precedence over SPARSE_J if more than one is supplied.
! If neither of the three flags is present, the nonstiff Adams
! option will be used. Similiarly, ABSERR_VECTOR takes
! precedence over ABSERR.
!
! Note on Jacobian Storage Formats:
!
! If you supply an analytic Jacobian PD, load the
! Jacobian elements DF(I)/DY(J), the partial
! derivative of F(I) with respect to Y(J), using
! the following formats. Here, Y is the solution,
! F is the derivative, and PD is the Jacobian.
!
! For a full Jacobian, load PD(I,J) with DF(I)/DY(J).
! Your code might look like this:
!    DO J = 1, NEQ
!       DO I = 1, NEQ
!          PD(I,J) = ... DF(I)/DY(J)
!       END DO
!    END DO
!
! For a banded Jacobian, load PD(I-J+MU+1,J) with
! DF(I)/DY(J) where ML is the lower bandwidth
! and MU is the upper bandwidth of the Jacobian.
! Your code might look like this:
!    DO J = 1, NEQ
!       I1 = MAX(1,J-ML)
!       I2 = MIN(N,J+MU)
!       DO I = I1, I2
!          K = I-J+MU+1
!          PD(K,J) = ... DF(I)/DY(J)
!       END DO
!    END DO
! ..
   IMPLICIT NONE
! ..
! .. Function Return Value ..
      TYPE (VODE_OPTS) :: OPTS
! ..
! .. Scalar Arguments ..
      REAL (WP), OPTIONAL, INTENT (IN) :: ABSERR,RELERR
      INTEGER, OPTIONAL, INTENT (IN) :: LOWER_BANDWIDTH,   &
        NEVENTS,UPPER_BANDWIDTH
      LOGICAL, OPTIONAL, INTENT (IN) :: BANDED_J, DENSE_J, &
        SPARSE_J, USER_SUPPLIED_JACOBIAN
! ..
! .. Array Arguments ..
      REAL (WP), OPTIONAL, INTENT (IN) :: ABSERR_VECTOR(:)
! ..
! .. Local Scalars ..
      INTEGER :: IER,IOPT,METH,MF,MFA,MFSIGN,MITER,ML,MOSS, &
        MU,NAE,NG,NRE
      LOGICAL :: BANDED,DENSE,SPARSE
      CHARACTER (80) :: MSG
! ..
! .. Intrinsic Functions ..
      INTRINSIC ALLOCATED, IABS, MAX, MINVAL, PRESENT, SIGN, SIZE
! ..
! .. FIRST EXECUTABLE STATEMENT SET_NORMAL_OPTS
! ..
      RUSER(1:LRWUSER) = ZERO
      IUSER(1:LIWUSER) = 0

!       Allow default error tolerances?
      ALLOW_DEFAULT_TOLS = .FALSE.

!       Maximum number of consecutive error test failures?
      CONSECUTIVE_EFAILS = KFH

!       Maximum number of consecutive corrector iteration failures?
      CONSECUTIVE_CFAILS = MXNCF

!       Use JACSP to approximate Jacobian?
      USE_JACSP = .FALSE.

!       Set the flag to indicate that SET_NORMAL_OPTS has been called.
      OPTS_CALLED = .TRUE.

!       Set the MA48 storage cleanup flag.
      MA48_WAS_USED = .FALSE.

!       Set the fast factor option for MA48,
      USE_FAST_FACTOR = .TRUE.

!       Set the constant Jacobian flags.
      J_IS_CONSTANT = .FALSE.
      J_HAS_BEEN_COMPUTED = .FALSE.

!       Determine the working precision and define the value for UMAX
!       expected by MA28. Note that it is different for single and
!       double precision.
      WPD = KIND(1.0D0)
      WPS = KIND(1.0E0)
      IF (WPD/=WP .AND. WPS/=WP) THEN
        MSG = 'Illegal working precision in SET_NORMAL_OPTS.'
        CALL XERRDV(MSG,80,2,0,0,0,0,ZERO,ZERO)
      END IF
      IF (WPD==WP) THEN
!       Working precision is double.
        UMAX = 0.999999999_WP
      ELSE
!       Working precision is single.
        UMAX = 0.9999_WP
      END IF

      MA28AD_CALLS = 0
      MA28BD_CALLS = 0
      MA28CD_CALLS = 0
      MC19AD_CALLS = 0
!_______________________________________________________________________
! *****MA48 build change point. Insert these statements.
!       MA48AD_CALLS = 0
!       MA48BD_CALLS = 0
!       MA48CD_CALLS = 0
!_______________________________________________________________________
      IRNCP = 0
      ICNCP = 0
      MINIRN = 0
      MINICN = 0
      MAX_MINIRN = 0
      MAX_MINICN = 0
      MAX_NNZ = 0

!       Set the flag to warn the user if |(y(t)| < abserr.
      YMAXWARN = .FALSE.

!       Load defaults for the optional input arrays for DVODE.
      IUSER(1:8) = 0
      RUSER(1:8) = ZERO

!       Set the method flag.
      MF = 10
      IF (PRESENT(SPARSE_J)) THEN
        IF (SPARSE_J) THEN
           MF = 227
           IF (PRESENT(USER_SUPPLIED_JACOBIAN)) THEN
             MSG = 'You have indicated you wish to supply an'
             CALL XERRDV(MSG,90,1,0,0,0,0,ZERO,ZERO)
             MSG = 'exact sparse Jacobian in function'
             CALL XERRDV(MSG,90,1,0,0,0,0,ZERO,ZERO)
             MSG = 'SET_NORMAL_OPTS. In order to do this,'
             CALL XERRDV(MSG,90,1,0,0,0,0,ZERO,ZERO)
             MSG = 'you must use SET_OPTS. Execution will'
             CALL XERRDV(MSG,90,1,0,0,0,0,ZERO,ZERO)
             MSG = 'continue.'
             CALL XERRDV(MSG,90,1,0,0,0,0,ZERO,ZERO)
           END IF
        END IF
      END IF

      IF (PRESENT(BANDED_J)) THEN
        IF (BANDED_J) THEN
          IF (PRESENT(USER_SUPPLIED_JACOBIAN)) THEN
            IF (USER_SUPPLIED_JACOBIAN) THEN
              MF = 24
            ELSE
              MF = 25
            END IF
          ELSE
            MF = 25
          END IF
        END IF
      END IF

      IF (PRESENT(DENSE_J)) THEN
        IF (DENSE_J) THEN
          IF (PRESENT(USER_SUPPLIED_JACOBIAN)) THEN
            IF (USER_SUPPLIED_JACOBIAN) THEN
              MF = 21
            ELSE
              MF = 22
            END IF
          ELSE
            MF = 22
          END IF
        END IF
      END IF

!       Check for errors in MF.
      MFA = IABS(MF)
      MOSS = MFA/100
      METH = (MFA-100*MOSS)/10
      MITER = MFA - 100*MOSS - 10*METH
      IF (METH<1 .OR. METH>2) THEN
        MSG = 'Illegal value of METH in SET_NORMAL_OPTS.'
        CALL XERRDV(MSG,100,2,0,0,0,0,ZERO,ZERO)
      END IF
      IF (MITER<0 .OR. MITER>7) THEN
        MSG = 'Illegal value of MITER in SET_NORMAL_OPTS.'
        CALL XERRDV(MSG,110,2,0,0,0,0,ZERO,ZERO)
      END IF
      IF (MOSS<0 .OR. MOSS>2) THEN
        MSG = 'Illegal value of MOSS in SET_NORMAL_OPTS.'
        CALL XERRDV(MSG,120,2,0,0,0,0,ZERO,ZERO)
      END IF

!       Reset MF, now that MOSS is known.
      MFSIGN = SIGN(1,MF)
      MF = MF - 100*MOSS*MFSIGN

      IF (MITER==0) THEN
        DENSE = .FALSE.
        BANDED = .FALSE.
        SPARSE = .FALSE.
      ELSE IF (MITER==1 .OR. MITER==2) THEN
        DENSE = .TRUE.
        BANDED = .FALSE.
        SPARSE = .FALSE.
      ELSE IF (MITER==3) THEN
        DENSE = .FALSE.
        BANDED = .FALSE.
        SPARSE = .FALSE.
      ELSE IF (MITER==4 .OR. MITER==5) THEN
        DENSE = .FALSE.
        BANDED = .TRUE.
        SPARSE = .FALSE.
      ELSE IF (MITER==6 .OR. MITER==7) THEN
        DENSE = .FALSE.
        BANDED = .FALSE.
        SPARSE = .TRUE.
      END IF

!       Define the banded Jacobian band widths.
      IF (BANDED) THEN
        IF (PRESENT(LOWER_BANDWIDTH)) THEN
          ML = LOWER_BANDWIDTH
          IUSER(1) = ML
        ELSE
          MSG = 'In SET_NORMAL_OPTS you have indicated a'
          CALL XERRDV(MSG,130,1,0,0,0,0,ZERO,ZERO)
          MSG = 'banded Jacobian but you have not supplied'
          CALL XERRDV(MSG,130,1,0,0,0,0,ZERO,ZERO)
          MSG = 'the lower bandwidth.'
          CALL XERRDV(MSG,130,2,0,0,0,0,ZERO,ZERO)
        END IF
        IF (PRESENT(UPPER_BANDWIDTH)) THEN
          MU = UPPER_BANDWIDTH
          IUSER(2) = MU
        ELSE
          MSG = 'In SET_NORMAL_OPTS you have indicated a'
          CALL XERRDV(MSG,140,1,0,0,0,0,ZERO,ZERO)
          MSG = 'banded Jacobian but you have not supplied'
          CALL XERRDV(MSG,140,1,0,0,0,0,ZERO,ZERO)
          MSG = 'the upper bandwidth.'
          CALL XERRDV(MSG,140,2,0,0,0,0,ZERO,ZERO)
        END IF
      END IF

!       Define the sparse Jacobian options.
      IF (SPARSE) THEN
!         Set the MA28 message flag.
        LP = 0
!         Set the MA28 pivot sequence frequency flag.
        REDO_PIVOT_SEQUENCE = .FALSE.
!         Set the MA28 singularity threshold.
        EPS = 1.0E-4_WP
!         Use scaling of the iteration matrix.
        SCALE_MATRIX = .TRUE.
!         Define the elbow room factor for the MA28 sparse work arrays.
        ELBOW_ROOM = 2
!         NZSWAG is a swag for the number of nonzeros in the Jacobian.
        NZ_SWAG = 0
!         Use partial pivoting.
        U_PIVOT = ONE
!         Indicate that SET_IAJA has not yet been called successfully.
        IAJA_CALLED = .FALSE.
!         Check for illegal method flags.
        IF (MOSS==2 .AND. MITER/=7) THEN
          MSG = 'In SET_NORMAL_OPTS MOSS=2 but MITER is not 7.'
          CALL XERRDV(MSG,150,2,0,0,0,0,ZERO,ZERO)
        END IF
        IF (MOSS==1 .AND. MITER/=6) THEN
          MSG = 'In SET_NORMAL_OPTS MOSS=1 but MITER is not 6.'
          CALL XERRDV(MSG,160,2,0,0,0,0,ZERO,ZERO)
        END IF
!         IF (MOSS==0 .AND. MITER/=7) THEN
!           MSG = 'In SET_NORMAL_OPTS MOSS=0 but MITER is not 7.'
!           CALL XERRDV(MSG,170,2,0,0,0,0,ZERO,ZERO)
!         END IF
      END IF

!       Define the number of event functions.
      IF (PRESENT(NEVENTS)) THEN
        IF (NEVENTS>0) THEN
          NG = NEVENTS
        ELSE
          NG = 0
        END IF
      ELSE
        NG = 0
      END IF

!       No solution bounds will be imposed.
      BOUNDS = .FALSE.

!       Load the user options into the solution structure.
      OPTS%MF = MF
      OPTS%METH = METH
      OPTS%MITER = MITER
      OPTS%MOSS = MOSS
      OPTS%DENSE = DENSE
      OPTS%BANDED = BANDED
      OPTS%SPARSE = SPARSE
      OPTS%NG = NG

      IOPT = 1
      OPTS%IOPT = IOPT

!       Process the error tolerances.

!       Relative error tolerance.
      NRE = 1
      ALLOCATE (OPTS%RTOL(NRE),STAT=IER)
      CALL CHECK_STAT(IER,10)
      IF (PRESENT(RELERR)) THEN
        IF (RELERR<ZERO) THEN
          MSG = 'RELERR must be nonnegative.'
          CALL XERRDV(MSG,180,2,0,0,0,0,ZERO,ZERO)
        END IF
        OPTS%RTOL = RELERR
      ELSE
        IF (ALLOW_DEFAULT_TOLS) THEN
           OPTS%RTOL = 1.0E-4_WP
           MSG = 'By not specifying RELERR, you have elected to use a default'
           CALL XERRDV(MSG,190,1,0,0,0,0,ZERO,ZERO)
           MSG = 'relative error tolerance equal to 1.0D-4. Please be aware a'
           CALL XERRDV(MSG,190,1,0,0,0,0,ZERO,ZERO)
           MSG = 'tolerance this large is not appropriate for all problems.'
           CALL XERRDV(MSG,190,1,0,0,0,0,ZERO,ZERO)
           MSG = 'Execution will continue'
           CALL XERRDV(MSG,190,1,0,0,0,0,ZERO,ZERO)
        ELSE
           MSG = 'You must specify a nonzero relative error tolerance.'
           CALL XERRDV(MSG,200,2,0,0,0,0,ZERO,ZERO)
        END IF
      END IF

!       Absolute error tolerance(s).
      IF (PRESENT(ABSERR_VECTOR)) THEN
        IF (MINVAL(ABSERR_VECTOR)<ZERO) THEN
          MSG = 'All components of ABSERR_VECTOR must'
          CALL XERRDV(MSG,210,1,0,0,0,0,ZERO,ZERO)
          MSG = 'be nonnegative.'
          CALL XERRDV(MSG,210,2,0,0,0,0,ZERO,ZERO)
        END IF
        NAE = SIZE(ABSERR_VECTOR)
      ELSE
        NAE = 1
      END IF
      ALLOCATE (OPTS%ATOL(NAE),STAT=IER)
      CALL CHECK_STAT(IER,20)
      IF (PRESENT(ABSERR_VECTOR)) THEN
        OPTS%ATOL = ABSERR_VECTOR
      ELSE IF (PRESENT(ABSERR)) THEN
        IF (ABSERR<ZERO) THEN
          MSG = 'ABSERR must be nonnegative.'
          CALL XERRDV(MSG,220,2,0,0,0,0,ZERO,ZERO)
        END IF
        OPTS%ATOL = ABSERR
      ELSE
        IF (ALLOW_DEFAULT_TOLS) THEN
           OPTS%ATOL = 1D-6
           MSG = 'By not specifying ABSERR, you have elected to use a default'
           CALL XERRDV(MSG,230,1,0,0,0,0,ZERO,ZERO)
           MSG = 'absolute error tolerance equal to 1.0D-6. Please be aware a'
           CALL XERRDV(MSG,230,1,0,0,0,0,ZERO,ZERO)
           MSG = 'tolerance this large is not appropriate for all problems.'
           CALL XERRDV(MSG,230,1,0,0,0,0,ZERO,ZERO)
           MSG = 'Execution will continue'
           CALL XERRDV(MSG,230,1,0,0,0,0,ZERO,ZERO)
        ELSE
           MSG = 'You must specify a vector of absolute error tolerances or'
           CALL XERRDV(MSG,240,1,0,0,0,0,ZERO,ZERO)
           MSG = 'a scalar error tolerance. It is recommended that you use'
           CALL XERRDV(MSG,240,1,0,0,0,0,ZERO,ZERO)
           MSG = 'a vector of absolute error tolerances.'
           CALL XERRDV(MSG,240,2,0,0,0,0,ZERO,ZERO)
        END IF
      END IF

!       ITOL error tolerance flag.
!          ITOL   RTOL     ATOL            EWT(i)
!            1   scalar   scalar  RTOL*ABS(Y(i)) + ATOL
!            2   scalar   array   RTOL*ABS(Y(i)) + ATOL(i)
      IF (PRESENT(ABSERR_VECTOR)) THEN
         OPTS%ITOL = 2
      ELSE
         OPTS%ITOL = 1
      END IF
      RETURN

END FUNCTION SET_NORMAL_OPTS
!_______________________________________________________________________

FUNCTION SET_INTERMEDIATE_OPTS(DENSE_J, BANDED_J, SPARSE_J,          &
  USER_SUPPLIED_JACOBIAN,                                            &
  LOWER_BANDWIDTH, UPPER_BANDWIDTH,                                  &
  RELERR, ABSERR, ABSERR_VECTOR,                                     &
  TCRIT, H0, HMAX, HMIN, MAXORD, MXSTEP, MXHNIL,                     &
  NZSWAG, USER_SUPPLIED_SPARSITY, MA28_RPS,                          &
  NEVENTS, CONSTRAINED, CLOWER, CUPPER, CHANGE_ONLY_f77_OPTIONS)     &
RESULT(OPTS)

! FUNCTION SET_INTERMEDIATE_OPTS:
!    Jacobian type:
!       DENSE_J, BANDED_J, SPARSE_J
!    Analytic Jacobian:
!       USER_SUPPLIED_JACOBIAN
!    If banded Jacobian:
!       LOWER_BANDWIDTH, UPPER_BANDWIDTH
!    Error tolerances:
!       RELERR, ABSERR, ABSERR_VECTOR
!    VODE.f77 optional parameters:
!       TCRIT, H0, HMAX, HMIN, MAXORD, MXSTEP, MXHNIL
!    Sparse flags:
!       NZSWAG, USER_SUPPLIED_SPARSITY, MA28_RPS
!    Rootfinding:
!       NEVENTS
!    Impose bounds on solution:
!       CONSTRAINED, CLOWER, CUPPER
!    Change one or more of HMAX, HMIN, TCRIT, MXSTEP, MXHNIL, MAXORD
!       CHANGE_ONLY_f77_OPTIONS
! RESULT(OPTS)

! SET_OPTIONS sets user parameters for DVODE via keywords. Values
! that are defined herein will be used internally by DVODE.
! All option keywords are OPTIONAL and order is not important.

! Note that DVODE_F90 requires that one of SET_NORMAL_OPTS or
! SET_INTERMEDIATE_OPTS or SET_OPTS is called before the first
! time DVODE_F90 is called.

!                     Quick Summary of Options

! DENSE_J                - Jacobian is sparse alternative to MF
! BANDED_J               - Jacobian is banded alternative to MF
! SPARSE_J               - Jacobian is sparse alternative to MF
! USER_SUPPLIED_JACOBIAN - User will supply Jacobian subroutine
!                          (user supplied subroutine JAC required)
! LOWER_BANDWIDTH        - Lower bandwidth ML in DVODE
! UPPER_BANDWIDTH        - Upper bandwidth MU in DVODE
! RELERR                 - Scalar relative error tolerance in DVODE
! ABSERR                 - Scalar absolute error tolerance in DVODE
! ABSERR_VECTOR          - Vector absolute error tolerance in DVODE
! TCRIT                  - Critical time TCRIT in DVODE
! H0                     - Starting step size in DVODE
! HMAX                   - Maximum step size in DVODE
! HMIN                   - Minimum step size in DVODE
! MAXORD                 - Maximum integration order in DVODE
! MXSTEP                 - Maximum number of integration steps
!                          in DVODE
! MXHNIL                 - Maximum number of T+H=T messages in DVODE
! NZSWAG                 - guess for the number of nonzeros in sparse
!                          Jacobian
! USER_SUPPLIED_SPARSITY - user will supply sparsity structure
!                          arrays by calling USERSETS_IAJA
! MA28_RPS               - Redo MA28AD pivot sequence if a singularity
!                          is encountered
! NEVENTS                - number of user defined root finding functions
!                          (user supplied subroutine G required)
! CONSTRAINED,           - array of solution component indices that
! CLOWER,                  are to be constrained by the lower and
! CUPPER                   upper bound arrays CLOWER and CUPPER so that
!                          CLOWER(I) <= Y(CONSTRAINED(I)) <= CUPPER(I)
!                     Options Types
! DENSE_J                - logical
! BANDED_J               - logical
! SPARSE_J               - logical
! USER_SUPPLIED_JACOBIAN - logical
! LOWER_BANDWIDTH        - integer
! UPPER_BANDWIDTH        - integer
! RELERR                 - real(wp) scalar
! ABSERR                 - real(wp) scalar
! ABSERR_VECTOR          - real(wp) vector
! TCRIT                  - real(wp) scalar
! H0                     - real(wp) scalar
! HMAX                   - real(wp) scalar
! HMIN                   - real(wp) scalar
! MAXORD                 - integer
! MXSTEP                 - integer
! MXHNIL                 - integer
! NZSWAG                 - integer
! USER_SUPPLIED_SPARSITY - logical
! MA28_RPS               - logical
! NEVENTS                - integer
! CONSTRAINED            - integer array
! CLOWER                 - real(wp) array
! CUPPER                 - real(wp) array
! CHANGE_ONLY_f77_OPTIONS- logical
! Argument list parameters:
! ABSERR                   = scalar absolute error tolerance
! ABSERR_VECTOR            = vector of absolute error tolerances
! RELERR                   = scalar relative error tolerance
! NEVENTS                  = Number of event functions (requires
!                            user-supplied GFUN)
! DENSE_J                  = use dense linear algebra if .TRUE.
! BANDED_J                 = use banded linear algebra if .TRUE.
!   LOWER_BANDWIDTH        = lower bandwidth of the Jacobian
!                            (required if BANDED_J = .TRUE.)
!   UPPER_BANDWIDTH        = upper bandwidth of the Jacobian
!                            (required if BANDED_J = .TRUE.)
! SPARSE_J                 = use sparse linear algebra if .TRUE.
!
!   NZSWAG                 = If you wish you may supply a guess,
!                            NZSWAG, at the number of nonzeros
!                            in the Jacobian matrix. In some cases
!                            this will speed up the determination
!                            of the necessary storage.
!   USER_SUPPLIED_SPARSITY = .TRUE. if you wish to supply the sparsity
!                            structure arrays directly by calling
!   MA28_RPS               = .TRUE. to force MA28AD to calculate a new
!                            pivot sequence for use in MA28BD if a
!                          
