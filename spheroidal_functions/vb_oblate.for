	  module vb_oblate
	  use regime
	  
	  private
	  
	  public coblfcn
	  
	  contains
      subroutine coblfcn(c,m,lnum,ioprad,x,r1c,ir1e,r1dc,ir1de,r2c,
     1                   ir2e,r2dc,ir2de,iopang,iopnorm,narg,arg,s1c,
     2                   is1e,s1dc,is1de, enr, maxd)
c
c  subroutine version of the fortran program coblfcn
c
c      version 1.03
c      September 2020
c
c  developed by arnie lee van buren
c  www.mathieuandspheroidalwavefunctions.com
c
c  purpose:     To calculate the first and second kind oblate radial
c               functions r1 and r2 and their first derivatives with
c               respect to the shape parameter x for a given value of
c               the complex size parameter c and the order m and for a
c               specified number lnum of values for the degree l = m,
c               m+1, ...,m+lnum-1.
c               To calculate the first kind oblate angular functions
c               s1 and their first derivatives with respect to the angle
c               coordinate eta for a given value of c and m and for a
c               specified number of degrees l = m, m+1,..., m+lnum-1
c
c  This program for complex c is a major advance in capability
c  from coblfcn version 1.01 and an improvement over version 1.02.
c  It is written in fortran 90. It was developed around oblfcn
c  version 1.24 which obtains spheroidal function values for real
c  values of the size parameter c = kd/2, where k is the wavenumber
c  and d is the interfoocal distance of the elliptical cross-section
c  of the oblate spheroid. A description of the methods used in
c  oblfcn is provided in the article 'Accurate calculation of oblate
c  spheroidal wave functions,' available at arXiv.org, identifier
c  1708.07929, August 2017, revised September 2019. Coblfcn provides
c  function values for c complex = real(c) + i aimag(c) = cr + i ci,
c  where the imaginary part ci often accounts for losses in wave
c  propagation. Ci is assumed positive in coblfcn. If the user has
c  a negative value for ci, just run coblfcn with ci positive instead
c  and take the complex conjugate of the results, including the function
c  values, eigenvalues, expansion coefficients, and normalization
c  factors.
c
c  A description of the methods used in coblfcn as well as a brief
c  description of this program is provided in the article 'Calculation
c  of oblate spheroidal wave functions with complex argument,' available
c  at arXiv.org, identifier 2009.01618, August 2020.
c
c  Coblfcn can be run in either double precision [64 bit; real(8) and
c  complex(8)], quadruple precision [128 bit; real(16) and complex(16)]
c  arithmetic, or a hybrid where everything except the Bouwkamp
c  procedure to refine the eigenvalue is run in double precision
c  arithmetic. In the latter case, coblfcn switches to quadruple
c  precision for the Bouwkamp procedure in those cases where double
c  precision fails to provide highly accurate eigenvalues. See the
c  discussion below about when to choose which arithmetic. The choice
c  of arithmetic is set in the module param located at the end of
c  coblfcn. Here, the kind parameters knd and knd1 are set by the
c  statements:
c      integer, parameter :: knd = selected_real_kind(8)
c      integer, parameter :: knd1 = selected_real_kind(16)
c  Set the value of knd in the parenthesis to either 8 for 64 bit
c  arithmetic or to 16 for 128 bit arithmetic. knd controls the
c  arithmetic for all but the Bouwkamp procedure
c  Set the value of knd1 to be used for the Bouwkamp procedure.
c  Note that if knd = 16, knd1 should also = 16.
c
c  Some computers may have more than 8 bytes for double precision
c  data and more than 16 bytes for quadruple precision data. In this
c  case just use the appropriate integers for the kind parameters in
c  module param. Also change the values of kindd and kindq set in
c  statement 5 below the comments section to the number of bytes for
c  double precision data and quadruple precision data, respectively.
c  Increasing the values of kindd and kindq will improve the accuracy
c  and can significantly extend the paramteter ranges over which
c  coblfcn provides useful results.
c
c  Following is a description of input and output. Included is a
c  discussion of the output file fort.60 that alerts the user when the
c  accuracy of the radial and/or angular functions fall below a
c  specified minimum value and when the eigenvalue routine fails. Next
c  is a discussion of the ranges of parameters for which coblfcn
c  provides useful results. Following that is a discussion of how to
c  obtain the expansion d coefficients and the various angular function
c  normalizations. Finally is a discussion of four optional output files
c  fort.20, fort.30, fort.40, and fort.50 that may be of interest.
c
c     Input and Output
c
c    Input and output parameters from the subroutine call statement are
c    defined below:
c
c          c      : desired complex value of the size parameter (= kd/2,
c                   where k is the complex wavenumber and d is the
c                   interfocal length) (either complex*16 or complex*32)
c
c          m      : desired value for the order m (integer)
c
c          lnum   : number of values desired for the degree l (integer)
c                   if lnum is less than 2*(real(c)+aimag(c))/pi it
c                   should be an even integer
c          ioprad : (integer)
c                 : =0 if radial functions are not computed
c                 : =2 if radial functions of both kinds and
c                      their first derivatives are computed
c
c          x      : value of the radial coordinate x (a nominal value of
c                   10.0d0 or 10.0q0 can be entered for x if ioprad
c                   = 0) (either real*8 or real*16)
c
c          r1c   :  either complex*16 or complex*32 vectors of length
c          r1dc     lnum containing the characteristics for the radial
c                   functions of the first kind r1 and their
c                   first derivatives
c
c          ir1e   : integer vectors of length lnum containing the
c          ir1de    exponents corresponding to r1c and r1dc
c
c          r2c    : complex*16 or complex*32 vectors of length lnum
c          r2dc     containing the characteristics for the radial
c                   functions of the second kind r2 and their first
c                   derivatives
c
c          ir2e   : integer vectors of length lnum containing the
c          ir2de    exponents corresponding to r2c and r2dc
c
c          iopang : (integer)
c                 : =0 if angular functions are not computed
c                 : =1 if angular functions of the first kind
c                      are computed
c                 : =2 if angular functions of the first kind and
c                      their first derivatives are computed
c
c          iopnorm: (integer)
c                 : =0 if not scaled. The angular functions have
c                      the same norm as the corresponding associated
c                      legendre function [i.e., we use the Meixner
c                      -Schafke normalization scheme.] This norm
c                      becomes very large as m becomes large. The
c                      angular functions are computed below as
c                      a characteristic and an exponent to avoid
c                      overflow.
c                 : =1 if angular functions of the first kind
c                      (and their first derivatives if computed)
c                      are scaled by the square root of the
c                      normalization of the corresponding
c                      associated Legendre function. The resulting
c                      scaled angular functions have unity norm.
c                      This is very useful since it removes the
c                      need to calculate a normalization factor
c                      when using the angular function values given
c                      here. It also eliminates any chance for
c                      overflow when the characteristics and exponents
c                      are combined to form the angular functions.
c
c          narg   : number of values of the angular coordinate eta for
c                   which angular functions are calculated (integer)
c
c          arg:     vector containing the values of eta for which
c                   angular functions are desired (real*8 or real*16)
c
c          s1c,   : two-dimensional arrays s1c(lnum,narg) and
c          s1dc     s1dc(lnum,narg) that contain narg calculated
c                   complex characteristics for the angular functions
c                   and their first derivatives for each of the lnum
c                   values of l (complex*16 or complex*32)
c                   For example, s1c(10,1) is the characteristic
c                   of the angular function for l = m + 10 -1 and
c                   the first value of eta given by arg(1)
c
c          is1e,  : integer arrays is1e(lnum,narg) and is1de(lnum,narg)
c          is1de    containing the exponents corresponding to s1c and
c                   s1dc
c
c
c     Output file fort.60
c
c  Coblfcn provides a file fort.60 that should be of interest to the
c  user, especially when using coblfcn outside the ranges where useful
c  results are expected. See the discussion of these ranges below.
c  Whenever the estimated accuracy of the radial functions falls below a
c  designated integer value, the associated values of x, c, m, and l are
c  written to fort.60. The integer is currently set equal to 6 in write
c  statements found after the line numbered 1405 below in subroutine
c  main. I expect that accuracies as low as 5 or 6 digits and possibly
c  even as low as 4 digits may be sufficient in many cases, especially
c  when they occur very rarely. The integer can be set higher if desired.
c
c  Coblfcn will also write to fort.60 whenever the estimated accuracy
c  of the Meixner and Schafke normalization for the angular functions
c  is less than the same integer used for alert about tne radial
c  function accuracy. It provides the values of m, l, and c where this
c  occurs together with the value ndec - jsubms -1, where ndec is the
c  precision and jsubms is the subtraction error in calculating the
c  normalization.
c
c  Fort.60 also alerts the user to the unlikely situation where the
c  eigenvalue procedure has failed for one or more values of l and the
c  same eigenvalue is obtained for two different values of l (of the
c  same parity). In this case one or more eigenvalues are missing and
c  the results are not usable. I do not expect this to be a possibility
c  unless coblfcn is used with ci and cr beyond the ranges stated below
c  for which useful results are expected.
c
c  If the user does not desire fort.60, it can be avoided by searching
c  for write(60 and placing a c in column 1 to commenting out the write
c  statement that are found. The open statement can also be commented
c  out.
c
c     Useful parameter ranges for coblfcn
c
c  The following discussion is provided to help the user choose which
c  arithmetic option to use. If the compiler does not support quadriple
c  arithmetic, then the only option is to use double precision. If
c  the compiler does support quadruple precision arithmetic, then it is
c  recommended that the hybrid version be used instead of the entirely
c  double precision version whenever double precision is expected to
c  provide adequate accuracy. Use of quadruple precision arithmetic
c  for the Bouwkamp procedure increases the run time somewhat but the
c  program is still reasonably fast and the results are more accurate.
c  The improvement in accuracy increases as ci increases. If the input
c  parameters are outside those for which the hybrid version is expected
c  to provide useful results, then the only choice is to use quadruple
c  precision if available. Note that this increases the run time by a
c  factor of 50 or more.
c
c  Coblfcn was tested extensively using a laptop pc and a Fortran
c  compiler that provides approximately 15 decimal digits in double
c  precision arithmetic and approximately 33 digits in quadruple
c  precision arithmetic. The estimated accuracy of the resulting
c  function values is given below in terms of decimal digits. It is
c  often the Wronskian comparison result but can be an estimate
c  based on subtraction errors in the calculation. See the discussion
c  below for the accuracy integer naccr in the section for fort.20.
c  If the user's computer provides a different number of digits, the
c  following estimates should be adjusted up or down depending on
c  whether more or fewer digits are provided. Testing included values
c  of x ranging from 0.000001 to 10 as well as the special case of
c  x = 0, values for cr up to 5000, and values of ci up to 200. Testing
c  for both the double precision and hybrid versions included all values
c  of the order m from 0 to 200 and from 210 to 1000 in steps of 10.
c  Testing for the quadruple precision version included values of m from
c  0 to 200 in steps of 10 and from 250 to 1000 in steps of 50. For all
c  versions the values of the degree l ranged from m to m + lnum -1,
c  where lnum was chosen sufficiently large that the magnitudes of r1
c  and r1d were less than 10**(-300).
c
c  An integer called minacc is used in coblfcn. This designates the
c  minimum number of accurate digits desired for the radial spheroidal
c  functions of the second kind. The value of minacc controls which
c  methods are used to calculate the radial functions. Minacc is set
c  equal to 8 for 64 bit arithmetic. It is recommended that this not
c  be changed. For 128 arithmetic minacc is set equal to 15 digits
c  for values of ci up to 20. This should provide 15 or more digits of
c  accuracy here. For ci > 20, minacc is set equal to 8 digits. Minacc
c  can be increased in this case but higher accuracy might not always
c  be achieved. Also the computation time will likely go up with larger
c  values of minacc. Minacc is set below following these introductory
c  comment statements.
c
c  In the following discussion the term useful results means that the
c  estimated accuracy for the radial functions observed during testing
c  never fell below 5 decimal digits unless otherwise stated. I expect
c  there are many applications where occasional 5 digit results are
c  acceptable. Possibly even an isolated 4 digit result is acceptable.
c  Note that there is no guarantee that the estimated accuracy for the
c  values for parameter values other than those that I tested will be
c  as high as I report below. The discussion below will focus on x
c  unequal to zero. It is expected that function values for x = 0 will
c  be useful for the same values of c and m that useful results are
c  obtained for at least one value of x. Here both the radial function
c  of the first kind r1 and its first derivative will be very accurate.
c  Whenever values of the radial function of the second kind r2 for
c  l - m even are less accurate than 5 digits, they are expected to be
c  proportionally smaller in magnitude than r1. Similarly for the first
c  derivatives of r2 and r1 when l - m is odd.
c
c  Double precision arithmetic
c
c  Using double precision arithmetic, coblfcn provides useful results
c  for ci (the imaginary part of c) as large as 10, for cr (the real
c  real part of c) as large as 5000, for m up to at least 1000 and for
c  all tested values of x down to 0.000001. When ci is less than about
c  5, coblfcn provides results with an accuracy similar to that
c  provided by the program oblfcn for real c. Extensive testing for
c  ci = 10 showed that the radial functions of the first kind r1 and
c  its first derivative r1d are almost always accurate to 10 or more
c  decimal digits. For all values of x except zero the radial functions
c  r2 and r2d are usually accurate to 8 or more decimal digits but
c  accuracies lower than this were seen, especially for larger cr and
c  smaller x. Nearly all accuracies less than 8 digits occurred near but
c  somewhat below the so-called breakpoint. [The breakpoint is defined
c  as the minimum value of l such that the magnitudes of r2 and r2d
c  begin to increase with increasing l while the magnitudes of r1 and
c  r1d begin to decrease. An approximate value for the breakpoint for
c  for small values of m is given by l = 2*(cr+ci)/pi]. No 5 digit
c  results were seen for cr up to 200 or so. Only a few 5 digit results
c  were seen for x >= 0.001, even for cr = 5000. The largest number of 5
c  digit results occurred for cr = 5000 and x = 0.000001. Even here,
c  however, there were no more than about 3 such results for each m. See
c  the discussions below about the estimated accuracy of the angular
c  functions.
c
c  Similar testing When ci was increased to 12 showed a few more 5 digit
c  results for small x. At least 5 digits are obtained for cr up to 5000
c  for all of the values of x down to 0.000001.
c
c  Testing with ci = 15 showed yet more 5 digits results and even 4
c  digit results for cr = 5000. However, at least 5 digits of accuracy
c  were obtained for cr up to 2000 for all values of x down to 0.000001.
c
c  When ci = 20, there are 5 or more digits of accuracy for all tested
c  values of x when cr <= 100. There are 5 or more digits of accuracy
c  at x >= 0.1 for cr = 150 and at x >= 0.2 for cr up to 2000. Testing
c  showed duplicated eigenvalues of the same parity for some values of
c  m when cr = 5000.
c
c  When ci = 25, there were at least 5 digits of accuracy for x >= 0.2
c  when cr <= 100, for x >= 0.3 when cr = 150 and for x >= 0.4 when cr
c  = 200 and 250.
c
c  Testing for yet higher values of ci showed a continued increase in
c  the minimum value of x and a decrease in the maximum value of cr for
c  which useful function values were obtained. If the user is interested
c  in values of x somewhat larger than 0.3 together with moderate to
c  small values of cr, then the double precision, or more likely the
c  hybrid version, of coblfcn may be useful when ci is greater than 25.
c  Otherwise, it will be necessary to use the quadruple precision
c  version. It is recommended that the file fort.60 described below be
c  used to assure that you are obtaining the accuracy you need and that
c  there are no repeated eigenvalues of the same parity. Testing showed
c  the appearance of duplicated eigenvalues for ci = 45 for some values
c  of m when cr was only 200.
c
c
c  Quadruple precision arithmetic
c
c  When lower accuracy occurs with 64 bit arithmetic, much higher
c  accuracy can be obtained using 128 bit arithmetic. However, coblfcn
c  runs faster by a factor up to 50 or more in 64 bit arithmetic
c  than it does for 128 bit arithmetic. Running coblfcn on an ordinary
c  laptop computer can take a long time using 128 bit arithmetic when
c  cr is very large.
c
c  Testing for values of ci up to 40 showed that coblfcn provides useful
c  results for cr up to at least 5000 and for m up to 1000 and x down to
c  0.000001. Estimated accuracies for the radial functions were 8 or
c  more digits except for a possible rare 6 or 7 digit result. See the
c  discussions of estimated accuracies given below for the special case
c  of x = 0 and for the angular functions.
c
c  Testing for ci = 50 shows useful results for cr up to at least 2000.
c  Accuracies for the radial functions were almost always at least 8
c  more digits but occasional accuracies as low as 5 digits were seen
c  near the breakpoint, primarily for cr >= 1000 and m > 200. When cr
c  = 2000, there were even a few 4 digit results for m >= 700 and a few
c  3 digit results for m >= 800. It is unlikely that values of m this
c  large will be required by the user.
c
c  Testing for ci = 60 shows useful results for cr <= 100 although a
c  few 4 digit results occurred at x = 0.00001 and 0.000001. Results
c  for cr = 200 were similar except that there were a few more 4 digit
c  results at x = 0.00001 and 0.00001 and some 4 digit results now at
c  x = 0.0001 and x = 0.001. Testing for cr = 500 shows that useful
c  results are only obtained for x >= 0.05 for all m and for x < 0.05
c  for m up to about 160. Testing for cr = 1000 and 2000 showed useful
c  results for x >= 0.2 for all m with a possible rare 4 digit result.
c
c  Testing for ci = 70 shows useful results for cr <= 100 when x >=
c  0.01, for cr = 150 when x >= 0.05, for cr = 200 when x >= 0.1, for
c  cr = 500 when x >= 0.3 and for cr = 1000 when x >= 0.5.
c
c  Testing for ci = 80 shows useful results for cr <= 20 when x >=
c  0.01, for cr = 50 when x >= 0.02, for cr = 100 when x >= 0.2,
c  for cr = 200 when x >= 0.3, for ci = 300 when x >= 0.5, and for
c  ci = 400 when x >= 1.1.
c
c  Testing for yet higher values of ci showed a continued increase in
c  the minimum value of x and a decrease in the maximum value of cr
c  for which useful function values were obtained. It is recommended
c  that the user, especially for ci greater than 80, use the file
c  fort.60 described below to assure that the desired accuracy is being
c  obtained and that there are no repeated eigenvalues of the same
c  parity.
c
c  Sometimes either r2 or r2d will have lower accuracy than the
c  estimate provided by coblfcn. This can happen for large values of
c  cr, values of x less than about 0.1 and a few values of l - m
c  somewhat below or near the so-called breakpoint. Here r2 for l - m
c  even or r2d for l - m odd may be less accurate. However, these
c  function values are usually somewhat smaller in magnitude than the
c  corresponding values for r1 or r1d and make a reduced contribution
c  to the numereical solution of physical problems utilizing oblate
c  spheroidal functions.
c
c     Estimated angular function accuracy
c
c  For both choices of arithmetic, the angular functions and their first
c  derivatives can lose accuracy for low l, high c, and eta near zero
c  due to subtraction errors in their series calculation. However,
c  their magnitude in this case is corresponding smaller than angular
c  functions for higher values of l and/or eta not near zero. The loss
c  in accuracy due to these subtraction errors should not adversely
c  effect numerical results for physical problems using these functions.
c
c  A second source of inaccuracy in the angular functions arises from
c  subtraction errors that can occur in the calculation of their
c  Meixner-Schafke normalization factor dmsnorm. Note that the loss in
c  accuracy here is not in addition to other losses in accuracy for the
c  angular functions but rather sets an upper limit to their accuracy.
c  These errors are largest for m = 0. They occur for values of l - m
c  somewhat less than the breakpoint and grow with increasing ci and to
c  some extent with increasing cr. They are zero for small ci and can
c  become as large as 6 digits for ci = 20 and 12 digits for ci = 50 as
c  cr increases to 5000. This loss in accuracy is not likely a problem
c  when using double precision arithmetic with 15 decimal digits since
c  as ci becomes larger than 20, the values of cr for which the radial
c  functions are accurate to 5 or more digits are progressively smaller,
c  being 100 for ci = 25 and less than this for higher ci. Using double
c  precision, the Meixner and Schafke normalization should be accurate
c  to at least 5 digits wherever all of the radial functions for a
c  given value of m are also accurate to at least 5 digits. The file
c  fort.60 mentioned above can also alert the user whenever the
c  estimated accuracy for the Meixner and Schafke normalization is less
c  than the same integer number of decimal digits selected for alerts
c  about the accuracy of the radial functions.
c
c  For higher values of ci when using quadruple precision, the loss of
c  accuracy in the normalization factor is even greater. For ci = 60,
c  the loss of accuracy can be as large as 24 digits for cr = 2000 and
c  25 digits for cr = 5000. For ci = 70 it is 26 digits for cr = 1000.
c  And for ci = 80 it is 25 digits for cr = 400. This should not be a
c  problem using 33 decimal digits since it still allows for accuracies
c  of at least 5 digits for the angular functions everywhere the radial
c  functions also have an accuracy of 5 or more digits.
c
c  A third source of inaccuracy in the angular functions arises from the
c  potential loss of accuracy in the eigenvalues at values of l near and
c  somewhat below the breakpoint when ci is not very small and ci is
c  moderate to large. This ia most likely to occur when using double
c  precision for the calculations including for the Bouwkamp procedure.
c  Use of quadruple precision, if available, for the Bouwkamp procedure
c  will help considerably here.
c
c     d Coefficients
c
c  The user may desire values for the d coefficients that appear in
c  the expression for the angular functions as well as in many of the
c  expressions used to calculate the radial functions. Ratios of
c  successive d coefficients are stored in the vector enr where enr(k)
c  = d(subscript 2k+ix) divided by d(subscript 2k-2+ix). The vector enr
c  is calculated in the subroutine dnorm in statement 20 and passed to
c  subroutine main. The number lim2 of d coefficients calculated for a
c  given l is chosen to be sufficient to compute radial and angulalar
c  functions for that l. The number of d coefficients needed to compute
c  r1 and r1d and s1 and s1d can range from less than 100 to somewhat
c  more than l/2 for large values of l. The number of d coefficients
c  required to compute r2 and r2d are comparable to this unless they are
c  computed using one of the Neumann function expansions. Here one can
c  require up to 200000 or so coefficients when x is near 0.01. Note
c  that the vector enr returned by the subroutine conver contains scaled
c  ratios where the scaling has been chosen to produce a symmetric
c  matrix for computing eigenvalues. The scaling factors are removed in
c  subroutine dnorm to obtain the desired d coefficient ratios.
c
c  The d coefficients themselves can be obtained starting with the value
c  for d with the subscript l - m. If iopnorm is set = 0, coblfcn uses
c  the Meixner_Schafke normalization scheme for the angular functions.
c  Here the angular functions have the same norm as the corresponding
c  associated Legendre functions. When c is real the calculation of the
c  normalizing factor is accurate with no associated subraction errors.
c  As ci increases and to a some extent as cr increases, the subtraction
c  errors increase for some values of l. They are near zero for ci up to
c  10, then increase to 6 digits at ci = 20 and to digits at ci = 50.
c  See the discussion of this in the section 'Estimated angular function
c  accuracy' given above. The subroutine dnorm computes d(subscript l-m)
c  for this normalization and returns it to subroutine main as a
c  characteristic dmlms and an exponent idmlmse. Use of an exponent
c  avoids possible overflow of d(subscript l-m) for extremely large c
c  and m. When the user sets iopnorm = 1 so that the angular functions
c  have unit norm, the corresponding characteristic and exponent for
c  d(subscript l-m) are calculated in subroutine s1 and returned to
c  subroutine main as dmlms1 and idmlms1e. Corresponding values for the
c  characteristic and exponent of d(subscript l-m) for the  Morse-
c  Feshbach and Flammer normalizations are computed in dnorm and
c  returned to main as dmlmf, idmlmfe and dmlf, idmlfe. Note that for
c  c complex, all three of these normalization calculations suffer
c  subtraction errors for lower values of l-m and non-small c. The
c  values for d(subscript l-m) will have reduced accuracy in this
c  case. When c is real only the Flammer normalization suffers
c  subtraction errors in its calculation.
c
c     Optional output files
c
c  Output data can be written to the following four files: fort.20,
c  fort.30, fort.40, fort.50. These are currently suppressed. If
c  desired, one or more of these files, say fort.40, can be obtained
c  by searching for write(40 statements and removing the c given in
c  column 1 of each of them as well as their associated format
c  statements. Also remove the c from the open(40 statement near the
c  beginning of the program. Note that the open statements allow
c  default values to apply for the specs. These can, of course, be
c  changed.
c
c   fort.20
c
c     This file contains values for all radial functions that have
c     been calculated.
c     The first line in the file contains the values for x, c, and
c     m, formatted as follows (see statements 260 and 265 in
c     subroutine main):
c
c                x      : e23.14 in 64 bit arithmetic; e38.30
c                       : in 128 bit arithmetic
c                c      : (e23.14,e23.14) in 64 bit arithmetic;
c                       : (e38.30,e38,30) in 128 bit arithmetic
c                m      : i5
c
c     Each subsequent line in fort.20 contains radial functions
c     for given values of l. The first line contains values for l = m,
c     the next for l=m+1 and continuing to l=m+lnum-1. The radial
c     functions are preceeded by the value of l and followed by the
c     accuracy, equal to the estimated number of accurate decimal digits
c     in the radial functions as measured using the Wronskian (except
c     for (1) the case where x is equal to zero, (2) the case where
c     the Wronskian is used to obtain a leading coefficient for the
c     Legendre function expansion method and (3) the case where both
c     the radial functions of the first kind R1 and those of the
c     second kind R2 are very large such that r2 and r2d are very nearly
c     equal to -i times r1 and -i times r1d, respectively. [see comments
c     below regarding naccr]
c
c       The output and corresponding format for each line is as follows
c       (see statements 1340, 1350, and 1380 in main).
c
c         l            : value for l (i5)
c         r1c(l-m+1)   : complex characteristic of the oblate radial
c                        function of the first kind r1 for
c                        the given value of l (f17.14,f17.14)
c         ir1e(l-m+1)  : exponent of r1 (i6)
c         r1dc(l-m+1)  : complex characteristic of the first derivative
c                        of the oblate radial function of the
c                        first kind r1d (f17.14,f17.14)
c         ir1de(l-m+1) : exponent of r1d (i6)
c         r2c(l-m+1)   : complex characteristic of the oblate radial
c                        function of the second kind r2 for
c                        the given value of l (f17.14,f14.7)
c         ir2e(l-m+1)  : exponent of the oblate radial function of
c                        second kind (i6). If the exponent for any
c                        of the functions is greater than 9999,
c                        the format can be increased to i7 or
c                        higher. Note that the procedures used in
c                        this program allow for exponents much
c                        larger than those allowed in 64 bit
c                        arithmetic on the users computer since
c                        the floating point function values
c                        provided are given as a characteristic
c                        and an integer exponent. Use of ratios in
c                        calculating the functions eliminates
c                        overflow and underflow during the
c                        calculations.
c         r2dc(l-m+1)  : complex characteristic of the first derivative
c                        of the oblate radial function of second kind
c                        r2d (f17.14,f17.4)
c         ir2de(l-m+1) : exponent of r2d (i6). [See comment above
c                        for ir2e.]
c         naccr(l-m+1) : estimated accuracy: usually equal to the
c                        number of decimal digits of agreement
c                        between the theoretical Wronskian and the
c                        calculated Wronskian (i2). This is
c                        indicated by the letter w following the
c                        integer naccr. When ci is small, r1 and r1d
c                        are very accurate so that naccr usually
c                        relects the accuracy of r2 and r2d.
c
c                        Several situations are described below where
c                        the Wronskian can not be used to estimate the
c                        accuracy. Other factors are used here to
c                        estimate naccr. This is indicated by using
c                        the letter e instead of w following the value
c                        for naccr in fort.20.
c
c                        Sometimes near the breakpoint value for l,
c                        the Legendre function expansions for r2 and r2d
c                        converge with acceptable subtraction error but
c                        the leading coefficient for the series with the
c                        P Legendre function is highly inaccurate. Here
c                        the Wronskian can often be used for small ci to
c                        obtain improved accuracy for the coefficient.
c                        The accuracy is estimated using subtraction
c                        errors involved in the calculations and the
c                        estimated accuracy of the eigenvalue in those
c                        cases where the Bouwkamp eigenvalue routine
c                        does not converge fully.
c
c                        When ci becomes large, there are often values
c                        of l where both r1 and r2 and their first
c                        derivatives are large in magnitude. Here values
c                        for r2 and r2d are given by -i times values of
c                        r1 and r1d. The accuracy is estimated using
c                        the magnitudes of r1 and r1d, the magnitude
c                        of the theoretical value of the wronskian, and
c                        the estimated accuracy of r1 and r1d. See the
c                        discussion of this below.
c
c                        When x is equal to zero, the Wronskian is used
c                        to obtain r2d when l-m is even and r2 when l-m
c                        is odd. Here the accuracy naccr is approximated
c                        using the estimated accuracy of the eigenvalue
c                        and the estimated accuracy of the joining
c                        factor used in using the Legendre function
c                        expansion to compute r2 and r2d. The accuracy
c                        estimate for the other functions is given by
c                        the estimated accuracy of the eigenvalue minus
c                        one digit. When r2 or r2d are less accurate due
c                        to large subtraction errors in the calculation
c                        of the joining factor, they are often somewhat
c                        small in magnitude, especially for small values
c                        of ci. Here the reduced accuracy is not likely
c                        to result in lower accuracy in solutions to
c                        problems utilizing these functions. When the
c                        estimated accuracy is zero, the functions are
c                        set equal to zero.
c
c                        There are two other situations where the value
c                        for naccr can be an estimated value: (1) When
c                        the function values for r2 and r2d are
c                        obtained using the integral expressions for
c                        small values of x, the accuracy given by the
c                        Wronskian value is adjusted downward to account
c                        for the likelihood that the accuracy for r2 for
c                        l - m even and r2d for l - m odd is lower than
c                        the Wronskian estimate. (2) The accuracy using
c                        the traditional Legendre expansion is also
c                        estimated and the accuracy estimate provided
c                        by coblfcn is taken to be the smaller of this
c                        and the Wronskian estimate. In both cases, the
c                        accuracy is designated by w rather than e since
c                        the Wronskian accuracy sets an upper bound.
c
c   fort.30
c
c     This file fort.30 contains values for all angular functions
c     that have been calculated. Its first line contains the values
c     for c and m, formatted as follows (see statements 50 and 55 in
c     subroutine main).
c
c                c      : (e23.14,e23.14) in 64 bit arithmetic;
c                       : (e38.30,e38,30) in 128 bit arithmetic
c                m      : i5
c
c     The second line in fort.30 contains the value for the first l
c     (=m), formatted as follows (see statement 270 in subroutine
c     main):
c
c                l      : i6
c
c     This is followed by a series of narg lines. Each line contains
c     a desired value of the angular coordinate eta followed by the
c     corresponding angular functions and accuracy. Specific output
c     and format for each line is as follows.
c
c        for iopang = 1:
c
c               barg   ; angular coordinate eta
c                        (f17.14; see statement 1460 in subroutine main)
c               s1c    : complex characteristic of the oblate angular
c                        function of first kind (f17.14,f17.14); see
c                        statement 1460 in subroutine main)
c               is1e   : exponent of the oblate angular function of
c                        first kind (i5; see statement 1460 in
c                        subroutine main)
c
c        for iopang = 2, each line also includes:
c               s1dc   : complex characteristic of the first derivative
c                        of the oblate angular function of first kind
c                        (f17.14,f17.14); see statement 1470 in
c                        subroutine main)
c               is1de  : exponent of the first derivative of the
c                        oblate angular function of first kind (i5;
c                        see statement 1470 in subroutine main)
c
c        for iopang = 1:
c               naccs  : accuracy: estimate of the number of decimal
c                        digits of accuracy in the angular function
c                        (i2). It is a conservative estimate based on
c                        the calculated subtraction error in the
c                        Legendre function series for the angular
c                        functions and the estimated accuracy of the
c                        normaliation factor. When the accuracy estimate
c                        is equal to 0, the corresponding angular
c                        functions are set equal to zero. (i2; see
c                        statements 1460 and 1470 in subroutine main).
c                        The calculated angular functions tend to be
c                        less accurate the larger the value of c, the
c                        smaller the value of l - m (for values less
c                        than approximately 2c divided by pi), and the
c                        closer eta is to zero (i.e., the closer theta
c                        is to 90 degrees).
c
c        for iopang = 2:
c              naccds  : accuracy: includes an estimate of the number
c                        of decimal digits of accuracy in the first
c                        derivative of the angular function (i2)
c
c   fort.40 and fort.50
c
c     These files are diagnostic files that contain information about
c     specific techniques used and numbers of terms required for the
c     radial function and angular function calculations, respectively.
c     They are annotated and should be self explanatory.
c
c
c        use param
c

		complex(knd), allocatable, dimension(:,:) :: enr
c  real scalars and arrays
        real(knd) arg(narg),ca,step1,step2,step3,x,xneu
c
c  complex scalars and arrays
        complex(knd) c,r1c(lnum),r1dc(lnum),r2c(lnum),
     1               r2dc(lnum),s1c(lnum,narg),s1dc(lnum,narg)
c
c  integer arrays
        dimension ir1e(lnum),ir1de(lnum),ir2e(lnum),ir2de(lnum)
        dimension is1e(lnum,narg),is1de(lnum,narg)
c
c  open output files
c        open(20, file='fort.20')
c        open(30, file='fort.30')
c        open(40, file='fort.40')
c        open(50, file='fort.50')
        open(60, file='fort.60')
c
c  Here is where the user sets kindd, the number of bytes available
c  in double precision data for the computer that coblfcn is run on.
c  Similarly kindq, the number of bytes available in quadruple
c  precision data is set here. This allows coblfcn to be run on
c  computers with values for kindd and kindq different than 8 and 16,
c  respectively.
c
5       kindd=8
        kindq=16
c
c  Here the minimum desired accuray minacc is set to 8 decimal
c  digits for double precision arithmetic (knd = kindd).
c  The minimum desired accuracy for quadruple precision arithmetic
c  (knd = kindq) is set to 15 decimal digits when aimag(c) <= 20
c  and to 8 decimal digits when aimag(c) > 20.
c  These can be changed if desired. See comments
c  below about changing minacc
c
        if(knd.eq.kindd) minacc=8
        if(knd.eq.kindq.and.aimag(c).le.20.0e0_knd) minacc=15
        if(knd.eq.kindq.and.aimag(c).gt.20.0e0_knd) minacc=8
c
c     ndec: the maximum number of decimal digits available in real
c           arithmetic.
c     nex:  the maximum exponent available in real arithmetic.
c
        ca=abs(c)
        ndec=precision(ca)
        nex=range(ca)-1
c
c  set array dimensions
        mnum=1
        mmin=m
        minc=0
        nbp=int(2*(real(c)+aimag(c))/3.1416)
        maxe=max(50,nbp+30)+10
        maxe2=maxe+maxe
        if(maxe2.lt.lnum) maxe2=lnum
        maxm=mmin+minc*(mnum-1)
        maxlp=lnum+maxm+1
        maxint=2*(nbp+33)+50*aimag(c)+6
        maxj=lnum+3*ndec+int(ca)+5+maxm
        maxp=max(lnum+3*ndec+int(ca)+5,maxlp+5)
        maxn=maxj
        maxpdr=2*int(ca)+4*ndec+int(100*x)+8
        neta=30
        if(ioprad.ne.2) go to 10
          step1=0.1e0_knd
          step2=0.1e0_knd
          step3=0.8e0_knd
          nstep1=1
          nstep2=1
          nstep3=3
          ngau=100*(2+int(ca/500.0e0_knd))
c
          if(x.lt.0.2e0_knd) then
          step1=x/4.0e0_knd
          step2=0.075e0_knd
          step3=1.0e0_knd-step1-step2
          nstep1=1
          nstep2=1
          nstep3=4
          ngau=100*(2+int(ca/500.0e0_knd))
          end if
c
          if(x.lt.0.1e0_knd) then
          step1=x/4.0e0_knd
          step2=0.025e0_knd
          step3=1.0e0_knd-step1-step2
          nstep1=1
          nstep2=1
          nstep3=4
          ngau=200
          if(ca.gt.500.0e0_knd) ngau=200*(2+int(ca/1000.0e0_knd))
          end if
c
          if(x.lt.0.05e0_knd) then
          step1=x/15.0e0_knd
          step2=0.002e0_knd
          step3=1.0e0_knd-step1-step2
          nstep1=1
          nstep2=1
          nstep3=2
          ngau=300
          if(ca.gt.500.0e0_knd) ngau=500
          if(ca.gt.1000.0e0_knd) ngau=800
          if(ca.gt.1500.0e0_knd) ngau=1000
          if(ca.gt.2000.0e0_knd) ngau=1200
          if(ca.gt.2500.0e0_knd) ngau=1500
          if(ca.gt.3000.0e0_knd) ngau=1700
          if(ca.gt.3500.0e0_knd) ngau=1900
          if(ca.gt.4000.0e0_knd) ngau=2200
          if(ca.gt.4500.0e0_knd) ngau=2500
          if(ca.gt.5000.0e0_knd) ngau=2500+300*int((ca-4500.0e0_knd)/
     1       500.0e0_knd)
          end if
c
          if(x.lt.0.01e0_knd) then
          step1=x/15.0e0_knd
          step2=0.002e0_knd
          step3=1.0e0_knd-step1-step2
          nstep1=1
          nstep2=1
          nstep3=2
          ngau=300
          if(ca.gt.500.0e0_knd) ngau=600
          if(ca.gt.1000.0e0_knd) ngau=800
          if(ca.gt.1500.0e0_knd) ngau=1000
          if(ca.gt.2000.0e0_knd) ngau=300*int(ca/500.0e0_knd)
          end if
c
          if(x.lt.0.001e0_knd) then
          step1=x/15.0e0_knd
          step2=0.002e0_knd
          step3=1.0e0_knd-step1-step2
          nstep1=3
          nstep2=1
          nstep3=2
          ngau=600
          if(ca.gt.300.0e0_knd) ngau=800
          if(ca.gt.1000.0e0_knd) ngau=900
          if(ca.gt.1500.0e0_knd) ngau=1000
          if(ca.gt.2000.0e0_knd) ngau=300*int(ca/500.0e0_knd)
            if(aimag(c).gt.5.0e0_knd) then
            ngau=ngau+(ngau/5)*min(5,int(aimag(c))-5)
            end if
            if(knd.eq.kindq.and.x.lt.0.00001e0_knd) ngau=2*ngau
            if(knd.eq.kindq.and.x.lt.0.000001e0_knd) ngau=2*ngau
            if(knd.eq.kindd.and.aimag(c).gt.4.0e0_knd.and.x.lt.
     1         0.00001e0_knd) ngau=2*ngau
            if(knd.eq.kindd.and.aimag(c).gt.4.0e0_knd.and.x.lt.
     1         0.000001e0_knd) ngau=2*ngau
          end if
c
          xneu=0.3e0_knd
          if(ca.gt.100.0e0_knd) xneu=0.04e0_knd
          if(ca.gt.600.0e0_knd) xneu=0.03e0_knd
          if(ca.gt.800.0e0_knd) xneu=0.01e0_knd
          if(aimag(c).gt.50.0e0_knd) xneu=0.01e0_knd
c
        if(x.lt.0.01e0_knd.or.ioprad.eq.0) go to 10
          if(knd.eq.kindd) then
          if(x.ge.0.01e0_knd) maxn=2*int(25/(x*x)+300/x+3*ca+
     1                               1250*knd)+5
          if(x.ge.0.1e0_knd) maxn=2*int((lnum+ca/5+0.5e0_knd*maxm+
     1                             200)*1.4e0_knd/x)+5
          if(x.ge.0.5e0_knd) maxn=2*int((lnum+ca/5+0.5e0_knd*maxm+
     1                             300)/x)+5
          if(x.ge.1.0e0_knd) maxn=2*int(lnum+ca/5+0.5e0_knd*maxm+
     1                             300)+5
          end if
          if(knd.eq.kindq) then
          if(x.ge.0.01e0_knd) maxn=2*int(25/(x*x)+400/x+3*ca+
     1                               1250*knd)+5
          if(x.ge.0.1e0_knd) maxn=2*int((lnum+ca/5+0.5e0_knd*maxm+
     1                             350)*1.4e0_knd/x)+5
          if(x.ge.0.5e0_knd) maxn=2*int((lnum+ca/5+0.5e0_knd*maxm+
     1                             400)/x)+5
          if(x.ge.1.0e0_knd) maxn=2*int(lnum+ca/5+0.5e0_knd*maxm+
     1                             400)+5
          end if
        maxn=maxn+maxm
10      maxp=max(maxn,maxp,maxpdr)
        maxq=lnum+3*ndec+int(ca)+maxm+maxm+4
        if(knd.eq.kindd.and.aimag(c).lt.10.0e0_knd.and.real(c).le.
     1     60.0e0_knd.and.real(c).ge.10.0e0_knd.and.mmin.le.40
     2     .and.x.le.0.99e0_knd.and.x.gt.0.1e0_knd)
     3         maxq=max(maxq,250-int(50*x)+4+maxm+maxm)
        maxdr=maxpdr/2+1
        maxd=maxn/2+1
        maxmp=maxm+maxm+5
        maxt=1
        jnebmax=30
        jnenmax=10
        if(x.lt.0.05e0_knd) jnenmax=1
        if(iopang.ne.0) maxt=narg
c
		if (allocated(enr)) then
			deallocate(enr)
		endif
		allocate(enr(lnum, 0:maxd))
		enr = 0
        call main (mmin,minc,mnum,lnum,c,ioprad,iopang,iopnorm,minacc,
     1             x,ngau,step1,nstep1,step2,nstep2,step3,nstep3,narg,
     2             arg,maxd,maxdr,maxe,maxe2,maxint,maxj,maxlp,maxm,
     3             maxmp,maxn,maxp,maxp1,maxpdr,maxq,maxt,neta,jnenmax,
     4             jnebmax,ndec,nex,xneu,r1c,ir1e,r1dc,ir1de,r2c,ir2e,
     5             r2dc,ir2de,s1c,is1e,s1dc,is1de,kindd,kindq, enr)
c
        end
c
c
        subroutine main (mmin,minc,mnum,lnum,cc,ioprad,iopang,iopnorm,
     1                   minacc,x,ngau,step1,nstep1,step2,nstep2,step3,
     2                   nstep3,narg,barg,maxd,maxdr,maxe,maxe2,maxint,
     3                   maxj,maxlp,maxm,maxmp,maxn,maxp,maxp1,maxpdr,
     4                   maxq,maxt,neta,jnenmax,jnebmax,ndec,nex,xneu,
     5                   r1c,ir1e,r1dc,ir1de,r2c,ir2e,r2dc,ir2de,s1,
     6                   is1,s1d,is1d,kindd,kindq, enr)
c
c  purpose:     To coordinate the calculation of both the oblate
c               spheroidal radial and angular functions and their
c               first derivatives using various algorithms.
c
c        use param
c
c  real(knd) and complex(knd) scalars
        real(knd) aj1,aj2,apcoef,apcoefn,api,arg1,c,coefn,coefn1,
     1            coefme,coefmo,coefme1,coefmo1,darg,dconp,dec,deta,
     2            em,etaval,factor,factor1,pcoef,pcoefn,pi,qdm0,
     3            qdm1,qm0,qm1,rm,rm2,step1,step2,step3,ten,term,
     4            termpq,teste,testeo,t1,t2,t3,t4,t5,t6,t7,wm,x,xb,
     5            xbninp,xl,xhigh,xlow,xneu
        real(knd) ang,apcoef1,etaval1,pcoefe,pcoefet,pcoefo,
     1            pdcoefe,pdcoefo,pdcoefet,pcoefe1,pcoefet1,pcoefo1,
     2            pdcoefe1,pdcoefo1,pdcoefet1
        real(knd1) t11,t21,t31,t41,t51,t61,t71
c
        complex(knd) cc,c2,c4,dfnorm,dmlf,dmfnorm,dmlmf,dmsnorm,dmlms,
     1               dmlms1,dneg,dc01,eigest,eigmat,eign,eigp,eigval,
     2               fac1,fac1d,fac2,fac2d,parg,r1cin,r1cm,r1dcin,r11c,
     3               r1dcm,r1d1c,r1ec,r1dec,r1dcest,r2sumest,r2ic,r2dic,
     4               r2lc,r2dlc,r2l1c,r2dl1c,r2nc,r2dnc,r2ec,r2dec,
     5               wronc,wronca,wroncb,wront
        complex(knd) c3
        complex(knd1) c21,c41
c
c  integer arrays with dimension lnum
        integer   ir1e(lnum),ir1de(lnum),ir2e(lnum),ir2de(lnum),
     1            match(lnum),iqdl(lnum),iql(lnum),ieigt(lnum),
     2            is1(lnum,narg),is1d(lnum,narg)
c
c  real(knd) arrays with dimension lnum
        real(knd) qdl(lnum),ql(lnum)
c
c  complex(knd) arrays with dimension lnum
        complex(knd) eig(lnum),r1c(lnum),r1dc(lnum),r2c(lnum),
     1               r2dc(lnum),s1(lnum,narg),s1d(lnum,narg)
c
c  integer and complex(knd) arrays with dimension lnum+1
        integer ifajo(lnum+1)
        complex(knd) fajo(lnum+1)
c
c  complex(knd) arrays with dimension maxd
        complex(knd) bliste(maxd),blisto(maxd),gliste(maxd),
     1               glisto(maxd),enr(lnum, maxd)
        complex(knd1) bliste1(maxd),blisto1(maxd),gliste1(maxd),
     1                glisto1(maxd)
c
c  complex(knd) arrays with dimension maxdr
        complex(knd) drhor(maxdr)
c
c  complex(knd) arrays with dimension maxint
        complex(knd) pint1(maxint),pint2(maxint),pint3(maxint),
     1               pint4(maxint),rpint1(maxint),rpint2(maxint)
c
c  complex(knd) array with dimension maxj
        complex(knd) sbesf(maxj),sbesdf(maxj),sbesfe(maxj),
     1               sbesdfe(maxj),sbesfsv(jnebmax,maxj),
     2               sbesdfsv(jnebmax,maxj)
c
c  integer and real(knd) arrays with dimension maxlp
        integer ibese(maxlp),ipnormint(maxlp),ipnormint1(maxlp),
     1          ibesee(maxlp),ibesesv(jnebmax,maxlp)
        real(knd) pnormint(maxlp),pnormint1(maxlp)
c
c  complex(knd) arrays with dimension maxlp
        complex(knd) sbesdr(maxlp),sbesn(maxlp),sneudr1(maxlp),
     1               sneudr2(maxlp),sneudre(maxlp),
     2               sneudrsv(jnenmax,maxlp),sbesne(maxlp),
     3               sbesdre(maxlp),sbesnsv(jnebmax,maxlp),
     4               sbesdrsv(jnebmax,maxlp)
c
c  real(knd) array with dimension maxmp
        real(knd) qdqr(maxmp)
c
c  complex(knc) array with dimension maxmp
        complex(knd) enrneg(maxmp)
c
c  real(knd) arrays with dimension maxn
        real(knd) prat1(maxn)
c
c  integer arrays with dimension maxn
        integer ineue1(maxn),ineue2(maxn),ineuee(maxn),
     1          ineuesv(jnenmax,maxn)
c
c  complex(knd) arrays with dimension maxn
        complex(knd) sneufe(maxn),sneudfe(maxn),sneufsv(jnenmax,maxn),
     1               sneudfsv(jnenmax,maxn),sneuf1(maxn),sneudf1(maxn),
     2               sneuf2(maxn),sneudf2(maxn),sneun1(maxn),
     3               sneun2(maxn),sneune(maxn),sneunsv(jnenmax,maxn)
c
c  real(knd) arrays with dimension given by maxp
        real(knd) alpha(maxp),beta(maxp),coefa(maxp),coefb(maxp),
     1            coefc(maxp),coefd(maxp),coefe(maxp),gamma(maxp),
     2            pdr(maxt,maxp),pdrat(maxt,maxp),pdratt(maxp),
     3            pr(maxt,maxp),prat(maxt,maxp),pratb(maxp),
     4            pratbsv(jnenmax,maxp),prattsv(jnenmax,maxp),
     5            pdrattsv(jnenmax,maxp),pratt(maxp)
c
c  real(knd) arrays with dimension given by maxp
        real(knd) pratb1(maxp),pratt1(maxp),pdratt1(maxp),
     1            pratbsv1(jnebmax,maxp),prattsv1(jnebmax,maxp),
     2            pdrattsv1(jnebmax,maxp)
c
c  real(knd) arrays with dimension maxpdr
        real(knd) prx(maxpdr),pdrx(maxpdr)
c
c  real(knd) arrays with dimension maxq
        real(knd) qrat(maxq),qdrat(maxq)
        real(knd) qrat1(maxq),qdrat1(maxq)
c
c  real(knd) and integer arrays with dimension maxt
        real(knd) arg(maxt),barg(maxt),etainp(maxt),pdnorm(maxt),
     1            pdnorma(maxt),pnorm(maxt),pnorma(maxt),pdtempe(maxt),
     2            pdtempo(maxt),ptempe(maxt),ptempo(maxt),xin(maxt),
     3            xlninp(maxt)
        integer ipdnorm(maxt),ipdnorma(maxt),ipnorm(maxt),
     1          ipnorma(maxt),ipdtempe(maxt),ipdtempo(maxt),
     2          iptempe(maxt),iptempo(maxt),is1e(maxt),is1de(maxt),
     3          naccs(maxt),naccds(maxt)
c
c  complex(knd) arrays with dimension maxt
        complex(knd) s1c(maxt),s1dc(maxt)
c
c  real(knd) arrays with dimension ngau
        real(knd) wr(ngau),xr(ngau)
c
c  complex(knd) arrays with dimension maxe2
        complex(knd) eigt(maxe2)
c
c  complex(knd) arrays with dimension maxe
        complex(knd) f(maxe),g(maxe)
        complex(knd) d(maxe),e(maxe)
c
c  real(knd) arrays with dimension neta
        real(knd) eta(neta),xbn(neta),xln(neta),wmeta2(neta)
c
c  miscellaneous integer arrays
        integer neeb(jnenmax),neeb1(jnebmax),limpsv(jnenmax),
     1          limp1sv(jnebmax),limnsv(jnenmax),jelimsv(jnenmax),
     2          jelim1sv(jnebmax),limjsv(jnebmax)
c
c
        c=abs(cc)
        ten=10.0e0_knd
        nfac=nex/3
        if(2*(nfac/2).ne.nfac) nfac=nfac-1
        teste=ten**(nfac)
        testeo=1.0e0_knd/teste
        dec=ten**(-ndec-1)
        factor=ten**(nex-20)
        factor1=1.0e0_knd/factor
        ifactor=nex-20
        jtest=ndec-minacc-2
        pi=acos(-1.0e0_knd)
        api=pi/180.0e0_knd
        c2=cc*cc
        c4=c2*c2
        c21=c2
        c41=c4

c
c  begin loops
          igau=0
          ibflag1=1
          ibflag2=1
            if(knd.eq.kindd) then
            nsubt=2
            nsub1mt=6
            end if
            if(knd.eq.kindq) then
            nsubt=4
            nsub1mt=10
            end if
c          if(knd.eq.kindd.and.ioprad.ne.0) write(40,20) x,cc
c20        format(1x,'x = ',e23.14,/,1x,'c = ',e23.14,e23.14)
c          if(knd.eq.kindq.and.ioprad.ne.0) write(40,25) x,cc
c25        format(1x,'x = ',e38.30,/,1x,'c = ',e38.30,e38.30)
          wront=1.0e0_knd/(cc*(x*x+1.0e0_knd))
            do 1540 mi=1,mnum
            m=mmin+minc*(mi-1)
            em=m
            m2=m+m
c            if(knd.eq.kindd.and.iopang.ne.0) write(50,30) cc,m
c30          format(1x,'c = ',e23.14,e23.14,'; m = ',i5)
c            if(knd.eq.kindq.and.iopang.ne.0) write(50,35) cc,m
c35          format(1x,'c = ',e38.30,e38.30,'; m = ',i5)
c            if(ioprad.ne.0) write(40,40) m
c40          format(1x,'m = ',i5)
c            if(knd.eq.kindd.and.iopang.ne.0) write(30,50) cc,m
c50          format(1x,'c = ',e23.14,e23.14,'; m = ',i5)
c            if(knd.eq.kindq.and.iopang.ne.0) write(30,55) cc,m
c55          format(1x,'c = ',e38.30,e38.30,'; m = ',i5)
            rm=m
            rm2=m+m
            icounter=0
            limcsav=0
            jbes=3*ndec+int(c)
            iopbes=1
            iopeta1=0
            jflageta1=0
            iopint=0
            iopleg=0
            iopleg1=0
            iopneu0=0
            iopeta=0
            if(ioprad.ne.2) go to 60
              if(knd.eq.kindd) then
              if(x.le.0.99e0_knd.and.c.le.20.0e0_knd) iopleg=1
              if(x.gt.0.99e0_knd.and.c.le.20.0e0_knd) iopneu0=1
              end if
              if(knd.eq.kindq) then
              if(x.le.0.99e0_knd.and.c.le.60.0e0_knd) iopleg=1
              if(x.gt.0.99e0_knd.and.c.le.40.0e0_knd) iopneu0=1
              end if
60          continue
            jneu1max=0
            jnen=0
            incnee=1
            if(x.gt.0.4e0_knd) incnee=2
            iplflag=0
            neelow=28
            if(x.gt.0.1e0_knd) neelow=26
            if(x.gt.0.2e0_knd) neelow=24
            if(x.gt.0.3e0_knd) neelow=22
            if(x.gt.0.4e0_knd) neelow=20
            if(x.gt.0.5e0_knd) neelow=18
            if(x.gt.0.6e0_knd) neelow=14
            if(x.gt.0.7e0_knd) neelow=12
            if(x.gt.0.8e0_knd) neelow=8
            if(x.gt.0.9e0_knd) neelow=2
            nee=neelow
            jnen1=0
            incnee1=1
            nee1=1
            idir=0
            iflageta1=0
            nsub1p=ndec
            nsubd1p=ndec
            if(iopang.eq.0) go to 70
            limps1=lnum+3*ndec+int(c)
            if((limps1+3).gt.maxp) limps1=maxp-3
            iopd=0
            if(iopang.eq.2) iopd=1
            call pleg(m,limps1,maxp,limcsav,iopd,ndec,barg,narg,maxt,
     1                pr,pdr,pdnorm,ipdnorm,pnorm,ipnorm,alpha,beta,
     2                gamma,coefa,coefb,coefc,coefd,coefe)
            limcsav=limps1
70          limj=lnum+3*ndec+int(c)+maxm
            if (x.ne.0.0e0_knd) go to 170
c
c  calculation of factors for radial functions when x = 0
            fac1=cmplx(1.0e0_knd,0.0e0_knd)
            ifac1=0
            if(m.eq.0) go to 100
              do 90 k=1,m
              fac1=fac1*cc/(k+k+1)
              if(abs(fac1).lt.factor) go to 80
              fac1=fac1/factor
              ifac1=ifac1+ifactor
              go to 90
80            if(abs(fac1).gt.factor1) go to 90
              fac1=fac1*factor
              ifac1=ifac1-ifactor
90            continue
100         iterm=int(log10(abs(fac1)))
            fac1=fac1*(ten**(-iterm))
            ifac1=ifac1+iterm
            fac1d=fac1*cc/(m+m+3)
            ifac1d=ifac1
            if(ioprad.eq.1) go to 170
            fac2=(m+m-1)*pi/(2.0e0_knd*cc)
            ifac2=0
            if(m.eq.0) go to 130
              do 120 k=1,m
              fac2=fac2*cc/(k*2.0e0_knd)
              if(abs(fac2).lt.factor) go to 110
              fac2=fac2/factor
              ifac2=ifac2+ifactor
              go to 120
110           if(abs(fac2).gt.factor1) go to 120
              fac2=fac2*factor
              ifac2=ifac2-ifactor
120           continue
130         iterm=int(log10(abs(fac2)))
            fac2=fac2*(ten**(-iterm))
            ifac2=ifac2+iterm
            fac2d=(m+m-1)*(m+m-3)*(m+m+1)*pi/(cc*cc*2.0e0_knd)
            ifac2d=0
            if(m.eq.0) go to 160
              do 150 k=1,m
              fac2d=fac2d*c/(k*2.0e0_knd)
              if(abs(fac2d).lt.factor) go to 140
              fac2d=fac2d/factor
              ifac2d=ifac2d+ifactor
              go to 150
140           if(abs(fac2d).gt.factor1) go to 150
              fac2d=fac2d*factor
              ifac2d=ifac2d-ifactor
150         continue
160         iterm=int(log10(abs(fac2d)))
            fac2d=fac2d*(ten**(-iterm))
            ifac2d=ifac2d+iterm
170         xb=sqrt(x*x+1.0e0_knd)
            if(ioprad.eq.0.or.x.eq.0.0e0_knd) go to 180
            prat1(1)=1.0e0_knd
            prat1(2)=rm2+1.0e0_knd
              do jp=3,limj
              aj1=jp-1
              aj2=jp-2
              prat1(jp)=(rm2+aj1)*(rm2+aj2)/(aj1*aj2)
              end do
            pcoefn=(x*x+1.0e0_knd)/(x*x)
            apcoefn=(rm/2.0e0_knd)*log10(pcoefn)
            ipcoefn=int(apcoefn)
            pcoefn=ten**(apcoefn-ipcoefn)
            if(mi.ne.1) go to 180
            call sphbes(cc,x,limj,maxj,maxlp,sbesf,sbesdf,sbesn,ibese,
     1                  sbesdr)
180         continue
c
c  obtain starting eigenvalues eigt(i) for the Bouwkamp procedure
            nbp=int(2*(real(cc)+aimag(cc))/3.1416)
            limeig=max(67,(4*nbp)/3)
            lime=max(50,nbp+30)
            lime2=lime+lime
              do 190 i=1,lime
              xl=m+i+i-2
              d(i)=xl*(xl+1.0e0_knd)/c2-(2.0e0_knd*xl*(xl+1.0e0_knd)-
     1             2.0e0_knd*em*em-1.0e0_knd)/((2.0e0_knd*xl-
     2             1.0e0_knd)*(2.0e0_knd*xl+3.0e0_knd))
190           continue
            nm1=lime-1
              do 200 i=1,nm1
              xl=m+i+i-2
              e(i)=(-1.0e0_knd/(2.0e0_knd*xl+3.0e0_knd))*
     1             sqrt(((xl+2.0e0_knd+em)*(xl+1.0e0_knd+em)*
     2             (xl+2.0e0_knd-em)*(xl+(1.0e0_knd)-em))/
     3             ((2.0e0_knd*xl+5.0e0_knd)*(2.0e0_knd*xl+
     4             1.0e0_knd)))
200           continue
            call cmtql1(lime,maxe,d,e)
              do 210 i=1,lime
              f(i)=c2*d(i)
210           continue
              do 220 i=1,lime
              xl=m+i+i-1
              d(i)=xl*(xl+1.0e0_knd)/c2-(2.0e0_knd*xl*(xl+1.0e0_knd)-
     1             2.0e0_knd*em*em-1.0e0_knd)/((2.0e0_knd*xl-
     2             1.0e0_knd)*(2.0e0_knd*xl+3.0e0_knd))
220           continue
            nm1=lime-1
              do 230 i=1,nm1
              xl=m+i+i-1
              e(i)=(-1.0e0_knd/(2.0e0_knd*xl+3.0e0_knd))*
     1             sqrt(((xl+2.0e0_knd+em)*(xl+1.0e0_knd+em)*
     2             (xl+2.0e0_knd-em)*(xl+(1.0e0_knd)-em))/
     3             ((2.0e0_knd*xl+5.0e0_knd)*(2.0e0_knd*xl+
     4             1.0e0_knd)))
230           continue
            call cmtql1(lime,maxe,d,e)
              do 240 i=1,lime
              g(i)=c2*d(i)
240           continue
            call eigorder(cc,lime,maxe,maxe2,f,g,m,lnum,limeig,eigt,
     1                    lipl,liplp,lips)
c
c  determine the number of leading decimal digits of agreement
c  for lower order paired eigenvalues
            listart=1
            matlim=min(lnum,nbp+3)
            matlim=max(matlim,lipl+1)
            if(matlim.gt.lnum) matlim=lnum
            limi=0
            match(2)=0
              do i=2,matlim,2
              match(i)=-int(log10(abs((eigt(i)-eigt(i-1))/eigt(i)
     1                  +ten*dec)))
              if(match(i).lt.0) match(i)=0
              if(match(i).gt.10) limi=i
              if(match(i).gt.ndec) match(i)=ndec
              if(match(i).gt.minacc+1) listart=i+1
              end do
              if(listart.gt.1) then
              iopleg=0
              iopneu0=0
              end if
255         continue
              if(knd.eq.kindd.and.c.gt.100.0e0_knd) then
               do i=2,limi,2
                  if(match(i).lt.11) then
                  match(i)=11
                  end if
                end do
              end if
c
c  prepare for do loop over index l
            iflag=0
            iflagint=0
            iflagneu=0
            iqlegflag=0
            ipint=0
            intlim=maxint-3
            istartint=0
            xlow=0.0005e0_knd
            if(aimag(cc).gt.2.0e0_knd) xlow=0.0000001e0_knd
            if(x.lt.xlow) istartint=2
            xhigh=0.2e0_knd
            if(x.gt.xhigh) istartint=2
            if(istartint.eq.2) iopint=0
            if(ioprad.eq.2.and.iopint.eq.0.and.x.le.0.99e0_knd) iopleg=1
            nflag=0
            iflageta=0
            ieta=0
            naccr=minacc+1
            naccrsav=minacc
            naccrsavp=minacc
            naccneu0p=0
            naccneu0p2=0
            nacclegp=minacc
            naccintp=minacc
            naccrplp=0
            naccetabp=0
            legflag=0
            jflagleg=0
            legstart=m
            ioppsum=1
            iopqnsum=1
            jneu1=0
            jneu1s=0
            jneu0=0
            jneu0s=0
            icflag1=0
            nacintsa=0
            jeta=0
            jeta1=0
            jbes=0
            jint=0
            jlegp=0
            jleg=0
            jleg1=0
            jtest0=0
            jtest1=0
            naccetas=0
            jetaflag=0
            naccleg1p=0
            ir1ep=0
            istarteta=1
            istartneu0=1
            istartleg=1
            naccflag=0
            jneu0x=0
            limax=m
            isteta1=1
            iopeige=0
            iopeigo=0
            incnflag=0
            ietacount=0
            ieigflage=1
            ieigflago=1
            iopeige=0
            iopeigo=0
            iflagl1=0
            max1e=0
            max1o=0
            kflag=0
              if(knd1.ne.knd) then
              if(real(cc).le.10.0e0_knd) kflag=2
              if(real(cc).gt.10.0e0_knd.and.real(cc).le.25.0e0_knd
     1           .and.aimag(cc).le.15.0e0_knd) kflag=2
              if(real(cc).gt.25.0e0_knd.and.real(cc).le.50.0e0_knd
     1           .and.aimag(cc).le.10.0e0_knd) kflag=2
              if(real(cc).gt.50.0e0_knd.and.real(cc).le.100.0e0_knd
     1           .and.aimag(cc).le.9.0e0_knd) kflag=2
              if(real(cc).gt.100.0e0_knd.and.real(cc).le.1000.0e0_knd
     1           .and.aimag(cc).le.7.0e0_knd) kflag=2
              if(real(cc).gt.1000.0e0_knd.and.real(cc).le.10000.0e0_knd
     1           .and.aimag(cc).le.6.0e0_knd) kflag=2
              else
              kflag=2
              end if
            if(knd.eq.kindd.and.aimag(cc).lt.10.0e0_knd.and.real(cc).le.
     1         60.0e0_knd.and.real(cc).ge.10.0e0_knd.and.m.le.40
     2         .and.x.le.0.99e0_knd.and.x.gt.0.1e0_knd) iflagl1=1
c            if(knd.eq.kindd.and.ioprad.ne.0) write(20,260) x,cc,m
c260         format(1x,e23.14,e23.14,e23.14,i5)
c            if(knd.eq.kindq.and.ioprad.ne.0) write(20,265) x,cc,m
c265         format(1x,e38.30,e38.30,e38.30,i5)
              do 1510 li=1,lnum
              l=m+(li-1)
              nacccor=0
              iopbesa=0
c              if(iopang.ne.0) write(30,270) l
c270           format(1x,i6)
c              if(iopang.ne.0) write(50,280) l
c280           format(1x,'l = ',i6)
              ix=l-m-2*((l-m)/2)
              if(ix.eq.0) naccflag=0
              if(li.eq.1) naccrsav=minacc
              if(li.eq.lips.or.li.eq.liplp) legstart=l
              naccr=-1
              naccint=0
              naccleg=0
              naccneu0=0
              naccneu1=0
              naccleg1=0
              naccetab=0
              naccrt=0
              jflagl=0
              neemax=0
                if(lips.ne.lipl+1.and.li.ne.1.and.li.eq.lips) then
                istartint=0
                if(x.lt.xlow) istartint=2
                if(x.gt.xhigh) istartint=2
                iflag=0
                iflagint=0
                istartleg=1
                istartneu0=1
                naccetas=0
                jetaflag=0
                iflageta=0
                iflagneu=0
                nee1=1
                isteta1=1
                iopeta=0
                iopint=0
                iopbes=1
                iopeta1=0
                end if
              if(x.lt.0.05e0_knd) istarteta=0
              if(li.eq.listart.and.iopneu0.eq.0.and.x.ge.xneu) iopneu0=1
              if(istartneu0.eq.0) iopneu0=0
              if(istarteta.eq.0) iopeta=0
              limdrad=3*ndec+int(c)+10
              if(ioprad.ne.0.and.li.ne.1) limdrad=jbes+jbes+20+
     1                                    int(sqrt(c))
              if(ioprad.eq.2.and.iopint.ne.0.and.li.gt.listart.and.
     1            jint.gt.jbes) limdrad=jint+jint+20+int(sqrt(c))
              limdang=3*ndec+int(c)
              if(iopang.ne.0.and.li.ne.1) limdang=jang+jang+20+
     1                                            int(sqrt(c))
              if(iopang.eq.0) limd=limdrad
              if(ioprad.eq.0) limd=limdang
              if(iopang.ne.0.and.ioprad.ne.0) limd=max(limdang,limdrad)
              if(li.eq.1) limnorm=limdrad
              if(li.gt.1) limnorm=jnorm+jnorm+20+int(sqrt(c))
              limd=max(limd,limnorm)
              if(ioprad.ne.2) go to 290
              if(iopleg.eq.1) limdleg=l-m+3*ndec+int(c)
              if(iopleg.eq.2) limdleg=jleg+jleg+20+int(sqrt(c))
              if(iopleg.ne.0.and.li.ge.listart) limd=max(limd,limdleg)
                if(knd.eq.kindd) then
                if(x.ge.0.01e0_knd) lneu0=2*int(25/(x*x)+300/x+3*c+
     1                                    10000)
                if(x.ge.0.1e0_knd) lneu0=2*int((l-m+c/5+0.5e0_knd*m+
     1                                   200)*1.4e0_knd/x)
                if(x.ge.0.5e0_knd) lneu0=2*int((l-m+c/5+0.5e0_knd*m+
     1                                   200)/x)
                if(x.ge.1.0e0_knd) lneu0=2*int(l-m+c/5+0.5e0_knd*m+
     1                                    200)
                end if
                if(knd.eq.kindq) then
                if(x.ge.0.01e0_knd) lneu0=2*int(25/(x*x)+400/x+3*c+
     1                                    14000)
                if(x.ge.0.1e0_knd) lneu0=2*int((l-m+c/5+0.5e0_knd*m+
     1                                   350)*1.4e0_knd/x)
                if(x.ge.0.5e0_knd) lneu0=2*int((l-m+c/5+0.5e0_knd*m+
     1                                   300)/x)
                if(x.ge.1.0e0_knd) lneu0=2*int(l-m+c/5+0.5e0_knd*m+
     1                                    300)
                end if
              if(iopneu0.eq.3.and.jtest0.gt.minacc.and.li.ne.listart+1)
     1               lneu0=jneu0max+jneu0max+40+int(sqrt(c)*
     2                     (2/min(1.0e0_knd,x))+100/x)
              if((iopneu0.ne.0.and.li.ge.listart).or.(x.ge.xneu
     1               .and.li.eq.listart)) limd=max(limd,lneu0)
                if(knd.eq.kindd) then
                if(x.ge.0.05e0_knd) leta=2*int(25/(x*x)+300/x+3*c+
     1                                   10000)
                if(x.ge.0.1e0_knd) leta=2*int((l-m+c/5+0.5e0_knd*m+
     1                                   200)*1.4e0_knd/x)
                if(x.ge.0.5e0_knd) leta=2*int((l-m+c/5+0.5e0_knd*m+
     1                                  300)/x)
                if(x.ge.1.0e0_knd) leta=2*int(l-m+c/5+0.5e0_knd*m+300)
                end if
                if(knd.eq.kindq) then
                if(x.ge.0.05e0_knd) leta=2*int(25/(x*x)+400/x+3*c+
     1                                   14000)
                if(x.ge.0.1e0_knd) leta=2*int((l-m+c/5+0.5e0_knd*m+
     1                                   350)*
     1                                  1.4e0_knd/x)
                if(x.ge.0.5e0_knd) leta=2*int((l-m+c/5+0.5e0_knd*m+
     1                                  400)/x)
                if(x.ge.1.0e0_knd) leta=2*int(l-m+c/5+0.5e0_knd*m+400)
                end if
              if(iopeta.eq.3.and.nacceta.gt.minacc)
     1            leta=jeta+jeta+40+
     2                 int(sqrt(c)*(2/min(1.0e0_knd,x))+5/x)
              if((iopeta.ne.0.and.li.ge.listart).or.(x.ge.0.05e0_knd
     1               .and.li.eq.listart)) limd=max(limd,leta)
290           continue
              if(ix.eq.0) limd=max(limd,max1e)
              if(ix.eq.0) limd=max(limd,max1o)
              if(limd.gt.maxn) limd=maxn-4
              if(2*(limd/2).ne.limd) limd=limd-1
              if(li.le.limeig)  eigmat=eigt(li)
              if(li.gt.limeig) eigmat=4.0e0_knd*eig(li-1)-6.0e0_knd*
     1                      eig(li-2)+4.0e0_knd*eig(li-3)-eig(li-4)
              if(li.le.8) eigest=(0.0e0_knd,0.0e0_knd)
              if(li.gt.8) eigest=4.0e0_knd*eig(li-2)-6.0e0_knd*
     1                      eig(li-4)+4.0e0_knd*eig(li-6)-eig(li-8)
c
c  use Bouwkamp procedure to obtain accurate eigenvalue
              ndec1=precision(bliste1(1))
                if(li.eq.1) then
                jlowe=1
                limdle=2
                min1e=2
                ienre=(3*ndec+int(c))/2
                if(kflag.ne.2) ienre=(3*ndec1+int(c))/2
                max1e=min(2*ienre+20,limd)
                jlow1e=1
                end if
                if(li.eq.2) then
                jlowo=1
                limdlo=3
                min1o=3
                ienro=(3*ndec+int(c))/2
                if(kflag.ne.2) ienro=(3*ndec1+int(c))/2
                max1o=min(2*ienro+20,limd)
                jlow1o=1
                end if
                if(li.gt.2.and.ix.eq.0) then
                max1e=min(limd,max1e)
                end if
                if(li.gt.2.and.ix.eq.1) then
                max1o=min(limd,max1o)
                end if
c
c  compute the coeficients in the bouwkamp method
              if(ix.eq.1) go to 310
c
c  beta coefficients (bliste) for l-m even
              if(limdle.gt.limd) go to 300
              j=jlowe
                do i=limdle,limd,2
                i2=i+i
                t1=i
                t2=i-1
                t3=m2+i
                t4=m2+i-1
                t5=m2+i2-1
                t6=m2+i2-3
                t7=m2+i2+1
                bliste(j)=c4*t1*t2*t3*t4/(t5*t5*t6*t7)
                j=j+1
                end do
                jsave=j
c
c  gamma coeficients (gliste) for l-m even
              j=jlowe
                do i=limdle-1,limd+1,2
                i2=i+i
                t1=m+i-1
                t2=m+i
                t3=m2*m2-1
                t4=m2+i2-3
                t5=m2+i2+1
                gliste(j)=t1*t2-(0.5e0_knd)*c2*((1.0e0_knd)-t3/
     1                    (t4*t5))
                j=j+1
                end do
                limdle=limd+2
                jlowe=jsave
300           continue
              go to 320
310           continue
c
c  beta coefficients (blisto) for l-m odd
              if(limdlo.gt.limd) go to 320
              j=jlowo
                do i=limdlo,limd,2
                i2=i+i
                t1=i
                t2=i-1
                t3=m2+i
                t4=m2+i-1
                t5=m2+i2-1
                t6=m2+i2-3
                t7=m2+i2+1
                blisto(j)=c4*t1*t2*t3*t4/(t5*t5*t6*t7)
                j=j+1
                end do
                jsave=j
c
c  gamma coeficient (glisto) for l-m odd
              j=jlowo
                do i=limdlo-1,limd+1,2
                i2=i+i
                t1=m+i-1
                t2=m+i
                t3=m2*m2-1
                t4=m2+i2-3
                t5=m2+i2+1
                glisto(j)=t1*t2-(0.5e0_knd)*c2*((1.0e0_knd)-t3/
     1                    (t4*t5))
                j=j+1
                end do
                jlowo=jsave
                limdlo=limd+1
320           continue
              if(kflag.eq.2) go to 350
              if(ix.eq.1) go to 340
c
c  beta coefficients (bliste1) for l-m even
              if(min1e.gt.max1e) go to 330
              j=jlow1e
                do i=min1e,max1e,2
                  if(knd1.ne.knd) then
                  i2=i+i
                  t11=i
                  t21=i-1
                  t31=m2+i
                  t41=m2+i-1
                  t51=m2+i2-1
                  t61=m2+i2-3
                  t71=m2+i2+1
                  bliste1(j)=c41*t11*t21*t31*t41/(t51*t51*t61*t71)
                  else
                  bliste1(j)=bliste(j)
                  end if
                j=j+1
                end do
                jsave=j
c
c  gamma coeficients (gliste1) for l-m even
              j=jlow1e
              if(li.gt.2) j=jlow1e+1
              imin=min1e-1
              if(li.gt.2) imin=min1e+1
                do i=imin,max1e+1,2
                  if(knd1.ne.knd) then
                  i2=i+i
                  t11=m+i-1
                  t21=m+i
                  t31=m2*m2-1
                  t41=m2+i2-3
                  t51=m2+i2+1
                  gliste1(j)=t11*t21-(0.5e0_knd1)*c21*((1.0e0_knd1)-t31/
     1                       (t41*t51))
                  else
                  gliste1(j)=gliste(j)
                  end if
                j=j+1
                end do
                jlow1e=jsave
                min1e=max1e+2
330           continue
              go to 350
340           continue
c
c  beta coefficients (blisto1) for l-m odd
              if(min1o.gt.max1o) go to 350
              j=jlow1o
                do i=min1o,max1o,2
                  if(knd1.ne.knd) then
                  i2=i+i
                  t11=i
                  t21=i-1
                  t31=m2+i
                  t41=m2+i-1
                  t51=m2+i2-1
                  t61=m2+i2-3
                  t71=m2+i2+1
                  blisto1(j)=c41*t11*t21*t31*t41/(t51*t51*t61*t71)
                  else
                  blisto1(j)=blisto(j)
                  end if
                j=j+1
                end do
                jsave=j
c
c  gamma coeficient (glisto1) for l-m odd
              j=jlow1o
              if(li.gt.2) j=jlow1o+1
              imin=min1o-1
              if(li.gt.2) imin=min1o+1
                do i=imin,max1o+1,2
                  if(knd1.ne.knd) then
                  i2=i+i
                  t11=m+i-1
                  t21=m+i
                  t31=m2*m2-1
                  t41=m2+i2-3
                  t51=m2+i2+1
                  glisto1(j)=t11*t21-(0.5e0_knd1)*c21*((1.0e0_knd1)-t31/
     1                    (t41*t51))
                  else
                  glisto1(j)=glisto(j)
                  end if
                j=j+1
                end do
              min1o=max1o+1
              jlow1o=jsave
350           continue
              itestm=ndec
              matest=0
              eigp=(0.0e0_knd,0.0e0_knd)
              eign=(0.0e0_knd,0.0e0_knd)
              if(li.gt.2.and.li.le.lipl) eigp=eig(li-2)
              if(li.gt.2.and.li.gt.lipl.and.li.le.limeig) eigp=eig(li-2)
                if(lips-2.gt.0.and.li-2.ge.lips.and.li-2.lt.liplp) then
                if(2*((lips-li+2)/2).ne.lips-li+2) eigp=eig(lips-1)
                if(2*((lips-li+2)/2).eq.lips-li+2) eigp=eig(lips-2)
                end if
              if(li.le.lipl) eign=eigt(li+2)
              if(li.gt.lipl.and.li.le.limeig-2) eign=eigt(li+2)
                if(li+2.ge.lips.and.li+2.lt.liplp) then
                if(2*((liplp-li-2)/2).ne.liplp-li-2) eign=eigt(liplp+1)
                if(2*((liplp-li-2)/2).eq.liplp-li-2) eign=eigt(liplp)
                end if
              if(li.le.lipl) mat=match(li+1-ix)
              if(li.gt.lipl) mat=0
              matest=mat
                if(li.le.lipl.and.li.gt.4) then
                nmatch=-int(log10(abs((eig(li-1-ix)-eig(li-2-ix))/
     1                  (eig(li-1-ix))+ten*dec)))
                int5=-int(log10(abs((eig(li-1)-eigt(li-1))/
     1                  (eig(li-1))+ten*dec)))
                int6=-int(log10(abs((eig(li-2)-eigt(li-2))/
     1                  (eig(li-1))+ten*dec)))
                  if(nmatch.gt.matest+2) then
                  if(ix.eq.1.and.int5.gt.int6.and.iflage.eq.0)
     1                eigmat=eig(li-1)
                  if(ix.eq.0.and.int5.gt.int6.and.iflago.eq.0)
     1                eigmat=eigt(li+1)
                  matest=nmatch-2
                  end if
                end if
              if(ix.eq.0.and.li.gt.2) iflagep=iflage
              if(ix.eq.0) call conver (l,m,lnum,cc,limd,bliste,gliste,
     1                                 bliste1,gliste1,ndec,maxd,ioprad,
     2                                 minacc,eigest,eigmat,lipl,lips,
     3                                 liplp,mat,eigp,eign,kindd,kindq,
     4                                 eigval,enr(li,:),ienre,itestme,
     5                                 naccre,ieigt,iopeige,iflage,
     6								   kflag)
              if(ix.eq.1) call conver (l,m,lnum,cc,limd,blisto,glisto,
     1                                 blisto1,glisto1,ndec,maxd,ioprad,
     2                                 minacc,eigest,eigmat,lipl,lips,
     3                                 liplp,mat,eigp,eign,kindd,kindq,
     4                                 eigval,enr(li,:),ienro,itestmo,
     5                                 naccre,ieigt,iopeigo,iflago,
     6								   kflag)
              if(ix.eq.0) itestm=itestme
              if(ix.eq.1) itestm=itestmo
              eig(li)=eigval
              if(ix.eq.0) max1e=2*ienre+20
              if(ix.eq.1) max1o=2*ienro+20
380           call dnorm (l,m,cc,ndec,nex,limd,maxd,enr(li,:),ioprad,
     1                    iopang,dc01,idc01,dfnorm,idfe,dmlf,idmlfe,
     2                    dmfnorm,idmfe,dmlmf,idmlmfe,dmsnorm,idmse,
     3                    dmlms,idmlmse,jnorm,jsubf,jsubmf,jsubms)
              jsub=max(jsubf,jsubmf)
              if(ioprad.eq.0) go to 1410
              if(li.lt.listart) go to 385
              if(ioprad.ne.2) go to 385
              if(li.eq.listart.and.jsub.le.ndec-naccrsav.and.x.le.
     1            0.99e0_knd.and.iopleg.eq.0.and.iopleg1.eq.0.and.
     2            iopint.eq.1) iopleg=1
              if(ndec-jsub.gt.naccrsav-2.and.x.le.0.99e0_knd.and.
     1        iopleg.eq.0.and.iopleg1.eq.0.and.l.ge.legstart) iopleg=1
                if(li.eq.listart.and.x.ge.xneu) then
                  if(jsubf.le.ndec-minacc) then
                  iopneu0=1
                  else
                  iopneu0=0
                  end if
                end if
              jsubtest=ndec-min(naccrsav,minacc)
              if(li.ne.listart.and.jsubf.le.jsubtest.and.x.ge.xneu
     1          .and.iopneu0.eq.0.and.(iopleg.eq.0.or.nacclegp.lt.
     2          minacc).and.(iopint.eq.0.or.naccintp.lt.minacc).and.
     3          (iopleg1.eq.0.or.naccleg1p.lt.minacc).and.iopeta.eq.0)
     4               iopneu0=4
              if(istartneu0.eq.0) iopneu0=0
              if((li.eq.listart.or.iopeta.ne.0).and.iopneu0.eq.4)
     1               iopneu0=1
              if(li.eq.listart.and.iopneu0.eq.1.and.iopint.eq.1)
     1            iopint=0
385           continue
              if(x.ne.0.0e0_knd) go to 550
c
c  determine oblate radial functions of both kinds when x = 0
              ioppsum=0
              limdr=int(c)+2*ndec
              if(ioprad.eq.2.and.m.ne.0) call dalt(l,m,cc,ndec,nex,
     1                 limdr,maxdr,maxmp,ioppsum,eigval,enrneg,drhor,
     2                 dneg,idneg,nsdneg,nsdrhor1, nsdrho)
	 
              ifsub=max(jsub,nsdneg)
              if(m.eq.0) dneg=(1.0e0_knd,0.0e0_knd)
              if(m.eq.0) idneg=0
              naccr=min(ndec-ifsub-1,naccre-1,itestm-1,ndec-2)
              if(naccr.lt.0) naccr=0
              naccr1=min(ndec-jsubmf-1,naccre-1,itestm-1,ndec-2)
              if(ix.eq.1) go to 420
              if(li.ne.1) fac1=-real(l-m,knd)*(l-m-1)*fac1/
     1                          (real(l+m,knd)*(l+m-1))
              iterm=int(log10(abs(fac1)))
              fac1=fac1*(ten**(-iterm))
              ifac1=ifac1+iterm
              r1c(li)=fac1*dc01/dmfnorm
              ir1e(li)=int(log10(abs(r1c(li))))
              r1c(li)=r1c(li)*(ten**(-ir1e(li)))
              ir1e(li)=ir1e(li)+ifac1+idc01-idmfe
              if(abs(r1c(li)).ge.1.0e0_knd) go to 390
              r1c(li)=r1c(li)*ten
              ir1e(li)=ir1e(li)-1
390           r1dc(li)=(0.0e0_knd,0.0e0_knd)
              ir1de(li)=0
              if(ioprad.eq.1) go to 450
              r2dc(li)=1.0e0_knd/(cc*r1c(li))
              ir2de(li)=int(log10(abs(r2dc(li))))
              r2dc(li)=r2dc(li)*(ten**(-ir2de(li)))
              ir2de(li)=ir2de(li)-ir1e(li)
              if(abs(r2dc(li)).ge.1.0e0_knd) go to 400
              r2dc(li)=r2dc(li)*ten
              ir2de(li)=ir2de(li)-1
400           if(naccr.eq.0) r2c(li)=(0.0e0_knd,0.0e0_knd)
              if(naccr.eq.0) ir2e(li)=0
              if(li.ne.1) fac2=-real(l+m-1,knd)*(l-m-1)*fac2/
     1                         (real(l-m,knd)*(l+m))
              if(naccr.eq.0) go to 450
              r2c(li)=fac2*dfnorm*dfnorm/(dneg*dc01*dmfnorm)
              ir2e(li)=int(log10(abs(r2c(li))))
              r2c(li)=r2c(li)*(ten**(-ir2e(li)))
              ir2e(li)=ir2e(li)+ifac2-idneg-idc01+idfe+idfe-idmfe
              if(abs(r2c(li)).ge.1.0e0_knd) go to 450
              r2c(li)=r2c(li)*ten
              ir2e(li)=ir2e(li)-1
410           go to 450
420           r1c(li)=(0.0e0_knd,0.0e0_knd)
              ir1e(li)=0
              if(li.ne.2) fac1d=-(l-m)*(l-m-1)*fac1d/((l+m)*(l+m-1))
              iterm=int(log10(abs(fac1d)))
              fac1d=fac1d*(ten**(-iterm))
              ifac1d=ifac1d+iterm
              r1dc(li)=fac1d*dc01/dmfnorm
              ir1de(li)=int(log10(abs(r1dc(li))))
              r1dc(li)=r1dc(li)*(ten**(-ir1de(li)))
              ir1de(li)=ir1de(li)+ifac1d+idc01-idmfe
              if(abs(r1dc(li)).ge.1.0e0_knd) go to 430
              r1dc(li)=r1dc(li)*ten
              ir1de(li)=ir1de(li)-1
430           if(ioprad.eq.1) go to 450
              r2c(li)=-1.0e0_knd/(c*r1dc(li))
              ir2e(li)=int(log10(abs(r2c(li))))
              r2c(li)=r2c(li)*(ten**(-ir2e(li)))
              ir2e(li)=ir2e(li)-ir1de(li)
              if(abs(r2c(li)).ge.1.0e0_knd) go to 440
              r2c(li)=r2c(li)*ten
              ir2e(li)=ir2e(li)-1
440           if(naccr.eq.0) r2dc(li)=(0.0e0_knd,0.0e0_knd)
              if(naccr.eq.0) ir2de(li)=0
              if(li.ne.2) fac2d=-real(l-m,knd)*(l+m)*fac2d/
     1                          (real(l+m-1,knd)*(l-m-1))
              if(naccr.eq.0) go to 450
              r2dc(li)=fac2d*dfnorm*dfnorm/(dneg*dc01*dmfnorm)
              ir2de(li)=int(log10(abs(r2dc(li))))
              r2dc(li)=r2dc(li)*(ten**(-ir2de(li)))
              ir2de(li)=ir2de(li)+ifac2d-idneg-idc01+idfe+idfe-idmfe
              if(abs(r2dc(li)).ge.1.0e0_knd) go to 450
              r2dc(li)=r2dc(li)*ten
              ir2de(li)=ir2de(li)-1
450           continue
c              write(40,460)
c460           format(5x,'calculated accurate values for r1 and r1d ',
c     1               'using nonzero term in traditional Bessel ',
c     2               'function expansion')
c              if(knd.eq.kindd) write(40,570) r1c(li),ir1e(li),r1dc(li),
c     1                                   ir1de(li)
c              if(knd.eq.kindq) write(40,575) r1c(li),ir1e(li),r1dc(li),
c     1                                   ir1de(li)
              if(ioprad.eq.1) go to 1330
c              if(ix.eq.0) write(40,500)
c500           format(5x,'calculated r2 using nonzero term in Legendre',
c     1              ' function expansion and r2d from Wronskian and r1')
c              if(ix.eq.1) write(40,510)
c510           format(5x,'calculated r2d using nonzero term in ',
c     1               'Legendre function expansion and r2 from ',
c     2               'Wronskian and r1d')
c              if(knd.eq.kindd) write(40,520) r2c(li),ir2e(li),r2dc(li),
c     1                                   ir2de(li)
c              if(knd.eq.kindq) write(40,525) r2c(li),ir2e(li),r2dc(li),
c     1                                    ir2de(li)
c520           format(10x,'r2 = ',f19.15,f19.15,i5,5x,'r2d = ',
c     1                f19.15,f19.15,i5)
c525           format(10x,'r2 = ',f35.31,f35.31,i5,/,12x,'r2d = ',
c     1                f35.31,f35.31,i5)
c              if(ix.eq.0) write(40,530) naccr
c530           format(12x,'r2 is accurate to ',I2,' decimal digits. r1,'
c     1               ' r1d, and r2d are nearly fully accurate.')
c              if(ix.eq.1) write(40,540) naccr
c540           format(12x,'r2d is accurate to ',I2,' decimal digits. r1,'
c     1               ' r1d, and r2 are nearly fully accurate.')
              go to 1330
550           continue
c
c  determine oblate radial functions of the first kind
c    r1 calculation using traditional Bessel functions series (eta=1)
c              write(40,560)
c560           format(4x,'r1 and r1d calculation')
              iopbesa=0
              naccr1=0
                if(iopbes.eq.0) then
                jbes=jnorm
                go to 585
                end if
              if(li.eq.1) limr1=3*ndec+int(c)
              if(li.ne.1) limr1=jbes+jbes+20+int(sqrt(c))
              limr1=min(limr1,limj-m-2)
              call r1bes(l,m,cc,x,limr1,ndec,maxd,enr(li,:),maxj,maxn,
     1                   maxlp,nex,iflag,sbesf,sbesdf,sbesn,ibese,
     2                   sbesdr,prat1,pcoefn,ipcoefn,dmfnorm,idmfe,
     3                   ir1ep,r11c,ir11e,r1d1c,ir1d1e,jbes1,nsub,
     4					 nsubd)
              jbes=jbes1
              iopbes=2
c              if(knd.eq.kindd) write(40,570) r11c,ir11e,r1d1c,ir1d1e
c              if(knd.eq.kindq) write(40,575) r11c,ir11e,r1d1c,ir1d1e
c570           format(10x,'r1 = ',f19.15,f19.15,i5,5x,'r1d = ',
c     1                f19.15,f19.15,i5)
c575           format(10x,'r1 = ',f35.31,f35.31,i5,/,10x,'r1d = ',
c     1                f35.31,f35.31,i5)
              r1c(li)=r11c
              ir1e(li)=ir11e
              r1dc(li)=r1d1c
              ir1de(li)=ir1d1e
              iopbesa=1
                if(l.eq.m.and.(nsub.gt.nsubt.or.nsubd.gt.nsubt.or.
     1              jsubmf.gt.nsubt)) then
                iopeta1=1
                iopbes=0
                idir=0
                nee1=1
                nsub1p=max(nsub,jsubmf)
                nsubd1p=max(nsubd,jsubmf)
                go to 580
                end if
                if(iopeta1.ne.0.or.isteta1.ne.0) then
                  if(nsub.le.nsubt.and.nsubd.le.nsubt.and.jsubmf.le.
     1               nsubt) then
                  iopeta1=0
                  if(li.gt.4*nbp/3) isteta1=0
                  else
                  iopbes=0
                  iopeta1=2
                  idir=0
                  nee1=1
                  iflageta1=0
                  nsub1p=max(nsub,jsubmf)
                  nsubd1p=max(nsubd,jsubmf)
                  end if
                go to 580
                end if
                if((li.le.lipl+1.or.(li.eq.lips.and.liplp.ne.lipl+1)
     1           .or.(li.eq.lips+1.and.liplp.ne.lipl+1)).and.l.ne.m
     2          .and.max(nsub,nsubd,jsubmf).gt.nsubt) then
                nee1=1
                idir=0
                iflageta1=0
                jflageta1=0
                iopeta1=1
                iopbes=0
                nsub1p=max(nsub,jsubmf)
                nsubd1p=max(nsubd,jsubmf)
                end if
580           continue
              naccr1=ndec-2-max(nsub,nsubd,jsubmf)
              naccr1=min(naccr1,naccre-1,itestm-1)
              naccr1c=min(naccre,itestm)-max(nsub,nsubd,jsubmf)
              if(naccr1.lt.0) naccr1=0
585           continue
c
c    r1 calculation using the variable eta expansion
              if(iopeta1.eq.0) go to 760
              liplopt=0
              nee1st=nee1
              nee1count=1
                if(li.eq.lips.and.liplp.ne.lipl+1.and.li.ne.1) then
                liplopt=1
                idir=0
                nsub1p=jsubmf
                nsubd1p=jsubmf
                nee1=1
                iflageta1=0
                jflageta1=0
                iopeta1=2
                end if
                if(ieta.eq.0) then
                  do jnet=1,neta
                  ang=jnet*pi*0.5e0_knd/(neta+1)
                  eta(jnet)=cos(ang)
                  wmeta2(jnet)=2.0e0_knd*(1.0e0_knd+eta(jnet))*
     1                         (sin(0.5e0_knd*ang)**2)
                  xbn(jnet)=sqrt(x*x+wmeta2(jnet))
                  xln(jnet)=eta(jnet)*x/xbn(jnet)
                  end do
                ieta=1
                end if
              if(iopeta1.eq.1) iopeta1=2
590           if(iopeta1.eq.3) go to 690
              etaval1=eta(nee1)
600           xbninp=xbn(nee1)
              netainp=1
              etainp(1)=eta(nee1)
              xlninp(1)=xln(nee1)
              limj=lnum+3*ndec+int(c)+m
              if(limj.lt.maxlp+1) limj=maxlp+1
              if(limj.gt.maxj) limj=maxj-2
              limp=limj-m
              if(jnen1.eq.0) go to 660
              jnenlim1=jnen1
              if(jnen1.gt.jnebmax) jnenlim1=jnebmax
              limplim=limp
              limjlim=limj
                do 650 jn=1,jnenlim1
                if(nee1.ne.neeb1(jn)) go to 650
                if(limplim.gt.limp1sv(jn)) limplim=limp1sv(jn)
                if(limjlim.gt.limjsv(jn)) limjlim=limjsv(jn)
                  do 610 je=1,limplim
                  pratb1(je)=pratbsv1(jn,je)
                  pratt1(je)=prattsv1(jn,je)
                  pdratt1(je)=pdrattsv1(jn,je)
610               continue
                  do 620 je=1,limjlim
                  sbesfe(je)=sbesfsv(jn,je)
                  sbesdfe(je)=sbesdfsv(jn,je)
620               continue
                  jelim=maxlp
                  if(maxlp.gt.limj+1) jelim=limj+1
                  if(jelim1.gt.jelim1sv(jn)) jelim1=jelim1sv(jn)
                  do 630 je=1,jelim1
                  sbesne(je)=sbesnsv(jn,je)
                  sbesdre(je)=sbesdrsv(jn,je)
                  ibesee(je)=ibesesv(jn,je)
630               continue
                go to 680
650             continue
660           continue
              jnen1=jnen1+1
              jnencur1=jnen1-(jnebmax*int((jnen1-1)/jnebmax))
              neeb1(jnencur1)=nee1
              call sphbes(cc,xbninp,limj,maxj,maxlp,sbesfe,sbesdfe,
     1                    sbesne,ibesee,sbesdre)
                do je=1,limj
                sbesfsv(jnencur1,je)=sbesfe(je)
                sbesdfsv(jnencur1,je)=sbesdfe(je)
                limjsv(jnencur1)=limj
                end do
                jelim1=maxlp
                if(maxlp.gt.limj+1) jelim1=limj+1
                  do 670 je=1,jelim1
                  sbesnsv(jnencur1,je)=sbesne(je)
                  sbesdrsv(jnencur1,je)=sbesdre(je)
670               ibesesv(jnencur1,je)=ibesee(je)
                  jelim1sv(jnencur1)=jelim1
              iopd=3
              call pleg(m,limp,maxp,limcsav,iopd,ndec,xlninp,netainp,
     1                  maxt,prat,pdrat,pdnorma,ipdnorma,pnorma,ipnorma,
     2                  alpha,beta,gamma,coefa,coefb,coefc,coefd,coefe)
              limcsav=max(limcsav,limp)
                do je=1,limp
                pratt1(je)=prat(1,je)
                pdratt1(je)=pdrat(1,je)
                prattsv1(jnencur1,je)=pratt1(je)
                pdrattsv1(jnencur1,je)=pdratt1(je)
                limp1sv(jnencur1)=limp
                end do
              limpd=2*(lnum+int(c)+ndec)
              if(limpd.gt.limp) limpd=limp
              iopd=2
              call pleg(m,limpd,maxp,limcsav,iopd,ndec,etainp,netainp,
     1                  maxt,prat,pdrat,pdnorma,ipdnorma,pnorma,
     2                  ipnorma,alpha,beta,gamma,coefa,coefb,coefc,
     3                  coefd,coefe)
                do je=1,limpd
                pratb1(je)=prat(1,je)
                pratbsv1(jnencur1,je)=pratb1(je)
                end do
              pratb1(limpd+1)=0.0e0_knd
              pratb1(limpd+2)=0.0e0_knd
              pratbsv1(jnencur1,limpd+1)=0.0e0_knd
              pratbsv1(jnencur1,limpd+2)=0.0e0_knd
680           continue
              pcoefe1=(x*x+1.0e0_knd)/(xbn(nee1)*xbn(nee1))
              apcoef1=(rm/2.0e0_knd)*log10(pcoefe1)
              ipcoefe1=int(apcoef1)
              pcoefe1=ten**(apcoef1-ipcoefe1)
              pcoefo1=pcoefe1*pratt1(2)/pratb1(2)
              ipcoefo1=ipcoefe1
              pdcoefe1=pcoefe1
              if(m.ne.0) pdcoefe1=-pcoefe1*rm*xln(nee1)*xbn(nee1)*
     1                    xbn(nee1)/((x*x+1.0e0_knd)*wmeta2(nee1))
              ipdcoefe1=ipcoefe1
              pdcoefo1=pdcoefe1*pdratt1(2)/pratb1(2)
              ipdcoefo1=ipdcoefe1
              if(li.eq.1) go to 690
                do jl=3,li+ix,2
                pcoefe1=pcoefe1*pratt1(jl)/pratb1(jl)
                iterm=log10(abs(pcoefe1))
                pcoefe1=pcoefe1*ten**(-iterm)
                ipcoefe1=ipcoefe1+iterm
                pdcoefe1=pdcoefe1*pdratt1(jl)/pratb1(jl)
                iterm=log10(abs(pdcoefe1))
                pdcoefe1=pdcoefe1*ten**(-iterm)
                ipdcoefe1=ipdcoefe1+iterm
                end do
              continue
              if(li.lt.3) go to 690
                do jl=4,li+1-ix,2
                pcoefo1=pcoefo1*pratt1(jl)/pratb1(jl)
                iterm=log10(abs(pcoefo1))
                pcoefo1=pcoefo1*ten**(-iterm)
                ipcoefo1=ipcoefo1+iterm
                pdcoefo1=pdcoefo1*pdratt1(jl)/pratb1(jl)
                iterm=log10(abs(pdcoefo1))
                pdcoefo1=pdcoefo1*ten**(-iterm)
                ipdcoefo1=ipdcoefo1+iterm
                end do
690           if(ix.eq.0) go to 700
              pcoefet1=pcoefo1
              ipcoefet1=ipcoefo1
              pcoefo1=pcoefo1*pratt1(li+2)/pratb1(li+2)
              iterm=int(log10(abs(pcoefo1)))
              pcoefo1=pcoefo1*ten**(-iterm)
              ipcoefo1=ipcoefo1+iterm
              pdcoefet1=pdcoefo1
              ipdcoefet1=ipdcoefo1
              pdcoefo1=pdcoefo1*pdratt1(li+2)/pratb1(li+2)
              iterm=int(log10(abs(pdcoefo1)))
              pdcoefo1=pdcoefo1*ten**(-iterm)
              ipdcoefo1=ipdcoefo1+iterm
              go to 710
700           pcoefet1=pcoefe1
              ipcoefet1=ipcoefe1
              pcoefe1=pcoefe1*pratt1(li+2)/pratb1(li+2)
              iterm=int(log10(abs(pcoefe1)))
              pcoefe1=pcoefe1*ten**(-iterm)
              ipcoefe1=ipcoefe1+iterm
              pdcoefet1=pdcoefe1
              ipdcoefet1=ipdcoefe1
              pdcoefe1=pdcoefe1*pdratt1(li+2)/pratb1(li+2)
              iterm=int(log10(abs(pdcoefe1)))
              pdcoefe1=pdcoefe1*ten**(-iterm)
              ipdcoefe1=ipdcoefe1+iterm
710           continue
              wm=wmeta2(nee1)
              limeta=l+3*ndec+int(c)
              if(iopeta1.eq.3) limeta=jeta1+jeta1+20
              if(limeta.gt.limp-2) limeta=limp-2
              if(limeta.gt.limd) limeta=limd
              limeta=min(limeta,maxj-m-2)
              call r1eta(l,m,cc,x,etaval1,nee1,limeta,ndec,nex,maxd,
     1                   maxlp,maxj,maxp,minacc,wm,enr(li,:),sbesfe,
     2                   sbesne,ibesee,sbesdfe,sbesdre,pdratt1,pratb1,
     3                   pratt1,pcoefet1,ipcoefet1,pdcoefet1,ipdcoefet1,
     4                   ir1ep,r1ec,ir1ee,r1dec,ir1dee,nsub1,nsubd1,
     5					 jeta1)
              if(iopbes.eq.0) jbes=jeta1
              if(iopbes.ne.0) jbes=max(jbes,jeta1)
720           continue
c              if(knd.eq.kindd) write(40,730) etaval1,nee1,r1ec,ir1ee,
c     1                                       r1dec,ir1dee
c              if(knd.eq.kindq) write(40,735) etaval1,nee1,r1ec,ir1ee,
c     1                                     r1dec,ir1dee
c730           format(15x,'eta = ',f12.9,'; nee = ',i4,/,10x,'r1 = ',
c     1                f19.15,f19.15,i5,5x,'r1d = ',f19.15,f19.15,i5)
c735           format(15x,'eta = ',f12.9,'; nee = ',i4,/,10x,'r1 = ',
c     1                f35.31,f35.31,i5,/,5x,'r1d = ',f35.31,f35.31,i5)
                if(nsub1.le.nsubt.or.nsubd1.le.nsubt) then
                if(idir.eq.0) idir=-1
                iflageta1=0
                nsub1p=nsub1
                nsubd1p=nsubd1
                r1c(li)=r1ec
                ir1e(li)=ir1ee
                r1dc(li)=r1dec
                ir1de(li)=ir1dee
                  if(nee1.eq.1.and.nsub1.gt.jsubmf+1.and.nsubd1.gt.
     1               jsubmf+1.and.li.gt.lipl.and.l.gt.nbp) then
                  iopbes=1
                  go to 750
                  end if
                if(jflageta1.ne.0) jflageta1=jflageta1-1
                if(jflageta1.eq.0) iopeta1=3
                go to 750
                end if
                if(idir.eq.0) then
                  if(nee1.eq.1) then
                  r1cm=r1ec
                  ir1em=ir1ee
                  r1dcm=r1dec
                  ir1dem=ir1dee
                  nsub1m=nsub1
                  nsubd1m=nsubd1
                  nee1m=1
                    if(iopbesa.eq.1.and.nsub1.lt.nsub1p.and.nsubd1.lt.
     1              nsubd1p) then
                    r1c(li)=r1ec
                    ir1e(li)=ir1ee
                    r1dc(li)=r1dec
                    ir1de(li)=ir1dee
                    end if
                    if(iopbesa.eq.1.and.nsub1.gt.nsub1p+1.and.nsubd1.gt.
     1              nsubd1p+1.and.li.gt.4*nbp/3) then
                    iopeta1=0
                    isteta1=0
                    go to 750
                    end if
                  end if
                  if(nee1.ne.1.and.((nsub1.lt.nsub1m.and.nsubd1.le.
     1               nsubd1m).or.(nsubd1.lt.nsubd1m.and.nsub1.le.
     2               nsub1m))) then
                  r1cm=r1ec
                  ir1em=ir1ee
                  r1dcm=r1dec
                  ir1dem=ir1dee
                  nsub1m=nsub1
                  nsubd1m=nsubd1
                  nee1m=nee1
                  end if
                  if(nee1.gt.neta-incnee1.or.((nsub1.gt.nsub1m+1.or.
     1              nsubd1.gt.nsubd1m+1).and.min(nsub1m,nsubd1m).lt.
     2              nsub1mt).or.(nsub1.gt.nsub1m+3.and.nsubd1.gt.
     3              nsubd1m+3)) then
                    r1c(li)=r1cm
                    ir1e(li)=ir1em
                    r1dc(li)=r1dcm
                    ir1de(li)=ir1dem
                    nsub1=nsub1m
                    nsubd1=nsubd1m
                    nsub1p=nsub1
                    nsubd1p=nsubd1
                    nee1=nee1m
                    iopeta1=2
                    idir=-1
                    if(nee1.eq.1) idir=0
                    go to 750
                  end if
                  if(nee1.gt.nee1m+3.and.((nsub1.le.nsub1m+1.and.
     1               nsubd1.le.nsubd1m+1).and.min(nsub1m,nsubd1m)
     2               .lt.nsub1mt)) then
                    r1c(li)=r1cm
                    ir1e(li)=ir1em
                    r1dc(li)=r1dcm
                    ir1de(li)=ir1dem
                    nsub1=nsub1m
                    nsubd1=nsubd1m
                    nsub1p=nsub1
                    nsubd1p=nsubd1
                    nee1=nee1m
                    iopeta1=2
                    idir=-1
                    go to 750
                  end if
                nee1=nee1+incnee1
                iopeta1=2
                nsub1p=nsub1
                nsubd1p=nsubd1
                go to 590
                end if
740             if((nsub1.gt.nsub1p.or.nsubd1.gt.nsubd1p).and.
     1               iflageta1.eq.1)
     2               then
                iflageta1=0
                nsub1=nsub1p
                nsubd1=nsubd1p
                iopeta1=2
                jflageta1=2
                nee1=nee1+incnee1*(1+nee1count)
                if(nee1.gt.neta) nee1=neta
                go to 750
                end if
                if(max(nsub1,nsubd1).eq.max(nsub1p,nsubd1p).and.nee1.ne.
     1              nee1st) then
                nee1count=nee1count+1
                else
                nee1count=1
                end if
              r1c(li)=r1ec
              ir1e(li)=ir1ee
              r1dc(li)=r1dec
              ir1de(li)=ir1dee
              nsub1p=nsub1
              nsubd1p=nsubd1
              iflageta1=1
              nee1=nee1-incnee1
              iopeta1=2
                if(nee1.eq.0.and.nee1st.eq.1) then
                nee1=2
                idir=0
                r1cm=r1ec
                ir1em=ir1ee
                r1dcm=r1dec
                ir1dem=ir1dee
                nsub1m=nsub1
                nsubd1m=nsubd1
                nee1m=1
                iflageta1=0
                go to 590
                end if
                if(nee1.eq.0.and.nee1st.ne.1) then
                  if(jsubmf.le.min(nsub1p,nsubd1p)) then
                  iopbes=1
                  iopeta1=2
                  else
                  nee1=1+nee1count*incnee1
                  if(nee1.gt.neta) nee1=neta
                  iopeta1=2
                  iflageta1=0
                  end if
                go to 750
                end if
              go to 590
750           continue
              nacetr1=ndec-2-max(nsub1,nsubd1)
              nacetr1=min(nacetr1,naccre-1,itestm-1)
              nacetr1c=min(naccre,itestm)-max(nsub1,nsubd1)
              if(nacetr1.lt.0) nacetr1=0
                if(iopbesa.eq.1.and.naccr1c.ge.nacetr1c) then
                iopbes=1
                r1c(li)=r11c
                ir1e(li)=ir11e
                r1dc(li)=r1d1c
                ir1de(li)=ir1d1e
                else
                naccr1=nacetr1
                naccr1c=nacetr1c
                end if
                if(nee1.eq.1.and.nsub1.gt.jsubmf+1.and.
     1             nsubd1.gt.jsubmf+1.and.li.gt.lipl.and.l.gt.nbp) then
                iopbes=1
                end if
                if(nsub1p.ge.nsub1mt-1.and.nsubd1p.ge.nsub1mt-1) then
                nee1=1
                idir=0
                iflageta1=0
                jflageta1=0
                iopeta1=1
                end if
760           continue
              if(ioprad.eq.1) go to 1330
c
c  determine oblate radial functions of the second kind
c
c              write(40,770)
c770           format(4x,'r2 and r2d calculation')
c
c  decide whether to use values of r1 and r1d for r2 and r2d
              naccrpl=ir1e(li)+ir1de(li)+int(log10(abs(r1c(li)*
     1                r1dc(li))*c*(x*x+1.0e0_knd)))
              if(naccrpl.lt.0) naccrpl=0
              if(naccrpl.gt.ndec) naccrpl=ndec
              if(min(naccrpl,naccr1).gt.0) then
                naccr=min(naccrpl,naccr1)
                r2c(li)=cmplx(-aimag(r1c(li)),real(r1c(li)),knd)
                ir2e(li)=ir1e(li)
                r2dc(li)=cmplx(-aimag(r1dc(li)),real(r1dc(li)),knd)
                ir2de(li)=ir1de(li)
                end if
              naccrt=0
              iopmatch=1
              nmatch=0
              iii=-int(log10(x)-0.699e0_knd)
              if(iii.lt.0) iii=0
                if(li.le.lipl.and.match(2).ne.0) then
                  if(ix.eq.0) then
                  nmatch=-int(log10(abs((eig(li)-eigt(li+1))/
     1                              (eig(li))+ten*dec)))
                  if(nmatch.gt.ndec) nmatch=ndec
                  if(min(nmatch-2,naccr1).ge.minacc.and.li.ge.listart)
     1               listart=li+2
                  end if
                  if(ix.eq.1) then
                  nmatch=-int(log10(abs((eig(li)-eig(li-1))/
     1                              (eig(li))+ten*dec)))
                  if(nmatch.gt.ndec) nmatch=ndec
                  if(min(nmatch-3,naccr1).ge.minacc.and.li.ge.
     1                 listart-1) listart=li+3
                  end if
                  if(li.ge.listart.or.(knd.eq.kindq.and.match(li+1-ix)
     1               -iii.lt.minacc).or.(knd.eq.kindd.and.match(li+1-ix)
     2               -iii.lt.4)) then
                    if(min(naccrpl,naccr1).ge.minacc) then
                    iopint=0
                    iopleg=0
                    iopleg1=0
                    iopneu0=0
                    iopeta=0
                    naccr=min(naccrpl,naccr1)
c                    write(40,775) l,naccr
c                    write(20,1350) l,r1c(li),ir1e(li),r1dc(li),
c     1                             ir1de(li),r2c(li),ir2e(li),r2dc(li),
c     2                             ir2de(li),naccr
                    go to 1400
                    end if
                    if(ix.eq.0) then
                    naccrt=min(match(li+1)-2,naccr1)
                    naccr=naccrt
                      if(min(naccrpl,naccr1).ge.max(naccrt,2)) then
                      naccr=min(naccrpl,naccr1)
                      iopmatch=0
c                      if(naccr.gt.0) write(40,775) l,naccr
c775                   format(8x,'values for r2 and r2dc for l = ',i6,
c     1                       ' accurate to ',i2,' digits are given by'
c     2                       ' ir1 and ir1d')
                      end if
                    end if
                    if(ix.eq.1) then
                    wroncb=r1c(li)*r1dc(li-1)*ten**(ir1e(li)+
     1                     ir1de(li-1))
                    wronca=r1c(li-1)*r1dc(li)*
     1                     ten**(ir1e(li-1)+ir1de(li))
                    wronc=wronca-wroncb
                    naccrt=-int(log10(abs((wronc-wront)/wront)+dec))
                    if(naccrt.gt.ndec) naccrt=ndec
                    if(naccrt.lt.0) naccrt=0
                    nacccor=-int(log10(abs((wronca-wroncb)/wronca)+
     1                      dec))
                    if(nacccor.lt.0) nacccor=0
                    naccrt=naccrt+nacccor
                    if(naccrt.gt.ndec) naccrt=ndec
                    naccrt=min(naccrt,nmatch,naccr1)
                    naccr=naccrt
                      if(min(naccrpl,naccr1).ge.max(naccrt,1)) then
                      naccr=min(naccrpl,naccr1)
                      iopmatch=0
c                      write(40,775) l,naccr
                      end if
                    end if
                  end if
                  if(li.lt.listart.and.((knd.eq.kindq.and.
     1               match(li+1-ix)-iii.ge.minacc).or.(knd.eq.kindd
     2               .and.match(li+1-ix)-iii.ge.4))) then
                    if(ix.eq.0) then
                    naccrt=min(match(li+1)-2,naccr1)
                    if(naccrt.lt.0) naccrt=0
                    naccr=naccrt
                    end if
                  go to 1390
                  end if
                if(x.le.0.99e0_knd.and.iopleg.eq.0.and.l.ge.legstart
     1             .and.ndec-jsub.ge.naccr) iopleg=1
                if(x.lt.0.01e0_knd.and.iopleg.eq.0.and.l.ge.legstart
     1             .and.ndec-jsub.ge.match(li+1-ix)-iii) iopleg=1
                if(x.ge.0.05e0_knd.and.x.le.0.9e0_knd.and.jsub.gt.
     1             ndec-minacc.and.iopleg.eq.0.and.iopeta.eq.0) iopeta=4
                if(li.eq.listart.and.iopeta.eq.4) iopeta=1
                if(x.gt.0.99e0_knd.and.iopneu0.eq.0.and.jsubf.lt.
     1             ndec-naccr) iopneu0=4
                if(li.eq.listart.and.iopneu0.eq.4) iopneu0=1
                if(x.le.xhigh.and.x.ge.xlow.and.iopint.eq.0.and.
     1              iopmatch.eq.1.and.match(li+1-ix)-iii.lt.minacc
     2                  .and.iopleg.eq.0) iopint=1
                if(nacclegp.gt.minacc) iopint=0
                end if
                if(li.gt.lipl.and.min(naccrpl,naccr1).gt.1) then
                naccr=min(naccrpl,naccr1)
                iopmatch=0
c                if(naccr.gt.0) write(40,775) l,naccr
                  if(naccr.ge.minacc) then
c                  write(20,1350) l,r1c(li),ir1e(li),r1dc(li),
c     1                           ir1de(li),r2c(li),ir2e(li),r2dc(li),
c     2                           ir2de(li),naccr
                  iopint=0
                  iopleg=0
                  iopleg1=0
                  if(min(naccrpl,naccr1).gt.minacc+2) iopneu0=0
                  iopeta=0
                  go to 1400
                  end if
                end if
                if(li.gt.lipl.or.match(2).eq.0) then
                if(x.le.0.99e0_knd.and.iopleg.eq.0.and.l.ge.legstart
     1             .and.ndec-jsub.gt.min(minacc,naccrsav-2)) iopleg=1
                if(x.gt.0.99e0_knd.and.ndec-jsubf.ge.naccrsav.and.
     1             iopneu0.eq.0.and.istartneu0.eq.1) iopneu0=4
                if(iopeta.ne.0.and.iopneu0.eq.4) iopneu0=1
                if(x.ge.0.05e0_knd.and.iopint.eq.0.and.iopleg.eq.0.and.
     1              iopeta.eq.0.and.iopneu0.eq.0) iopeta=4
                if(x.le.xhigh.and.x.ge.xlow.and.iopint.eq.0) iopint=1
                if(nacclegp.gt.minacc) iopint=0
                end if
              r1cin=r1c(li)
              ir1ein=ir1e(li)
              r1dcin=r1dc(li)
              ir1dein=ir1de(li)
780           continue
c
c     r2 calculation using integration technique
              if(iopint.eq.0.or.istartint.eq.2) go to 830
                if(li.gt.intlim-5) then
                istartint=2
                iopint=0
                go to 830
                end if
              if(iopint.gt.1) go to 800
              limint=2*(nbp+33)+50*aimag(cc)+2
              if(igau.eq.1) go to 790
              call gauss(ngau,xr,wr)
              igau=1
790           if(iopint.eq.1.and.ipint.eq.1) iopint=2
              if(ipint.eq.1) go to 800
              ngqs=nstep1+nstep2+nstep3
              parg=cc*sqrt(x*x+1.0e0_knd)
              lim1max=4*int(abs(real(parg))+abs(aimag(parg)))+25
              call pint(cc,m,lnum,x,limint,maxint,maxlp,maxmp,lim1max,
     1                  ndec,nex,ngqs,ngau,wr,xr,step1,nstep1,step2,
     2                  nstep2,step3,nstep3,intlim,rpint1,rpint2,pint1,
     3                  pint2,pint3,pint4,norme,pnormint,ipnormint,
     4                  coefme,coefmo)
              iopint=2
800           continue
              if(iopint.eq.2) limint=2*(nbp+33)+50*aimag(cc)
              if(iopint.eq.3) limint=jint+jint+20+int(sqrt(c))
              if(limint.gt.intlim) limint=intlim-1
              if(2*(limint/2).ne.limint) limint=limint-1
              call r2int(l,m,cc,x,limint,ndec,nex,maxd,enr(li,:),dc01,
     1                   idc01,maxint,maxmp,maxlp,intlim,rpint1,rpint2,
     2                   pint1,pint2,pint3,pint4,norme,pnormint,
     3                   ipnormint,coefme,coefmo,ipint,r2ic,ir2ie,
     4                   r2dic,ir2die,jint,coefn,icoefn,isub,isubd)
              if(ipint.eq.0) ipint=1
              naccrs=naccr
              wronca=r1c(li)*r2dic*ten**(ir1e(li)+ir2die)
              wroncb=r2ic*r1dc(li)*ten**(ir2ie+ir1de(li))
              wronc=wronca-wroncb
              naccint=-int(log10(abs((wronc-wront)/wront)+dec))
              if(naccint.lt.0) naccint=0
              if(naccint.gt.ndec) naccint=ndec
              nacccor=-int(log10(abs((wronca-wroncb)/wronca)+dec))
              if(nacccor.lt.0) nacccor=0
              if(nacccor.gt.naccrpl) nacccor=naccrpl
              if(nacccor.gt.0) naccint=min(naccint+nacccor,naccr1)
              if(x.gt.0.1e0_knd) go to 805
              naccintw=naccint
              jacc=0
              if((wronca*wroncb).ne.0.0e0_knd)
     1              jacc=int(log10(abs(wronca/wroncb)))
              if(jacc.gt.0) naccint=min(naccint,ndec-isub-1)
              if(jacc.lt.0) naccint=min(naccint,ndec-isubd-1)
              if(jacc.eq.0) naccint=min(naccint,ndec-max(isub,isubd)-1)
              naccc=0
                if(jacc.gt.0.and.(ir2ie-ir1e(li)).lt.nex-10) then
                naccc=int(log10(abs(r2ic/(r1c(li)*
     1                 (10.0e0_knd**(ir1e(li)-ir2ie))+
     2                 (0.0e0_knd,1.0e0_knd)*r2ic)))-0.5e0_knd)
                if(naccc.gt.0) naccc=0
                end if
                if(jacc.lt.0.and.(ir2die-ir1de(li)).lt.nex-10) then
                naccc=int(log10(abs(r2dic/(r1dc(li)*
     1                 (10.0e0_knd**(ir1de(li)-ir2die))+
     2                 (0.0e0_knd,1.0e0_knd)*r2dic)))-0.5e0_knd)
                if(naccc.gt.0) naccc=0
                end if
              naccint=min(naccintw,naccint-naccc)
              iopint=3
              if(naccint.lt.0) naccint=0
805           continue
                if(naccint.gt.naccr.or.(naccint.eq.naccr.and.naccr
     1             .eq.naccrt.and.x.lt.0.1e0_knd)) then
                naccr=naccint
                jflagl=0
                r2c(li)=r2ic
                ir2e(li)=ir2ie
                r2dc(li)=r2dic
                ir2de(li)=ir2die
                iopmatch=0
                end if
810           continue
                if(naccint.ge.minacc) iflagint=1
                if(naccint.gt.minacc+1.and.naccintp.gt.minacc+1) then
                iopleg1=0
                iopleg=0
                iopneu0=0
                iopeta=0
                end if
                if(naccint.le.minacc) then
                if(x.le.0.99e0_knd.and.ndec-jsub.gt.min(naccint,
     1             naccrs,naccrsav-2).and.l.ge.legstart.and.iopleg.eq.0)
     2             iopleg=1
                if(x.ge.xneu) iflagneu=1
                if(iflagneu.eq.1.and.(ndec-jsubf).ge.minacc.and.
     1             iopneu0.eq.0.and.iopeta.eq.0.and.istartneu0.eq.1)
     2               iopneu0=4
                if(iopneu0.eq.4.and.(l.eq.listart.or.iopeta.ne.0))
     1               iopneu0=1
                if(iflagneu.eq.1.and.iopneu0.eq.0.and.iopeta.eq.0
     1             .and.iflageta.eq.0.and.x.ge.0.05e0_knd) iopeta=4
                if(iopeta.eq.4.and.l.eq.listart) iopeta=1
                end if
              naccintp=naccint
c              if(knd.eq.kindd) write(40,820) naccint,r2ic,ir2ie,r2dic,
c     1                                       ir2die
c              if(knd.eq.kindq) write(40,825) naccint,r2ic,ir2ie,r2dic,
c     1                                       ir2die
c820           format(15x,'accuracy in decimal digits = ',i2,/,10x,
c     1               'r2 = ',f19.15,f19.15,i5,5x,'r2d = ',
c     2               f19.15,f19.15,i5)
c825           format(15x,'accuracy in decimal digits = ',i2,/,10x,
c     1               'r2 = ',f35.31,f35.31,i5,/10x,'r2d = ',f35.31,
c     2               f35.31,i5)
                if(istartint.eq.1) then
                istartint=0
                if(naccint.eq.0.and.li.gt.max(nbp,liplp)) istartint=2
                go to 1320
                end if
830           continue
c
c     r2 calculation using a Legendre expansion and joining factor
              if(l.lt.legstart) iopleg=0
              if(ndec-jsub.lt.3) iopleg=0
              if(ndec-jsub.lt.max(naccr,nmatch-iii))
     1             iopleg=0
                if(iopleg.eq.0.or.(naccr.ge.minacc.and.naccr.eq.naccint)
     1             .or.(naccr.eq.naccrt.and.nmatch-2.ge.minacc.and.
     2             x.ge.0.01e0_knd)) then
                if(iopleg.eq.2) iopleg=1
                go to 940
                end if
              if(jflagleg.eq.1) go to 900
              jflagleg=1
              limdr=c+2*ndec+50*x+2
              if(limdr.gt.maxdr) limdr=maxdr
              if(ioppsum.eq.0) go to 840
              xin(1)=x
              limpleg=limdr+limdr
              if(limpleg.gt.maxp-2) limpleg=maxp-2
              iopd=4
              call pleg(m,limpleg,maxp,limcsav,iopd,ndec,xin,1,maxt,
     1                  prat,pdrat,pdnorma,ipdnorma,pnorma,ipnorma,
     2                  alpha,beta,gamma,coefa,coefb,coefc,coefd,coefe)
              limcsav=max(limcsav,limpleg)
                do jj=1,limpleg
                prx(jj)=prat(1,jj)
                pdrx(jj)=pdrat(1,jj)
                end do
840             if(iqlegflag.eq.0) then
                limq=lnum+3*ndec+int(c)+2
                if(knd.eq.kindd.and.aimag(cc).lt.10.0e0_knd.and.
     1             real(cc).le.60.0e0_knd.and.real(cc).ge.10.0e0_knd
     2             .and.m.le.40.and.x.le.0.99e0_knd.and.x.gt.0.1e0_knd)
     3            limq=max(limq,250-int(50*x)+2)
                call qleg(m,lnum,limq,maxmp,maxq,x,ndec,nex,iflagl1,
     1               qdrat,qdqr,qdm1,iqdm1,qdl,iqdl,qrat,qm1,iqm1,ql,
     2               iql,termpq,itermpq,qrat1,qdrat1,qm0,qdm0)
                iqlegflag=1
                end if
              fajo(1)=cc/(rm2-1.0e0_knd)
              ifajo(1)=0
              if(m.eq.0) go to 870
                do im=1,m
                fajo(1)=fajo(1)*(im+im)/cc
                if(abs(fajo(1)).lt.1.0e+10_knd) go to 850
                fajo(1)=fajo(1)*(1.0e-10_knd)
                ifajo(1)=ifajo(1)+10
850             continue
                if(abs(fajo(1)).gt.1.0e-10_knd) go to 860
                fajo(1)=fajo(1)*(1.0e+10_knd)
                ifajo(1)=ifajo(1)-10
860             continue
                end do
870           continue
              if(lnum.eq.1) go to 900
              fajo(2)=-cc*fajo(1)/(rm2-3.0e0_knd)
              ifajo(2)=ifajo(1)
              if(lnum.eq.2) go to 900
                do jl=3,lnum-1,2
                fajo(jl)=fajo(jl-2)*(jl+m+m-1)/(jl-2)
                ifajo(jl)=ifajo(jl-2)
                  if(abs(fajo(jl)).gt.1.0e10_knd) then
                  fajo(jl)=fajo(jl)*1.0e-10_knd
                  ifajo(jl)=ifajo(jl)+10
                  end if
                fajo(jl+1)=fajo(jl-1)*(jl+m+m-1)/(jl)
                ifajo(jl+1)=ifajo(jl-1)
                  if(abs(fajo(jl+1)).gt.1.0e10_knd) then
                  fajo(jl+1)=fajo(jl+1)*1.0e-10_knd
                  ifajo(jl+1)=ifajo(jl+1)+10
                  end if
                end do
                if(2*(lnum/2).ne.lnum.and.lnum.ge.4) then
                fajo(lnum)=fajo(lnum-2)*(lnum+m+m-1)/(lnum-2)
                ifajo(lnum)=ifajo(lnum-2)
                end if
900           continue
              limleg=l-m+3*ndec+int(c)
              limdr=int(c)+2*ndec+int(50*x)
              if(iopleg.eq.2) limleg=jleg+jleg+20+int(sqrt(c))
              if(iopleg.eq.2) limdr=jlegp+10+int(sqrt(c)/2)
              if(limdr.gt.maxdr-4) limdr=maxdr-4
              if(limleg.gt.limq-4) limleg=limq-4
              nsdrhor1=0
              nsdneg=0
              nsdrho=0
              call dalt(l,m,cc,ndec,nex,limdr,maxdr,maxmp,ioppsum,
     1                  eigval,enrneg,drhor,dneg,idneg,nsdneg,nsdrhor1,
     2                  nsdrho)
              kflagl=0
              ifsub=max(jsub,nsdneg)
              call r2leg(l,m,cc,x,lnum,minacc,limleg,limdr,iflagp,ndec,
     1                   nex,maxd,maxmp,maxpdr,maxdr,maxq,enr(li,:),
     2                   enrneg,drhor,nsdrhor1,nsdrho,dc01,idc01,dneg,
     3                   idneg,nsdneg,dfnorm,idfe,dmfnorm,idmfe,prx,
     4                   pdrx,qdrat,qdqr,qdm1,iqdm1,qdl,iqdl,qrat,qm1,
     5                   iqm1,ql,iql,fajo,ifajo,ifsub,jsub,termpq,
     6                   itermpq,ioppsum,iopqnsum,r1cin,ir1ein,r1dcin,
     7                   ir1dein,naccr1,naccrpl,itestm,r2lc,ir2le,r2dlc,
     8                   ir2dle,jleg,jlegp,jflagl,naccleg,kflagl,isub,
     9                   isubd,nacccor)
              iopleg=2
                if(naccleg.gt.naccr.or.(naccleg.eq.naccr.and.x
     1              .lt.0.1e0_knd)) then
                naccr=naccleg
                r2c(li)=r2lc
                ir2e(li)=ir2le
                r2dc(li)=r2dlc
                ir2de(li)=ir2dle
                iopmatch=0
                end if
910           continue
              nacclegs=naccleg
              if(kflagl.eq.1) nacclegs=0
                if(naccleg.gt.minacc+1.and.nacclegp
     1            .gt.minacc+1.and.li.gt.max(liplp+10,nbp)) then
                istartint=2
                iopint=0
                end if
                if(naccleg.lt.minacc) then
                if(knd.eq.kindd.and.aimag(cc).lt.10.0e0_knd.and.
     1             real(cc).le.60.0e0_knd.and.real(cc).ge.10.0e0_knd
     2             .and.m.le.40.and.x.le.0.99e0_knd.and.x.gt.0.1e0_knd
     3             .and.li.lt.nbp.and.naccr.lt.minacc) iopleg1=1
                if(iopleg1.eq.0.and.x.ge.xneu.and.iopneu0.eq.0.and.
     1             istartneu0.eq.1.and.iopeta.eq.0.and.ndec-jsubf.gt.
     2             naccr) iopneu0=4
                if(iopneu0.eq.4.and.(li.eq.listart.or.iopeta.ne.0))
     1               iopneu0=1
                end if
                if(naccleg.lt.naccrsav-2.and.naccleg.lt.naccrsavp-2.and.
     1             naccleg.lt.minacc.and.li.gt.2) then
                if(iopint.ne.0.and.naccint.gt.naccleg+2.and.
     1             naccint.gt.nacclegp+2.and.istartleg.eq.1) iopleg=0
                if(iopint.eq.0) itest=naccrsav
                if(iopint.ne.0) itest=min(naccint,naccrsav)
                  if(iopleg1.eq.0.and.istartleg.eq.1) then
                  legstinc=min(naccrsav,naccrsavp,minacc)-
     1               max(naccleg,nacclegp)
                  if(legstinc.lt.4) legstinc=1
                  legstart=l+legstinc
                  if(abs(aimag(cc)).ge.5.0e0_knd.and.((ir1e(li)+
     1               ir1de(li)).gt.2.or.li.lt.nbp)) legstart=l+1
                  if(legstart.lt.l+1) legstart=l+1
                  if(legstart.eq.l+1.and.iopleg.eq.0) iopleg=1
                  end if
                if(iopleg1.ne.0.and.istartleg.eq.1) legstart=
     1             l+minacc-nacclegs
                end if
                if(naccleg.ge.minacc.and.nacclegp.ge.minacc) then
                iopleg1=0
                iopneu0=0
                iopeta=0
                irtest=max(ir1e(li),ir1de(li))
                  if(li.gt.max(nbp,liplp).and.irtest.lt.-10) then
                  iopint=0
                  istartleg=0
                  istartneu0=0
                  end if
                end if
                if(naccleg.eq.minacc.and.nacclegp.lt.minacc) then
                  if(iopneu0.ne.0.and.istartneu0.eq.1)
     1               then
                  iopneu0=4
                  if(x.ge.0.05e0_knd) iopeta=4
                  end if
                end if
c920           if(knd.eq.kindd) write(40,820) naccleg,r2lc,ir2le,r2dlc,
c     1                                       ir2dle
c              if(knd.eq.kindq) write(40,825) naccleg,r2lc,ir2le,r2dlc,
c     1                                       ir2dle
              nacclegp=naccleg
              if(naccleg.le.0) iopleg=0
940           continue
c
c     r2 calculation using the Baber and Hasse Legendre expansion
              if(iopleg1.eq.0) go to 960
                if(iqlegflag.eq.0) then
                limq=max(lnum+3*ndec+int(c)+2,250-int(50*x)+2)
                call qleg(m,lnum,limq,maxmp,maxq,x,ndec,nex,iflagl1,
     1                    qdrat,qdqr,qdm1,iqdm1,qdl,iqdl,qrat,qm1,iqm1,
     2                    ql,iql,termpq,itermpq,qrat1,qdrat1,qm0,qdm0)
                iqlegflag=1
                end if
              if(iopleg1.eq.1) limleg1=250-int(50*x)
              if(iopleg1.ne.1) limleg1=jleg1+20
              if(limleg1.gt.250-int(50*x)) limleg1=250-
     1                                   int(50*x)
              call r2leg1(l,m,cc,x,limleg1,maxq,ndec,eigval,qrat1,
     1                    qdrat1,qm0,qdm0,r1cin,ir1ein,r2l1c,ir2l1e,
     2                    r2dl1c,ir2dl1e,jleg1)
              wronca=r1c(li)*r2dl1c*ten**(ir1e(li)+ir2dl1e)
              wroncb=r2l1c*r1dc(li)*ten**(ir2l1e+ir1de(li))
              wronc=wronca-wroncb
              naccleg1=-int(log10(abs((wronc-wront)/wront)+dec))
              if(naccleg1.lt.0) naccleg1=0
              if(naccleg1.gt.ndec) naccleg1=ndec
                if(naccleg1.gt.0) then
                nacccor=-int(log10(abs((wronca-wroncb)/wronca)+dec))
                if(nacccor.lt.0) nacccor=0
                if(nacccor.gt.naccrpl) nacccor=naccrpl
                naccleg1=min(naccleg1+nacccor,naccr1)
                end if
              if(iopleg.eq.2.and.naccleg.ge.naccleg1) iopleg1=0
                if(naccleg1.lt.minacc+1.and.iopleg.eq.0) then
                iopleg=1
                end if
              if(naccleg1.ge.minacc.and.iopneu0.ne.0.and.iopeta.eq.0)
     1            iopneu0=4
                naccleg1p=naccleg1
                if(naccleg1.gt.naccr) then
                r2c(li)=r2l1c
                ir2e(li)=ir2l1e
                r2dc(li)=r2dl1c
                ir2de(li)=ir2dl1e
                naccr=naccleg1
                iopleg1=2
                iopmatch=0
                end if
c              if(knd.eq.kindd) write(40,820) naccleg1,r2l1c,ir2l1e,
c     1                                       r2dl1c,ir2dl1e
c              if(knd.eq.kindq) write(40,825) naccleg1,r2l1c,ir2l1e,
c     1                                       r2dl1c,ir2dl1e
960           continue
c
c     r2 calculation using Neumann function expansion with eta set
c     equal to 0
              if(iopneu0.eq.0.or.istartneu0.eq.0.or.ndec-jsubf.le.naccr
     1           .or.(naccr.ge.minacc.and.naccr.ne.naccrt.and.x.lt.
     2           0.1e0_knd)) go to 1080
              if(iopneu0.gt.1) go to 1050
              if(ibflag2.eq.0) go to 1040
                if(knd.eq.kindd) then
                if(x.ge.0.01e0_knd) limn=2*int(25/(x*x)+300/x+
     1                                   3*c+10000)+2
                if(x.ge.0.1e0_knd) limn=2*int((lnum+c/5+0.5e0_knd*
     1                                   maxm+200)*1.4e0_knd/x)+2
                if(x.ge.0.5e0_knd) limn=2*int((lnum+c/5+0.5e0_knd*
     1                                   maxm+200)/x)+2
                if(x.ge.1.0e0_knd) limn=2*int(lnum+c/5+0.5e0_knd*maxm+
     1                                  200)+2
                end if
                if(knd.eq.kindq) then
                if(x.ge.0.01e0_knd) limn=2*int(25/(x*x)+400/x+
     1                                   3*c+14000)+2
                if(x.ge.0.1e0_knd) limn=2*int((lnum+c/5+0.5e0_knd*
     1                                   maxm+350)*1.4e0_knd/x)+2
                if(x.ge.0.5e0_knd) limn=2*int((lnum+c/5+0.5e0_knd*
     1                                   maxm+300)/x)+2
                if(x.ge.1.0e0_knd) limn=2*int(lnum+c/5+0.5e0_knd*maxm+
     1                                  300)+2
                end if
              limn=limn+maxm
              if(limn.gt.maxn) limn=maxn-2
              limbesf=4*int(real(cc*xb)+abs(aimag(cc*xb)))+25
              call sphneu(cc,xb,limn,maxn,maxlp,limbesf,sneuf2,sneun2,
     1                    ineue2,sneudf2,sneudr2)
              ibflag2=0
1040          continue
              iopneu0=iopneu0+1
1050          continue
              if(iopneu0.eq.4) go to 1080
                if(knd.eq.kindd) then
                if(x.ge.0.01e0_knd) limneu0=2*int(25/(x*x)+300/x+3*c+
     1                                     10000)
                if(x.ge.0.1e0_knd) limneu0=2*int((l-m+c/5+0.5e0_knd*m+
     1                                     200)*1.4e0_knd/x)
                if(x.ge.0.5e0_knd) limneu0=2*int((l-m+c/5+0.5e0_knd*m+
     1                                     200)/x)
                if(x.ge.1.0e0_knd) limneu0=2*int(l-m+c/5+0.5e0_knd*m+
     1                                     200)
                end if
                if(knd.eq.kindq) then
                if(x.ge.0.01e0_knd) limneu0=2*int(25/(x*x)+400/x+3*c+
     1                                      14000)
                if(x.ge.0.1e0_knd) limneu0=2*int((l-m+c/5+0.5e0_knd*m+
     1                                     350)*1.4e0_knd/x)
                if(x.ge.0.5e0_knd) limneu0=2*int((l-m+c/5+0.5e0_knd*m+
     1                                     300)/x)
                if(x.ge.1.0e0_knd) limneu0=2*int(l-m+c/5+0.5e0_knd*m+
     1                                     300)
                end if
              if(iopneu0.eq.3.and.jtest0.ge.minacc)
     1            limneu0=jneu0max+jneu0max+40+
     2                    int(sqrt(c)*(2/min(1.0e0_knd,x))+100/x)
              if(limneu0.gt.limd-2) limneu0=limd-2
              limneu0=min(limneu0,limn-m-2)
              iopneu0=3
              call r2neu0(l,m,cc,x,limneu0,ndec,nex,maxd,maxlp,maxn,
     1                    minacc,enr(li,:),sneuf2,sneun2,ineue2,sneudf2,
     2                    sneudr2,dfnorm,idfe,r1dcin,ir1dein,r2nc,ir2ne,
     3                    r2dnc,ir2dne,jneu0,jtest0,jsub0)
              jneu0max=max(jneu0s,jneu0)
              jneu0s=jneu0
              naccrs=naccr
              wronca=r1c(li)*r2dnc*ten**(ir1e(li)+ir2dne)
              wroncb=r2nc*r1dc(li)*ten**(ir2ne+ir1de(li))
              wronc=wronca-wroncb
              naccneu0=-int(log10(abs((wronc-wront)/wront)+dec))
              if(naccneu0.lt.0) naccneu0=0
              if(naccneu0.gt.ndec) naccneu0=ndec
              naccneu0w=naccneu0
              nacccor=-int(log10(abs((wronca-wroncb)/wronca)+dec))
              if(nacccor.lt.0) nacccor=0
              if(nacccor.gt.naccrpl) nacccor=naccrpl
              if(nacccor.gt.0) naccneu0=min(naccneu0+nacccor,ndec-
     1                         jsubf,ndec-jsub0,naccr1)
              if(naccneu0.lt.0) naccneu0=0
              if(naccneu0.lt.naccneu0w) naccneu0=naccneu0w
              if(iopeta.ne.0.and.naccneu0.ge.naccetas) iopeta=0
                if(naccneu0.ge.minacc.and.naccneu0p.ge.minacc.and.
     1             naccneu0p2.ge.minacc) then
                if(l.gt.max(liplp,nbp).and.x.gt.0.1e0_knd) iopint=0
                iopleg=0
                iopleg1=0
                iopeta=0
                end if
              if(naccneu0.lt.minacc.and.iopeta.eq.0.and.x.ge.0.05e0_knd)
     1              iopeta=1
                if(naccneu0.gt.naccr.or.(naccneu0.eq.naccr.and.naccleg
     1             .eq.naccr.and.jflagl.eq.1)) then
                naccr=naccneu0
                jflagl=0
                r2c(li)=r2nc
                ir2e(li)=ir2ne
                r2dc(li)=r2dnc
                ir2de(li)=ir2dne
                iopmatch=0
                end if
                if(naccneu0.gt.minacc+1.and.naccneu0p.gt.minacc+1.and.
     1             li.gt.liplp+10.and.x.ge.0.1e0_knd) then
                 istartint=2
                 iopint=0
                 end if
c              if(knd.eq.kindd) write(40,820) naccneu0,r2nc,ir2ne,r2dnc,
c     1                                       ir2dne
c              if(knd.eq.kindq) write(40,825) naccneu0,r2nc,ir2ne,r2dnc,
c     1                                       ir2dne
              naccneu0p2=naccneu0p
              if(naccneu0.eq.0) iopneu0=0
1080          continue
c
c      r2 calculation using the variable eta expansion
              if(iopneu0.ne.0.and.iopneu0.ne.4.and.iopeta.eq.0.and.x.ge.
     1           0.05e0_knd.and.naccr.lt.minacc.and.istarteta.ne.0)
     2            iopeta=1
              if(iopneu0.eq.4) iopneu0=1
              if(iopeta.eq.0.or.naccr.ge.minacc.or.istarteta.eq.0)
     1             go to 1310
              naccetab=0
              iopnee=0
              if(nee.lt.neelow) nee=neelow
              naccetamax=0
              neemax=nee
              naccnmax=0
              netatry=1
              naccd=0
              jetam=0
                if(ieta.eq.0) then
                  do jnet=1,neta
                  ang=jnet*pi*0.5e0_knd/(neta+1)
                  eta(jnet)=cos(ang)
                  wmeta2(jnet)=2.0e0_knd*(1.0e0_knd+eta(jnet))*
     1                         (sin(0.5e0_knd*ang)**2)
                  xbn(jnet)=sqrt(x*x+wmeta2(jnet))
                  xln(jnet)=eta(jnet)*x/xbn(jnet)
                  end do
                  ieta=1
                end if
              if(iopeta.eq.1) iopeta=2
1090          if(iopeta.eq.4) go to 1320
              if(iopeta.eq.3) go to 1180
              etaval=eta(nee)
1100          xbninp=xbn(nee)
              netainp=1
              etainp(1)=eta(nee)
              xlninp(1)=xln(nee)
                if(knd.eq.kindd) then
                if(x.ge.0.05e0_knd) limn=2*int(25/(x*x)+300/x+
     1                                   3*c+10000)+2
                if(x.ge.0.1e0_knd) limn=2*((lnum+c/5+0.5e0_knd*maxm+
     1                                  200)*1.4e0_knd/x)+2
                if(x.ge.0.5e0_knd) limn=2*((lnum+c/5+0.5e0_knd*maxm+
     1                                  300)/x)+2
                if(x.ge.1.0e0_knd) limn=2*(lnum+c/5+0.5e0_knd*maxm+
     1                                  300)+2
                end if
                if(knd.eq.kindq) then
                if(x.ge.0.05e0_knd) limn=2*int(25/(x*x)+400/x+
     1                                   3*c+14000)+2
                if(x.ge.0.1e0_knd) limn=2*((lnum+c/5+0.5e0_knd*maxm+
     1                                  350)*1.4e0_knd/x)+2
                if(x.ge.0.5e0_knd) limn=2*((lnum+c/5+0.5e0_knd*maxm+
     1                                  400)/x)+2
                if(x.ge.1.0e0_knd) limn=2*(lnum+c/5+0.5e0_knd*maxm+
     1                                  400)+2
                end if
              limn=limn+m
              limp=limn-m
              limpd=2*(lnum+int(c)+ndec)
              if(limpd.gt.limp) limpd=limp
              if(jnen.eq.0) go to 1160
              jnenlim=jnen
              if(jnen.gt.jnenmax) jnenlim=jnenmax
              limplim=limp
              limnlim=limn
                do 1150 jn=1,jnenlim
                if(nee.ne.neeb(jn)) go to 1150
                if(limplim.gt.limpsv(jn)) limplim=limpsv(jn)
                if(limnlim.gt.limnsv(jn)) limnlim=limnsv(jn)
                  do je=1,limpd
                  pratb(je)=pratbsv(jn,je)
                  end do
                  do je=1,limplim
                  pratt(je)=prattsv(jn,je)
                  pdratt(je)=pdrattsv(jn,je)
                  end do
                  do je=1,limnlim
                  sneufe(je)=sneufsv(jn,je)
                  sneudfe(je)=sneudfsv(jn,je)
                  end do
                  jelim=maxlp
                  if(maxlp.gt.limn+1) jelim=limn+1
                  if(jelim.gt.jelimsv(jn)) jelim=jelimsv(jn)
                  do 1130 je=1,jelim
                  sneune(je)=sneunsv(jn,je)
                  sneudre(je)=sneudrsv(jn,je)
                  ineuee(je)=ineuesv(jn,je)
1130              continue
c                write(40,1140) etaval
c1140            format(8x,'r2eta: reused expansion functions for eta ='
c     1                 ,f13.9,'.')
                go to 1180
1150           continue
1160          continue
              jnen=jnen+1
              jnencur=jnen-(jnenmax*int((jnen-1)/jnenmax))
              neeb(jnencur)=nee
              limbesf=4*int(real(cc*xbninp)+abs(aimag(cc*xbninp)))+25
              call sphneu(cc,xbninp,limn,maxn,maxlp,limbesf,sneufe,
     1                    sneune,ineuee,sneudfe,sneudre)
                do je=1,limn
                sneufsv(jnencur,je)=sneufe(je)
                sneudfsv(jnencur,je)=sneudfe(je)
                limnsv(jnencur)=limn
                end do
                jelim=maxlp
                if(maxlp.gt.limn+1) jelim=limn+1
                  do 1170 je=1,jelim
                  sneunsv(jnencur,je)=sneune(je)
                  sneudrsv(jnencur,je)=sneudre(je)
                  ineuesv(jnencur,je)=ineuee(je)
1170              continue
                  jelimsv(jnencur)=jelim
              iopd=3
              call pleg(m,limp,maxp,limcsav,iopd,ndec,xlninp,netainp,
     1                  maxt,prat,pdrat,pdnorma,ipdnorma,pnorma,ipnorma,
     2                  alpha,beta,gamma,coefa,coefb,coefc,coefd,coefe)
              limcsav=max(limcsav,limp)
                do je=1,limp
                pratt(je)=prat(1,je)
                pdratt(je)=pdrat(1,je)
                prattsv(jnencur,je)=pratt(je)
                pdrattsv(jnencur,je)=pdratt(je)
                limpsv(jnencur)=limp
                end do
              iopd=2
              call pleg(m,limpd,maxp,limcsav,iopd,ndec,etainp,netainp,
     1                  maxt,prat,pdrat,pdnorma,ipdnorma,pnorma,ipnorma,
     2                  alpha,beta,gamma,coefa,coefb,coefc,coefd,coefe)
                do je=1,limpd
                pratb(je)=prat(1,je)
                pratbsv(jnencur,je)=pratb(je)
                end do
              pratb(limpd+1)=0.0e0_knd
              pratb(limpd+2)=0.0e0_knd
              pratbsv(jnencur,limpd+1)=0.0e0_knd
              pratbsv(jnencur,limpd+2)=0.0e0_knd
1180          continue
              pcoefe=(x*x+1.0e0_knd)/(xbn(nee)*xbn(nee))
              apcoef=(rm/2.0e0_knd)*log10(pcoefe)
              ipcoefe=int(apcoef)
              pcoefe=ten**(apcoef-ipcoefe)
              pcoefo=pcoefe*pratt(2)/pratb(2)
              ipcoefo=ipcoefe
              pdcoefe=pcoefe
              if(m.ne.0) pdcoefe=-pcoefe*rm*xln(nee)*xbn(nee)*
     1                    xbn(nee)/((x*x+1.0e0_knd)*wmeta2(nee))
              ipdcoefe=ipcoefe
              pdcoefo=pdcoefe*pdratt(2)/pratb(2)
              ipdcoefo=ipdcoefe
              if(li.lt.3) go to 1190
                do jl=3,li+ix,2
                pcoefe=pcoefe*pratt(jl)/pratb(jl)
                iterm=log10(abs(pcoefe))
                pcoefe=pcoefe*ten**(-iterm)
                ipcoefe=ipcoefe+iterm
                pdcoefe=pdcoefe*pdratt(jl)/pratb(jl)
                iterm=log10(abs(pdcoefe))
                pdcoefe=pdcoefe*ten**(-iterm)
                ipdcoefe=ipdcoefe+iterm
                end do
              continue
              if(li.lt.4) go to 1190
                do jl=4,li+1-ix,2
                pcoefo=pcoefo*pratt(jl)/pratb(jl)
                iterm=log10(abs(pcoefo))
                pcoefo=pcoefo*ten**(-iterm)
                ipcoefo=ipcoefo+iterm
                pdcoefo=pdcoefo*pdratt(jl)/pratb(jl)
                iterm=log10(abs(pdcoefo))
                pdcoefo=pdcoefo*ten**(-iterm)
                ipdcoefo=ipdcoefo+iterm
                end do
1190          if(ix.eq.0) go to 1200
              pcoefet=pcoefo
              ipcoefet=ipcoefo
              pcoefo=pcoefo*pratt(li+2)/pratb(li+2)
              iterm=int(log10(abs(pcoefo)))
              pcoefo=pcoefo*ten**(-iterm)
              ipcoefo=ipcoefo+iterm
              pdcoefet=pdcoefo
              ipdcoefet=ipdcoefo
              pdcoefo=pdcoefo*pdratt(li+2)/pratb(li+2)
              iterm=int(log10(abs(pdcoefo)))
              pdcoefo=pdcoefo*ten**(-iterm)
              ipdcoefo=ipdcoefo+iterm
              go to 1210
1200          pcoefet=pcoefe
              ipcoefet=ipcoefe
              pcoefe=pcoefe*pratt(li+2)/pratb(li+2)
              iterm=int(log10(abs(pcoefe)))
              pcoefe=pcoefe*ten**(-iterm)
              ipcoefe=ipcoefe+iterm
              pdcoefet=pdcoefe
              ipdcoefet=ipdcoefe
              pdcoefe=pdcoefe*pdratt(li+2)/pratb(li+2)
              iterm=int(log10(abs(pdcoefe)))
              pdcoefe=pdcoefe*ten**(-iterm)
              ipdcoefe=ipdcoefe+iterm
1210          continue
                if(knd.eq.kindd) then
                if(x.ge.0.05e0_knd) limeta=2*int(25/(x*x)+300/x+3*c+
     1                                     1250*knd)
                if(x.ge.0.1e0_knd) limeta=2*int((l-m+c/5+0.5e0_knd*m+
     1                                    200)*1.4e0_knd/x)
                if(x.ge.0.5e0_knd) limeta=2*int((l-m+c/5+0.5e0_knd*m+
     1                                    300)/x)
                if(x.ge.1.0e0_knd) limeta=2*int(l-m+c/5+0.5e0_knd*m+
     1                                     300)
                end if
                if(knd.eq.kindq) then
                if(x.ge.0.05e0_knd) limeta=2*int(25/(x*x)+400/x+
     1                                     3*c+1250*knd)
                if(x.ge.0.1e0_knd) limeta=2*int((l-m+c/5+0.5e0_knd*m+
     1                                    350)*1.4e0_knd/x)
                if(x.ge.0.5e0_knd) limeta=2*int((l-m+c/5+0.5e0_knd*m+
     1                                    400)/x)
                if(x.ge.1.0e0_knd) limeta=2*int(l-m+c/5+0.5e0_knd*m+
     1                                     400)
                end if
              if(iopeta.eq.3.and.naccrsav.gt.minacc)
     1             limeta=jeta+jeta+40+
     2                    int(sqrt(c)*(2/min(1.0e0_knd,x))+5/x)
              if(iopeta.eq.3.and.naccrsav.le.minacc)
     1                        limeta=jeta+jeta+500+c
              if(iopeta.eq.2) limeta=max(limeta,jeta+jeta+500+int(c))
              if(limeta.gt.limp-2) limeta=limp-2
              if(limeta.gt.limd) limeta=limd
              wm=wmeta2(nee)
              call r2eta(l,m,cc,x,etaval,nee,limeta,ndec,maxd,
     1                   maxlp,maxn,maxp,minacc,wm,enr(li,:),sneufe,
     2                   sneune,ineuee,sneudfe,sneudre,pdratt,pratb,
     3                   pratt,pcoefet,ipcoefet,pdcoefet,ipdcoefet,
     4                   r1cin,ir1ein,r1dcin,ir1dein,naccr1,naccrpl,
     5                   naccnmax,naccr,r2ec,ir2ee,r2dec,ir2dee,nacceta,
     6                   jeta,naccd)
              netatry=netatry+1
              naccetas=nacceta
              naccetasc=min(nacceta,naccr1)
              naccrs=naccr
                if(naccetasc.gt.naccrs.or.(naccetasc.eq.naccrs.and.
     1             naccleg.eq.naccrs.and.jflagl.eq.1))
     2             then
                naccr=naccetasc
                naccetab=naccetasc
                jflagl=0
                r2c(li)=r2ec
                ir2e(li)=ir2ee
                r2dc(li)=r2dec
                ir2de(li)=ir2dee
                iopmatch=0
                end if
                if(naccetas.gt.naccetamax) then
                neemax=nee
                naccetamax=naccetas
                jetam=jeta
                end if
c1270          if(knd.eq.kindd) write(40,1280) naccetasc,etaval,nee,r2ec,
c     1                                    ir2ee,r2dec,ir2dee
c              if(knd.eq.kindq) write(40,1285) naccetasc,etaval,nee,r2ec,
c     1                                    ir2ee,r2dec,ir2dee
c1280          format(15x,'r2eta accuracy = ',i2,' decimal digits; eta',
c     1               ' = ',f12.9,'; nee = ',i4,/,10x,'r2 = ', f19.15,
c     2               f19.15,i5,5x,'r2d = ',f19.15,f19.15,i5)
c1285          format(15x,'r2eta accuracy = ',i2,' decimal digits; eta',
c     1               ' = ',f12.9,'; nee = ',i4,/,10x,'r2 = ', f35.31,
c     2               f35.31,i5,/10x,'r2d = ',f35.31,f35.31,i5)
              iopeta=3
                if(naccetas.lt.naccetamax-2.or.nee.eq.30) then
                nee=neemax-incnee
                iopeta=2
                iplflag=0
                go to 1310
                end if
              jetaflag=0
                if(naccetas.ge.minacc) then
                jetaflag=1
                ietacount=ietacount+1
                if(ietacount.ge.5) incnflag=1
                  if(iplflag.eq.0.and.nee.gt.incnee+incnee) then
                  nee=nee-incnee
                  iopeta=2
                  end if
                iplflag=1
                go to 1320
                end if
              iopeta=2
              if(iplflag.eq.1.and.incnflag.eq.1.and.netatry.eq.2)
     1               iopnee=0
              ietacount=0
              incnflag=0
              if(iopnee.eq.0) go to 1290
              nee=neemax-incnee
              if(nee.eq.0) nee=incnee
              go to 1310
1290          if(nee.eq.neta) go to 1300
              nee=nee+incnee
              go to 1090
1300          continue
              iopeta=3
              if(naccetas.lt.minacc.and.nee.eq.neta)
     1               nee=nee-incnee
              if(naccetas.lt.minacc) iopeta=2
              if(naccetas.lt.minacc) iplflag=0
              if(naccneu0.ne.0.and.naccetas.ge.naccneu0) iopneu0=0
1310          continue
                if(naccr.lt.minacc.and.iopint.eq.0.and.istartint.ne.2
     1              .and.li.le.intlim-5) then
                istartint=1
                iopint=1
                go to 780
                end if
1320          continue
                if(iopint.ne.0.and.naccint.lt.naccr-5.and.naccint.lt.
     1          naccrsav-5.and.x.gt.0.1e0_knd) then
                iopint=0
                if(li.gt.max(nbp,liplp)) istartint=2
                end if
              if(iopeta.eq.4) iopeta=1
              if(iopint.ne.0.and.naccint.gt.nacintsa) nacintsa=naccint
                if(l.eq.m.and.iopint.ne.0.and.naccint.eq.naccr) then
                if(iopeta.ne.0) iflageta=1
                iopeta=0
                end if
                if(neemax.eq.30.and.iopneu0.eq.0) then
                if(istartneu0.eq.1) iopneu0=1
                end if
              if(naccr.eq.ndec) naccr=ndec-1
              if(naccr.gt.0) go to 1330
              naccr=0
              r2c(li)=(1.0e0_knd,0.0e0_knd)
              ir2e(li)=-ndec
              r2dc(li)=(1.0e0_knd,0.0e0_knd)
              ir2de(li)=-ndec
1330          continue
                if(ioprad.eq.1) then
c                write(20,1340) l,r1c(li),ir1e(li),
c     1               r1dc(li),ir1de(li),naccr1
c1340            format(1x,i5,2x,2(f17.14,f17.14,i6,2x),i2)
                go to 1400
                end if
                if(x.eq.0.0e0_knd) then
c                write(20,1350) l,r1c(li),ir1e(li),r1dc(li),ir1de(li),
c     1                   r2c(li),ir2e(li),r2dc(li),ir2de(li),naccr
c1350            format(1x,i5,2x,2(f17.14,f17.14,i6,2x),/,8x,
c     1                   2(f17.14,f17.14,i6,2x),i2,'e')
                go to 1400
                end if
                if(ix.eq.0.and.li.lt.lipl.and.match(2).ne.0) then
                naccrps=naccr
                naccflag=1
                iopmatchp=iopmatch
                jjflagl=0
                if(jflagl.eq.1.and.naccr.eq.naccleg.and.naccint.ne.
     1                 naccr) jjflagl=1
                  if(iopmatch.eq.1) then
c                  write(40,1360) li+m-1,li+m
                  end if
                naccrsav=naccr
                go to 1400
                end if
                if(ix.eq.1.and.naccflag.eq.1) then
                naccflag=0
                wronca=r1c(li-1)*r1dc(li)*ten**(ir1e(li-1)+
     1                 ir1de(li))
                wroncb=r1c(li)*r1dc(li-1)*
     1                 ten**(ir1e(li)+ir1de(li-1))
                wronc=wronca-wroncb
                naccrp=-int(log10(abs((wronc-wront)/wront)+dec))
                if(naccrp.lt.0) naccrp=0
                if(naccrp.gt.ndec) naccrp=ndec
                nacccor=-int(log10(abs((wronca-wroncb)/wronca)+dec))
                if(nacccor.lt.0) nacccor=0
                if(nacccor.gt.naccrpl) nacccor=naccrpl
                naccrp=naccrp+nacccor
                naccrp=min(naccrp,naccr1p,nmatch)
                  if(iopmatchp.eq.1.or.(iopmatchp.eq.0.and.
     1                 naccrp.gt.naccrps)) then
                  r2c(li-1)=r1c(li)
                  ir2e(li-1)=ir1e(li)
                  r2dc(li-1)=r1dc(li)
                  ir2de(li-1)=ir1de(li)
                  jjflagl=0
c                  write(40,1360) li+m-2,li+m-1
c                  write(40,1395) naccrp,l-1
                  end if
                  if(iopmatchp.eq.0.and.naccrps.ge.naccrp)
     1                 naccrp=naccrps
                  if(jjflagl.eq.1.or.naccrp.eq.min(naccrplp,naccr1p))
     1                then
                  continue
c                  write(20,1350) l-1,r1c(li-1),ir1e(li-1),
c     1                r1dc(li-1),ir1de(li-1),r2c(li-1),ir2e(li-1),
c     2                r2dc(li-1),ir2de(li-1),naccrp
                  else
                  continue
c                  write(20,1380) l-1,r1c(li-1),ir1e(li-1),
c     1                r1dc(li-1),ir1de(li-1),r2c(li-1),ir2e(li-1),
c     2                r2dc(li-1),ir2de(li-1),naccrp
                  end if
                  if(iopmatch.eq.1) then
                  r2c(li)=-r1c(li-1)
                  ir2e(li)=ir1e(li-1)
                  r2dc(li)=-r1dc(li-1)
                  ir2de(li)=ir1de(li-1)
c                  write(40,1370) li+m-1,li+m-2
c                  write(40,1395) naccr,l
                  end if
                  jjflagl=0
                  if(iopmatch.eq.0.and.naccr.eq.naccleg.and.
     1              naccint.ne.naccr.and.jflagl.eq.1) jjflagl=1
                  if(jjflagl.eq.1.or.naccr.eq.min(naccrpl,naccr1))
     1                then
                  continue
c                  write(20,1350) l,r1c(li),ir1e(li),r1dc(li),ir1de(li),
c     1                         r2c(li),ir2e(li),r2dc(li),ir2de(li),
c     2                         naccr
                  else
                  continue
c                  write(20,1380) l,r1c(li),ir1e(li),r1dc(li),ir1de(li),
c     1                         r2c(li),ir2e(li),r2dc(li),ir2de(li),
c     2                         naccr
                  end if
                naccrsav=naccr
                go to 1400
                end if
c1360          format(8x,'Values for r2 and r2d for ','l = ',i5,
c     1               ' are given by r1 and r1d for l = ',i5)
c1370          format(8x,'Values for r2 and r2d for ','l = ',i5,
c     1               ' are given by -r1 and -r1d for l = ',i5)
                if((jflagl.eq.1.and.naccr.eq.naccleg.and.
     1              naccr.ne.naccint).or.naccr.eq.min(naccrpl,naccr1))
     2               then
                continue
c                write(20,1350) l,r1c(li),ir1e(li),r1dc(li),ir1de(li),
c     1                   r2c(li),ir2e(li),r2dc(li),ir2de(li),naccr
                else
                continue
c                write(20,1380) l,r1c(li),ir1e(li),r1dc(li),ir1de(li),
c     1                   r2c(li),ir2e(li),r2dc(li),ir2de(li),naccr
                end if
              go to 1400
c1380          format(1x,i5,2x,2(f17.14,f17.14,i6,2x),/,8x,
c     1                   2(f17.14,f17.14,i6,2x),i2,'w')
1390          continue
                if(ix.eq.1) then
                  if(min(naccrplp,naccr1p).ge.nmatch-iii) then
                  naccrp=min(naccrplp,naccr1p)
c                  write(40,775) l-1,naccrp
c                  write(20,1350) l-1,r1c(li-1),ir1e(li-1),r1dc(li-1),
c     1                     ir1de(li-1),r2c(li-1),ir2e(li-1),r2dc(li-1),
c     2                     ir2de(li-1),naccrp
                  else
                  r2c(li-1)=r1c(li)
                  ir2e(li-1)=ir1e(li)
                  r2dc(li-1)=r1dc(li)
                  ir2de(li-1)=ir1de(li)
                  wronca=r1c(li-1)*r2dc(li-1)*ten**(ir1e(li-1)+
     1                   ir2de(li-1))
                  wroncb=r2c(li-1)*r1dc(li-1)*
     1                   ten**(ir2e(li-1)+ir1de(li-1))
                  wronc=wronca-wroncb
                  naccrp=-int(log10(abs((wronc-wront)/wront)+dec))
                  if(naccrp.lt.0) naccrp=0
                  nacccor=-int(log10(abs((wronca-wroncb)/wronca)+dec))
                  if(nacccor.lt.0) nacccor=0
                  if(nacccor.gt.naccrpl) nacccor=naccrpl
                  naccrp=naccrp+nacccor
                  if(naccrp.gt.ndec-1) naccrp=ndec-1
                  naccrp=min(naccrp,naccr1p,nmatch)
c                  write(40,1395) naccrp,l-1
c                  write(20,1380) l-1,r1c(li-1),ir1e(li-1),r1dc(li-1),
c     1                     ir1de(li-1),r2c(li-1),ir2e(li-1),r2dc(li-1),
c     2                     ir2de(li-1),naccrp
                  end if
                  if(min(naccrpl,naccr1).ge.nmatch-iii) then
                  naccr=min(naccrpl,naccr1)
c                  write(40,775) l,naccr
c                  write(20,1350) l,r1c(li),ir1e(li),r1dc(li),ir1de(li),
c     1                           r2c(li),ir2e(li),r2dc(li),ir2de(li),
c     2                           naccr
                  else
                  r2c(li)=-r1c(li-1)
                  ir2e(li)=ir1e(li-1)
                  r2dc(li)=-r1dc(li-1)
                  ir2de(li)=ir1de(li-1)
                  wronca=r1c(li)*r2dc(li)*ten**(ir1e(li)+
     1                   ir2de(li))
                  wroncb=r2c(li)*r1dc(li)*
     1                   ten**(ir2e(li)+ir1de(li))
                  wronc=wronca-wroncb
                  naccr=-int(log10(abs((wronc-wront)/wront)+dec))
                  if(naccr.lt.0) naccr=0
                  nacccor=-int(log10(abs((wronca-wroncb)/wronca)+dec))
                  if(nacccor.lt.0) nacccor=0
                  if(nacccor.gt.naccrpl) nacccor=naccrpl
                  naccr=naccr+nacccor
                  if(naccr.gt.ndec-1) naccr=ndec-1
                  naccr=min(naccr,naccr1,nmatch)
c                  write(40,1395) naccr,l
c                  write(20,1380) l,r1c(li),ir1e(li),r1dc(li),ir1de(li),
c     1                           r2c(li),ir2e(li),r2dc(li),ir2de(li),
c     2                           naccr
                  end if
c1395            format(10x,'Accuracy using eigenvalue match is ',i3,
c     1                 ' digits for l = ',i5)
                end if
1400          continue
              irtest=max(ir1e(li),ir1de(li))
                if(naccr.eq.naccleg.and.naccrsav.eq.nacclegp.and.li.gt.
     1             max(liplp,nbp).and.irtest.lt.-10) then
                iopneu0=0
                end if
              ir1ep=ir1e(li)
              naccrsavp=naccrsav
              naccrsav=naccr
              naccr1p=naccr1
              naccrtp=naccrt
              naccetabp=naccetab
1405          continue
                if(ioprad.eq.2.and.li.le.lipl.and.
     1               match(2).ne.0) then
                continue
                if(ix.eq.1.and.naccrp.lt.6) write(60,*)
     1             ' est. acc. = ',naccrp, ' digits for m = ',m,
     2             ' l = ',l-1,' x = ',x,' c = ',cc
                if(ix.eq.1.and.naccr.lt.6) write(60,*)
     1            ' est. acc. = ',naccr, ' digits for m = ',
     2             m,' l = ',l,' x = ',x,' c = ',cc
                end if                                 
                if(ioprad.eq.2.and.li.gt.lipl.and.naccr.lt.6) then
                continue
                write(60,*) ' est. acc. = ',naccr,' digits for m = ',m,
     1           ' l = ', l,' x = ',x,' c = ',cc
                end if
              if(ioprad.eq.1.and.naccr1.lt.6) write(60,*)
     1           'est. r1 acc. = ',naccr1,' digits for m = ',m,' l = ',
     2           l,' x = ',x,' c = ',cc
              naccetamax=0
              naccrp=naccr
              if(ioprad.eq.2) naccrplp=naccrpl
              if(ioprad.eq.2) naccneu0p=naccneu0
              naccrep=naccre
1410            if(ndec-jsubms-1.lt.6) then
                write(60,*) ' est. MS norm acc. = ',ndec-jsubms-1,
     1            ' digits for m = ',m,' l = ', l,' c = ',cc
                end if
              if(iopang.eq.0) go to 1510
c
c  determine first kind oblate angular functions
              if(l.eq.m) lims1=3*ndec+int(c)
              if(l.ne.m) lims1=jang+jang+20+c/25
              if(lims1.gt.maxp) lims1=maxp
              call s1leg(l,m,cc,iopang,iopnorm,barg,narg,lims1,ndec,nex,
     1                   maxt,maxd,maxp,enr(li,:),dmlms,idmlmse,jsubms,
     2                   pr,pdr,pdnorm,ipdnorm,pnorm,ipnorm,pdtempe,
     3                   ipdtempe,pdtempo,ipdtempo,ptempe,iptempe,
     4                   ptempo,iptempo,itestm,naccre,kindd,kindq,s1c,
     5                   is1e,s1dc,is1de,naccs,naccds,jang,dmlms1,
     6					 idmlms1e)
                do 1500 jarg=1,narg
                s1(li,jarg)=s1c(jarg)
                s1d(li,jarg)=s1dc(jarg)
                is1(li,jarg)=is1e(jarg)
                is1d(li,jarg)=is1de(jarg)
c                if(iopang.eq.1) write(50,1420)
c     1                barg(jarg),naccs(jarg)
c                if(iopang.eq.2) write(50,1430)
c     1                barg(jarg),naccs(jarg),naccds(jarg)
c1420            format(1x,'eta = ',e24.15,'   accuracy = ',i2,
c     1                 ' digits.')
c1430            format(1x,'eta = ',e24.15,'   s1 and s1d accuracy = ',
c     1                 i2,' and ',i2,' digits.')
c                if(iopang.eq.1) write(30,1460)
c     1                barg(jarg),s1c(jarg),is1e(jarg),naccs(jarg)
c                if(iopang.eq.2) write(30,1470)
c     1                barg(jarg),s1c(jarg),is1e(jarg),s1dc(jarg),
c     2                is1de(jarg),naccs(jarg),naccds(jarg)
c                if(knd.eq.kindd.and.iopang.eq.1) write(50,1480)
c     1                s1c(jarg),is1e(jarg)
c                if(knd.eq.kindd.and.iopang.eq.2) write(50,1490)
c     1                s1c(jarg),is1e(jarg),s1dc(jarg),is1de(jarg)
c                if(knd.eq.kindq.and.iopang.eq.1) write(50,1485)
c     1                s1c(jarg),is1e(jarg)
c                if(knd.eq.kindq.and.iopang.eq.2) write(50,1495)
c     1                s1c(jarg),is1e(jarg),s1dc(jarg),is1de(jarg)
c1460            format(1x,f17.14,2x,f17.14,1x,f17.14,2x,i5,2x,', ',i2)
c1470            format(1x,f17.14,2x,f17.14,1x,f17.14,2x,i5,2x,f17.14,
c     1                 1x,f17.14,2x,i5,2x,i2,', ',i2)
c1480            format(12x,'s1 = ',f17.14,f17.14,2x,i5)
c1485            format(12x,'s1 = ',f35.31,f35.31,2x,i5)
c1490            format(12x,'s1 = ',f17.14,1x,f17.14,2x,i5,5x,'s1d = ',
c     1                 f17.14,1x,f17.14,2x,i5)
c1495            format(12x,'s1 = ',f35.31,1x,f35.31,2x,i5,/12x,'s1d = ',
c     1                 f35.31,1x,f35.31,2x,i5)
1500            continue
1510          continue
1515      continue
          ijmax=min(lnum,limeig)
            do i=1,ijmax-2
              do j=i+2,ijmax,2
              iegch=-int(log10(abs((eig(i)-eig(j))/eig(i))+dec))
                if(iegch.gt.min(ieigt(i),ieigt(j))) then
                write(60,*) m,cc,i+m-1,j+m-1,' duplicated eigenvalue'
                end if
              end do
            end do
1540      continue
        continue
        return
        end
c
c
        subroutine s1leg (l,m,cc,iopang,iopnorm,barg,narg,lims1,ndec,
     1                    nex,maxt,maxd,maxp,enr,dmlms,idmlmse,
     2                    jsubms,pr,pdr,pdnorm,ipdnorm,pnorm,ipnorm,
     3                    pdtempe,ipdtempe,pdtempo,ipdtempo,ptempe,
     4                    iptempe,ptempo,iptempo,itestm,naccre,kindd,
     5                    kindq,s1c,is1e,s1dc,is1de,naccs,naccds,jang,
     6                    dmlms1,idmlms1e)
c
c  purpose:     To calculate the oblate angular functions of the first
c               kind and their first derivatives with respect to eta.
c
c  parameters:
c
c     input :   l       : l
c               m       : m
c               cc      : complex c
c               iopang  : index = 1 when angular functions of the
c                         first kind are calculated; = 2 when the
c                         first derivatives with respect to eta are
c                         also calculated
c               iopnorm : = 1 when the angular functions (and
c                         first derivatives) are scaled by the
c                         square root of the normalization of the
c                         corresponding Legendre function, giving them
c                         unity norm; iopnorm = 0 otherwise
c               barg    : array of eta values for which angular
c                         functions are desired
c               narg    : number of eta values
c               lims1   : approximately twice the maximum number
c                         of terms available to be taken in the
c                         sums
c               ndec    : number of decimal digits available in
c                         real arithmetic
c               nex     : maximum exponent in real arithmetic
c               maxt    : dimension of barg, pdnorm, ipdnorm, pnorm,
c                         ipnorm, pdtempe, ipdtempe, pdtempo, ipdtempo,
c                         ptempe, iptempe, ptempo, iptempo, s1c, is1e,
c                         s1dc, is1de, and naccs arrays. first dimension
c                         of the doubly dimensioned arrays pr and pdr
c               maxd    : dimension of enr array
c               maxp    : second dimension of pr and pdr arrays
c               enr     : array of d coefficient ratios
c               dmlmse  : characteristic of the d coefficient with
c                         index l - m when using Meixner-Schafke
c                         normalization for the angular functions
c               idmlmse : exponent associated with dmsnorm
c               jsubms  : effective subtraction error in calculating
c                         the Meixner-Schafke normalization
c               pr      : array of ratios of successive first kind
c                         associated Legendre functions of the same
c                         parity
c               pdr     : array of ratios of successive derivatives of
c                         first kind associated Legendre functions of
c                         the same parity
c               pdnorm  : array of characteristics of the first
c                         derivatives of associated Legendre functions
c                         of the first kind of order m and degree m
c               ipdnorm : array of exponents corresponding to pdnorm
c               pnorm   : array of characteristics of the associated
c                         Legendre functions of the first kind of order
c                         m and degree m
c               ipnorm  : array of exponents corresponding to pnorm
c               pdtempe : storage array of characteristics of the ratio
c                         of the first derivative of the associated
c                         Legendre function of order m and degree l - 2
c                         or l - 1, depending on whether l - m is even
c                         or odd, to the first derivative of the
c                         function of order m and degree m
c               ipdtempe: array of exponents corresponding to pdtempe
c               pdtempo : storage array of characteristics of the ratio
c                         of the first derivative of the associated
c                         Legendre function of order m and degree l - 2
c                         or l - 1, depending on whether l - m is odd
c                         or even, to the first derivtive of the
c                         function of order m and degree m
c               ipdtempo: array of exponents corresponding to pdtempo
c               ptempe  : storage array of characteristics of the ratio
c                         of the associated Legendre function of order
c                         m and degree l - 2 or l - 1, depending on
c                         whether l - m is even or odd, to the function
c                         of order m and degree m
c               iptempe : array of exponents corresponding to ptempe
c               ptempo  : storage array of characteristics of the ratio
c                         of the associated Legendre function of order
c                         m and degree l - 2 or l - 1, depending on
c                         whether l - m is odd or even, to the
c                         function of order m and degree m
c               iptempo : array of exponents corresponding to ptempo
c               itestm  : number of leading decimal digits of agreement
c                         between the forward and the reverse recursion
c                         involved in calculation of the d coefficients
c               naccre  : estimated accuracy of the eigenvalue
c               kindd   : number of bytes for real data in double
c                         precision
c               kindq   : number of bytes for real data in quadruple
c                         precision
c
c
c     output:   s1c    : array of characteristics of oblate
c                        angular functions of the first kind
c               is1e   : array of exponents of oblate angular
c                        functions of the first kind
c               s1dc   : array of characteristics of derivative with
c                        respect to eta of oblate angular functions
c                        of the first kind
c               is1de  : array of exponents of derivatives with respect
c                        to eta of oblate angular functions of first
c                        kind
c               naccs  : array of integer estimates of the number of
c                        accurate decimal digits in the values obtained
c                        for s1
c               naccds : array of integer estimates of the number of
c                        accurate decimal digits in the values obtained
c                        for s1d
c               jang   : maximum value of the index j in the forward
c                        sum for r1 and r1d, i.e., the highest enr(j)
c                        used
c               dmlms1  : characteristic of the d coefficient with index
c                         l - m when the angular functions have unity
c                         norm
c               idmlms1e: exponent associated with dmlms1
c
c        use param
c
c  real(knd) scalars and arrays
        real(knd) adec,aj,aj2,dcon,dec,factor,fterm,rm2,rm2m1,rm2m3,
     1            rm2p1,ten,teste,testeo
        real(knd) barg(maxt),pdr(maxt,maxp),pdnorm(maxt),
     1            pnorm(maxt),pr(maxt,maxp),pdtemp(maxt),ptemp(maxt),
     2            pdtempe(maxt),ptempe(maxt),pdtempo(maxt),ptempo(maxt)
c
c  complex(knd) scalars and arrays
        complex(knd) cc,coef,dnew,dnewd,dmlms,dmlms1,dold,doldd,s1,s1d
        complex(knd) enr(maxd),s1c(maxt),s1dc(maxt)
c
c  integer arrays
        integer ipdnorm(maxt),ipnorm(maxt),ipdtemp(maxt),iptemp(maxt),
     1          ipdtempe(maxt),iptempe(maxt),ipdtempo(maxt),
     2          iptempo(maxt),is1de(maxt),is1e(maxt),naccs(maxt),
     3          naccds(maxt)
c
        ten=10.0e0_knd
        dec=ten**(-ndec-1)
        dcon=dec
        adec=1000.0e0_knd*dec
        nfac=nex/3
        if(2*(nfac/2).ne.nfac) nfac=nfac-1
        teste=ten**nfac
        testeo=1.0e0_knd/teste
        kflag=0
        rm2=m+m
        rm2m1=m+m-1
        rm2p1=m+m+1
        rm2m3=m+m-3
        if(l.gt.(m+1)) go to 30
          do 20 k=1,narg
          if(pnorm(k).eq.0.0e0_knd) go to 20
          if(l.eq.(m+1)) go to 10
          ptempe(k)=pr(k,1)
          iptempe(k)=0
          pdtempe(k)=pdr(k,1)
          ipdtempe(k)=0
          ptemp(k)=ptempe(k)
          pdtemp(k)=pdtempe(k)
          iptemp(k)=0
          ipdtemp(k)=0
          go to 20
10        ptempo(k)=pr(k,2)
          iptempo(k)=0
          pdtempo(k)=pdr(k,2)
          ipdtempo(k)=0
          ptemp(k)=ptempo(k)
          pdtemp(k)=pdtempo(k)
          iptemp(k)=0
          ipdtemp(k)=0
20        continue
30      continue
        lm2=(l-m)/2
        ix=l-m-2*lm2
        ixx=ix-1
        ixx2=ixx+2
        if(l.lt.(m+2)) go to 110
          do 100 k=1,narg
          if(pnorm(k).eq.0.0e0_knd) go to 100
          if(ix.ne.0) go to 60
          ptempe(k)=ptempe(k)*pr(k,l-m+1)
            if(abs(ptempe(k)).gt.1.0e+10_knd) then
            ptempe(k)=ptempe(k)*(1.0e-10_knd)
            iptempe(k)=iptempe(k)+10
            end if
          ptemp(k)=ptempe(k)
          iptemp(k)=iptempe(k)
          if(abs(barg(k)).lt.adec) go to 100
          pdtempe(k)=pdtempe(k)*pdr(k,l-m+1)
            if(abs(pdtempe(k)).gt.1.0e+10_knd) then
            pdtempe(k)=pdtempe(k)*(1.0e-10_knd)
            ipdtempe(k)=ipdtempe(k)+10
            end if
          pdtemp(k)=pdtempe(k)
          ipdtemp(k)=ipdtempe(k)
          go to 100
60        if(abs(barg(k)).lt.adec) go to 80
          ptempo(k)=ptempo(k)*pr(k,l-m+1)
            if(abs(ptempo(k)).gt.1.0e+10_knd) then
            ptempo(k)=ptempo(k)*(1.0e-10_knd)
            iptempo(k)=iptempo(k)+10
            end if
          ptemp(k)=ptempo(k)
          iptemp(k)=iptempo(k)
80        pdtempo(k)=pdtempo(k)*pdr(k,l-m+1)
            if(abs(pdtempo(k)).gt.1.0e+10_knd) then
            pdtempo(k)=pdtempo(k)*(1.0e-10_knd)
            ipdtempo(k)=ipdtempo(k)+10
            end if
          pdtemp(k)=pdtempo(k)
          ipdtemp(k)=ipdtempo(k)
100        continue
110     continue
        lim=lims1/2-ix
        jlow=lm2+1
        jang=0
c
c  compute the associated Legendre function normalization factor
        factor=1.0e0_knd
        ifactor=0
        if(iopnorm.eq.0) go to 210
        if(m.eq.0) go to 170
          do 160 j=1,m
          aj=j
          factor=factor*(aj+aj)*(aj+aj-1.0e0_knd)
            if(factor.gt.teste) then
            factor=factor*testeo
            ifactor=ifactor+nfac
            end if
160       continue
170     if(l.eq.m) go to 190
          do 180 j=1,l-m
          aj=j
          factor=factor*(rm2+aj)/(aj)
            if(factor.gt.teste) then
            factor=factor*testeo
            ifactor=ifactor+nfac
            end if
180       continue
190     continue
        factor=factor*2.0e0_knd/(l+l+1.0e0_knd)
          if(2*(ifactor/2).ne.ifactor) then
          factor=factor*ten
          ifactor=ifactor-1
          end if
        factor=sqrt(factor)
        ifactor=ifactor/2
        iterm=int(log10(factor))
        factor=factor*(ten**(-iterm))
        ifactor=ifactor+iterm
        dmlms1=dmlms/factor
        idmlms1e=idmlmse-ifactor
c        if(knd.eq.kindd) write(50,200) factor,ifactor
c        if(knd.eq.kindq) write(50,205) factor,ifactor
c200     format(1x,'square root of Legendre norm = ',f19.15,2x,i5)
c205     format(1x,'square root of Legendre norm = ',f35.31,2x,i5)
210     continue
c
c  compute the angular function s1
          do 380 k=1,narg
          if(pnorm(k).eq.0.0e0_knd) go to 220
          if((ix.eq.1).and.(abs(barg(k)).lt.adec)) go to 220
          if(((abs(abs(barg(k))-1.0e0_knd)).lt.adec).
     1         and.(m.ne.0)) go to 220
          go to 230
220       s1c(k)=(0.0e0_knd,0.0e0_knd)
          is1e(k)=0
          naccs(k)=ndec
          go to 300
230       dold=(1.0e0_knd,0.0e0_knd)
          s1=dold
          is1=0
          fterm=s1
          lflag=0
            do 240 j=jlow,lim
            dnew=dold*enr(j)*pr(k,j+j+ixx2)
            s1=s1+dnew
            if(abs(s1).gt.fterm) fterm=abs(s1)
            if(abs(dnew/s1).lt.dcon) go to 250
              if(abs(s1).gt.teste) then
              s1=s1*testeo
              dnew=dnew*testeo
              fterm=fterm*testeo
              is1=is1+nfac
              kflag=1
              end if
            dold=dnew
240         continue
250       if(j.gt.jang) jang=j
c          write(50,260) barg(k),j
c260       format(8x,'s1 calculation for eta = ',f13.8,' converged in ',
c     1           i6,' terms.')
          if(lm2.lt.1.or.kflag.eq.1) go to 280
          dold=(1.0e0_knd,0.0e0_knd)
          j=lm2
            do 270 jj=1,lm2
            dnew=dold/(pr(k,j+j+ixx2)*enr(j))
            s1=s1+dnew
            if(abs(s1).gt.fterm) fterm=abs(s1)
            if(abs(dnew/s1).lt.dcon) go to 280
            dold=dnew
            j=j-1
270         continue
280       s1c(k)=s1*dmlms*ptemp(k)*pnorm(k)/factor
          if(abs(s1c(k)).ne.0.0e0_knd) iterm=int(log10(abs(s1c(k))))
          if(abs(s1c(k)).eq.0.0e0_knd) iterm=0
          s1c(k)=s1c(k)*(ten**(-iterm))
          is1e(k)=is1+iptemp(k)+ipnorm(k)+iterm-ifactor+idmlmse
          if(abs(s1c(k)).ge.1.0e0_knd) go to 290
          s1c(k)=s1c(k)*ten
          is1e(k)=is1e(k)-1
290       continue
          if(abs(s1).eq.0.0e0_knd) naccs(k)=0
            if(abs(s1).ne.0.0e0_knd) then
            iacc=int(log10(abs((fterm)/(s1))))
            if(iacc.lt.0) iacc=0
            if(iacc.gt.ndec) iacc=ndec
            naccs(k)=min(ndec-2-iacc,naccre-1,itestm-1,
     1               ndec-1-jsubms)
            end if
          if(naccs(k).gt.0) go to 300
          naccs(k)=0
          naccds(k)=0
          s1c(k)=(0.0e0_knd,0.0e0_knd)
          is1e(k)=0
          s1dc(k)=(0.0e0_knd,0.0e0_knd)
          is1de(k)=0
          go to 380
c
c       compute the first derivative of the angular function when
c       iopang equals 2
300       if(iopang.ne.2) go to 380
          if(pnorm(k).eq.0.0e0_knd) go to 310
          if((ix.eq.0).and.(abs(barg(k)).lt.adec)) go to 310
          if(((abs(abs(barg(k))-1.0e0_knd)).lt.adec).and.(m.ne.0)
     1        .and.(m.ne.2)) go to 310
          go to 320
310       s1dc(k)=(0.0e0_knd,0.0e0_knd)
          is1de(k)=0
          naccds(k)=ndec
          go to 370
320       doldd=(1.0e0_knd,0.0e0_knd)
          s1d=doldd
          is1d=0
          if(l.eq.0) s1d=(0.0e0_knd,0.0e0_knd)
          fterm=s1d
            do 330 j=jlow,lim
            dnewd=doldd*enr(j)*pdr(k,j+j+ixx2)
            s1d=s1d+dnewd
            if(abs(s1d).gt.fterm) fterm=abs(s1d)
            if(abs(dnewd/s1d).lt.dcon) go to 340
              if(abs(s1d).gt.teste) then
              s1d=s1d*testeo
              dnewd=dnewd*testeo
              fterm=fterm*testeo
              is1d=is1d+nfac
              end if
            doldd=dnewd
330         continue
340       if(lm2.lt.1.or.kflag.eq.1) go to 360
          doldd=(1.0e0_knd,0.0e0_knd)
          j=lm2
          ja=lm2
          if(m.eq.0.and.ix.eq.0) ja=lm2-1
          if(ja.eq.0) go to 360
            do 350 jj=1,ja
            dnewd=doldd/(pdr(k,j+j+ixx2)*enr(j))
            s1d=s1d+dnewd
            if(abs(s1d).gt.fterm) fterm=abs(s1d)
            if(abs(dnewd/s1d).lt.dcon) go to 360
            doldd=dnewd
            j=j-1
350         continue
360       s1dc(k)=s1d*dmlms*pdtemp(k)*pdnorm(k)/factor
          if(abs(s1dc(k)).ne.0.0e0_knd) iterm=int(log10(abs(s1dc(k))))
          if(abs(s1dc(k)).eq.0.0e0_knd) iterm=0
          s1dc(k)=s1dc(k)*ten**(-iterm)
          is1de(k)=is1d+ipdtemp(k)+ipdnorm(k)+iterm-ifactor+idmlmse
          naccds(k)=0
            if(abs(s1d).ne.0.0e0_knd) then
            iacc=int(log10(abs((fterm)/(s1d))))
            if(iacc.lt.0) iacc=0
            if(iacc.gt.ndec) iacc=ndec
            naccds(k)=min(ndec-2-iacc,naccre-1,itestm-1,
     1                ndec-1-jsubms)
            end if
          if(naccds(k).lt.0) naccds(k)=0
          if(abs(s1dc(k)).ge.1.0e0_knd) go to 370
          s1dc(k)=s1dc(k)*ten
          is1de(k)=is1de(k)-1
370       continue
          if(naccds(k).eq.0) s1dc(k)=(0.0e0_knd,0.0e0_knd)
          if(naccds(k).eq.0) is1de(k)=0
380       continue
        return
        end
c
c
        subroutine r1bes(l,m,cc,x,limr1,ndec,maxd,enr,maxj,maxn,maxlp,
     1                   nex,iflag,sbesf,sbesdf,sbesn,ibese,sbesdr,
     2                   prat1,pcoefn,ipcoefn,dmfnorm,idmfe,ir1ep,r1c,
     3                   ir1e,r1dc,ir1de,jbes,nsub,ndsub)
c
c  purpose:     To calculate the oblate radial function of the
c               first kind and its first derivative with respect
c               to x, using the traditional expansion of spherical
c               Bessel functions of the first kind with argument
c               c*x, i.e., with eta = 1.
c
c  parameters:
c
c     input:    l      : l
c               m      : m
c               cc     : omplex c
c               x      : x
c               limr1  : approximately twice the maximum number of
c                        terms available to be taken in the series
c               ndec   : number of decimal digits available in
c                        real arithmetic
c               maxd   : dimension of enr array
c               enr    : d coefficient ratios
c               maxj   : dimension of sbesf and sbesdf arrays
c               maxn   : dimension of prat1 array
c               maxlp  : dimension of the sbesdr, sbesn, and ibese
c                        arrays
c               nex    : maximum exponent available in real(knd)
c                        arithmetic
c               iflag  : integer = 1 if forward series not needed;
c                        =0 if the forward series is computed
c               sbesf  : array of ratios of successive first kind
c                        spherical Bessel functions of the same parity
c               sbesdf : array of ratios of successive derivatives of
c                        first kind spherical Bessel functions of the
c                        same parity
c               sbesn  : array of characteristics for Bessel functions
c               ibese  : array of exponents corresponding to sbesn
c               sbesdr : value of ratio of first derivative of
c                        spherical Bessel function to the corresponding
c                        Bessel function
c               prat1  : array of ratios of successive coefficients in
c                        r1 and r1d sum
c               pcoefn : characteristic of coefficient for term in both
c                        r1 and r1d sums that contains Bessel function
c                        of order l
c               ipcoefn: exponent (to the base 10) corresponding to
c                        pcoefn
c               dmfnorm: characteristic of Morse-Feshbach normalization
c                        sum of the d coefficients. equal to the
c                        reciprocal of the value of the d coefficient
c                        d(n = l - m) using this normalization for the
c                        angular functions
c               idmfe  : exponent associated with dmfnorm
c               ir1ep  : exponent for the value of r1 for l-1
c
c     output  : r1c    : characteristic of oblate radial function
c                        of the first kind
c               ir1e   : exponent of oblate radial function of the
c                        first kind
c               r1dc   : characteristic of derivative with respect
c                        to x of oblate radial function of the first
c                        kind
c               ir1de  : exponent of derivative with respect to x of
c                        oblate radial function of the first kind
c               jbes   : maximum value of the index j in the forward
c                        sum for r1 and r1d, i.e., the highest enr(j)
c                        used
c               nsub   : subtraction error in calculating r1
c               ndsub  : subtraction error in calculating r1d
c
c        use param
c
c  real(knd) scalars and arrays
        real(knd) c,dec,em,pcoefn,r1dcoef,sposr,sposi,sposar,sposai,
     1            sdposr,sdposi,ten,teste,testeo,x,x2
        real(knd) prat1(maxn)
c  complex(knd) scalars and arrays
        complex(knd) cc,coef,dmfnorm,dnew,dnewd,dold,doldd,r1c,r1ca,
     1               r1d,r1dc,r1temp,r1tempa,r1dtemp,r1dtempa,
     2               r1top,r1topa,r1dtop,term,termd
        complex(knd) enr(maxd),sbesdf(maxj),sbesdr(maxlp),sbesf(maxj),
     1               sbesn(maxlp)
c
c  integer array
        integer ibese(maxlp)
c
c  convergence ratio dec is set according to the requested accuracy
        ten=10.0e0_knd
        dec=ten**(-ndec-1)
        lm2=(l-m)/2
c
c  ix=0 for l-m even, ix=1 for l-m odd
        ix=l-m-2*lm2
        lim=limr1/2-ix
        nfac=nex/2
        if(2*(nfac/2).ne.nfac) nfac=nfac-1
        teste=ten**nfac
        testeo=1.0e0_knd/teste
        ir1tope=0
        mml=m+m-1+ix
        em=m
        jflag=0
        kflag=0
        nsuba=0
        x2=x*x
        term=(0.0e0_knd,0.0e0_knd)
        r1dcoef=-em/(x*(x*x+1.0e0_knd))
        coef=-(sbesn(m+2)/(sbesn(m+1)*sbesdr(m+1)))*
     1       ten**(ibese(m+2)-ibese(m+1))
        r1topa=(0.0e0_knd,0.0e0_knd)
        sposar=0.0e0_knd
        sposai=0.0e0_knd
c
c  forward summation of numerator series for both r1 and r1d
        dold=(1.0e0_knd,0.0e0_knd)
        doldd=dold
        sposr=1.0e0_knd
        sdposr=1.0e0_knd
        sposi=0.0e0_knd
        sdposi=0.0e0_knd
        r1top=dold
        r1dtop=doldd
          if(lm2.eq.0.and.ix.eq.0) then
          jflag=1
          r1topa=-x2
          sposar=0.0e0_knd
          sposai=0.0e0_knd
          r1dtop=coef
          sdposr=0.0e0_knd
          sdposi=0.0e0_knd
          if(real(coef).gt.0.0e0_knd) sdposr=real(coef)
          if(aimag(coef).gt.0.0e0_knd) sdposi=aimag(coef)
          end if
          if(iflag.eq.1) then
c          write(40,10)
c10        format(8x,'r1bes: forward series not used.')
          jtop=lm2
          go to 50
          end if
          do 20 j=lm2+1,lim
          jj=j+j+ix
          dnew=-dold*enr(j)*sbesf(jj+m)*prat1(jj+1)
          dnewd=-doldd*enr(j)*sbesdf(jj+m)*prat1(jj+1)
          r1top=r1top+dnew
          r1dtop=r1dtop+dnewd
            if(jflag.eq.1) then
            r1topa=r1topa+dnew
            if(real(dnew).gt.0.0e0_knd) sposar=sposar+real(dnew)
            if(aimag(dnew).gt.0.0e0_knd) sposai=sposai+aimag(dnew)
            end if
          if(real(dnew).gt.0.0e0_knd) sposr=sposr+real(dnew)
          if(aimag(dnew).gt.0.0e0_knd) sposi=sposi+aimag(dnew)
          if(real(dnewd).gt.0.0e0_knd) sdposr=sdposr+real(dnewd)
          if(aimag(dnewd).gt.0.0e0_knd) sdposi=sdposi+aimag(dnewd)
          if((abs(dnew/r1top)+abs(dnewd/r1dtop)).lt.dec) go to 30
            if(abs(r1top).gt.teste) then
            r1top=r1top*testeo
            dnew=dnew*testeo
            sposr=sposr*testeo
            sposi=sposi*testeo
            r1dtop=r1dtop*testeo
            dnewd=dnewd*testeo
            sdposr=sdposr*testeo
            sdposi=sdposi*testeo
            ir1tope=ir1tope+nfac
            kflag=1
              if(jflag.eq.1) then
              r1topa=r1topa*testeo
              sposar=sposar*testeo
              sposai=sposai*testeo
              end if
            end if
          dold=dnew
          doldd=dnewd
20        continue
30      continue
        jtop=min(j,lim)
        jterms=jtop-lm2
c        write(40,40) lim,jterms
c40      format(8x,'r1bes: ',i6,' total terms ',
c     1         'available; forward series converged in ',i6,' terms.')
50	continue
c
c  backward summation of numerator series for r1 and r1d
        if(lm2.lt.1.or.kflag.eq.1) go to 80
        dold=(1.0e0_knd,0.0e0_knd)
        doldd=(1.0e0_knd,0.0e0_knd)
          do 70 j=lm2,1,-1
          jj=j+j+ix
          dnew=-dold/(sbesf(jj+m)*prat1(jj+1)*enr(j))
          dnewd=-doldd/(sbesdf(jj+m)*prat1(jj+1)*enr(j))
            if(j.eq.1.and.ix.eq.0) then
            jflag=1
            term=-x2*dnew
            r1topa=r1top+term
            dnewd=coef*dnewd
            if(real(term).gt.0.0e0_knd) sposar=sposar+real(term)
            if(aimag(term).gt.0.0e0_knd) sposai=sposai+aimag(term)
            end if
          r1top=r1top+dnew
          r1dtop=r1dtop+dnewd
          if(real(dnew).gt.0.0e0_knd) sposr=sposr+real(dnew)
          if(aimag(dnew).gt.0.0e0_knd) sposi=sposi+aimag(dnew)
          if(real(dnewd).gt.0.0e0_knd) sdposr=sdposr+real(dnewd)
          if(aimag(dnewd).gt.0.0e0_knd) sdposi=sdposi+aimag(dnewd)
          if(j.eq.1) go to 80
          if((abs(dnew/r1top)+abs(dnewd/r1dtop)).lt.dec) go to 80
            if(abs(r1top).gt.teste) then
            r1top=r1top*testeo
            dnew=dnew*testeo
            ir1tope=ir1tope+nfac
            r1dtop=r1dtop*testeo
            dnewd=dnewd*testeo
            sposr=sposr*testeo
            sposi=sposi*testeo
            sdposr=sdposr*testeo
            sdposi=sdposi*testeo
            iflag=1
            if(sposr.lt.abs(real(r1top)*dec))
     1           sposr=abs(real(r1top))*dec
            if(sposi.lt.abs(aimag(r1top)*dec))
     1           sposi=abs(aimag(r1top))*dec
            if(sdposr.lt.abs(real(r1dtop)*dec))
     1           sdposr=abs(real(r1dtop))*dec
            if(sdposi.lt.abs(aimag(r1dtop)*dec))
     1           sdposi=abs(aimag(r1dtop))*dec
            end if
60        dold=dnew
          doldd=dnewd
70        continue
80        continue
        nsubr=0
        if(sposr.ne.0.0e0_knd.and.real(r1top).ne.0.0e0_knd) nsubr=
     1         int(log10(sposr/abs(real(r1top))))
        if(nsubr.lt.0) nsubr=0
        if(nsubr.gt.ndec) nsubr=ndec
        nsubi=0
        if(sposi.ne.0.0e0_knd.and.aimag(r1top).ne.0.0e0_knd) nsubi=
     1         int(log10(sposi/abs(aimag(r1top))))
        if(nsubi.lt.0) nsubi=0
        if(nsubi.gt.ndec) nsubi=ndec
        nsubar=0
        nsubai=0
          if(jflag.eq.1) then
          if(real(term).gt.0.0e0_knd) sposar=sposar+real(term)
          nsubar=0
          if(sposar.ne.0.0e0_knd.and.real(r1topa).ne.0.0e0_knd)
     1       nsubar=int(log10(sposar/abs(real(r1topa))))
          if(nsubar.lt.0) nsubar=0
          if(nsubar.gt.ndec) nsubar=ndec
          nsubai=0
          if(aimag(term).gt.0.0e0_knd) sposai=sposai+aimag(term)
          if(sposai.ne.0.0e0_knd.and.aimag(r1topa).ne.0.0e0_knd)
     1       nsubai=int(log10(sposai/abs(aimag(r1topa))))
          if(nsubai.lt.0) nsubai=0
          if(nsubai.gt.ndec) nsubai=ndec
          end if
        ndsubr=0
        if(sdposr.ne.0.0e0_knd.and.real(r1dtop).ne.0.0e0_knd)
     1        ndsubr=int(log10(sdposr/abs(real(r1dtop))))
        if(ndsubr.lt.0) ndsubr=0
        if(ndsubr.gt.ndec) ndsubr=ndec
        ndsubi=0
        if(sdposi.ne.0.0e0_knd.and.aimag(r1dtop).ne.0.0e0_knd)
     1        ndsubi=int(log10(sdposi/abs(aimag(r1dtop))))
        if(ndsubi.lt.0) ndsubi=0
        if(ndsubi.gt.ndec) ndsubi=ndec
c
c  compute r1 and r1d
        r1temp=r1top*sbesn(l+1)*pcoefn/dmfnorm
        iterm=0
        if(abs(r1temp).ne.0.0e0_knd) iterm=int(log10(abs(r1temp)))
        ir1e=ir1tope+ibese(l+1)+ipcoefn-idmfe+iterm
        r1c=r1temp*(ten**(-iterm))
        if(abs(r1c).ge.1.0e0_knd) go to 90
        r1c=r1c*ten
        ir1e=ir1e-1
90      continue
          if(jflag.eq.0) then
          r1ca=r1c
          ir1ea=ir1e
          else
          r1tempa=r1topa*sbesn(l+1)*pcoefn/dmfnorm
          iterm=0
          if(abs(r1tempa).ne.0.0e0_knd) iterm=int(log10(abs(r1tempa)))
          ir1ea=ir1tope+ibese(l+1)+ipcoefn-idmfe+iterm
          r1ca=r1tempa*(ten**(-iterm))
          end if
        continue
        r1dtemp=r1dcoef*r1ca
        r1dtempa=(cc*r1dtop*sbesn(l+1)*sbesdr(l+1)*pcoefn/
     1            dmfnorm)*ten**(ibese(l+1)+ipcoefn+ir1tope-idmfe-
     2            ir1ea)
        r1dc=r1dtemp+r1dtempa
        ndsub1r=0
        if(real(r1dtemp).ne.0.0e0_knd.and.real(r1dc).ne.0.0e0_knd)
     1    ndsub1r=int(log10(abs(real(r1dtemp)/real(r1dc))))
        if(ndsub1r.lt.0) ndsub1r=0
        if(ndsub1r.gt.ndec) ndsub1r=ndec
        ndsub1i=0
        if(aimag(r1dtemp).ne.0.0e0_knd.and.aimag(r1dc).ne.0.0e0_knd)
     1   ndsub1i=int(log10(abs(aimag(r1dtemp)/aimag(r1dc))))
        if(ndsub1i.lt.0) ndsub1i=0
        if(ndsub1i.gt.ndec) ndsub1i=ndec
        n1r=0
        if(real(r1dtemp).ne.0.0e0_knd.and.real(r1dtempa).ne.0.0e0_knd)
     1   n1r=int(log10(abs(real(r1dtemp)/real(r1dtempa))))
        if(n1r.gt.0.and.jflag.eq.1) ndsubr=max(nsubar,ndsubr-n1r)+
     1          ndsub1r
        if(n1r.le.0.and.jflag.eq.1) ndsubr=max(ndsubr,nsubar+n1r)+
     1          ndsub1r
        if(n1r.gt.0.and.jflag.eq.0) ndsubr=max(nsubr,ndsubr-n1r)+ndsub1r
        if(n1r.le.0.and.jflag.eq.0) ndsubr=max(ndsubr,nsubr+n1r)+ndsub1r
        if(ndsubr.gt.ndec) ndsubr=ndec
        n1i=0
        if(aimag(r1dtemp).ne.0.0e0_knd.and.aimag(r1dtempa).ne.
     1     0.0e0_knd)
     2   n1i=int(log10(abs(aimag(r1dtemp)/aimag(r1dtempa))))
        if(n1i.gt.0.and.jflag.eq.1) ndsubi=max(nsubai,ndsubi-n1i)+
     1          ndsub1i
        if(n1i.le.0.and.jflag.eq.1) ndsubi=max(ndsubi,nsubai+n1i)+
     1          ndsub1i
        if(n1i.gt.0.and.jflag.eq.0) ndsubi=max(nsubi,ndsubi-n1i)+ndsub1i
        if(n1i.le.0.and.jflag.eq.0) ndsubi=max(ndsubi,nsubi+n1i)+ndsub1i
        if(ndsubi.gt.ndec) ndsubi=ndec
100     iterm=0
        if(abs(r1dc).ne.0.0e0_knd) iterm=int(log10(abs(r1dc)))
        ir1de=ir1ea+iterm
 	r1dc=r1dc*ten**(-iterm)
        if(abs(r1dc).ge.1.0e0_knd) go to 110
        r1dc=r1dc*ten
        ir1de=ir1de-1
110	continue
c        if(nsubr+ndsubr.gt.0) write(40,120) nsubr,ndsubr
c120     format(15x,'sub. errors in the real parts of the num.',
c     1         ' of r1 and r1d are ',i2,' and ',i2,' digits.')
c        if(nsubi+ndsubi.gt.0) write(40,125) nsubi,ndsubi
c125     format(15x,'sub. errors in the imag. parts of the num.',
c     1         ' of r1 and r1d are ',i2,' and ',i2,' digits.')
        jbes=jtop
        nsub=max(nsubr,nsubi)
        ndsub=max(ndsubr,ndsubi)
130     return
        end
c
c
        subroutine r1eta (l,m,cc,x,eta,nee,limeta,ndec,nex,maxd,maxlp,
     1                    maxj,maxp,minacc,wm,enr,sbesf,sbesn,ibese,
     2                    sbesdf,sbesdr,pdratt,pratb,pratt,pcoefn,
     3                    ipcoefn,pdcoefn,ipdcoefn,ir1ep,r1c,ir1e,
     4                    r1dc,ir1de,naccs1,naccs2,jeta)
c
c  purpose:     To calculate the oblate radial function of the
c               first kind and its first derivative with respect
c               to x, using an expansion of spherical Bessel
c               functions.
c
c  parameters:
c
c     input:    l       : l
c               m       : m
c               cc      : complex c
c               x       : x
c               eta     : value for eta used in calculation
c               nee     : index in the array of eta values in the main
c                         program that corresponds to the value of eta
c                         used in r1eta calculations
c               limeta  : maximum number of terms available in the sums
c                         for r1 and r1d
c               ndec    : number of decimal digits available in
c                         real arithmetic
c               nex     : maximum exponent available in real(knd)
c                         arithmetic
c               maxd    : dimension of enr array
c               maxlp   : maximum  l value desired; dimension
c                         of the sbesn, sbesdr, and ibese arrays
c               maxj    : dimension of sbesf and sbesdf arrays
c               maxp    : dimension of pdratt, pratb, and pratt arrays
c               minacc  : minimum number of accurate decimal digits
c                         that are requested
c               wm      : value of 1 - eta*eta computed in a way that
c                         avoids the subtraction error that would occur
c                         if it were computed directly when eta is near
c                         unity
c                         subtraction error that would occur
c               enr     : array of ratios of successive d coefficients
c               sbesf   : array of ratios of successive spherical
c                         Bessel functions of the same parity
c               sbesn   : array of characteristics for Bessel functions
c               ibese   : array of exponents corresponding to sbesn
c               sbesdf  : array of ratios of successive first
c                         derivatives of spherical Bessel functions of
c                         the same parity
c               sbesdr  : array of ratios of first derivatives of the
c                         spherical Bessel functions to the
c                         corresponding functions
c               pdratt  : array of ratios of successive first
c                         derivatives of the associated Legendre
c                         functions of the first kind of the same parity
c                         (used in numerator series)
c               pratb   : array of ratios of successive associated
c                         Legendre functions of the first kind of the
c                         same parity (used in denominator series)
c               pratt   : array of ratios of successive associated
c                         Legendre functions of the first kind of the
c                         same parity (used in numerator series)
c               pcoefn  : characteristic of the ratio of the numerator
c                         and denominator associated Legendre functions
c                         of the first kind of order m and degree l
c               ipcoefn : exponent corresponding to pcoefn
c               pdcoefn : characteristic of the ratio of the first
c                         derivative of the associated Legendre function
c                         of the first kind in the numerator and the
c                         associated Legendre function of the first kind
c                         in the denominator, both of order m and
c                         degree l
c               ipdcoefn: exponent corresponding to pdcoefn
c               ir1ep   : exponent for the value of r1 for l-1
c
c     output:   r1c     : characteristic of the radial function of the
c                         first kind
c               irie    : exponent corresponding to r1c
c               r1dc    : characteristic of the first derivative with
c                         respect to x of the radial function of the
c                         first kind
c               ir1de   : exponent corresponding to r1dc
c               naccs1  : larger of (1) the subtraction error for either
c                         the real or imaginary part of the numerator
c                         series for r1, whichever one has the larger
c                         magnitude and (2) the subtraction error for
c                         either the real or imaginary part of the
c                         denominator series, whichever one has the
c                         larger magnitude
c               naccs2  : larger of (1) the subtraction error for either
c                         the real or imaginary part of the numerator
c                         series for r1d, whichever one has the larger
c                         magnitude and (2) the subtraction error for
c                         either the real or imaginary part of the
c                         denominator series, whichever one has the
c                         larger magnitude
c               jeta    : maximum number of terms taken in the numerator
c                         and denominator sums for r1 and r1d
c
c        use param
c
c  real(knd) scalars and arrays
        real(knd) c,dec,eta,etas,pcoefn,pdcoefn,r1dcoef1,rm,rm2,
     1             sumdnpr,sumdnpi,sumdpr,sumdpi,sumnpr,sumnpi,ten,
     2             test,testd,teste,testeo,wm,xet,xets,x
        real(knd) pratb(maxp),pratt(maxp),pdratt(maxp)
c  complex(knd) scalars and arrays
        complex(knd) cc,denom,dnew,dnewd,dnewd1,dnewd2,dold,doldd1,
     1               doldd2,reld12,r1c,r1dc,r1dcoef2,r1dtemp,
     2               r1temp,r1test
        complex(knd) enr(maxd),sbesdr(maxlp),sbesn(maxlp),sbesf(maxj),
     1               sbesdf(maxj)
c
c  integer arrays
        integer ibese(maxlp)
c
        ten=10.0e0_knd
        dec=ten**(-ndec-2)
        nfac=nex/3
        if(2*(nfac/2).ne.nfac) nfac=nfac-1
        teste=ten**nfac
        testeo=1.0e0_knd/teste
        ir1tempe=0
        iflag=1
        rm=m
        etas=eta*eta
        xet=sqrt(x*x+wm)
        xets=xet*xet
        r1dcoef1=eta*wm/(xets*xet)
        r1dcoef2=cc*x/xet
        reld12=(r1dcoef2/r1dcoef1)*sbesdr(l+1)*(pcoefn/
     1         pdcoefn)*ten**(ipcoefn-ipdcoefn)
        rm2=rm*2.0e0_knd
        lm2=(l-m)/2
c
c  ix = 0 for l-m even; ix = 1 for l-m odd
        ix=l-m-2*lm2
        lim=limeta/2-ix
c
c  compute radial function of the first kind and its first derivative
c
c  backward series for denominator
        idenom=0
        denom=(1.0e0_knd,0.0e0_knd)
        sumdpr=1.0e0_knd
        sumdpi=0.0e0_knd
        if (lm2.eq.0) go to 20
        dold=(1.0e0_knd,0.0e0_knd)
          do 10 j=lm2,1,-1
          jj=j+j+ix
          dnew=dold/(pratb(jj+1)*enr(j))
          denom=denom+dnew
          if(real(dnew).gt.0.0e0_knd) sumdpr=sumdpr+real(dnew)
          if(aimag(dnew).gt.0.0e0_knd) sumdpi=sumdpi+aimag(dnew)
          if(abs(dnew/denom).lt.dec) go to 20
          dold=dnew
10        continue
20      continue
c
c  forward series for denominator
        dold=(1.0e0_knd,0.0e0_knd)
          do 30 j=lm2+1,lim
          jj=j+j+ix
          dnew=dold*enr(j)*pratb(jj+1)
          denom=denom+dnew
          if(real(dnew).gt.0.0e0_knd) sumdpr=sumdpr+real(dnew)
          if(aimag(dnew).gt.0.0e0_knd) sumdpi=sumdpi+aimag(dnew)
          if(abs(dnew/denom).lt.dec) go to 40
            if(abs(denom).gt.teste) then
            denom=denom*testeo
            dnew=dnew*testeo
            sumdpr=sumdpr*testeo
            sumdpi=sumdpi*testeo
            idenom=idenom+nfac
            end if
25        dold=dnew
30        continue
40      continue
        jden=j
        ndensr=0
        if(sumdpr/real(denom).ne.0.0e0_knd) ndensr=
     1                              int(log10(abs(sumdpr/real(denom))))
        if(ndensr.lt.0) ndensr=0
        if(ndensr.gt.ndec) ndensr=ndec
        ndensi=0
        if(sumdpi.ne.0.0e0_knd) ndensi=int(log10(abs(sumdpi/
     1           aimag(denom))))
        if(ndensi.lt.0) ndensi=0
        if(ndensi.gt.ndec) ndensi=ndec
        iterm=int(log10(abs(denom)))
        idenom=idenom+iterm
        denom=denom*ten**(-iterm)
          if(abs(real(denom)).gt.abs(aimag(denom))) then
          ndens=ndensr
          else
          ndens=ndensi
          end if
c
c  backward series for numerator
        dold=(1.0e0_knd,0.0e0_knd)
        doldd1=dold
        doldd2=reld12
        r1temp=dold
        sumnpr=1.0e0_knd
        sumnpi=0.0e0_knd
        r1dtemp=doldd2
        sumdnpr=0.0e0_knd
        sumdnpi=0.0e0_knd
        if(real(doldd2).gt.0.0e0_knd) sumdnpr=real(doldd2)
        if(aimag(doldd2).gt.0.0e0_knd) sumdnpi=aimag(doldd2)
        if(l.ne.0) r1dtemp=r1dtemp+(1.0e0_knd,0.0e0_knd)
        if(l.ne.0) sumdnpr=sumdnpr+(1.0e0_knd,0.0e0_knd)
        if(lm2.eq.0) go to 60
          do 50 j=lm2,1,-1
          jj=j+j+ix
          dnew=-dold/(sbesf(jj+m)*pratt(jj+1)*enr(j))
          dnewd1=-doldd1/(sbesf(jj+m)*pdratt(jj+1)*enr(j))
          dnewd2=-doldd2/(sbesdf(jj+m)*pratt(jj+1)*enr(j))
          r1temp=r1temp+dnew
          if(real(dnew).gt.0.0e0_knd) sumnpr=sumnpr+real(dnew)
          if(aimag(dnew).gt.0.0e0_knd) sumnpi=sumnpi+aimag(dnew)
          dnewd=dnewd1+dnewd2
          r1dtemp=r1dtemp+dnewd
          if(real(dnewd).gt.0.0e0_knd) sumdnpr=sumdnpr+real(dnewd)
          if(aimag(dnewd).gt.0.0e0_knd) sumdnpi=sumdnpi+aimag(dnewd)
          if(abs(dnew/r1temp)+abs(dnewd/r1dtemp).lt.dec) go to 60
          if(abs(r1temp).lt.teste) go to 45
          r1temp=r1temp*testeo
          dnew=dnew*testeo
          ir1tempe=ir1tempe+nfac
          r1dtemp=r1dtemp*testeo
          dnewd1=dnewd1*testeo
          dnewd2=dnewd2*testeo
          sumnpr=sumnpr*testeo
          sumnpi=sumnpi*testeo
          iflag=0
          if(sumnpr.lt.abs(real(r1temp)*dec))
     1       sumnpr=abs(real(r1temp))*dec
          if(sumnpi.lt.abs(aimag(r1temp)*dec))
     1       sumnpi=abs(aimag(r1temp))*dec
          sumdnpr=sumdnpr*testeo
          sumdnpi=sumdnpi*testeo
          if(sumdnpr.lt.abs(real(r1dtemp)*dec))
     1       sumdnpr=abs(real(r1dtemp))*dec
          if(sumdnpi.lt.abs(aimag(r1dtemp)*dec))
     1       sumdnpi=abs(aimag(r1dtemp))*dec
45        dold=dnew
          doldd1=dnewd1
          doldd2=dnewd2
50      continue
60      continue
        if(m.eq.0.and.jj.eq.2) r1dtemp=r1dtemp-dnewd1
c
c  forward series for numerator
          if(iflag.eq.0) then
          j=lm2
          go to 130
          end if
        dold=(1.0e0_knd,0.0e0_knd)
        doldd1=dold
        doldd2=reld12
        dnewsum=(0.0e0_knd,0.0e0_knd)
        dnewdsum=(0.0e0_knd,0.0e0_knd)
        doldd=reld12
        if(l.ne.0) doldd=reld12+(1.0e0_knd,0.0e0_knd)
          do 110 j=lm2+1,lim-1
          jj=j+j+ix
          dnew=-dold*enr(j)*sbesf(jj+m)*pratt(jj+1)
          dnewd1=-doldd1*enr(j)*sbesf(jj+m)*pdratt(jj+1)
          dnewd2=-doldd2*enr(j)*sbesdf(jj+m)*pratt(jj+1)
          r1temp=r1temp+dnew
          if(real(dnew).gt.0.0e0_knd) sumnpr=sumnpr+real(dnew)
          if(aimag(dnew).gt.0.0e0_knd) sumnpi=sumnpi+aimag(dnew)
          if(abs(dnew).ne.0.0e0_knd) test=abs(dnew/r1temp)
          dnewd=dnewd1+dnewd2
          r1dtemp=r1dtemp+dnewd
          if(real(dnewd).gt.0.0e0_knd) sumdnpr=sumdnpr+real(dnewd)
          if(aimag(dnewd).gt.0.0e0_knd) sumdnpi=sumdnpi+aimag(dnewd)
          if(abs(dnewd).ne.0.0e0_knd) testd=abs(dnewd/r1dtemp)
          if(test.lt.dec.and.testd.lt.dec) go to 130
            if(abs(r1temp).gt.teste) then
            r1temp=r1temp*testeo
            dnew=dnew*testeo
            sumnpr=sumnpr*testeo
            sumnpi=sumnpi*testeo
            ir1tempe=ir1tempe+nfac
            r1dtemp=r1dtemp*testeo
            sumdnpr=sumdnpr*testeo
            sumdnpi=sumdnpi*testeo
            dnewd1=dnewd1*testeo
            dnewd2=dnewd2*testeo
            end if
          dold=dnew
          doldd1=dnewd1
          doldd2=dnewd2
110       continue
130     jnum=j
        naccns1r=0
        if(sumnpr.ne.0.0e0_knd) naccns1r=int(log10(abs(sumnpr/
     1                                   real(r1temp))))
        if(naccns1r.lt.0) naccns1r=0
        if(naccns1r.gt.ndec) naccns1r=ndec
        naccns1i=0
        if(sumnpi.ne.0.0e0_knd) naccns1i=int(log10(abs(sumnpi/
     1                                   aimag(r1temp))))
        if(naccns1i.lt.0) naccns1i=0
        if(naccns1i.gt.ndec) naccns1i=ndec
        naccns2r=0
        if(sumdnpr.ne.0.0e0_knd) naccns2r=int(log10(abs(sumdnpr/
     1                                   real(r1dtemp))))
        if(naccns2r.lt.0) naccns2r=0
        if(naccns2r.gt.ndec) naccns2r=ndec
        naccns2i=0
        if(sumdnpi.ne.0.0e0_knd) naccns2i=int(log10(abs(sumdnpi/
     1                                   aimag(r1dtemp))))
        if(naccns2i.lt.0) naccns2i=0
        if(naccns2i.gt.ndec) naccns2i=ndec
c
c  combining results to form the radial function characteristics
c  r1c and r1dc and corresponding exponents ir1e and ir1de
        r1c=r1temp*sbesn(l+1)*pcoefn/denom
          if(abs(real(r1temp)).gt.abs(aimag(r1temp))) then
          naccs1=max(naccns1r,ndens)
          else
          naccs1=max(naccns1i,ndens)
          end if
        iterm=0
        if(abs(r1c).ne.0.0e0_knd) iterm=int(log10(abs(r1c)))
        ir1e=ir1tempe+ibese(l+1)+ipcoefn-idenom+iterm
        r1c=r1c*ten**(-iterm)
        r1dc=(r1dcoef1*r1dtemp*sbesn(l+1)*pdcoefn/denom)*
     1       ten**(ibese(l+1)+ipdcoefn-idenom-ir1e+ir1tempe)
          if(abs(real(r1dtemp)).gt.abs(aimag(r1dtemp))) then
          naccs2=max(naccns2r,ndens)
          else
          naccs2=max(naccns2i,ndens)
          end if
        iterm=0
        if(abs(r1dc).ne.0.0e0_knd) iterm=int(log10(abs(r1dc)))
        ir1de=ir1e+iterm
 	r1dc=r1dc*ten**(-iterm)
        ir1e=ir1e
        ir1de=ir1de
c        write(40,140) jnum,jden,lim
c140     format(8x,'r1eta: numerator, denominator converged in ',
c     1         i6,' ,',i6,' terms; ',i6,' terms available.')
c        if(naccns1r+naccns1i.gt.0) write(40,150) naccns1r,naccns1i
c150     format(15x,'subt. errors in the real and imag. parts of',
c     1         ' r1 numer. are ',i2,' and ',i2,' digits.')
c        if(naccns2r+naccns2i.gt.0) write(40,155) naccns2r,naccns2i
c155     format(15x,'subt. errors in the real and imag. parts of',
c     1         ' r1d numer. are ',i2,' and ',i2,' digits.')
c        if(ndensr+ndensi.gt.0) write(40,160) ndensr,ndensi
c160     format(15x,'subt. errors in the real and imag. parts of',
c     1         ' the denom. are ',i2,' and ',i2,' digits')
        if(abs(r1c).ge.1.0e0_knd) go to 170
        r1c=r1c*ten
        ir1e=ir1e-1
170     continue
        if(abs(r1dc).ge.1.0e0_knd) go to 180
        r1dc=r1dc*ten
        ir1de=ir1de-1
180     continue
190     jeta=max(jden,jnum)
        return
        end
c
c
        subroutine r2int (l,m,cc,x,limint,ndec,nex,maxd,enr,dc01,idc01,
     1                    maxint,maxmp,maxlp,intlim,rpint1,rpint2,
     2                    pint1,pint2,pint3,pint4,norme,pnorm,ipnorm,
     3                    coefme,coefmo,ipint,r2c,ir2e,r2dc,ir2de,jint,
     4                    coefn,icoefn,isub,isubd)
c
c
c  purpose:     To calculate values of the radial function of the
c               second kind and its first derivative using an integral
c               representation of the radial functions in terms of the
c               angular function of the first kind together with a
c               Neumann function kernal. The angular function is
c               expanded in a series of associated Legendre functions.
c               Gaussian quadrature is used (in subroutine pint) to
c               evaluate the resulting integrals involving associated
c               Legendre functions times the Neumann function kernel.
c               This subroutine performs the summation of the
c               integrals times d coefficients to obtain r2 and r2d.
c
c  parameters:
c
c     input:    l      : l
c               m      : m
c               cc     : complex c
c               x      : x
c               limint : approximately twice the maximum number of
c                        terms available to be taken in the series
c               ndec   : number of decimal digits available in
c                        real arithmetic
c               nex    : maximum exponent in real(knd) arithmetic
c               maxd   : dimension of enr array
c               enr    : d coefficient ratios
c               dc01   : characteristic of the first d coefficient,
c                        either d0 or d1, depending on whether l-m
c                        is even or odd
c               idc01  : exponent (base 10) of the first d coefficient
c               maxint : dimension of pint and rpint arrays
c               maxmp  : dimension of norme array
c               maxlp  : dimension of the pnorm and ipnorm arrays
c               intlim : highest order of Legendre function for which
c                        integrals are not totally inaccurate due to
c                        subtraction errors in their calculation.
c                        series for calculating r2 and r2d will not
c                        include contributions from terms involving
c                        integrals above order intlim
c               rpint1 : arrays of ratios of successive integrals of
c                        either the first or the third kind, depending
c                        on whether l-m is even or odd
c               rpint2 : array of ratios of successive integrals of
c                        either the second or the fourth kind,
c                        depending on whether l-m is even or odd
c               pint1  : array of scaled values for the integrals of
c                        the first kind
c               pint2  : array of scaled values for the integrals of
c                        the second kind
c               pint3  : array of scaled values for the integrals of
c                        the third kind
c               pint4  : array of scaled values for the integrals of
c                        the fourth kind
c               norme  : exponent used to scale the Neumann function
c                        of order m involved in the integrals
c               pnorm  : array of characteristics of the scaling factors
c                        used for the associated Legendre functions in
c                        the integrals to avoid overflow
c               ipnorm : array of exponents (base 10) corresponding to
c                        pnorm
c               coefme : coefficient used to multiply r2 to get one of
c                        the two contributions to r2d when l-m is even
c               coefmo : coefficient used to multiply r2 to get one of
c                        the two contributions to r2d when l-m is odd
c               ipint  : equal to zero the first time r2int is called;
c                        equal to unity otherwise
c
c     output:   r2c    : characteristic of oblate radial function
c                        of the second kind
c               ir2e   : exponent of oblate radial function of the
c                        second kind
c               r2dc   : characteristic of derivative with respect
c                        to x of oblate radial function of the second
c                        kind
c               ir2de  : exponent of derivative with respect to x of
c                        oblate radial function of the second kind
c               jint   : maximum value of the index j in the forward
c                        sum for r2 and r2d, i.e., the highest enr(j)
c                        used
c               coefn  : characteristic of coefficient that is only
c                        calculated once (for l = m) and is then
c                        used for all values of l
c               icoefn : exponent for coefn
c               isub   : larger of the subtraction error in the real
c                        and imginary parts of r2
c               isubd  : larger of the subtraction error in the real
c                        and imaginary parts of r2d
c
c        use param
c
c  real(knd) scalars and arrays
        real(knd) aj,arr,c,coefa,coefn,coefme,coefmo,dec,dcon,ri,rm,
     1             rm2,r2dposr,r2dposi,r2posr,r2posi,ten,teste,
     2             testeo,x
        real(knd) pnorm(maxlp)
c
c  complex(knd) scalars and arrays
        complex(knd) cc,coefl,dnew,dnewd,dold,doldd,dc01,r2c,r2dc,
     1                r2dtemp,r2temp,rs
        complex(knd) enr(maxd),pint1(maxint),pint2(maxint),
     1                pint3(maxint),pint4(maxint),rpint1(maxint),
     2                rpint2(maxint)
c
c  integer arrays
        integer ipnorm(maxlp)
c
        ten=10.0e0_knd
        nfac=nex/3
        if(2*(nfac/2).ne.nfac) nfac=nfac-1
        teste=ten**(nfac)
        testeo=1.0e0_knd/teste
        rm=m
        rm2=rm+rm
        lm2=(l-m)/2
        ix=l-m-2*lm2
        ixx=ix-1
        ixx2=ixx+2
        lim=limint/2-ix-1
        if(limint.gt.intlim-2) lim=intlim/2-ix-1
c
c  compute the leading coefficient
        if(ipint.ne.0) go to 20
        icoefn=norme
        coefn=0.5e0_knd
        if(m.eq.0) go to 20
          do 10 i=1,m
          ri=i
	  coefn=coefn/(ri+ri)
            if(coefn.lt.testeo) then
            coefn=coefn*teste
            icoefn=icoefn-nfac
            end if
10      continue
          iterm=int(log10(abs(coefn)))
          coefn=coefn*ten**(-iterm)
          icoefn=icoefn+iterm
20      continue
        if(ix.eq.0) coefa=(rm2+1.0e0_knd)*coefn
        if(ix.eq.1) coefa=(rm2+3.0e0_knd)*coefn
        if((ix.eq.0).and.(2*(lm2/2).ne.lm2)) coefa=-coefa
        if((ix.eq.1).and.(2*((l-m-1)/4).ne.(l-m-1)/2)) coefa=-coefa
	coefl=coefa/dc01
        icoefl=-idc01+icoefn
        dec=ten**(-ndec-1)
        dcon=dec
        jlow=lm2+1
c
c  compute the integrals involving the angular functions by summing
c  d coefficients times corresponding integrals of Legendre
c  functions
c
c  forward summation of series for r2 and r2d
        iflag=0
        jint=lim
        r2temp=(0.0e0_knd,0.0e0_knd)
        r2dtemp=(0.0e0_knd,0.0e0_knd)
        r2posr=0.0e0_knd
        r2posi=0.0e0_knd
        r2dposr=0.0e0_knd
        r2dposi=0.0e0_knd
        if(jlow.gt.lim) go to 40
        dold=(1.0e0_knd,0.0e0_knd)
        doldd=(1.0e0_knd,0.0e0_knd)
        r2dtemp=doldd
        r2temp=dold
        r2posr=1.0e0_knd
        r2posi=0.0e0_knd
        r2dposr=1.0e0_knd
        r2dposi=0.0e0_knd
          do 30 j=jlow,lim
          dnew=dold*enr(j)*rpint1(j+j+ixx2)
          dnewd=doldd*enr(j)*rpint2(j+j+ixx2)
          r2temp=r2temp+dnew
          r2dtemp=r2dtemp+dnewd
          if(real(dnew).gt.0.0e0_knd) r2posr=r2posr+real(dnew)
          if(real(dnewd).gt.0.0e0_knd) r2dposr=r2dposr+real(dnewd)
          if(aimag(dnew).gt.0.0e0_knd) r2posi=r2posi+aimag(dnew)
          if(aimag(dnewd).gt.0.0e0_knd) r2dposi=r2dposi+aimag(dnewd)
          if((abs(dnew/r2temp)+abs(dnewd/r2dtemp)).lt.dcon) go to 40
          dold=dnew
          doldd=dnewd
30        continue
c
c  backward summation of series for r2 and r2d
40      jint=min(j,lim)
        if(j.eq.0) jint=lim
        if(lm2.lt.1.or.iflag.eq.1) go to 70
        dold=(1.0e0_knd,0.0e0_knd)
        doldd=(1.0e0_knd,0.0e0_knd)
        j=lm2
          do 60 jj=1,lm2
          dnew=dold/(rpint1(j+j+ixx2)*enr(j))
          dnewd=doldd/(rpint2(j+j+ixx2)*enr(j))
          if(j.gt.lim) go to 50
          r2temp=r2temp+dnew
          r2dtemp=r2dtemp+dnewd
          if(real(dnew).gt.0.0e0_knd) r2posr=r2posr+real(dnew)
          if(real(dnewd).gt.0.0e0_knd) r2dposr=r2dposr+real(dnewd)
          if(aimag(dnew).gt.0.0e0_knd) r2posi=r2posi+aimag(dnew)
          if(aimag(dnewd).gt.0.0e0_knd) r2dposi=r2dposi+aimag(dnewd)
          if((abs(dnew/r2temp)+abs(dnewd/r2dtemp)).lt.dcon) go to 70
50        dold=dnew
          doldd=dnewd
          j=j-1
60        continue
70      continue
        isubr=0
        if(r2posr.ne.0.0e0_knd) isubr=int(log10(abs(r2posr/
     1                                real(r2temp))+dec))
        if(isubr.lt.0) isubr=0
        if(isubr.gt.ndec) isubr=ndec
        isubdr=0
        if(r2dposr.ne.0.0e0_knd) isubdr=int(log10(abs(r2dposr/
     1                                  real(r2dtemp))+dec))
        if(isubdr.lt.0) isubdr=0
        if(isubdr.gt.ndec) isubdr=ndec
        isubi=0
        if(r2posi.ne.0.0e0_knd) isubi=int(log10(abs(r2posi/
     1                                aimag(r2temp))+dec))
        if(isubi.lt.0) isubi=0
        if(isubi.gt.ndec) isubi=ndec
        isubdi=0
        if(r2dposi.ne.0.0e0_knd) isubdi=int(log10(abs(r2dposi/
     1                                  aimag(r2dtemp))+dec))
        if(isubdi.lt.0) isubdi=0
        if(isubdi.gt.ndec) isubdi=ndec
        isub=isubr
          if(aimag(r2temp).ne.0.0e0_knd) then
          iterm=int(log10(abs(aimag(r2temp)/real(r2temp))))
          if(iterm.lt.0) isub=max(isubr,isubi+iterm)
          if(iterm.gt.0) isub=max(isubi,isubr-iterm)
          end if
        isubd=isubdr
          if(aimag(r2dtemp).ne.0.0e0_knd) then
          iterm=int(log10(abs(aimag(r2dtemp)/real(r2dtemp))))
          if(iterm.lt.0) isubd=max(isubdr,isubdi+iterm)
          if(iterm.gt.0) isubd=max(isubdi,isubdr-iterm)
          end if
	r2temp=r2temp*coefl*pnorm(l-m+1)
        if(ix.eq.0) r2temp=r2temp*pint1(l-m+1)
        if(ix.eq.1) r2temp=x*r2temp*pint3(l-m+1)
        iterm=0
        if(abs(r2temp).ne.0.0e0_knd) iterm=int(log10(abs(r2temp)))
        ir2e=iterm+ipnorm(l-m+1)+icoefl
        r2c=r2temp*ten**(-iterm)
        if(abs(r2c).ge.1.0e0_knd) go to 80
        r2c=r2c*ten
        ir2e=ir2e-1
80	r2dtemp=-r2dtemp*coefl*pnorm(l-m+1)*cc*x
        rs=r2dtemp
        if(ix.eq.0) r2dtemp=r2dtemp*pint2(l-m+1)+r2temp*coefme
        if(ix.eq.1) r2dtemp=r2dtemp*pint4(l-m+1)*x+r2temp*coefmo
        jsuba=0
        if(ix.eq.0.and.m.ne.0) jsuba=int(log10(abs(r2temp*coefme/
     1                               r2dtemp)+dec))
        if(ix.eq.1) jsuba=int(log10(abs(r2temp*coefmo/r2dtemp)+dec))
        jsubb=0
        if(ix.eq.0.and.m.ne.0) jsubb=int(log10(abs(rs*pint2(l-m+1)/
     1                               r2dtemp)+dec))
        if(ix.eq.1) jsubb=int(log10(abs(rs*x*pint4(l-m+1)/r2dtemp)+
     1                     dec))
        if(m.ne.0.or.ix.eq.1) isubd=max(isub+jsuba,isubd+jsubb,0)
        if(isubd.gt.ndec) isubd=ndec
c        write(40,90) jint,lim,isub,isubd
c90      format(8x,'r2int: converged in ',i5,' terms; 'i5,
c     1         ' available; ',i2,' and ',i2,' digits of sub. error'
c     2         ' in r2 and r2d')
        jterm=0
        if(abs(r2dtemp).ne.0.0e0_knd) jterm=int(log10(abs(r2dtemp)))
        ir2de=jterm+ipnorm(l-m+1)+icoefl
        r2dc=r2dtemp*ten**(-jterm)
        if(abs(r2dc).ge.1.0e0_knd) go to 100
        r2dc=r2dc*ten
        ir2de=ir2de-1
100     continue
        return
        end
c
c
        subroutine r2leg (l,m,cc,x,lnum,minacc,limleg,limdr,iflagp,ndec,
     1                    nex,maxd,maxmp,maxpdr,maxdr,maxq,enr,enrneg,
     2                    drhor,nsdrhor1,nsdrho,dc01,idc01,dneg,idneg,
     3                    nsdneg,dfnorm,idfe,dmfnorm,idmfe,prx,pdrx,qdr,
     4                    qdqr,qdm1,iqdm1,qdl,iqdl,qr,qm1,iqm1,ql,iql,
     5                    fajo,ifajo,ifsub,jsub,termpq,itermpq,ioppsum,
     6                    iopqnsum,r1c,ir1e,r1dc,ir1de,naccr1,naccrpl,
     7                    itestm,r2c,ir2e,r2dc,ir2de,jleg,jlegp,jflagl,
     8                    naccleg,kflagl,nsubleg,nsubdleg,nacccor)
c
c  purpose:     To evaluate the oblate radial function of the
c               second kind and its first derivative with respect
c               to x using the traditional expansion in associated
c               Legendre functions.
c
c  parameters:
c
c     input :   l       : l
c               m       : m
c               c       : complex c
c               x       : radial coordinate x
c               lnum    : number of l values desired
c               minacc  : desired accuracy in decimal digits
c               limleg  : approximately twice the maximum number
c                         of terms available to be taken in qsum,
c                         (sum involving q's time d coefficients)
c               limdr   : maximum number of terms available to be
c                         taken in psum (sum involving p's time
c                         d rho coefficients)
c               iflagp  : integer flag set = 1 if psum series converges
c                         fully; set = 0 otherwise
c               ndec    : number of decimal digits available in
c                         real arithmetic
c               nex     : maximum exponent is real(knd) arithmetic
c               maxd    : dimension of enr array
c               maxmp   : dimension of enrneg array
c               maxpdr  : dimension of prx and pdrx arrays
c               maxdr   : dimension of drhor array
c               maxq    : dimension of qr and qdr arrays
c               enr     : array of d coefficient ratios
c               enrneg  : array of d coefficient ratios with
c                         negative subscripts
c               drhor   : array of d rho coefficient ratios
c               nsdrhor1: subtraction error in the step calculating
c                         drhor(1) from drhor(2)
c               nsdrho  : subtraction error in calculating drhor(1)
c               dc01    : characteristic of the ratio of the first d
c                         coefficient with nonnegative subscript, either
c                         d0 or d1 depending on whether l-m is even or
c                         odd, to the d coefficient with subscript l - m
c               idc01   : exponent (base 10) corresponding to dc01
c               dneg    : characteristic of the ratio of the d
c                         coefficient with subscript -2m+ix to the
c                         d coefficient with subscript ix, where
c                         ix = 0 or 1 depending on whether l - m
c                         is even or odd
c               idneg   : exponent corresponding to dneg
c               nsdneg  : subtraction error in calculating dneg
c               dfnorm  : characteristic of Flammer normalization sum of
c                         d coefficients. equal to the reciprocal of
c                         the value of the d coefficient d(n = l - m)
c                         using this normalization for the angular
c                         functions
c               idfe    : exponent associated with dfnorm
c               dmfnorm : characteristic of Morse-Feshbach normalization
c                         sum of the d coefficients. equal to the
c                         reciprocal of the value of the d coefficient
c                         d(n = l - m) using this normalization for the
c                         angular functions
c               idmfe   : exponent associated with dmfnorm
c               prx     : ratios of successive Legendre functions of
c                         the first kind of the same parity
c               pdrx    : ratios of successive first derivatives of
c                         Legendre functions of the first kind of the
c                         same parity
c               qdr     : ratios of first derivatives of successive
c                         Legendre functions of the second kind
c               qdqr    : array of ratios of derivatives of associated
c                         Legendre functions of the second kind to the
c                         corresponding Legendre function for degrees
c                         from -m to m-1c
c               qdm1    : characteristic of the first derivative of
c                         the associated Legendre function of the second
c                         kind with order m and degree m-1
c               iqdm1   : exponent corresponding to qdm1
c               qdl     : array of characteristics of the first
c                         derivatives of the associated Legendre
c                         functions of the second kind with order m
c                         and degrees from m to m+lnum-1, scaled by
c                                        -m/2
c                         (2m-1)!!(x*x+1)
c               iqdl    : array of exponents corresponding to qdl
c               qr      : array of ratios of successive associated
c                         Legendre functions of the second kind
c               qm1     : characteristic of the associated Legendre
c                         function of the second kind with order m
c                         and degree m-1
c               iqm1    : exponent corresponding to qm1
c               ql      : array of characteristics of the associated
c                         Legendre function of the second kind with
c                         order m and degrees from m to m+lnum-1
c                                                  -m/2
c                         scaled by (2m-1)!!(x*x+1)
c               iql     : array of exponents corresponding to ql
c               fajo    : characteristic of the joining factor of the
c                         second kind
c               ifajo   : exponent corresponding to fajo
c               ifsub   : subtraction error in forming fajo coming from
c                         dfnorm, dmfnorm, and dneg
c               jsub    : larger of the subtraction errors for the
c                         Flammer and the Morse and Feshbach
c                         normalizations
c               termpq  : characteristic of the relative size of the
c                         maximum terms in the positive degree q series
c                         and the p series used to calculate r2 and r2d
c               itermpq : exponent corresponding to termpq
c               ioppsum : integer flag = 0 if psum need not be computed
c                         since its contribution to r2 and r2d is
c                         negligible; = 1 if psum is computed
c               iopqnsum: integer flag = 0 if qnsum need not be computed
c                         since its contribution to r2 and r2d is
c                         negligible; = 1 if qnsum is computed
c               r1c     : characteristic of the radial function of the
c                         first kind
c               ir1e    : exponent corresponding to r1c
c               r1dc    : characteristic of the first derivative of the
c                         radial function of the first kind
c               ir1de   : exponent corresponding to r1dc
c               naccr1  : estimated accuracy of r1c and r1dc
c               naccrpl : degree in decimal digits that r2 = ir1 and
c                         r2d = ir1d
c               itestm  : number of digits of match between the forward
c                         and backward recursions in calculating enr
c
c     output:   r2c     : characteristic of oblate
c                         radial function of the second kind
c               ir2e    : exponent of oblate radial function of the
c                         second kind
c               r2dc    : characteristic of derivative with
c                         respect to x of oblate radial function
c                         of the second kind
c               ir2de   : exponent of derivative with respect to x of
c                         oblate radial function of second kind
c               jleg    : maximum number of terms taken in qsum
c               jlegp   : maximum number of terms taken in psum
c               jflagl  : equal to 1 if a more accurate value
c                         for the leading coefficient drhor(1)
c                         in psum is obtained using the Wronskian;
c                         equal to 0 otherwise.
c               naccleg : Wronskian estimate if jflagl = 0; estimate
c                         of accuracy when jflag1 = 1
c               kflagl  : equal to one if either qsum or psum
c                         becomes so large that the summation
c                         is exited and r2 and r2d are set equal
c                         to 10.0e0_knd**nex; equal to 0 otherwise
c               nsubleg : subtraction error in decimal digits
c                         in the calculation of r2c
c               nsubdleg: subtraction error in decimal digits
c                         in the calculation of r2dc
c               nacccor : subtraction error in forming Wronskian
c                         using r2 and r2d obtained in r2leg
c
c        use param
c
c  real(knd) scalars and arrays
        real(knd) dconp,dconq,dconqn,dec,pdsumpi,pdsumpr,psumpi,psumpr,
     1            qdm1,qm1,qndsumpi,qndsumpr,qnsumpi,qnsumpr,qsumpi,
     2            qsumpr,qdsumpi,qdsumpr,spdsumpi,spdsumpr,spsumpi,
     3            spsumpr,ten,termpq,test,testd,testm,testdm,testp,tm,x
        real(knd) prx(maxpdr),pdrx(maxpdr),qdl(lnum),qdr(maxq),
     1            qdqr(maxmp),ql(lnum),qr(maxq)
c
c  complex(knd) scalars and arrays
        complex(knd) cc,dfnorm,dmfnorm,dneg,dnegjf,dnew,dnewd,dold,
     1               doldd,dc01,psum,pdsum,qndsum,qdsum,qdsump,qnsum,
     2               qsum,qterm,r1c,r1dc,r2c,r2dc,spsum,spdsum,wronc,
     3               wronca,wroncb,wront,xden,xdrhor,xrhs
        complex(knd) drhor(maxdr),enr(maxd),enrneg(maxmp),fajo(lnum+1)
     1
c
c  integer arrays
        integer ifajo(lnum+1),iqdl(lnum),iql(lnum)
c
        ten=10.0e0_knd
        dec=ten**(-ndec-1)
        dconp=dec
        testp=ten**(nex-10)
        lm2=(l-m)/2
        ix=l-m-2*lm2
        imxp=m+m+ix
        ixx=1-ix
        lim1=limleg/2-ix
        lim2=limdr-1
        wront=1.0e0_knd/(cc*(x*x+1.0e0_knd))
        if(ioppsum.eq.0) lim2=0
        rm=m
        tm=rm+rm
        dconq=dec
        dconqn=dec
        dnegjf=dneg*dc01
        if(m.eq.0) dnegjf=dc01
        iterm=int(log10(abs(dnegjf)))
        dnegjf=dnegjf*ten**(-iterm)
        idnegjf=idneg+idc01+iterm
        if(m.eq.0) idnegjf=idc01+iterm
        fajo(l-m+1)=fajo(l-m+1)*dmfnorm*dnegjf/dfnorm
        iterm=int(log10(abs(fajo(l-m+1))))
        fajo(l-m+1)=fajo(l-m+1)*ten**(-iterm)
        ifajo(l-m+1)=ifajo(l-m+1)+idnegjf+iterm+idmfe-idfe
c
c  begin calculation of series for r2
c
c  calculate d*q sum over positive n using pyramid summation
c
c  backward summation
        qsum=(1.0e0_knd,0.0e0_knd)
        qdsum=(1.0e0_knd,0.0e0_knd)
        qsumpr=1.0e0_knd
        qsumpi=0.0e0_knd
        qdsumpr=1.0e0_knd
        qdsumpi=0.0e0_knd
        if(lm2.eq.0) go to 20
        dold=(1.0e0_knd,0.0e0_knd)
        doldd=(1.0e0_knd,0.0e0_knd)
        j=lm2
          do 10 jj=1,lm2
          dnew=-dold/(qr(j+j+imxp)*qr(j+j+imxp-1)*enr(j))
          qsum=qsum+dnew
          if(real(dnew).gt.0.0e0_knd) qsumpr=qsumpr+real(dnew)
          if(aimag(dnew).gt.0.0e0_knd) qsumpi=qsumpi+aimag(dnew)
          dnewd=-doldd/(qdr(j+j+imxp)*qdr(j+j+imxp-1)*enr(j))
          qdsum=qdsum+dnewd
          if(real(dnewd).gt.0.0e0_knd) qdsumpr=qdsumpr+real(dnewd)
          if(aimag(dnewd).gt.0.0e0_knd) qdsumpi=qdsumpi+
     1                                               aimag(dnewd)
            if(abs(qsum).gt.testp.or.abs(qdsum).gt.testp) then
            r2c=(1.0e0_knd,0.0e0_knd)
            r2dc=(1.0e0_knd,0.0e0_knd)
            ir2e=nex
            ir2de=nex
            nsubleg=ndec
            nsubdleg=ndec
            jleg=j
            jlegp=0
            jnleg=0
            kflagl=1
            go to 180
            end if
          if((abs(dnew/qsum)+abs(dnewd/qdsum)).lt.dconq) go to 20
          dold=dnew
          doldd=dnewd
          j=j-1
10        continue
20      continue
c
c  forward summation
        jlow=lm2+1
        dold=(1.0e0_knd,0.0e0_knd)
        doldd=(1.0e0_knd,0.0e0_knd)
          do 30 j=jlow,lim1
          dnew=-dold*enr(j)*qr(j+j+imxp)*qr(j+j+imxp-1)
          qsum=qsum+dnew
          if(real(dnew).gt.0.0e0_knd) qsumpr=qsumpr+real(dnew)
          if(aimag(dnew).gt.0.0e0_knd) qsumpi=qsumpi+aimag(dnew)
          dnewd=-doldd*enr(j)*qdr(j+j+imxp)*qdr(j+j+imxp-1)
          qdsum=qdsum+dnewd
          if(real(dnewd).gt.0.0e0_knd) qdsumpr=qdsumpr+real(dnewd)
          if(aimag(dnewd).gt.0.0e0_knd) qdsumpi=qdsumpi+aimag(dnewd)
          if((abs(dnew/qsum)+abs(dnewd/qdsum)).lt.dconq) go to 40
          dold=dnew
          doldd=dnewd
30        continue
40      continue
        jleg=j
        if(jleg.gt.lim1) jleg=lim1
        nsqsumr=0
        if(real(qsum)*qsumpr.ne.0.0e0_knd) nsqsumr=int(log10(abs(qsumpr
     1                                            /real(qsum))))
        if(nsqsumr.gt.ndec) nsqsumr=ndec
        if(nsqsumr.lt.0) nsqsumr=0
        nsqsumi=0
        if(aimag(qsum)*qsumpi.ne.0.0e0_knd) nsqsumi=
     1                              int(log10(abs(qsumpi/aimag(qsum))))
        if(nsqsumi.gt.ndec) nsqsumi=ndec
        if(nsqsumi.lt.0) nsqsumi=0
        nsqsum=nsqsumr
          if(aimag(qsum).ne.0.0e0_knd) then
          iterm=int(log10(abs(aimag(qsum)/real(qsum))))
          if(iterm.lt.0) nsqsum=max(nsqsumr,nsqsumi+iterm)
          if(iterm.gt.0) nsqsum=max(nsqsumi,nsqsumr-iterm)
          end if
        nsqdsumr=0
        if(real(qdsum)*qdsumpr.ne.0.0e0_knd) nsqdsumr=int(log10(abs
     1                                           (qdsumpr/real(qdsum))))
        if(nsqdsumr.gt.ndec) nsqdsumr=ndec
        if(nsqdsumr.lt.0) nsqdsumr=0
        nsqdsumi=0
        if(aimag(qdsum)*qdsumpi.ne.0.0e0_knd) nsqdsumi=int(log10(abs
     1                                         (qdsumpi/aimag(qdsum))))
        if(nsqdsumi.gt.ndec) nsqdsumi=ndec
        if(nsqdsumi.lt.0) nsqdsumi=0
        nsqdsum=nsqdsumr
          if(aimag(qdsum).ne.0.0e0_knd) then
          iterm=int(log10(abs(aimag(qdsum)/real(qdsum))))
          if(iterm.lt.0) nsqdsum=max(nsqdsumr,nsqdsumi+iterm)
          if(iterm.gt.0) nsqdsum=max(nsqdsumi,nsqdsumr-iterm)
          end if
        qsum=qsum*ql(l-m+1)/(fajo(l-m+1)*termpq)
        iterm=0
        if(abs(qsum).ne.0.0e0_knd) iterm=int(log10(abs(qsum)))
        qsum=qsum*(ten**(-iterm))
        iqsum=iql(l-m+1)-ifajo(l-m+1)-itermpq+iterm
        qdsum=qdsum*qdl(l-m+1)/(fajo(l-m+1)*termpq)
        iterm=0
        if(abs(qdsum).ne.0.0e0_knd) iterm=int(log10(abs(qdsum)))
        qdsum=qdsum*(ten**(-iterm))
        iqdsum=iqdl(l-m+1)-ifajo(l-m+1)-itermpq+iterm
          if(2*(m/2).eq.m) then
          qsum=-qsum
          qdsum=-qdsum
          end if
          if(2*((l-m+1)/4).ne.(l-m+1)/2) then
          qsum=-qsum
          qdsum=-qdsum
          end if
45      continue
c
c  calculate d*q sum over negative n
        qnsum=(0.0e0_knd,0.0e0_knd)
        qndsum=(0.0e0_knd,0.0e0_knd)
        qnsumpr=0.0e0_knd
        qnsumpi=0.0e0_knd
        qndsumpr=0.0e0_knd
        qndsumpi=0.0e0_knd
        iqnsum=0
        iqndsum=0
        j2=0
        nmterm=0
        nsqnsum=0
        nsqndsum=0
        if(iopqnsum.eq.0.or.m.eq.0) go to 90
        nmterm=m
        qnsum=enrneg(m)
        qndsum=qnsum*qdqr(m+m)
        j2=1
        if(ix.eq.1) go to 50
        qnsum=qnsum*qr(m+m-1)
        qndsum=qnsum*qdqr(m+m-1)
50      continue
        if(real(qnsum).gt.0.0e0_knd) qnsumpr=real(qnsum)
        if(aimag(qnsum).gt.0.0e0_knd) qnsumpi=aimag(qnsum)
        if(real(qndsum).gt.0.0e0_knd) qndsumpr=real(qndsum)
        if(aimag(qndsum).gt.0.0e0_knd) qndsumpi=aimag(qndsum)
          if(m.eq.1) then
          jnleg=1
          go to 80
          end if
        dold=qnsum
          do 60 j=2,m
          dnew=-dold*enrneg(m-j+1)*qr(imxp-j-j+1)*qr(imxp-j-j+2)
          qnsum=qnsum+dnew
          if(real(dnew).gt.0.0e0_knd) qnsumpr=qnsumpr+real(dnew)
          if(aimag(dnew).gt.0.0e0_knd) qnsumpi=qnsumpi+aimag(dnew)
          dnewd=dnew*qdqr(imxp-j-j+1)
          qndsum=qndsum+dnewd
          if(real(dnewd).gt.0.0e0_knd) qndsumpr=qndsumpr+real(dnewd)
          if(aimag(dnewd).gt.0.0e0_knd) qndsumpi=qndsumpi+aimag(dnewd)
          if(abs(dnew/qnsum)+abs(dnewd/qndsum).lt.dec) go to 70
          dold=dnew
60        continue
70      jnleg=j
        if(jnleg.gt.m) jnleg=m
80      nsqnsumr=0
        if(real(qnsum)*qnsumpr.ne.0.0e0_knd) nsqnsumr=int(log10(abs
     1                                        (qnsumpr/real(qnsum))))
        if(nsqnsumr.gt.ndec) nsqnsumr=ndec
        if(nsqnsumr.lt.0) nsqnsumr=0
        nsqnsumi=0
        if(aimag(qnsum)*qnsumpi.ne.0.0e0_knd) nsqnsumi=int(log10(abs
     1                                        (qnsumpi/aimag(qnsum))))
        if(nsqnsumi.gt.ndec) nsqnsumi=ndec
        if(nsqnsumi.lt.0) nsqnsumi=0
        nsqnsum=nsqnsumr
          if(aimag(qnsum).ne.0.0e0_knd) then
          iterm=int(log10(abs(aimag(qnsum)/real(qnsum))))
          if(iterm.lt.0) nsqnsum=max(nsqnsumr,nsqnsumi+iterm)
          if(iterm.gt.0) nsqnsum=max(nsqnsumi,nsqnsumr-iterm)
          end if
        nsqnsum=max(nsqnsum,nsdneg)
        nsqndsumr=0
        if(real(qndsum)*qndsumpr.ne.0.0e0_knd) nsqndsumr=int(log10(abs
     1                                         (qndsumpr/real(qndsum))))
        if(nsqndsumr.gt.ndec) nsqndsumr=ndec
        if(nsqndsumr.lt.0) nsqndsumr=0
        nsqndsumi=0
        if(aimag(qndsum)*qndsumpi.ne.0.0e0_knd) nsqndsumi=int(log10(abs
     1                                        (qndsumpi/aimag(qndsum))))
        if(nsqndsumi.gt.ndec) nsqndsumi=ndec
        if(nsqndsumi.lt.0) nsqndsumi=0
        nsqndsum=nsqndsumr
          if(aimag(qndsum).ne.0.0e0_knd) then
          iterm=int(log10(abs(aimag(qndsum)/real(qndsum))))
          if(iterm.lt.0) nsqndsum=max(nsqndsumr,nsqndsumi+iterm)
          if(iterm.gt.0) nsqndsum=max(nsqndsumi,nsqndsumr-iterm)
          end if
        nsqndsun=max(nsqndsum,nsdneg)
        qnsum=qnsum*qm1*dc01/(fajo(l-m+1)*termpq)
        iterm=int(log10(abs(qnsum)))
        qnsum=qnsum*(ten**(-iterm))
        iqnsum=iqm1+idc01-ifajo(l-m+1)-itermpq+iterm
        qnsum=qnsum*(ten**(iqnsum-iqsum))
        qndsum=qndsum*qm1*dc01/(fajo(l-m+1)*termpq)
        iterm=int(log10(abs(qndsum)))
        qndsum=qndsum*(ten**(-iterm))
        iqndsum=iqm1+idc01-ifajo(l-m+1)-itermpq+iterm
        qndsum=qndsum*(ten**(iqndsum-iqdsum))
          if(2*(m/2).ne.m) then
          qnsum=-qnsum
          qndsum=-qndsum
          end if
90      continue
c
c       calculate d(rho|n)*p summation
        psum=(0.0e0_knd,0.0e0_knd)
        pdsum=(0.0e0_knd,0.0e0_knd)
        ipsum=0
        ipdsum=0
        jlegp=0
        nspsumr=0
        nspsumi=0
        nspdsumr=0
        nspdsumi=0
        nspsum=0
        nspdsum=0
        if(ioppsum.eq.0) go to 160
        psum=prx(ixx+1)*drhor(1)
        pdsum=pdrx(ixx+1)*drhor(1)
        dold=psum
        doldd=pdsum
        if(m.ne.0.or.ix.ne.1) go to 100
        pdsum=(0.0e0_knd,0.0e0_knd)
        doldd=drhor(1)
100     continue
        spsum=psum
        spdsum=pdsum
        psumpr=0.0e0_knd
        if(real(psum).gt.0.0e0_knd) psumpr=real(psum)
        psumpi=0.0e0_knd
        if(aimag(psum).gt.0.0e0_knd) psumpi=aimag(psum)
        pdsumpr=0.0e0_knd
        if(real(pdsum).gt.0.0e0_knd) pdsumpr=real(pdsum)
        pdsumpi=0.0e0_knd
        if(aimag(pdsum).gt.0.0e0_knd) pdsumpi=aimag(pdsum)
        spsumpr=psumpr
        spsumpi=psumpi
        spdsumpr=pdsumpr
        spdsumpi=pdsumpi
        testm=1.0e0_knd
        testdm=1.0e0_knd
        iflagp=0
        jlegpf=1
        jlegpd=1
          do 130 j=2,lim2
          dnew=dold*drhor(j)*prx(j+j-ix)
          psum=psum+dnew
          if(real(dnew).gt.0.0e0_knd) psumpr=psumpr+real(dnew)
          if(aimag(dnew).gt.0.0e0_knd) psumpi=psumpi+aimag(dnew)
          dnewd=doldd*drhor(j)*pdrx(j+j-ix)
          pdsum=pdsum+dnewd
          if(real(dnewd).gt.0.0e0_knd) pdsumpr=pdsumpr+real(dnewd)
          if(aimag(dnewd).gt.0.0e0_knd) pdsumpi=pdsumpi+aimag(dnewd)
            if(abs(psum).gt.testp.or.abs(pdsum).gt.testp) then
            r2c=(1.0e0_knd,0.0e0_knd)
            r2dc=(1.0e0_knd,0.0e0_knd)
            ir2e=nex
            ir2de=nex
            nsubleg=ndec
            nsubdleg=ndec
            jleg=j
            jlegp=0
            jnleg=0
            kflagl=1
            go to 180
            end if
          test=abs(dnew/psum)
          testd=abs(dnewd/pdsum)
          if(test.gt.testm.or.test.eq.0.0e0_knd) go to 110
          testm=test
          spsum=psum
          spsumpr=psumpr
          spsumpi=psumpi
          jlegpf=j
110       if(testd.gt.testdm.or.testd.eq.0.0e0_knd) go to 120
          testdm=testd
          spdsum=pdsum
          spdsumpr=pdsumpr
          spdsumpi=pdsumpi
          jlegpd=j
120       if(test+testd.lt.dconp) go to 140
          dold=dnew
          doldd=dnewd
130       continue
        go to 150
140     continue
        iflagp=1
150     jlegp=j
        if(jlegp.gt.lim2) jlegp=lim2
        jlegp=max(jlegpf,jlegpd)
        psum=spsum
        pdsum=spdsum
        psumpr=spsumpr
        psumpi=spsumpi
        pdsumpr=spdsumpr
        pdsumpi=spdsumpi
        nspsumr=0
        if(real(psum)*psumpr.ne.0.0e0_knd) nspsumr=int(log10(abs(psumpr
     1                                            /real(psum))))
        if(nspsumr.gt.ndec) nspsumr=ndec
        if(nspsumr.lt.0) nspsumr=0
        nspsumi=0
        if(aimag(psum)*psumpi.ne.0.0e0_knd) nspsumi=
     1                           int(log10(abs(psumpi/aimag(psum))))
        if(nspsumi.gt.ndec) nspsumi=ndec
        if(nspsumi.lt.0) nspsumi=0
        nspsum=nspsumr
          if(aimag(psum).ne.0.0e0_knd) then
          iterm=int(log10(abs(aimag(psum)/real(psum))))
          if(iterm.lt.0) nspsum=max(nspsumr,nspsumi+iterm)
          if(iterm.gt.0) nspsum=max(nspsumi,nspsumr-iterm)
          end if
        nspsum=max(nspsum,nsdrhor1-nsdrho)
        nspdsumr=0
        if(real(pdsum)*pdsumpr.ne.0.0e0_knd) nspdsumr=int(log10(abs
     1                                           (pdsumpr/real(pdsum))))
        if(nspdsumr.gt.ndec) nspdsumr=ndec
        if(nspdsumr.lt.0) nspdsumr=0
        nspdsumi=0
        if(aimag(pdsum)*pdsumpi.ne.0.0e0_knd) nspdsumi=int(log10(abs
     1                                         (pdsumpi/aimag(pdsum))))
        if(nspdsumi.gt.ndec) nspdsumi=ndec
        if(nspdsumi.lt.0) nspdsumi=0
        nspdsum=nspdsumr
          if(aimag(pdsum).ne.0.0e0_knd) then
          iterm=int(log10(abs(aimag(pdsum)/real(pdsum))))
          if(iterm.lt.0) nspdsum=max(nspdsumr,nspdsumi+iterm)
          if(iterm.gt.0) nspdsum=max(nspdsumi,nspdsumr-iterm)
          end if
        nspdsum=max(nspdsum,nsdrho-nsdrhor1)
        psum=-psum*dnegjf*termpq/fajo(l-m+1)
        iterm=0
        if(psum.ne.0.0e0_knd) iterm=int(log10(abs(psum)))
        psum=psum*(ten**(-iterm))
        ipsum=idnegjf+itermpq-ifajo(l-m+1)+iterm
        psum=psum*(ten**(ipsum-iqsum))
        pdsum=pdsum*dnegjf*termpq/fajo(l-m+1)
        if(m.ne.0) pdsum=-pdsum*rm*x/(x*x+1.0e0_knd)
        iterm=0
        if(pdsum.ne.0.0e0_knd) iterm=int(log10(abs(pdsum)))
        pdsum=pdsum*(ten**(-iterm))
        ipdsum=idnegjf+itermpq-ifajo(l-m+1)+iterm
        pdsum=pdsum*(ten**(ipdsum-iqdsum))
        if(2*((l-m)/2).eq.(l-m)) pdsum=-pdsum
160     continue
        r2c=qsum+qnsum+psum
        r2dc=qdsum+qndsum+pdsum
        wronca=r1c*r2dc*(ten**(ir1e+iqdsum))
        wroncb=r1dc*r2c*(ten**(ir1de+iqsum))
        wronc=wronca-wroncb
        naccleg=-int(log10(abs((wronc-wront)/wront)+dec))
        if(naccleg.lt.0) naccleg=0
        if(naccleg.gt.ndec) naccleg=ndec
        nacccor=-int(log10(abs((wronca-wroncb)/wronca)+dec))
        if(nacccor.lt.0) nacccor=0
        if(nacccor.gt.naccrpl) nacccor=naccrpl
        nacclega=naccleg
        if(naccleg.gt.0) naccleg=min(naccleg+nacccor,ndec-jsub,naccr1)
        nstest=max(nspsum,nspdsum)
        iflag2=0
        if(nsdrhor1.ne.0.and.naccleg.lt.minacc.and.nacclega.gt.1
     1     .and.x.le.0.01e0_knd) iflag2=1
          if(iflag2.eq.1) then
          ncw=int(log10(abs(wronca/wroncb)))
          if(ncw.gt.ndec) ncw=ndec
          if(ncw.lt.-ndec) ncw=-ndec
          nsqc=int(log10(abs(psum/qsum)))
          nsqdc=int(log10(abs(pdsum/qdsum)))
          nsqnc=0
          if(qnsum.ne.0.0e0_knd) nsqdc=int(log10(abs(psum/qnsum)))
          nsqndc=0
          if(qndsum.ne.0.0e0_knd) nsqndc=int(log10(abs(pdsum/qndsum)))
          if(ix.eq.0) nsc=max(nsqsum-nsqc,nsqnsum-nsqnc,nspsum,0)
          if(ix.eq.1) nsc=max(nsqdsum-nsqdc,nsqndsum-nsqndc,nspdsum,0)
          if(ix.eq.0) nacclest=ndec-max(nsc-ncw,0)
          if(ix.eq.1) nacclest=ndec-max(nsc+ncw,0)
          nacclest=min(nacclest,naccr1,ndec-nacccor-2,ndec-ifsub)
          if(nacclest.lt.0) nacclest=0
            if(nacclest.gt.naccleg) then
            xrhs=wront-(qdsum+qndsum)*r1c*ten**(ir1e+iqdsum)+
     1              (qsum+qnsum)*r1dc*ten**(ir1de+iqsum)
            xden=(r1c*pdsum*ten**(ir1e+iqdsum)-r1dc*psum*
     1              ten**(ir1de+iqsum))/drhor(1)
            xdrhor=xrhs/xden
            psum=psum*xdrhor/drhor(1)
            pdsum=pdsum*xdrhor/drhor(1)
            r2c=qsum+qnsum+psum
            r2dc=qdsum+qndsum+pdsum
            jflagl=1
            naccleg=nacclest
            end if
          end if
        nqs=0
        if(qsum/r2c.eq.0.0e0_knd) nqs=-ndec
        if(qsum/r2c.ne.0.0e0_knd) nqs=int(log10(abs(qsum/r2c)))
        nqns=0
        if(qnsum/r2c.eq.0.0e0_knd) nqns=-ndec
        if(m.ne.0.and.iopqnsum.ne.0.and.qnsum/r2c.ne.0.0e0_knd)
     1                nqns=int(log10(abs(qnsum/r2c)))
        nps=0
        if(psum/r2c.eq.0.0e0_knd) nps=-ndec
        if(ioppsum.ne.0.and.psum/r2c.ne.0.0e0_knd)
     1                nps=int(log10(abs(psum/r2c)))
        nsqsum=max(nsqsum,nsdneg)+nqs
        if(nsqsum.lt.0) nsqsum=0
        if(nsqsum.gt.ndec) nsqsum=ndec
        if(jflagl.eq.0) nspsum=max(nspsum,nsdrho)+nps
        if(jflagl.eq.1) nspsum=max(nspsum,nsdrho-nsdrhor1)+nps
        if(nspsum.lt.0) nspsum=0
        if(nspsum.gt.ndec) nspsum=ndec
        nsqnsum=max(nsqnsum,nsdneg)+nqns
        if(nsqnsum.lt.0) nsqnsum=0
        if(nsqnsum.gt.ndec) nsqnsum=ndec
        nsubleg=max(nsqsum,nsqnsum,nspsum,jsub)
        r2dc=qdsum+qndsum+pdsum
        nqds=0
        if(qdsum/r2dc.eq.0.0e0_knd) nqds=-ndec
        if(qdsum/r2dc.ne.0.0e0_knd) nqds=int(log10(abs(qdsum/r2dc)))
        nqnds=0
        if(qndsum/r2dc.eq.0.0e0_knd) nqnds=-ndec
        if(m.ne.0.and.iopqnsum.ne.0.and.qndsum/r2dc.ne.0.0e0_knd)
     1                nqnds=int(log10(abs(qndsum/r2dc)))
        if(qnsum/r2c.eq.0.0e0_knd.and.qndsum/r2dc.eq.0.0e0_knd)
     1     iopqnsum=0
        npds=0
        if(pdsum/r2c.eq.0.0e0_knd) npds=-ndec
        if(ioppsum.ne.0.and.pdsum/r2dc.ne.0.0e0_knd)
     1                npds=int(log10(abs(pdsum/r2dc)))
        if(psum/r2c.eq.0.0e0_knd.and.pdsum/r2dc.eq.0.0e0_knd) ioppsum=0
        nsqdsum=max(nsqdsum,nsdneg)+nqds
        if(nsqdsum.lt.0) nsqdsum=0
        if(nsqdsum.gt.ndec) nsqdsum=ndec
        if(jflagl.eq.0) nspdsum=max(nspdsum,nsdrho)+npds
        if(jflagl.eq.1) nspdsum=max(nspdsum,nsdrho-nsdrhor1)+npds
        if(nspdsum.lt.0) nspdsum=0
        if(nspdsum.gt.ndec) nspdsum=ndec
        nsqndsum=max(nsqndsum,nsdneg)+nqnds
        if(nsqndsum.lt.0) nsqndsum=0
        if(nsqndsum.gt.ndec) nsqndsum=ndec
        nsubdleg=max(nsqdsum,nsqndsum,nspdsum,jsub)
        naccleg=min(naccleg,ndec-max(nsubleg,nsubdleg))
        if(naccleg.lt.0) naccleg=0
        if(ioppsum.ne.0.and.naccleg.gt.2.and.nps.lt.(-ndec-ndec).and.
     1      npds.lt.(-ndec-ndec)) ioppsum=0
        if(iopqnsum.ne.0.and.naccleg.ne.0.and.nqns.lt.(-ndec-ndec).and.
     1      nqnds.lt.(-ndec-ndec)) iopqnsum=0
        iterm=int(log10(abs(r2c)))
        r2c=r2c*(ten**(-iterm))
        ir2e=iqsum+iterm
        if(abs(r2c).ge.1.0e0_knd) go to 170
        r2c=r2c*ten
        ir2e=ir2e-1
170     continue
        iterm=int(log10(abs(r2dc)))
        r2dc=r2dc*(ten**(-iterm))
        ir2de=iqdsum+iterm
        if(abs(r2dc).ge.1.0e0_knd) go to 180
        r2dc=r2dc*ten
        ir2de=ir2de-1
180     continue
c        if(ioppsum.eq.1.and.iopqnsum.eq.1) write(40,190) jleg,jlegp,
c     1                            jnleg,lim1,lim2,m,nsubleg,nsubdleg
c190     format(8x,'r2leg: qsum, psum and qnsum series converged in ',i6,
c     1        ',' i6,' and ',i4,' terms; ',i6,',' i6,' and ' i4,
c     2        ' terms avail.',/,15x,i2,' and ',i2,' digits of sub.',
c     3        ' error in r2 and r2d.')
c        if(ioppsum.eq.1.and.iopqnsum.eq.0) write(40,200) jleg,jlegp,
c     1                                     lim1,lim2,nsubleg,nsubdleg
c200     format(8x,'r2leg: qsum and psum series converged in ',i6,
c     1        ' and ',i6,' terms; ',i6,' and ',i6,' terms avail.',/,
c     2        15x,i2,' and ',i2,' digits of sub. error in r2 and r2d;',
c     3        ' qnsum is negligible.')
c        if(ioppsum.eq.0.and.iopqnsum.eq.1) write(40,210) jleg,jnleg,
c     1                                     lim1,m,nsubleg,nsubdleg
c210     format(8x,'r2leg: qsum and qnsum series converged in ',i6,
c     1        ' and ',i4,' terms; ',i6,' and ',i4,' terms avail.',/,
c     2         15x,i2,' and ',i2,' digits of sub. error in r2 and r2d;'
c     3         ' psum is negligible.')
c        if(ioppsum.eq.0.and.iopqnsum.eq.0) write(40,220) jleg,lim1,
c     1                                            nsubleg,nsubdleg
c220     format(8x,'r2leg: qsum series converged in ',i6,' terms with ',
c     1         i6,' terms avail.; 'i2,' and ',i2,' digits of',/,15x,
c     2         'sub. error in r2 and r2d; psum and qnsum are ',
c     3         'negligible.')
c        if(jflagl.eq.1) WRITE(40,230)
c230     format(15x,'Wronskian used to improve accuracy of the',
c     1          ' the leading psum coefficient drhor(1).')
        return
        end
c
c
        subroutine r2leg1(l,m,cc,x,limq,maxq,ndec,eigval,qr,qdr,
     1                    qm0,qdm0,r1c,ir1e,r2c,ir2e,r2dc,ir2de,jleg1)
c
c  purpose:     To evaluate the oblate radial function of the
c               second kind and its first derivative with respect
c               to x using an expansion in associated Legendre
c               functions of the second kind. This expansion is
c               due to Baber and Hasse.
c
c  parameters:
c
c     input :   l       : l
c               m       : m
c               cc      : complex c
c               x       : radial coordinate x
c               limq    : the maximum number of terms available
c               maxq    : dimension of qr,qdr,aratio,coef1,
c                         coef2,coef3
c               ndec    : number of decimal digits available in real
c                         arithemetic
c               eigval  : eigenvalue
c               qr      : array of ratios of successive associated
c                         Legendre functions of the second kind
c                                                 m    m
c                         beginning with qr(1) = Q  / Q . Note that
c                                                 1    0
c                         this differs from the choice used in r2leg
c                         where functions with negative degree are
c                         also used and where ratios involving degrees
c                         less than m are inverted.
c               qdr     : ratios of first derivatives of successive
c                         Legendre functions of the second kind
c                                                   m     m
c                         beginning with qdr(1) = Q'  / Q' . Note that
c                                                   1     0
c                         this differs from the choice used in r2leg
c                         where functions with negative degree are
c                         also used and where ratios involving a degree
c                         less than m are inverted.
c               qm0     : associated Legendre function of the second
c                         kind with order m and degree 0.
c               qdm0    : first derivative of the associated Legendre
c                         function of the second kind with order m and
c                         degree 0.
c               r1c     : charcteristic of the corresponding radial
c                         function of the first kind
c               ir1e    : exponent of the corresponding radial function
c                         of the first kind
c
c     output:   r2c     : characteristic of oblate
c                         radial function of the second kind
c               ir2e    : exponent of oblate radial function of the
c                         second kind
c               r2dc    : characteristic of derivative with
c                         respect to x of the oblate radial function
c                         of the second kind
c               ir2de   : exponent of derivative with respect to x of
c                         the oblate radial function of second kind
c               jleg1   : number of terms taken in the series
c
c        use param
c
c  real(knd) scalars and arrays
        real(knd) aj,dec,em,qm0,qdm0,ten,term,x,xx,qr(maxq),qdr(maxq)
c
c  complex(knd) scalars
        complex(knd) arg,cc,eigval,ea,eb,e3,ra,rb,r1c,r1,r2,r2c,r2d,
     1               r2dc,r3,sumi,sumdi,sumr,sumdr,termi,termdi,termr,
     2               termdr
c  complex(knd) arays
        complex(knd) aratio(maxq),coef1(maxq),coef2(maxq),coef3(maxq)
c
        ten=10.0e0_knd
        em=m
        m2=m+m
        dec=ten**(-ndec-1)
        arg=cc*x
        xx=x*x+1.0e0_knd
        r1=r1c*(ten**(ir1e))
        ea=sin(arg)/cc
        ra=cos(arg)/cc
          if(m.ne.0) then
          ea=ea/em
          ra=ra/em
          end if
        la=l-(l/4)*4+1
          if(la.eq.1) then
          r3=ra
          e3=-ea
          end if
          if(la.eq.2) then
          e3=ra
          r3=ea
          end if
          if(la.eq.3) then
          r3=-ra
          e3=ea
          end if
          if(la.eq.4) then
          e3=-ra
          r3=-ea
          end if
        m2=m+m
          do 10 j=1,limq
          n=j-m-1
          coef1(j)=cc*2*(n+m+1)*(n+m2+1)/(n+n+m2+3)
          term=(n+m)*(n+m+1)
          coef2(j)=term-eigval-cc*cc
          coef3(j)=cc*(n+n)*(n+m)/(n+n+m2-1)
10        continue
        aratio(1)=coef2(1)/coef1(1)
        aratio(limq)=0.0e0_knd
          do 20 j=limq-1,m+2,-1
          aratio(j-1)=coef3(j)/(coef1(j)*aratio(j)-coef2(j))
          if(abs(aratio(j-1)).gt.1.0e0_knd) go to 30
20        continue
30      continue
          do 40 n=2,j-1
          aratio(n)=(coef2(n)+coef3(n)/aratio(n-1))/coef1(n)
40        continue
        sumi=(1.0e0_knd,0.0e0_knd)
        sumdi=(1.0e0_knd,0.0e0_knd)
        termi=(1.0e0_knd,0.0e0_knd)
        termdi=(1.0e0_knd,0.0e0_knd)
        sumr=aratio(1)*qr(1)
        sumdr=aratio(1)*qdr(1)
        termr=sumr
        termdr=sumdr
          do 50 j=2,limq-2,2
          termi=-termi*aratio(j-1)*aratio(j)*qr(j-1)*qr(j)
          termdi=-termdi*aratio(j-1)*aratio(j)*qdr(j-1)*qdr(j)
          termr=-termr*aratio(j)*aratio(j+1)*qr(j)*qr(j+1)
          termdr=-termdr*aratio(j)*aratio(j+1)*qdr(j)*qdr(j+1)
          sumi=sumi+termi
          sumdi=sumdi+termdi
          sumr=sumr+termr
          sumdr=sumdr+termdr
          if(abs(termi/sumi)+abs(termdi/sumdi).gt.dec) go to 50
          if(abs(termr/sumr)+abs(termdr/sumdr).gt.dec) go to 50
          go to 60
50        continue
60      continue
        jleg1=min(j+1,limq-1)
        sumi=sumi*qm0
        sumdi=sumdi*qdm0
        sumr=sumr*qm0
        sumdr=sumdr*qdm0
        r2=r3*sumi+e3*sumr
        r2d=r3*sumdi+e3*sumdr
          if(4*((m+3)/4).eq.m+3) then
          r2d=-r2d
          end if
          if(2*(m/2).ne.m) then
          r2=-r2
          end if
        ir2e=int(log10(abs(r2)))
        r2c=r2*(ten**(-ir2e))
        if(abs(r2c).ge.1.0e0_knd) go to 70
        r2c=r2c*ten
        ir2e=ir2e-1
        r2d=r2d+cc*r1
70      continue
        ir2de=int(log10(abs(r2d)))
        r2dc=r2d*(ten**(-ir2de))
        if(abs(r2dc).ge.1.0e0_knd) go to 80
        r2dc=r2dc*ten
        ir2de=ir2de-1
80      continue
c        write(40,90) jleg1,limq
c90      format(8x,'r2leg1: series converged in ',i6,
c     1         ' terms; ',i6,' terms avail.')
        return
        end
c
c
        subroutine r2neu0 (l,m,cc,x,limneu,ndec,nex,maxd,maxlp,maxn,
     1                     minacc,enr,sneuf,sneun,ineue,sneudf,
     2                     sneudr,dfnorm,idfe,r1dc,ir1de,r2c,ir2e,r2dc,
     3                     ir2de,jneu,jtest,nsub0)
c
c  purpose:     To calculate the oblate radial function of the
c               second kind and its first derivative with respect
c               to x, using a series expansion of spherical
c               Neumann functions with argument c*sqrt(x*x+1);
c               (eta = 0)
c
c  parameters:
c
c     input:    l      : l
c               m      : m
c               cc     : complex c
c               x      : x
c               limneu : maximum number of terms to be taken in the
c                        series summations for r2 and r2d
c               ndec   : number of decimal digits available in
c                        real arithmetic
c               nex    : maximum exponent if real arithmetic
c               maxd   : dimension of enr array
c               maxlp  : dimension of the sneun, sneudn, ineue, and
c                        ineude arrays
c               maxn   : dimension of sneuf and sneudf arrays
c               minacc : number of decimal digits of desired accuracy
c                        of the resulting radial functions
c               enr    : array of ratios of successive d coefficients
c               sneuf  : array of ratios of successive spherical Neumann
c                        functions of the same parity
c               sneun  : array of characteristics for Neumann functions
c               ineue  : array of exponents for Neumann functions
c               sneudf : array of ratios of successive first derivatives
c                        of spherical Neumann functions of same parity
c               sneudr : array of ratios of first derivatives of Neumann
c                        functions to the corresponding functions
c               dfnorm : characteristic of Flammer normalization sum of
c                        d coefficients. equal to the reciprocal of
c                        the value of the d coefficient d(n = l - m)
c                        using this normalization for the angular
c                        functions
c               idfe   : exponent associated with dfnorm
c               r1dc   : charcteristic of corresponding first derivative
c                        of the radial function of the first kind
c               ir1de  : exponent of corresponding first derivative of
c                        the radial function of the first kind
c
c     output:   r2c    : characteristic of oblate radial function
c                        of the second kind
c               ir2e   : exponent of oblate radial function of the
c                        second kind
c               r2dc   : characteristic of derivative with respect
c                        to x of oblate radial function of the second
c                        kind
c               ir2de  : exponent of derivative with respect to x of
c                        oblate radial function of the second kind
c               jneu   : index of term where best convergence is
c                        achieved for r2 or for r2d, whichever term is
c                        larger
c               jtest  : smaller of the number of digits of convergence
c                        of the forward sums for r2 and r2d
c               nsub0  : larger of the subtraction errors in the series
c                        for r2 and r2d
c
c        use param
c
c  real(knd) scalars
        real(knd) con,dconb,dconf,dconi,dec,rj1,rj2,rm,r2est,sumpi,
     1            sumpr,sumdpi,sumdpr,ten,test,testd,testdm,teste,
     2            testeo,testm,txi,txr,txdi,txdr,x
c  complex(knd) scalars and arrays
        complex(knd) cc,dfnorm,dnew,dnewd,dold,doldd,r1dc,r2,r2c,r2d,
     1               r2dc,r2dstore,r2dtemp,r2temp,sr2temp,sr2dtemp,
     2               sumcoef
        complex(knd) enr(maxd),sneudr(maxlp),sneun(maxn),sneuf(maxn),
     1               sneudf(maxn)
c
c  integer arrays
        integer ineue(maxn)
c
        ten=10.0e0_knd
        nfac=nex/3
        if(2*(nfac/2).ne.nfac) nfac=nfac-1
        teste=ten**(nfac)
        testeo=1.0e0_knd/teste
        iscale=0
        rm=m
        dec=ten**(-ndec-1)
        dconf=ten**(-ndec-2)
        dconi=ten**(ndec+5)
        iexp=-ir1de-ineue(l+1)
          if(iexp.lt.nfac) then
          sumcoef=(ten**(iexp))/
     1            (cc*(x*x+1.0e0_knd)*r1dc*sneun(l+1))
          r2est=abs(sumcoef*dfnorm)*ten**(idfe)
          else
          r2est=ten**nfac
          end if
        dconb=r2est/dconi
        con=x/sqrt(x*x+1.0e0_knd)
        lm2=(l-m)/2
c
c  ix = 0 for l-m even; ix = 1 for l-m odd
        ix=l-m-2*lm2
        mml=m+m-1+ix
        lim=limneu/2-ix
c
c  compute radial function of the second kind
c
c  backward series
        r2temp=(1.0e0_knd,0.0e0_knd)
        sumpr=1.0e0_knd
        sumpi=0.0e0_knd
        r2dtemp=(1.0e0_knd,0.0e0_knd)
        sumdpr=1.0e0_knd
        sumdpi=0.0e0_knd
        if(r2est.gt.dconi) go to 20
        if (lm2.eq.0) go to 20
        dold=(1.0e0_knd,0.0e0_knd)
        doldd=(1.0e0_knd,0.0e0_knd)
          do 10 j=lm2,1,-1
          jj=j+j+ix
          rj1=jj-ix
          rj2=jj+mml
          dnew=dold*rj1/(rj2*sneuf(jj+m)*enr(j))
          dnewd=doldd*rj1/(rj2*sneudf(jj+m)*enr(j))
          r2temp=r2temp+dnew
          r2dtemp=r2dtemp+dnewd
          if(real(dnew).gt.0.0e0_knd) sumpr=sumpr+real(dnew)
          if(aimag(dnew).gt.0.0e0_knd) sumpi=sumpi+aimag(dnew)
          if(real(dnewd).gt.0.0e0_knd) sumdpr=sumdpr+real(dnewd)
          if(aimag(dnewd).gt.0.0e0_knd) sumdpi=sumdpi+aimag(dnewd)
          dold=dnew
          doldd=dnewd
10      continue
20      continue
c
c  forward series
        dold=(1.0e0_knd,0.0e0_knd)
        doldd=(1.0e0_knd,0.0e0_knd)
        testm=1.0e0_knd
        testdm=1.0e0_knd
        sr2temp=r2temp
        sr2dtemp=r2dtemp
        jneu=lim
        itest=0
          do 70 j=lm2+1,lim
          jj=j+j+ix
          rj1=jj-ix
          rj2=jj+mml
          dnew=dold*enr(j)*sneuf(jj+m)*rj2/rj1
          dnewd=doldd*enr(j)*sneudf(jj+m)*rj2/rj1
          r2temp=r2temp+dnew
          r2dtemp=r2dtemp+dnewd
          if(real(dnew).gt.0.0e0_knd) sumpr=sumpr+real(dnew)
          if(aimag(dnew).gt.0.0e0_knd) sumpi=sumpi+aimag(dnew)
          if(real(dnewd).gt.0.0e0_knd) sumdpr=sumdpr+real(dnewd)
          if(aimag(dnewd).gt.0.0e0_knd) sumdpi=sumdpi+aimag(dnewd)
          test=abs(dnew/r2temp)
          testd=abs(dnewd/r2dtemp)
          if(test+testd.lt.testm+testdm) go to 30
          go to 40
30        if(test.ne.0.0e0_knd) testm=test
          if(testd.ne.0.0e0_knd) testdm=testd
          sr2temp=r2temp
          txr=sumpr
          txi=sumpi
          sr2dtemp=r2dtemp
          txdr=sumdpr
          txdi=sumdpi
          jneu=j
40        continue
            if(test+testd.lt.dconf) then
              if(itest.eq.1) then
              go to 90
              else
              itest=1
              end if
            else
            itest=0
            end if
            if(abs(r2temp).gt.teste) then
            r2temp=r2temp*testeo
            r2dtemp=r2dtemp*testeo
            sr2temp=sr2temp*testeo
            sr2dtemp=sr2dtemp*testeo
            sumpr=sumpr*testeo
            sumpi=sumpi*testeo
            sumdpr=sumdpr*testeo
            sumdpi=sumdpi*testeo
            txr=txr*testeo
            txi=txi*testeo
            txdr=txdr*testeo
            txdi=txdi*testeo
            dnew=dnew*testeo
            dnewd=dnewd*testeo
            iscale=iscale+nfac
            end if
          dold=dnew
          doldd=dnewd
70        continue
        r2temp=sr2temp
        r2dtemp=sr2dtemp
90      continue
        nterms=min(j,lim)
        jtestm=-int(log10(testm+dec))
        if(jtestm.lt.0) jtestm=0
        if(jtestm.gt.ndec) jtestm=ndec
        jtestdm=-int(log10(testdm+dec))
        if(jtestdm.lt.0) jtestdm=0
        if(jtestdm.gt.ndec) jtestdm=ndec
        jtest=min(jtestm,jtestdm)
        naccs1r=0
        if(txr.ne.0.0e0_knd) naccs1r=int(log10(abs(txr/real(r2temp))))
        if(naccs1r.lt.0) naccs1r=0
        if(naccs1r.gt.ndec) naccs1r=ndec
        naccs1i=0
        if(txi.ne.0.0e0_knd) naccs1i=int(log10(abs(txi/aimag(r2temp))))
        if(naccs1i.lt.0) naccs1i=0
        if(naccs1i.gt.ndec) naccs1i=ndec
        naccs2r=0
        if(txdr.ne.0.0e0_knd) naccs2r=int(log10(abs(txdr/
     1                                real(r2dtemp))))
        if(naccs2r.lt.0) naccs2r=0
        if(naccs2r.gt.ndec) naccs2r=ndec
        naccs2i=0
        if(txdi.ne.0.0e0_knd) naccs2i=int(log10(abs(txdi/
     1                                aimag(r2dtemp))))
        if(naccs2i.lt.0) naccs2i=0
        if(naccs2i.gt.ndec) naccs2i=ndec
        naccs1=max(naccs1r,naccs1i)
        naccs2=max(naccs2r,naccs2i)
c
c  combining results to form the radial function characteristics
c  r2c and r2dc and corresponding exponents ir2e and ir2de
        r2=r2temp*sneun(l+1)/dfnorm
        if(ix.eq.1) r2=r2*con
        iterm=int(log10(abs(r2)))
        ir2e=ineue(l+1)-idfe+iterm+iscale
        r2c=r2*ten**(-iterm)
        if(abs(r2c).ge.1.0e0_knd) go to 100
        r2c=r2c*ten
        ir2e=ir2e-1
100	continue
        r2d=r2dtemp*sneun(l+1)*cc*con*sneudr(l+1)/dfnorm
        r2dstore=r2d*con
        ndsub=0
        if(ix.eq.1) r2d=r2dstore+r2/(x*(x*x+1.0e0_knd))
        if(ix.eq.1) ndsub=-log10(abs(r2d/r2dstore))
        if(ndsub.lt.0) ndsub=0
        naccs2=naccs2+ndsub
        nsub0=max(naccs1,naccs2)
c        write(40,110) nterms,lim,jtestm,naccs1,jtestdm,naccs2
c110     format(8x,'r2neu0 (eta=0) : numerator converged in ',i6,
c     1         ' terms; ',i6,' terms available.',/,15x,'r2 converged ',
c     2         'to ',i3,' digits with ',i3,' digits sub. error.',/,
c     3         15x,'r2d converged to ',i3,' digits with ',i3,' digits ',
c     4         'sub. error.')
c        if(ix.eq.1.and.ndsub.gt.0) write(40,120) ndsub
c120     format(15x,'subtraction error in forming r2d = ',i2,' digits.')
        iterm=int(log10(abs(r2d)))
        ir2de=ineue(l+1)-idfe+iterm+iscale
 	r2dc=r2d*ten**(-iterm)
        if(abs(r2dc).ge.1.0e0_knd) go to 130
        r2dc=r2dc*ten
        ir2de=ir2de-1
130	continue
        return
        end
c
c
        subroutine r2eta (l,m,cc,x,eta,nee,limeta,ndec,maxd,maxlp,maxn,
     1                    maxp,minacc,wm,enr,sneuf,sneun,ineue,
     2                    sneudf,sneudr,pdratt,pratb,pratt,pcoefn,
     3                    ipcoefn,pdcoefn,ipdcoefn,r1c,ir1e,r1dc,ir1de,
     4                    naccr1,naccrpl,naccnmax,naccr,r2c,ir2e,r2dc,
     5                    ir2de,nacceta,jeta,naccd)
c
c  purpose:     To calculate the oblate radial function of the
c               second kind and its first derivative with respect
c               to x, using an expansion of spherical Neumann
c               functions.
c
c  parameters:
c
c     input:    l       : l
c               m       : m
c               cc      : complex c
c               x       : x
c               eta     : value for eta used in calculation
c               nee     : index in the array of eta values in the main
c                         program that corresponds to the value of eta
c                         used in r2eta calculations
c               limeta  : maximum number of terms available in the sums
c                         for r2 and r2d
c               ndec    : number of decimal digits available in
c                         real arithmetic
c               maxd    : dimension of enr array
c               maxlp   : maximum  l value desired; dimension
c                         of the sneun, sneudr, and ineue arrays
c               maxn    : dimension of sneuf and sneudf arrays
c               maxp    : dimension of pdratt, pratb, and pratt arrays
c               minacc  : minimum number of accurate decimal digits
c                         that are requested
c               wm      : value of 1 - eta*eta computed in a way that
c                         avoids the subtraction error that would occur
c                         if it were computed directly when eta is near
c                         unity
c               enr     : array of ratios of successive d coefficients
c               sneuf   : array of ratios of successive spherical
c                         Neumann functions of the same parity
c               sneun   : array of characteristics for Neumann functions
c               ineue   : array of exponents corresponding to sneun
c               sneudf  : array of ratios of successive first
c                         derivatives of spherical Neumann functions of
c                         the same parity
c               sneudr  : array of ratios of first derivatives of the
c                         spherical Neumann functions to the
c                         corresponding functions
c               pdratt  : array of ratios of successive first
c                         derivatives of the associated Legendre
c                         functions of the first kind of the same parity
c                         (used in numerator series)
c               pratb   : array of ratios of successive associated
c                         Legendre functions of the first kind of the
c                         same parity (used in denominator series)
c               pratt   : array of ratios of successive associated
c                         Legendre functions of the first kind of the
c                         same parity (used in numerator series)
c               pcoefn  : characteristic of the ratio of the numerator
c                         and denominator associated Legendre functions
c                         of the first kind of order m and degree l
c               ipcoefn : exponent corresponding to pcoefn
c               pdcoefn : characteristic of the ratio of the first
c                         derivative of the associated Legendre function
c                         of the first kind in the numerator and the
c                         associated Legendre function of the first kind
c                         in the denominator, both of order m and
c                         degree l
c               ipdcoefn: exponent corresponding to pdcoefn
c               r1c     : characteristic of the radial function of the
c                         first kind (calculated in r1bes)
c               irie    : exponent corresponding to r1c
c               r1dc    : characteristic of the first derivative with
c                         respect to x of the radial function of the
c                         first kind (calculated in r1bes)
c               ir1de   : exponent corresponding to r1dc
c               naccr1  : estimated accuracy of r1 and r1d
c               naccrpl : degree in decimal digits to which r2 = ir1
c                         and r2d = ir1d
c               naccnmax: maximum accuracy (in decimal digits) obtained
c                         for the current value of l from previous
c                         r2eta calculations
c               naccr   : accuracy of radial functions calculated for
c                         this value of l earlier using other methods.
c                         if this is the only method used, naccr is
c                         given either by the default value of -1 or
c                         the estimate based on the number of matching
c                         digits of paired eigenvalues
c
c     output:   r2c     : characteristic of the oblate radial function
c                         of the second kind
c               ir2e    : exponent of the oblate radial function of the
c                         second kind
c               r2dc    : characteristic of the first derivative with
c                         respect to x of the oblate radial function
c                         of the second kind
c               ir2de   : exponent corresponding to r2dc
c               nacceta : estimated number of accurate decimal digits in
c                         r2 and r2d, computed from the Wronskian
c               jeta    : maximum number of terms taken in the numerator
c                         sums for r2 and r2d
c               naccd:    estimated accuracy of the denominator series
c
c     input/
c     output:   naccnmax: maximum accuracy (in decimal digits) obtained
c                         for the current value of l from all previous
c                         r2eta calculations (input) and including the
c                         curent r2eta calculation (output)
c
c        use param
c
c  real(knd) scalars and arrays
        real(knd) c,dcon,dconi,dec,eta,etas,pcoefn,pdcoefn,rm,rm2,
     1            r2dcoef1,r2est,r2test,sumdnpr,sumdnpi,sumdpi,
     2            sumdpr,sumnpi,sumnpr,ten,test,testd,testdm,testm,
     3            txi,txr,txdi,txdr,wm,xet,xets,x,pratb(maxp),
     4            pratt(maxp),pdratt(maxp)
c  complex(knd) scalars and arrays
        complex(knd) cc,denom,dnew,dnewd,dnewd1,dnewd2,dold,doldd1,
     1               doldd2,reld12,r1c,r1dc,r2c,r2dc,r2dcoef2,
     2               r2dtemp,r2temp,sr2temp,sr2dtemp,sumcoef,wronc,
     3               wronca,wroncb,wront
        complex(knd) enr(maxd),sneudr(maxlp),sneun(maxn),sneuf(maxn),
     1               sneudf(maxn)
c
c  integer arrays
        integer ineue(maxn)
c
        ten=10.0e0_knd
        dec=ten**(-ndec-2)
        rm=m
        c=abs(cc)
        dcon=ten**(-ndec)
        dconi=ten**(ndec+5)
        etas=eta*eta
        xet=sqrt(x*x+wm)
        xets=xet*xet
        naccrpl=ir1e+ir1de+int(log10(abs(r1c*r1dc)*c*(x*x+1.0e0_knd)))
        if(naccrpl.lt.0) naccrpl=0
        if(naccrpl.gt.ndec) naccrpl=ndec
        sumcoef=(ten**(-ir1de-ineue(l+1)-ipcoefn))/
     1          (cc*(x*x+1.0e0_knd)*r1dc*sneun(l+1)*pcoefn)
          if(naccrpl.gt.1) then
          sumcoef=sumcoef*(ten**(ir1de+ir1de))*cc*(x*x+1.0e0_knd)
          end if
        r2dcoef1=eta*wm/(xets*xet)
        r2dcoef2=cc*x/xet
        reld12=(r2dcoef2/r2dcoef1)*sneudr(l+1)*(pcoefn/
     1         pdcoefn)*ten**(ipcoefn-ipdcoefn)
        rm2=2*m
        lm2=(l-m)/2
        limb=4*int(abs(xet*real(cc))+abs(xet*aimag(cc)))
c
c  ix = 0 for l-m even; ix = 1 for l-m odd
        ix=l-m-2*lm2
        lim=limeta/2-ix
c
c  compute radial function of the second kind and its first derivative
c
c  backward series for denominator
        denom=(1.0e0_knd,0.0e0_knd)
        sumdpr=1.0e0_knd
        sumdpi=0.0e0_knd
        if (lm2.lt.1) go to 20
        dold=(1.0e0_knd,0.0e0_knd)
          do 10 j=lm2,1,-1
          jj=j+j+ix
          dnew=dold/(pratb(jj+1)*enr(j))
          denom=denom+dnew
          if(real(dnew).gt.0.0e0_knd) sumdpr=sumdpr+real(dnew)
          if(aimag(dnew).gt.0.0e0_knd) sumdpi=sumdpi+aimag(dnew)
          if(abs(dnew/denom).lt.dec) go to 20
          dold=dnew
10        continue
20      continue
c
c  forward series for denominator
        dold=1.0e0_knd
          do 30 j=lm2+1,lim
          jj=j+j+ix
          dnew=dold*enr(j)*pratb(jj+1)
          denom=denom+dnew
          if(real(dnew).gt.0.0e0_knd) sumdpr=sumdpr+real(dnew)
          if(aimag(dnew).gt.0.0e0_knd) sumdpi=sumdpi+aimag(dnew)
          if(abs(dnew/denom).lt.dec) go to 40
          dold=dnew
30        continue
40      continue
        jden=j
        numsubr=0
        if(sumdpr.ne.0.0e0_knd) numsubr=int(log10(abs(sumdpr/
     1                                  real(denom))))
        if(numsubr.lt.0) numsubr=0
        if(numsubr.gt.ndec) numsubr=ndec
        numsubi=0
        if(sumdpi.ne.0.0e0_knd) numsubi=int(log10(abs(sumdpi/
     1                                  aimag(denom))))
        if(numsubi.lt.0) numsubi=0
        if(numsubi.gt.ndec) numsubi=ndec
        numsub=max(numsubi,numsubr)
        naccd=ndec-max(2,int(log10(abs(cc))))-numsub
        if(naccd.lt.0) naccd=0
        r2est=abs(sumcoef*denom)
        r2test=r2est*dconi
c
c  backward series for numerator
        dold=(1.0e0_knd,0.0e0_knd)
        doldd1=(1.0e0_knd,0.0e0_knd)
        doldd2=reld12
        r2temp=(1.0e0_knd,0.0e0_knd)
        sumnpr=1.0e0_knd
        sumnpi=0.0e0_knd
        r2dtemp=doldd2
        sumdnpr=0.0e0_knd
        sumdnpi=0.0e0_knd
        if(real(doldd2).gt.0.0e0_knd) sumdnpr=real(doldd2)
        if(aimag(doldd2).gt.0.0e0_knd) sumdnpi=aimag(doldd2)
        if(l.ne.0) r2dtemp=r2dtemp+(1.0e0_knd,0.0e0_knd)
        if(l.ne.0) sumdnpr=sumdnpr+1.0e0_knd
        if(lm2.lt.1) go to 60
          do 50 j=lm2,1,-1
          jj=j+j+ix
          dnew=-dold/(sneuf(jj+m)*pratt(jj+1)*enr(j))
          dnewd1=-doldd1/(sneuf(jj+m)*pdratt(jj+1)*enr(j))
          dnewd2=-doldd2/(sneudf(jj+m)*pratt(jj+1)*enr(j))
          r2temp=r2temp+dnew
          dnewd=dnewd1+dnewd2
          r2dtemp=r2dtemp+dnewd
          if(real(dnew).gt.0.0e0_knd) sumnpr=sumnpr+real(dnew)
          if(aimag(dnew).gt.0.0e0_knd) sumnpi=sumnpi+aimag(dnew)
          if(real(dnewd).gt.0.0e0_knd) sumdnpr=sumdnpr+real(dnewd)
          if(aimag(dnewd).gt.0.0e0_knd) sumdnpi=sumdnpi+aimag(dnewd)
          if(abs(dnew/r2temp)+abs(dnewd/r2dtemp).lt.dec) go to 60
          dold=dnew
          doldd1=dnewd1
          doldd2=dnewd2
50      continue
60      continue
        if(m.eq.0.and.jj.eq.2) r2dtemp=r2dtemp-dnewd1
c
c  forward series for numerator
        dold=(1.0e0_knd,0.0e0_knd)
        doldd1=(1.0e0_knd,0.0e0_knd)
        doldd2=reld12
        test=1.0e0_knd
        testd=1.0e0_knd
        testm=1.0e0_knd
        testdm=1.0e0_knd
        js=lm2
        jds=lm2
        sr2temp=r2temp
        sr2dtemp=r2dtemp
        txr=sumnpr
        txi=sumnpi
        txdr=sumdnpr
        txdi=sumdnpi
        dnewsum=(0.0e0_knd,0.0e0_knd)
        dnewdsum=(0.0e0_knd,0.0e0_knd)
        doldd=reld12
        if(l.ne.0) doldd=reld12+(1.0e0_knd,0.0e0_knd)
          do 110 j=lm2+1,lim-1
          jj=j+j+ix
          dnew=-dold*enr(j)*sneuf(jj+m)*pratt(jj+1)
          dnewd1=-doldd1*enr(j)*sneuf(jj+m)*pdratt(jj+1)
          dnewd2=-doldd2*enr(j)*sneudf(jj+m)*pratt(jj+1)
          r2temp=r2temp+dnew
          if(real(dnew).gt.0.0e0_knd) sumnpr=sumnpr+real(dnew)
          if(aimag(dnew).gt.0.0e0_knd) sumnpi=sumnpi+aimag(dnew)
          if(abs(dnew).ne.0.0e0_knd) test=abs(dnew/r2temp)
          if(test.ge.testm) go to 80
          testm=test
          sr2temp=r2temp
          js=j
          txr=sumnpr
          txi=sumnpi
80        dnewd=dnewd1+dnewd2
          r2dtemp=r2dtemp+dnewd
          if(real(dnewd).gt.0.0e0_knd) sumdnpr=sumdnpr+real(dnewd)
          if(aimag(dnewd).gt.0.0e0_knd) sumdnpi=sumdnpi+aimag(dnewd)
          if(abs(dnewd).ne.0.0e0_knd) testd=abs(dnewd/r2dtemp)
          if(testd.ge.testdm) go to 100
          testdm=testd
          sr2dtemp=r2dtemp
          jds=j
          txdr=sumdnpr
          txdi=sumdnpi
100       if ((test+testd).lt.dcon.and.jj+m.gt.limb) go to 130
          if(abs(r2temp)+abs(r2dtemp).gt.r2test) go to 120
          if((abs(dnew).eq.0.0e0_knd).and.(abs(dnewd).eq.0.0e0_knd))
     1            go to 120
          dold=dnew
          doldd1=dnewd1
          doldd2=dnewd2
110       continue
120     r2temp=sr2temp
        r2dtemp=sr2dtemp
130     jmax=j
        jeta=js
        if(jds.gt.js) jeta=jds
        jtestm=-int(log10(abs(testm)))
        if(jtestm.lt.0) jtestm=0
        if(jtestm.gt.ndec) jtestm=ndec
        jtestdm=-int(log10(abs(testdm)))
        if(jtestdm.lt.0) jtestdm=0
        if(jtestdm.gt.ndec) jtestdm=ndec
        naccns1r=0
        if(txr.ne.0.0e0_knd) naccns1r=int(log10(abs(txr/real(r2temp))))
        if(naccns1r.lt.0) naccns1r=0
        if(naccns1r.gt.ndec) naccns1r=ndec
        naccns1i=0
        if(txi.ne.0.0e0_knd) naccns1i=int(log10(abs(txi/
     1                                  aimag(r2temp))))
        if(naccns1i.lt.0) naccns1i=0
        if(naccns1i.gt.ndec) naccns1i=ndec
        naccns2r=0
        if(txdr.ne.0.0e0_knd) naccns2r=int(log10(abs(txdr/
     1                                  real(r2dtemp))))
        if(naccns2r.lt.0) naccns2r=0
        if(naccns2r.gt.ndec) naccns2r=ndec
        naccns2i=0
        if(txdi.ne.0.0e0_knd) naccns2i=int(log10(abs(txdi/
     1                                 aimag(r2dtemp))))
        if(naccns2i.lt.0) naccns2i=0
        if(naccns2i.gt.ndec) naccns2i=ndec
        naccns1=max(naccns1r,naccns1i)
        naccns2=max(naccns2r,naccns2i)
        naccn1=min(jtestm-2,ndec-2-naccns1)
        naccn2=min(jtestdm-2,ndec-2-naccns2)
        naccn=min(naccn1,naccn2)
        if(naccn.lt.0) naccn=0
        naccnmaxp=naccnmax
        if(naccn.gt.naccnmax) naccnmax=naccn
c
c  combining results to form the radial function characteristics
c  r2c and r2dc and corresponding exponents ir2e and ir2de
        r2c=r2temp*sneun(l+1)*pcoefn/denom
        iterm=0
        if(abs(r2c).ne.0.0e0_knd) iterm=int(log10(abs(r2c)))
        ir2e=ineue(l+1)+ipcoefn+iterm
        r2c=r2c*ten**(-iterm)
        r2dc=(r2dcoef1*r2dtemp*sneun(l+1)*pdcoefn/denom)*
     1       ten**(ineue(l+1)+ipdcoefn-ir2e)
        iterm=0
        if(abs(r2dc).ne.0.0e0_knd) iterm=int(log10(abs(r2dc)))
        ir2de=ir2e+iterm
 	r2dc=r2dc*ten**(-iterm)
c        write(40,140) jmax,jden,lim,js,jtestm,naccns1,jds,jtestdm,
c     1               naccns2,naccn,naccd
c140     format(8x,'r2eta: numerator, denominator converged in ',
c     1         i6,' ,',i6,' terms; ',i6,' terms available.',/,
c     2         15x,'best r2 at ',i6,' terms with convergence to ',i2,
c     3         ' digits; ',i2,' digits subtr. error.',/,15x,
c     4         'best r2d at ',i6,' terms with convergence to ',i2,
c     5         ' digits; ',i2,' digits subtr. error.',/,15x,
c     6         'estimated numerator and denominator accuracy is ',i2,
c     7         ' and ',i2,' digits.')
        wronca=r1c*r2dc*(ten**(ir1e+ir2de))
        wroncb=r2c*r1dc*(ten**(ir2e+ir1de))
        wronc=wronca-wroncb
        wront=(1.0e0_knd,0.0e0_knd)/(cc*(x*x+1.0e0_knd))
        nacceta=-int(log10(abs((wronc-wront)/wront)+dec))
        if(nacceta.gt.ndec) nacceta=ndec
        if(nacceta.lt.0) nacceta=0
        naccetaw=nacceta
        nacccor=-int(log10(abs((wronca-wroncb)/wronca)+dec))
        if(nacccor.lt.0) nacccor=0
        if(nacccor.gt.naccrpl) nacccor=naccrpl
        if(nacccor.gt.0) nacceta=min(nacceta+nacccor,ndec-numsub)
        if(nacceta.lt.0) nacceta=0
        if(nacceta.lt.naccetaw) nacceta=naccetaw
160     if(abs(r2c).ge.1.0e0_knd) go to 170
        r2c=r2c*ten
        ir2e=ir2e-1
170     continue
        if(abs(r2dc).ge.1.0e0_knd) go to 180
        r2dc=r2dc*ten
        ir2de=ir2de-1
180     continue
        if(naccn.gt.nacceta) naccnmax=max(naccnmaxp,nacceta)
        return
        end
c
c
c
        subroutine cmtql1(n,maxe,d,e)
c
c  purpose:     to compute the eigenvalues of a complex symmetric
c               tridiagonal matrix
c
c  subroutine published by Jane Cullum and Ralph Willoughby in 'Lanczos
c  algorithsm for Large symmetric eigenvalue computations,' 2002;
c  I have removed the parameter ierr from the call statement and set
c  it initially equal to zero. I also relabeled the subdiagonal elements
c  e prior to input rather than after input.
c
c  parameters:
c
c     input:    n:    order of the matrix
c               maxe: dimension of both vectors d and e
c               d:    diagonal elements of the tridiagonal matrix
c               e:    subdiagonal elements of the matrix,
c                     beginning with e(1); e(n) is set equal to zero
c
c     output:   d  : eigenvalues
c

c        use param
c
        real(knd) machep,eps,temp,t0,t1,zero,half,one,two
        complex(knd) d(maxe),e(maxe),b,c,f,g,p,r,s,w,czero,cone
c
        ndec=precision(temp)
        machep=10.0e0_knd**(-ndec)
        eps=100.0e0_knd*machep
        zero=0.0e0_knd
        half=0.5e0_knd
        one=1.0e0_knd
        two=2.0e0_knd
        czero=cmplx(zero,zero,knd)
        cone=cmplx(one,zero,knd)
        ierr=0
        e(n)=czero
          do 140 l=1,n
          j=0
20          do 30 m=l,n
            if(m.eq.n) go to 40
            temp=abs(d(m))+abs(d(m+1))
            if(abs(e(m)).le.temp*machep) go to 40
30          continue
40        p=d(l)
          if(m.eq.l) go to 100
          if(j.eq.100) go to 150
          j=j+1
          g=(d(l+1)-p)*half
          t0=abs(g)
          t1=abs(e(l))
          if(t0.gt.t1) go to 50
          w=g/e(l)
          r=sqrt(cone+w*w)
          t0=abs(w+r)
          t1=abs(w-r)
          temp=one
          if(t1.gt.t0) temp=-one
          g=d(m)-p+e(l)/(w+temp*r)
          go to 60
50        continue
          w=e(l)/g
          r=sqrt(cone+w*w)
          t0=sqrt(cone+r)
          t1=sqrt(cone-r)
          temp=one
          if(t1.gt.t0) temp=-one
          g=d(m)-p+w*e(l)/(cone+temp*r)
60        continue
          s=cone
          c=-cone
          p=czero
          mml=m-l
            do 90 i1=1,mml
            i=m-i1
            f=s*e(i)
            b=-c*e(i)
            t0=abs(g)
            t1=abs(f)
            if(t1.gt.t0) go to 70
            w=f/g
            r=sqrt(cone+w*w)
            e(i+1)=g*r
            c=cone/r
            s=w*c
            go to 80
70          continue
            w=g/f
            r=sqrt(cone+w*w)
            e(i+1)=f*r
            s=cone/r
            c=w*s
80          continue
            temp=abs(w)**2+one
            t0=sqrt(temp)
            t1=abs(r)
            ierr=-l
            if(t1.le.eps*t0) go to 160
            ierr=0
            g=d(i+1)-p
            r=(d(i)-g)*s+two*c*b
            p=s*r
            d(i+1)=g+p
            g=b-c*r
90          continue
          d(l)=d(l)-p
          e(l)=g
          e(m)=zero
          go to 20
100       if(l.eq.1) go to 120
            do 110 i1=2,l
            i=l+2-i1
            if(abs(p).ge.abs(d(i-1))) go to 130
            d(i)=d(i-1)
110         continue
120       i=1
130       d(i)=p
140       continue
        go to 160
150     ierr=l
160     return
        end
c
c
        subroutine eigorder(c,n,np,np2,f,g,m,lnum,limeig,eigst,lipl,
     1                      liplp,lips)
c
c  purpose:     to order the starting values of the eigenvalues so that
c               the magnitudes of the resulting radial functions change
c               smoothly with increasing order.
c
c  parameters:
c
c     input:
c               c     : complex size parameter
c               n     : number of even starting values and the
c                       number of odd starting values obtained
c                       from the tridiagonal matrices
c               np    : dimension of f and g
c               np2   : dimension of eigst
c               f     : vector of unordered even starting values
c               g     : vector of unordered odd starting values
c               m     : value of m
c               lnum  : number of values of l-m for which spheroidal
c                       functions are desired, dimension of eigst
c               limeig: number of eigenvalue estimates that are desired
c
c     output:   eigst : array of ordered eigenvalue estimates,
c                       either lnum values or 2*imax values,
c                       whichever is smaller
c               lipl  : maximum value of l-m+1 for which eigenvalue
c                       estimates are either negative or paired
c               liplp : value of l-m+1 following prolate-like
c                       eigenvalues
c               lips  : value of l-m+1 for the first prolate-like
c                       eigenvalue or the first non-paired eigenvalue
c                       if there are no prolate-like eigenvalues
c
c        use param
        real(knd) cm,testpr,testpi
        complex(knd) c,cp,eigst(np2),f(np),g(np),fm(np),gm(np),fp(np),
     1               gp(np),fs(np),gs(np),p,testp,temp
c
        imax=limeig/2
        if(2*(limeig/2).ne.limeig) imax=imax+1
        iupp=min(imax,lnum/2+1)
        cm2=abs(c)**2
        lipl=0
c
c  order even eigenvalues in ascending real part and retain imax of them
          do i = 1,n-1
          k = i
          p=f(i)
            do j = i + 1, n
              if(real(f(j)).lt.real(p)) then
              k = j
              p = f(j)
              end if
            end do
            if(k.ne.i) then
            f(k)=f(i)
            f(i)=p
            end if
          end do
c
c  order odd eigenvalues in ascending real part and retain imax of them
          do i = 1,n-1
          k = i
          p=g(i)
            do j = i + 1, n
              if(real(g(j)).lt.real(p)) then
              k = j
              p = g(j)
              end if
            end do
            if(k.ne.i) then
            g(k)=g(i)
            g(i)=p
            end if
          end do
c
c  determine eigenvalues with negative real parts
        limp=0
        limpl=0
          do i = 1,imax
            if(real(f(i)).lt.0.0e0_knd.and.real(g(i)).lt.0.0e0_knd)
     1          limpl=limpl+1
            if(real(f(i)).lt.0.0e0_knd.or.real(g(i)).lt.0.0e0_knd) then
            limp=limp+1
            else
            go to 10
            end if
          end do
10      continue
c
c  Identify any paired eigenvalues other than those
c  with negative real parts
c
        ilimit=limp+1
        jci=aimag(c)
        iflag=0
          do i=ilimit,imax
          jmin=max(i-jci,limp+1)
          limpsav=limp
            do j=jmin,imax
            if(abs(real(f(i))-real(g(j)))/
     1           abs(real(f(i))).gt.0.001e0_knd) go to 20
            if(aimag(c).ne.0.0e0_knd.and.abs(aimag(f(i))-aimag(g(j)))/
     1           abs(aimag(f(i))).gt.0.001e0_knd) go to 20
            limp=limp+1
            limpl=limpl+1
            temp=f(i)
              do k=i,limp+1,-1
              f(k)=f(k-1)
              end do
            f(limp)=temp
            temp=g(j)
              do k=j,limp+1,-1
              g(k)=g(k-1)
              end do
            g(limp)=temp
            go to 30
20          continue
            end do
          if(limp.eq.limpsav.and.iflag.eq.1.and.i-is.gt.jci) go to 40
            if(limp.eq.limpsav.and.iflag.eq.0) then
            iflag=1
            is=i
            end if
30        continue
          end do
40      continue
        lipl=limp+limp
        lips=limp+limp+1
c
c  locate prolate-like eigenvalues
        nmax=int(abs(aimag(c)))/2
        nf=0
        ng=0
        if(nmax.eq.0) go to 130
        cp=(0.0e0_knd,-1.0e0_knd)*c
          do j=1,nmax
          ka=j-1
          n=2*(j-1)+1
          testp=n*cp+m*m-(n*n+5)/8.0e0_knd-
     1           n*(n*n+11-32*m*m)/(64.0e0_knd*cp)
          testpr=real(testp)
          testpi=aimag(testp)
          if(2*(ka/2).ne.ka) go to 70
            do i=limpl+1,imax
            if(abs((real(f(i))-testpr)/real(f(i))).gt.0.02e0_knd)
     1         go to 60
            if(aimag(c).ne.0.0e0_knd.and.abs((aimag(f(i))-testpi)/
     1         aimag(f(i))).gt.0.02e0_knd) go to 60
            nf=nf+1
            fp(nf)=f(i)
            isav=i
              do k=i,imax-nf
              f(k)=f(k+1)
              end do
            go to 90
60          continue
            end do
          go to 100
70          do i=limpl+1,imax
            if(abs((real(g(i))-testpr)/real(g(i))).gt.0.02e0_knd)
     1         go to 80
            if(aimag(c).ne.0.0e0_knd.and.abs((aimag(g(i))-testpi)/
     1         aimag(g(i))).gt.0.02e0_knd) go to 80
            ng=ng+1
            gp(ng)=g(i)
              do k=i,imax-ng
              g(k)=g(k+1)
              end do
            go to 90
80          continue
            end do
          go to 100
90        continue
          end do
100      continue
c
c  place prolate-like eigenvalues in series after last
c  paired eigenvalue or last eigenvalue with negative real
c  part if pairing ends before this
         if(nf.eq.0) go to 120
           i=limp+1
             do j=imax-nf,i,-1
             f(j+nf)=f(j)
             end do
             do j=1,nf
             f(i+j-1)=fp(j)
             end do
             do j=imax-ng,i,-1
             g(j+ng)=g(j)
             end do
             do j=1,ng
             g(i+j-1)=gp(j)
             end do
             lips=i+i-1
120      continue
130        do i=1,iupp
           eigst(i+i-1)=f(i)
             if(i+i.le.limeig) then
             eigst(i+i)=g(i)
             end if
           end do
         liplp=lips+nf+ng
         return
         end
c
c
        subroutine conver (l,m,lnum,cc,limd,blist,glist,blist1,glist1,
     1                     ndec,maxd,ioprad,minacc,eigest,eigmat,lipl,
     2                     lips,liplp,match,eigp,eign,kindd,kindq,
     3                     eigval,enr,ienr,itestm,naccre,ieigt,iopeig,
     4                     iflag,kflag)
c
c  purpose:     To determine a converged eigenvalue using the
c               boukwamp method.
c  parameters:
c
c     input:    l     : l
c               m     : m
c               lnum  : number of l values desired
c               cc    : complex c
c               limd  : number of enr values computed
c               blist : array of coefficients used in recursion relation
c                       (for knd arithmetic)
c               glist : array of coefficients used in recursion relation
c                       (for knd arithmetic)
c               blist1: array of coefficients used in recursion relation
c                       (for knd1 arithmetic)
c               glist1: array of coefficients used in recursion relation
c                       (for knd1 arithmetic)
c               ndec  : number of decimal digits available in real
c                       arithmetic
c               maxd  : dimension of enr,blist,glist arrays
c               ioprad: integer input equal to 0 if no radial functions
c                       are desired, equal to 1 if only radial functions
c                       of the first kind are desired, or equal to 2
c                       if radial functions of both the first and second
c                       kinds are desired
c               minacc: minacc desired accuracy
c               eigest: estimate of eigenvalue using extrapolation of
c                       the four previous eigenvalues of the same parity
c                       for l-m+1 greater than 8 and less than lipl+1.
c               eigmat: estimate of the eigenvalue obtained from matrix
c                       output values for l-m+1 up to max(67,4*nbp/3).
c                       for higher values of l-m+1, this estimate is
c                       given by extrapolation of the four previous
c                       eigenvalues.
c               lipl  : maximum value of l-m+1 for which eigenvalue
c                       estimates are either negative or paired
c               lips  : value of l-m+1 for first prolate-like
c                       eigenvalue, if there are any
c               liplp : value of l-m+1 immediately following the
c                       prolate-like eigenvalues
c               match : number of leading decimal digits of agreement
c                       between the matrix estimate for this eigenvalue
c                       and the matrix estimate for the eigenvalue of
c                       the neighboring l, either l + 1 if l - m is even
c                       or l - 1 if l - m is odd. Set equal to zero if
c                       l - m + 1 is greater than lipl
c               eigp  : previous eigenvalue of same parity
c               eign  : estimate of next eigenvalue of same parity
c               kindd : number of bytes for real data in double
c                       precision
c               kindq : number of bytes for real data in quadruple
c                       precision
c
c     output:   eigval: converged eigenvalue
c               enr   : array of scaled ratios of successive d
c                       coefficients
c               ienr  : index n of last d coefficient ratio used
c                       in computing first term in the denominator
c                       of the eigenvalue correction
c               itestm: number of matching digits for the forward and
c                       backward recursion for d coefficient ratios
c               naccre: estimated accuracy of d coefficients based on
c                       degree of convergence of eigenvalue or
c                       on the value of match (see above for the
c                       definition of match)
c               ieigt : integer vector used in a test in subroutine
c                       main to make sure that none of the lnum
c                       eigenvalues for a given value of m is a
c                       duplicate of another one of the same parity.
c                       If this happens, an error message is written
c                       to file 60 giving the values of l where the
c                       duplication occurred. Such duplication shows
c                       that an eigenvalue is missing and the function
c                       values for that m should not be used.
c
c     input/
c     output    iopeig: equal to unity for values of l-m+1 up to lipl
c                       when eigest is more accurate than eigmat;
c                       otherwise equal to zero. separate value for
c                       l-m even and for l-m odd.
c               iflag : equal to zero when the converged eigenvalue
c                       is used for the eigenvalue; equal to unity if
c                       the estimate for the eigenvalue is used for
c                       the eigenvalue
c               kflag : flag used to control the use of knd = 16
c                       arithmetic when running coblfcn in knd = 8
c                       mode with option to use knd = 16 arithmetic
c                       in the Bouwkamp procedure.
c                       equal to one when Bouwkamp procedure is run
c                       in knd = 16 arithmetic.
c                       equal to zero or two when it is run in knd = 8
c                       arithmetic.
c
c        use param
c
c  real(knd) scalars
        real(knd) c,dec,eigdec,ten
        real(knd1) dec1,eigdec1,eigdec2
c
c  complex(knd) scalars and arrays
        complex(knd) cc,cora,corb,de,dl,eign,eigp,eigstart,eigs,
     1               eigest,eigmat,eigval,enrc,blist(maxd),enr(maxd),
     2               enrf(maxd),glist(maxd)
        complex(knd1) cora1,corb1,de1,dl1,eigval1,enrc1,blist1(maxd),
     1                enr1(maxd),glist1(maxd)
c  integer array
        integer ieigt(lnum)
c
        c=abs(cc)
        ten=10.0e0_knd
        nbp=int(2*(real(cc)+aimag(cc))/3.1416)
        dec=ten**(-ndec-1)
        eigdec=dec*100.0e0_knd
        eigstart=eigmat
        if(l-m+1.gt.lipl) iopeig=0
        if(iopeig.eq.1) eigstart=eigest
        eigval=eigstart
        ndec1=precision(eigval1)
        dec1=10.0e0_knd1**(-ndec1-1)
        eigdec1=dec1*100.0e0_knd1
        eigdec2=dec*100.0e0_knd
        eigval1=eigval
        lm2=(l-m)/2
        lmtest=max(50,nbp)
        limdb=2*ienr+20
        if(l.eq.m.or.l.eq.m+1) limdb=2*ienr
        if(limdb.gt.limd) limdb=limd
c
c  begin bouwkamp procedure
        iflag=0
        ix=l-m-2*lm2
        ifc=1
        lim2=limdb/2-ix
        iglim=lim2+1
        irio=lm2+1
        iw1=lm2+2
        itry=1
        if(kflag.eq.1) go to 80
40      enr(1)=eigval-glist(1)
        if(lm2.lt.1) go to 45
c
c  evaluate the continued fraction
          do i=1,lm2
          enr(i+1)=-blist(i)/enr(i)-glist(i+1)+eigval
          end do
45      enr(lim2)=-blist(lim2)/(glist(iglim)-eigval)
        iw15=lim2-1
        ip=iw1+iw15
        if(iw15.lt.iw1) go to 50
c
c  evaluate the continued fraction
          do i=iw1,iw15
          ipi=ip-i
          enr(ipi)=-blist(ipi)/(glist(ipi+1)-eigval+enr(ipi+1))
          end do
50      enrc=-blist(irio)/(glist(irio+1)-eigval+enr(irio+1))
        de=enrc*enrc/blist(irio)
        corb=de
        if(lim2.lt.iw1) go to 55
c
c  compute first sum in the denominator of the correction
          do i=iw1,lim2
          de=enr(i)*enr(i)/blist(i)*de
          corb=corb+de
          if(abs(de/corb).lt.dec) go to 55
          end do
55      ienr=i
        if(ienr.lt.lim2-10.and.l.gt.m+1) ienr=lim2-12
        de=1.0e0_knd
        cora=de
        if(lm2.eq.0) go to 60
c
c  compute second term in the denominator of the correction
          do i=1,lm2
          de=blist(irio-i)/(enr(irio-i)*enr(irio-i))*de
          cora=cora+de
          if(abs(de/cora).lt.dec) go to 60
          end do
c
c  compute the correction to the eigenvalue
60      dl=(enrc-enr(irio)+dec)/(cora+corb+dec)
        eigval=dl+eigval
c
c  eigenvalue accurate enough?
        if(abs(dl/eigval).lt.eigdec) go to 70
        if(abs((enrc-enr(irio))/enrc).lt.eigdec) go to 70
        ifc=ifc+1
        if(ifc.le.20) go to 40
70      continue
        int1=-int(log10(abs(dl/eigval)+dec))
        int5=-int(log10(abs((cora+corb)*eigval/enrc)+dec))
        int2=-int(log10(abs((eigval-eigstart)/eigval)+dec))
        if(kflag.eq.0.and.int5.gt.0.and.int2.gt.2) kflag=1
          if(kflag.eq.1) then
          eigval1=eigstart
          ifc=1
          go to 80
          else
          go to 120
          end if
75      continue
80      enr1(1)=eigval1-glist1(1)
        if(lm2.lt.1) go to 85
c
c  evaluate the continued fraction
          do i=1,lm2
          enr1(i+1)=-blist1(i)/enr1(i)-glist1(i+1)+eigval1
          end do
85      enr1(lim2)=-blist1(lim2)/(glist1(iglim)-eigval1)
        iw15=lim2-1
        ip=iw1+iw15
        if(iw15.lt.iw1) go to 90
c
c  evaluate the continued fraction
          do i=iw1,iw15
          ipi=ip-i
          enr1(ipi)=-blist1(ipi)/(glist1(ipi+1)-eigval1+enr1(ipi+1))
          end do
90      enrc1=-blist1(irio)/(glist1(irio+1)-eigval1+enr1(irio+1))
        de1=enrc1*enrc1/blist1(irio)
        corb1=de1
        if(lim2.lt.iw1) go to 95
c
c  compute first sum in the denominator of the correction
          do i=iw1,lim2
          de1=enr1(i)*enr1(i)/blist1(i)*de1
          corb1=corb1+de1
          if(abs(de1/corb1).lt.dec1) go to 95
          end do
95      ienr=i
        ienrs=ienr
        if(ienr.lt.lim2-10.and.l.gt.m+1) ienr=lim2-12
        de1=1.0e0_knd1
        cora1=de1
        if(lm2.eq.0) go to 100
c
c  compute second term in the denominator of the correction
          do i=1,lm2
          de1=blist1(irio-i)/(enr1(irio-i)*enr1(irio-i))*de1
          cora1=cora1+de1
          if(abs(de1/cora1).lt.dec1) go to 100
          end do
c
c  compute the correction to the eigenvalue
100     dl1=(enrc1-enr1(irio)+dec1)/(cora1+corb1+dec1)
        eigval1=dl1+eigval1
c
c  eigenvalue accurate enough?
        if(abs(dl1/eigval1).lt.eigdec2) go to 110
        if(abs((enrc1-enr1(irio))/enrc1).lt.eigdec1) go to 110
        ifc=ifc+1
        if(ifc.le.20) go to 80
110     continue
        int1=-int(log10(abs(dl1/eigval1)+dec1))
        eigval=eigval1
        int5=-int(log10(abs((cora1+corb1)*eigval1/enrc1)+dec1))
        if(kflag.eq.1.and.l.gt.nbp.and.l.gt.liplp.and.int5.lt.1) kflag=2
120     continue
        iflag=0
        if(int1.gt.ndec) int1=ndec
        int2=-int(log10(abs((eigval-eigstart)/eigval)+dec))
        if(int2.gt.ndec) int2=ndec
        if(int2.lt.3.and.(l-m+1.lt.lips.or.l-m+1.ge.liplp)) iflag=1
        if(int2.lt.1.and.l-m+1.ge.lips.and.l-m+1.lt.liplp) iflag=1
          if(iflag.eq.1.and.int1.gt.5) then
          if(real(eigval).gt.real(eigp).and.real(eigval).lt.real(eign)
     1       .and.int1.gt.5) iflag=0
          end if
            if(iflag.eq.1.and.itry.eq.2) then
            icheck=-int(log10(abs((eigval-eigs)/eigval)+dec))
            if(icheck.gt.5) iflag=0
            end if
            if(iflag.eq.1.and.itry.eq.2) then
            icheck1=-int(log10(abs((eigs-eigstart)/eigs)+dec))
            icheck2=-int(log10(abs((eigval-eigstart)/eigval)+dec))
              if(icheck1.gt.icheck2) then
              eigval=eigs
              int1=int1s
              int2=int2s
              ifc=ifcs
              if(int1s.gt.5) iflag=0
              end if
            end if
          if(iflag.eq.1.and.itry.eq.1.and.(int1.gt.5.or.ifc.eq.21).and.
     1          match.lt.6) then
          eigs=eigval
          eigval=eigp+0.5e0_knd*(eign-eigp)
          ifcs=ifc
          ifc=1
          itry=2
          int1s=int1
          int2s=int2
          if(kflag.ne.1) go to 40
            if(kflag.eq.1) then
            eigval1=eigval
            go to 80
            end if
          end if
        itry=1
        if(int2.lt.match-2.and.iopeig.eq.0) iflag=1
        if(int2.gt.int1) iflag=1
        iopeigs=iopeig
        iopeig=0
          if(l-m+1.lt.lipl-1.and.match.gt.3) then
            if(iopeig.eq.0) then
            int4=-int(log10(abs((eigval-eigest)/eigval)+
     1           dec))
            int3=int2
            end if
            if(iopeig.eq.1) then
            int3=-int(log10(abs((eigval-eigmat)/eigval)+
     1           dec))
            int4=int2
            end if
          if(int4.gt.max(4,int3)) iopeig=1
          if(int3.ge.int4) iopeig=0
          end if
        if(iflag.eq.1) eigval=eigstart
140     continue
c        if(knd.eq.kindd.and.ioprad.ne.0.and.iflag.eq.0) write(40,150) l,
c     1           eigval,eigstart
c        if(knd.eq.kindd.and.ioprad.eq.0.and.iflag.eq.0) write(50,150) l,
c     1           eigval,eigstart
c        if(knd.eq.kindq.and.ioprad.ne.0.and.iflag.eq.0) write(40,155) l,
c     1           eigval,eigstart
c        if(knd.eq.kindq.and.ioprad.eq.0.and.iflag.eq.0) write(50,155) l,
c     1           eigval,eigstart
c150     format(1x,'l =',i5,6x,'eigenvalue =',e24.15,e24.15,/,16x,
c     1               ' estimate =',e24.15,e24.15)
c155     format(1x,'l =',i5,6x,'eigenvalue =',e39.31,e39.31,/,16x,
c     1               ' estimate =',e39.31,e39.31)
c        if(knd.eq.kindd.and.ioprad.ne.0.and.iflag.eq.1) write(40,160) l,
c     1          eigstart
c        if(knd.eq.kindd.and.ioprad.eq.0.and.iflag.eq.1) write(50,160) l,
c     1          eigstart
c        if(knd.eq.kindq.and.ioprad.ne.0.and.iflag.eq.1) write(40,165) l,
c     1          eigstart
c        if(knd.eq.kindq.and.ioprad.eq.0.and.iflag.eq.1) write(50,165) l,
c     1          eigstart
c160     format(1x,'l =',i5,6x,'eigenvalue =',e24.15,e24.15' obtained'
c     1               ' from tridiagonal matrix')
c165     format(1x,'l =',i5,6x,'eigenvalue =',e39.31,e39.31,/,30x,
c     1               ' obtained from tridiagonal matrix')
        if(iflag.eq.0) naccre=min(int1,ndec)
        if(iflag.eq.1) naccre=max(int2,match)
        if(naccre.gt.ndec) naccre=ndec
        if(iflag.eq.0) ieigt(l-m+1)=max(5,int1-2)
        if(iflag.eq.1) ieigt(l-m+1)=max(5,match-1)
        ieigt(l-m+1)=min(8,ieigt(l-m+1))
c
c  calculate the d coefficient ratios (enr)
        lim2=limd/2-ix
        iesub=0
          if(lm2.eq.0.and.real(cc).gt.50.0e0_knd) then
          iesub=-int(log10(abs((eigval-glist(1))/eigval)+dec))
          end if
        if(real(cc).gt.50.0e0_knd.and.iesub.lt.4) go to 200
        enr(lim2)=-blist(lim2)/(glist(lim2+1)-eigval)
          do i=lim2-1,lm2+1,-1
          enr(i)=-blist(i)/(glist(i+1)-eigval+enr(i+1))
          end do
        if(lm2.eq.0) go to 190
        enr(1)=eigval-glist(1)
        if(lm2.eq.1) go to 190
          do n=2,lm2
          enr(n)=-blist(n-1)/enr(n-1)-glist(n)+eigval
          end do
190     continue
        itestm=ndec
c        if(ioprad.ne.0) write(40,195) lim2,naccre
c        if(ioprad.eq.0) write(50,195) lim2,naccre
c195     format(7x,i7,' d coefs. estimated eigenvalue accuracy is ',
c     1         i3,' digits')
        go to 260
200     enr(lim2)=-blist(lim2)/(glist(lim2+1)-eigval)
          do 210 i=lim2-1,1,-1
          enr(i)=-blist(i)/(glist(i+1)-eigval+enr(i+1))
210       continue
        enrf(1)=eigval-glist(1)
        nlim=1
        itestm=-int(log10(abs((enrf(1)-enr(1))/enrf(1))+
     1           dec))
        if(itestm.lt.0) itestm=0
        if(itestm.ge.ndec) go to 230
          do 220 n=2,lim2
          enrf(n)=-blist(n-1)/enrf(n-1)-glist(n)+eigval
          itest=-int(log10(abs((enrf(n)-enr(n))/enrf(n))+
     1           dec))
          if(itest.lt.0) itest=0
          if(itest.lt.itestm-4) go to 230
          if(itest.le.itestm) go to 220
          itestm=itest
          nlim=n
          if(itestm.ge.ndec) go to 230
220       continue
230     nlimp=2*(nlim-1)+ix
c        if(ioprad.ne.0) write(40,240) lim2,itestm,nlimp,naccre
c        if(ioprad.eq.0) write(50,240) lim2,itestm,nlimp,naccre
c240     format(7x,i7,' d coefs. Forward and backward recursion d ',
c     1         'ratios match to ',i2,' digits at n = ',i6,/,24x,
c     2         'estimated eigenvalue accuracy is ',i3,' digits')
          do 250 n=1,nlim
          enr(n)=enrf(n)
250       continue
260     continue
        return
        end
c
c
        subroutine dnorm (l,m,cc,ndec,nex,limd,maxd,enr,ioprad,iopang,
     1                    dc01,idc01,dfnorm,idfe,dmlf,idmlfe,dmfnorm,
     2                    idmfe,dmlmf,idmlmfe,dmsnorm,idmse,dmlms,
     3                    idmlmse,jnorm,jsubf,jsubmf,jsubms)
c
c  purpose:     To compute d coefficient ratios from n values and to
c               calculate the normalization of the d coefficients.
c
c  parameters:
c
c     input:    l       : l
c               m       : m
c               cc      : complex c
c               ndec    : number of decimal digits available in
c                         real arithmetic
c               nex     : maximum exponent in real arithmetic
c               limd    : approximately twice the maximum number
c                         of terms available to be taken in the sum
c               maxd    : dimension of enr array
c               enr     : array of ratios of scaled d coefficients
c               ioprad  : set equal to zero if no radial functions are
c                         desired, set equal to two otherwise
c               iopang  : set equal to zero if no angular functions
c                         are desired; set equal to 1 if angular
c                         functions are desired; set equal to 2 if
c                         their first derivatives are also desired
c
c     output:   enr     : array of ratios of d coefficients.
c                         enr(i) = ratio of the d coefficient with
c                         subscript 2*i+ix to the d coefficient with
c                         subscript 2*(i-1)+ix. Here ix =0 when l-m is
c                         even and ix=1 when l-m is odd.
c                         If the user needs the d coefficent ratios,
c                         they are available below right before
c                         statement 20.
c               dc01    : characteristic of ratio of first d
c                         coefficient (either d0 or d1, depending on
c                         whether l-m is even or odd) to the d
c                         coefficient of order equal to l-m
c               idc01   : exponent associated with dc01
c               dfnorm  : characteristic of Flammer normalization sum of
c                         d coefficients. equal to the reciprocal of
c                         the value of the d coefficient d(n = l - m)
c                         using this normalization for the angular
c                         functions
c               idfe    : exponent associated with dfnorm
c               dmlf    : characteristic of the d coefficient with index
c                         l-m in the Flammer normalization
c               idmlfe  : exponent associated with dmlf
c               dmfnorm : characteristic of Morse-Feshbach normalization
c                         sum of the d coefficients. equal to the
c                         reciprocal of the value of the d coefficient
c                         d(n = l - m) using this normalization for the
c                         angular functions
c               idmfe   : exponent associated with dmfnorm
c               dmlmf   : characteristic of the d coefficient with index
c                         l-m in the Morse-Feshbach normalization
c               idmlmfe : exponent associated with dmlmf
c               dmsnorm : characteristic of Meixner-Schafke normalization
c                         sum of the d coefficients. equal to the
c                         reciprocal of the value of the d coefficient
c                         d(n = l - m) using this normalization for the
c                         angular functions
c               idmse   : exponent associated with dmsnorm
c               dmlms   : characteristic of the d coefficient with index
c                         l-m in the Meixner-Schafke normalization
c               idmlmse : exponent associated with dmlms
c               jnorm   : maximum index of enr required for convergence
c                         of dfnorm and dmfnorm
c               jsubf   : effective number of decimal digits of subtraction
c                         error incurred in calculating dfnorm
c               jsubmf  : effective number of decimal digits of subtraction
c                         error incurred in calculating dmfnorm
c               jsubms  : effective number of decimal digits of subtraction
c                         error incurred in calculating dmsnorm
c
c        use param
c
c  real(knd) scalars
        real(knd) aj,arr,c,dec,ea,rm2,rm2m1,rm2m3,rm2p1,sgn,sumpr,
     1            sumpi,ten,teste,testeo
c
c  complex(knd) scalars and array
        complex(knd) cc,coef,csq,dfnorm,dmlf,dmfnorm,dmlmf,dmsnorm,
     1               dmlms,dc01,enr(maxd),term
c
        c=abs(cc)
        ten=10.0e0_knd
        rm2=m+m
        rm2m1=m+m-1
        rm2p1=m+m+1
        rm2m3=m+m-3
        dec=ten**(-ndec-1)
        nfac=nex/3
        if(2*(nfac/2).ne.nfac) nfac=nfac-1
        teste=ten**nfac
        testeo=1.0e0_knd/teste
        ir1tempe=0
        nbp=int(2*(real(cc)+aimag(cc))/3.1416)
        csq=cc*cc
        lm2=(l-m)/2
        ix=l-m-2*lm2
        mml=m+m-1+ix
        lim2=limd/2-ix
        sgn=1.0e0_knd
        do 20 i=1,lim2
          arr=ix+i+i
          ea=arr+arr+rm2
          enr(i)=-(ea-1.0e0_knd)*(ea+1.0e0_knd)*enr(i)/((arr+rm2)*
     1             (arr+rm2-1.0e0_knd)*csq)
            if(i.le.lm2) then
            if(real(enr(i)).lt.(0.0e0_knd))  sgn=sgn*(-1.0e0_knd)
            end if
20        continue
c
c  compute the Morse-Feshbach normalizing factor
        term=(1.0e0_knd,0.0e0_knd)
        dmfnorm=term
        sumpr=1.0e0_knd
        sumpi=0.0e0_knd
        jlow=l-m+2
        jterm=lm2
        iflag=0
        idmfe=0
        fterm=dmfnorm
          do 30 j=jlow,limd,2
          aj=j
          jterm=jterm+1
          term=term*(aj+rm2)*enr(jterm)*(aj+rm2-1.0e0_knd)/
     1         (aj*(aj-1.0e0_knd))
          if(real(term).gt.0.0e0_knd) sumpr=sumpr+real(term)
          if(aimag(term).gt.0.0e0_knd) sumpi=sumpi+aimag(term)
          dmfnorm=dmfnorm+term
          if(abs(term/dmfnorm).lt.dec) go to 40
          if(abs(dmfnorm).lt.teste) go to 30
          dmfnorm=dmfnorm*testeo
          term=term*testeo
          sumpr=sumpr*testeo
          sumpi=sumpi*testeo
          idmfe=idmfe+nfac
          iflag=1
30        continue
40      jlow=l-m
        jmf=jterm
        if(jlow.lt.2.or.iflag.eq.1) go to 60
        term=(1.0e0_knd,0.0e0_knd)
        jterm=lm2
          do 50 j=jlow,2,-2
          aj=j
          term=term*aj*(aj-1.0e0_knd)/((aj+rm2
     1         -1.0e0_knd)*(aj+rm2)*enr(jterm))
          if(real(term).gt.0.0e0_knd) sumpr=sumpr+real(term)
          if(aimag(term).gt.0.0e0_knd) sumpi=sumpi+aimag(term)
          jterm=jterm-1
          dmfnorm=dmfnorm+term
          if(abs(term/dmfnorm).lt.dec) go to 60
50        continue
60      continue
        jsubr=ndec
        if(real(dmfnorm).ne.0.0e0_knd) jsubr=
     1              int(log10(abs(sumpr/real(dmfnorm))+dec))
        if(jsubr.lt.0) jsubr=0
        if(jsubr.gt.ndec) jsubr=ndec
        jsubi=0
        if(aimag(dmfnorm).ne.0.0e0_knd) jsubi=
     1              int(log10(abs(sumpi/aimag(dmfnorm))+dec))
        if(jsubi.lt.0) jsubi=0
        if(jsubi.gt.ndec) jsubi=ndec
        iterm=0
        if(abs(dmfnorm).ne.0.0e0_knd) iterm=int(log10(abs(dmfnorm)))
        idmfe=idmfe+iterm
        dmfnorm=dmfnorm*ten**(-iterm)
        dmlmf=1.0e0_knd/dmfnorm
        idmlmfe=-idmfe
c        if(ioprad.ne.0) write(40,70) jmf,jsubr,jsubi
c        if(ioprad.eq.0) write(50,70) jmf,jsubr,jsubi
c70      format(15x,'Morse-Feshbach norm. converged',
c     1         ' in ',i6,' terms with ',i2,' and ',i2,' digits of',
c     2         ' subtr.',/,23x,'error in the real and imag. parts')
        jsubmf=jsubr
          if(aimag(cc).ne.0.0e0_knd.and.aimag(dmfnorm).ne.0.0e0_knd)
     1        then
          jcor=int(log10(abs(real(dmfnorm)/aimag(dmfnorm))+dec))
          if(jcor.ge.0) jsubmf=max(jsubr,jsubi-jcor)
          if(jcor.lt.0) jsubmf=max(jsubr+jcor,jsubi)
          end if
          if(jsubmf.gt.ndec) jsubmf=ndec
c
c  compute the Flammer normalizing factor
        term=(1.0e0_knd,0.0e0_knd)
        sumpr=1.0e0_knd
        sumpi=0.0e0_knd
        dfnorm=term
        idfe=0
        iflag=0
          do 80 j=lm2+1,lim2
          jj=j+j+ix
          term=-term*enr(j)*(jj+mml)/(jj-ix)
          dfnorm=dfnorm+term
          if(real(term).gt.0.0e0+knd) sumpr=sumpr+real(term)
          if(aimag(term).gt.0.0e0_knd) sumpi=sumpi+aimag(term)
          if(abs(term/dfnorm).lt.dec) go to 90
          if(abs(dfnorm).lt.teste) go to 80
          dfnorm=dfnorm*testeo
          term=term*testeo
          sumpr=sumpr*testeo
          sumpi=sumpi*testeo
          idfe=idfe+nfac
          iflag=1
80        continue
90      continue
        jf=min(j,lim2)
        if(lm2.lt.1.or.iflag.eq.1) go to 110
        term=(1.0e0_knd,0.0e0_knd)
          do 100 j=lm2,1,-1
          jj=j+j+ix
          term=-term*(jj-ix)/((jj+mml)*enr(j))
          dfnorm=dfnorm+term
          if(real(term).gt.0.0e0_knd) sumpr=sumpr+real(term)
          if(aimag(term).gt.0.0e0_knd) sumpi=sumpi+aimag(term)
          if(abs(term/dfnorm).lt.dec) go to 110
100       continue
110     continue
        jnorm=max(jf,jmf)
        jsubr=ndec
        if(real(dfnorm).ne.0.0e0_knd) jsubr=
     1            int(log10(abs(sumpr/real(dfnorm))+dec))
        if(jsubr.lt.0) jsubr=0
        if(jsubr.gt.ndec) jsubr=ndec
        jsubi=0
        if(aimag(dfnorm).ne.0.0e0_knd) jsubi=
     1            int(log10(abs(sumpi/aimag(dfnorm))+dec))
        if(jsubi.lt.0) jsubi=0
        if(jsubi.gt.ndec) jsubi=ndec
        iterm=0
        if(abs(dfnorm).ne.0.0e0_knd) iterm=int(log10(abs(dfnorm)))
        idfe=idfe+iterm
        dfnorm=dfnorm*ten**(-iterm)
        dmlf=1.0e0_knd/dfnorm
        idmlfe=-idfe
c        if(ioprad.ne.0) write(40,120) jf,jsubr,jsubi
c        if(ioprad.eq.0) write(50,120) jf,jsubr,jsubi
c120     format(15x,'Flammer norm. converged in ',
c     1         i6,' terms with ',i2,' and ',i2,' digits of ',
c     2         'subt.',/,22x,' error in real and imag. parts.')
        jsubf=jsubr
          if(aimag(cc).ne.0.0e0_knd.and.aimag(dfnorm).ne.0.0e0_knd) then
          jcor=int(log10(abs(real(dfnorm)/aimag(dfnorm))+dec))
          if(jcor.ge.0) jsubf=max(jsubr,jsubi-jcor)
          if(jcor.lt.0) jsubf=max(jsubr+jcor,jsubi)
          end if
          if(jsubf.gt.ndec) jsubf=ndec
c
c  compute the d0(c|ml) or d1(c|ml)
        idc01=0
   	dc01=(1.0e0_knd,0.0e0_knd)
        if(lm2.eq.0) go to 140
          do 130 kjl=1,lm2
          kkjl=lm2-kjl+1
    	  dc01=dc01/enr(kkjl)
            if(abs(dc01).gt.teste) then
            dc01=dc01*testeo
            idc01=idc01+nfac
            end if
            if(abs(dc01).lt.testeo) then
            dc01=dc01*teste
            idc01=idc01-nfac
            end if
130       continue
        iterm=int(log10(abs(dc01)))
        dc01=dc01*(ten**(-iterm))
        idc01=idc01+iterm
140     continue
c
c  compute the Meixner-Schafke normalizing factor
        jflag=0
        idmse=0
        coef=(1.0e0_knd,0.0e0_knd)
        dmsnorm=coef
        fterm=1.0e0_knd
        jlow=l-m+2
        jterm=lm2
          do 150 j=jlow,limd,2
          aj=j
          aj2=aj+aj
          jterm=jterm+1
          coef=coef*(aj+rm2)*enr(jterm)*(aj+rm2m1)*enr(jterm)
     1         *(aj2+rm2m3)/(aj*(aj-1.0e0_knd)*(aj2+rm2p1))
          dmsnorm=dmsnorm+coef
          if(abs(dmsnorm).gt.fterm) fterm=abs(dmsnorm)
          if(abs(coef/dmsnorm).lt.dec) go to 160
            if(abs(dmsnorm).gt.teste) then
            dmsnorm=dmsnorm*testeo
            coef=coef*testeo
            idmse=idmse+nfac
            fterm=fterm*testeo
            jflag=1
            end if
150       continue
160     jlow=l-m
        jn=jterm
        jnorm=max(jnorm,jn)
        if(jlow.lt.2.or.jflag.eq.1) go to 180
        coef=(1.0e0_knd,0.0e0_knd)
        jterm=lm2
        j=jlow
          do 170 jj=2,jlow,2
          aj=j
          aj2=aj+aj
          coef=coef*aj*(aj-1.0e0_knd)*(aj2+rm2p1)/((aj2+rm2m3)*
     1            enr(jterm)*enr(jterm)*(aj+rm2)*(aj+rm2m1))
          jterm=jterm-1
          j=j-2
          dmsnorm=dmsnorm+coef
          if(abs(dmsnorm).gt.fterm) fterm=abs(dmsnorm)
          if(abs(coef/dmsnorm).lt.dec) go to 180
170       continue
180     jsubms=int(log10((fterm/abs(dmsnorm))+dec))
        if(jsubms.lt.0) jsubms=0
        if(jsubms.gt.ndec) jsubms=ndec
        iterm=int(log10(abs(dmsnorm)))
        dmsnorm=dmsnorm*ten**(-iterm)
        idmse=idmse+iterm
          if(2*(idmse/2).ne.idmse) then
          idmse=idmse-1
          dmsnorm=ten*dmsnorm
          end if
        dmlms=sgn/sqrt(dmsnorm)
        idmlmse=-idmse/2
c        if(iopang.ne.0) write(50,190) jn,jsubms
c190     format(5x,' Meixner-Schafke normalization converged in ',
c     1         i6,' terms with ',i2,' digits of subt. error')
200     continue
        return
        end
c
c
        subroutine dalt (l,m,cc,ndec,nex,limdr,maxdr,maxmp,ioppsum,
     1                   eigval,enrneg,drhor,dneg,idneg,nsdneg,nsdrhor1,
     2                   nsdrho)
c
c  purpose:     To calculate d ratios with negative subscripts
c               and d-rho ratios.
c  parameters:
c
c     input:    l       : l
c               m       : m
c               cc      : complex c
c               ndec    : number of decimal digits of precision in
c                         real(knd) arithmetic
c               nex     : maximum exponent in real(knd) arithmetic
c               limdr   : number of ratios of successive d-rho
c                         coefficients calculated
c               maxdr   : dimension of drhor array
c               maxmp   : dimension of enrneg array
c               ioppsum : integer index = 0 if no d rho coefficients
c                         are calculated (psum not needed for r2leg)
c               eigval  : eigenvalue
c
c     output:   enrneg  : array of d coefficient ratios with
c                         negative subscripts
c               drhor   : array of d rho coefficient ratios
c               dneg    : characteristic of the ratio of the d
c                         coefficient with index -2m+ix to the
c                         d coefficient with index ix, where
c                         ix = 0 if l-m is even and ix = 1 if
c                         l-m is odd
c               idneg   : exponent (base 10) of dneg
c               nsdneg  : subtaction error in calculating dneg
c               nsdrhor1: subtraction error in the step calculating
c                         drhor(1) from drhor(2)
c               nsdrho  : subtraction error in calculating drhor(1)
c
c        use param
c
c  real(knd) scalars and arrays
        real(knd) r,rm,rn,t,ten,teste,testeo,uterm,wterm
        real(knd) amnsdrho,ansdneg,ansdrho,asub,bsub
c
c  complex(knd) scalars and arrays
        complex(knd) cc,dneg,eigval
        complex(knd) enrneg(maxmp),drhor(maxdr),enrnega(maxmp)
        complex(knd) vterm
c
        ten=10.0e0_knd
        nfac=nex/3
        if(2*(nfac/2).ne.nfac) nfac=nfac-1
        teste=ten**(nfac)
        testeo=1.0e0_knd/teste
c  if l-m is even, ix=0; if l-m is odd, ix=1
        ix=(l-m)-2*((l-m)/2)
        t=m+m-ix-ix
c
c  calculate ratios of d coefficients with negative subscripts
c
c       enrneg(k) = { d(-2m+2k-2)/d(-2m+2k), l-m even    }
c                   { d(-2m+1+2k-2)/d(-2m+1+2k), l-m odd }
c
        rm=m
          if(m.eq.0) then
          dneg=(1.0e0_knd,0.0e0_knd)
          idneg=0
          go to 30
          end if
          do 10 i=1,m+1
          enrneg(i)=(0.0e0_knd,0.0e0_knd)
10        continue
c
c  first calculate enrneg(1)
        n=2-2*m+ix
        rn=n
        r=n+m+m
        uterm=r*(r-1.0e0_knd)/((rn+r+1.0e0_knd)*(rn+r-1.0e0_knd))
        r=n+m-2
        vterm=(2.0e0_knd*r*(r+1.0e0_knd)-2.0e0_knd*rm*rm-
     1        1.0e0_knd)/((2.0e0_knd*r+3.0e0_knd)*(2.0e0_knd*r-
     2        1.0e0_knd))-(r*(r+1.0e0_knd)-eigval)/(cc*cc)
c
c       calculations continue up to and including
c       enrneg(k=m) = { d(-2)/d(0), l-m even }
c                     { d(-1)/d(1), l-m odd  }
c
        enrneg(1)=-uterm/vterm
        dneg=enrneg(1)
        idneg=0
        ansdneg=0.0e0_knd
        nsdneg=0
c
c  backward recursion beginning with enrneg(1) and
c  ending with enrneg(m)
c
          do i=2,2*m-2,2
          ii=i-2*m+ix
          n=ii+2
          j=i/2
          rn=n
          r=n+m+m

          uterm=r*(r-1.0e0_knd)/((rn+r+1.0e0_knd)*(rn+r-
     1               1.0e0_knd))
          r=n+m-2

          vterm=(2.0e0_knd*r*(r+1.0e0_knd)-2.0e0_knd*rm*rm-
     1          1.0e0_knd)/((2.0e0_knd*r+3.0e0_knd)*(2.0e0_knd*r-
     2          1.0e0_knd))-(r*(r+1.0e0_knd)-eigval)/(cc*cc)
          r=n-4
          wterm=(r+2.0e0_knd)*(r+1.0e0_knd)/((r+r+rm+rm+
     1                3.0e0_knd)*(r+r+rm+rm+1.0e0_knd))
          enrneg(j+1)=-uterm/(wterm*enrneg(j)+vterm)
          dneg=dneg*enrneg(j+1)
            if(wterm*enrneg(j)*vterm.ne.(0.0e0_knd,0.0e0))
     1         then
            asub=log10(abs(vterm/(vterm+
     1                       wterm*enrneg(j))))
            if(asub.gt.0.0e0_knd) ansdneg=ansdneg+asub
            bsub=log10(abs(vterm/(wterm*enrneg(j))))
            if(bsub.gt.0.0e0_knd) ansdneg=max(0.0e0_knd,ansdneg-bsub)
            if(int(ansdneg).gt.nsdneg) nsdneg=ansdneg
            end if
            if(abs(dneg).gt.teste) then
            dneg=dneg*testeo
            idneg=idneg+nfac
            end if
            if(abs(dneg).lt.testeo) then
            dneg=dneg*teste
            idneg=idneg-nfac
            end if
          end do
          if(nsdneg.gt.ndec) nsdneg=ndec
        iterm=int(log10(abs(dneg)))
        dneg=dneg*(ten**(-iterm))
        idneg=idneg+iterm
c
c  calculate ratios of d rho coefficients
c
c       drhor(k-m) = { d(rho|2k)/d(rh0|2k-2), l-m even  }
c                    { d(rho|2k-1)/d(rho|2k-3), l-m odd }
c
30      if(ioppsum.eq.0) go to 60
          ansdrho=0.0e0_knd
          amnsdrho=0.0e0_knd
          nsdrho=0
          mnsdrho=0
          do 40 i=1,limdr
          drhor(i)=(0.0e0_knd,0.0e0_knd)
40        continue
          do 50 i=2*limdr,6,-2
          n=4-i+ix-m-m
          ii=(i-2)/2
          rn=n
          r=n+m+m
          uterm=r*(r-1.0e0_knd)/((rn+r+1.0e0_knd)*(rn+r-1.0e0_knd))
          r=n+m-2
          vterm=(2.0e0_knd*r*(r+1.0e0_knd)-2.0e0_knd*rm*rm-
     1           1.0e0_knd)/((2.0e0_knd*r+3.0e0_knd)*(2.0e0_knd*r-
     2           1.0e0_knd))-(r*(r+1.0e0_knd)-eigval)/(cc*cc)
          r=n-4
          wterm=(r+2.0e0_knd)*(r+1.0e0_knd)/((r+r+rm+rm+3.0e0_knd)*
     1          (r+r+rm+rm+1.0e0_knd))
          drhor(ii)=-uterm/(wterm*drhor(ii+1)+vterm)
            if(wterm*drhor(ii+1)*vterm.ne.(0.0e0_knd,0.0e0))
     1         then
            asub=log10(abs(vterm/(vterm+wterm*drhor(ii+1))))
            if(asub.gt.0.0e0_knd) ansdrho=ansdrho+asub
            bsub=log10(abs(vterm/(wterm*drhor(ii+1))))
            if(bsub.gt.0.0e0_knd) ansdrho=max(0.0e0_knd,ansdrho-bsub)
            if(ansdrho.gt.amnsdrho) amnsdrho=ansdrho
            end if
50        continue
        n=-2*m+ix
        r=n+m-2
        vterm=(2.0e0_knd*r*(r+1.0e0_knd)-2.0e0_knd*rm*rm-1.0e0_knd)/
     1        ((2.0e0_knd*r+3.0e0_knd)*(2.0e0_knd*r-1.0e0_knd))-
     2        (r*(r+1.0e0_knd)-eigval)/(cc*cc)
        r=n-4
        wterm=(r+2.0e0_knd)*(r+1.0e0_knd)/((r+r+rm+rm+3.0e0_knd)*
     1        (r+r+rm+rm+1.0e0_knd))
c
c       the final value of ii is 1;
c       drhor(1) has a special value:
c       drhor(1) = { d(rho|2m+2)/d(-2m), l-m even  }
c                  { d(rho|2m+1)/d(-2m+1), l-m odd }
c
        drhor(1)=1.0e0_knd/((t-1.0e0_knd)*(t+1.0e0_knd)*
     1           (wterm*drhor(2)+vterm))
          if(wterm*drhor(2)*vterm.ne.(0.0e0_knd,0.0e0))
     1         then
          asub=log10(abs(vterm/(vterm+wterm*drhor(2))))
          if(asub.gt.0.0e0_knd) ansdrho=ansdrho+asub
          bsub=log10(abs(vterm/(wterm*drhor(2))))
          if(bsub.gt.0.0e0_knd) ansdrho=max(0.0e0_knd,ansdrho-bsub)
          if(ansdrho.gt.amnsdrho) amnsdrho=ansdrho
          end if
        nsdrho=int(amnsdrho)+1
        nsdrhor1=asub
        if(nsdrhor1.lt.0) nsdrhor1=0
        if(ix.eq.1) drhor(1)=-drhor(1)
60      continue
        return
        end
c
c
	subroutine gauss (n,x,w)
c
c  purpose:     To evaluate the coordinates and weighting factors
c               for an nth order Gaussian quadrature
c
c  parameters:
c
c     input:    n  : order of quadrature
c
c     output:   x  : coordinate values for quadrature
c               w  : weighting factors
c
c        use param
c
c  real(knd) scalars and arrays
        real(knd) delta,der,pi,ri,s,t,ten,test,u,v,z
        real(knd) x(n),w(n)
c
        ndec=precision(pi)
        ten=10.0e0_knd
        test=ten**(-ndec)
        imax=(n+1)/2
        pi=acos(-1.0_knd)
          do 40 i=1,imax
          ri=i
	  z=cos(pi*(ri-0.25e0_knd)/(n+0.5e0_knd))
            do 20 j=1,30
            u=0.0e0_knd
	    v=1.0e0_knd
	      do 10 k=1,n
	      t=u
              u=v
	      v=((k+k-1)*z*u-(k-1)*t)/k
10   	      continue
            s=z*z-1.0e0_knd
	    der=n*(z*v-u)/s
	    delta=-v/der-0.5e0_knd*v*v*((n*n*s-n*z*z-n)*v+
     1            2.0e0_knd*n*z*u)/(der*der*der*s*s)
            z=z+delta
	    if(abs(delta/z).lt.test) go to 30
20          continue
30        continue
	  x(i)=-z
	  x(n+1-i)=z
	  w(i)=2.0e0_knd/((1.0e0_knd-z*z)*der*der)
	  w(n+1-i)=w(i)
40	  continue
	return
	end
c
c
        subroutine pleg (m,lim,maxp,limcsav,iopd,ndec,barg,narg,maxt,pr,
     1                   pdr,pdnorm,ipdnorm,pnorm,ipnorm,alpha,beta,
     2                   gamma,coefa,coefb,coefc,coefd,coefe)
c
c  purpose:     To calculate ratios of successive associated Legendre
c               functions of the first kind for given arguments barg.
c               to calculate corresponding ratios of their first
c               derivatives. To calculate the characteristics and
c               exponents of both the Legendre functions of the first
c               kind and their first derivatives.
c
c  parameters:
c
c     input:    m      : m
c               lim    : two less than the number of associated Legendre
c                        function ratios calculated for given arguments
c               maxp   : dimension of alpha, beta, gamma, coefa, coefb,
c                        coefc, coefd, and coefe arrays and second
c                        dimension of pr and pdr arrays
c               limcsav: integer equal to the number of coefficients in
c                        each of the arrays alpha, beta, gamma, coefa,
c                        coefb, coefc, coefd, and coefe arrays that
c                        have already been calculated in earlier calls
c                        to pleg for this value of m and will not be
c                        calculated again. [Note that the minimum
c                        array index for the coefficients is 3 and
c                        the maximum array index is limcsav+2]
c               iopd   : integer that is set = 0 if derivatives of
c                        Legendre functions (i.e., their ratios)
c                        are not required when iopang = 1 and the
c                        first derivatives of the angular functions
c                        are not requested.
c                        iopd is set = 1 when iopang = 2 and pleg is
c                        also being used to obtain ratios of first
c                        derivatives of Legendre functions for use in
c                        computing the first derivatives of the angular
c                        functions.
c                        iopd is set = 2 when pleg is being used to
c                        compute ratios of Legendre functions for use in
c                        the calculation of the denominator term used in
c                        calculating the radial functions of the second
c                        kind and their first derivatives in r2eta. Also
c                        used in r1eta in the same way when calculating
c                        radial functions of the first kind and their
c                        first derivatives.
c                        iopd is set = 3 when pleg is being used to
c                        compute ratios of both the Legendre functions
c                        and their first derivatives for use in the
c                        calculation of the numerator terms used
c                        in r2eta to calculate the radial functions of
c                        the second kind and their first deriatives. Also
c                        used in r1eta in the same way when calculating
c                        radial functions of the first kind and their
c                        first derivatives.
c                        iopd is set = 4 when pleg is being used to
c                        compute ratios of both the legendre functions
c                        and their first derivatives for use is the
c                        calculation of r2 and r2d in r2leg.
c               ndec   : number of decimal digits in real(knd)
c                        arithmetic
c               barg   : array of narg values of eta for which Legendre
c                        functions are to be calculated
c               narg   : number of specified values of eta in barg array
c               maxt   : dimension of barg array
c
c     output:   pr     : array of ratios of successive first kind
c                        associated Legendre functions of the same
c                        parity
c               pdr    : array of ratios of successive derivatives of
c                        first kind associated Legendre functions of
c                        the same parity
c               pdnorm : array of characteristics of the first
c                        derivatives of associated Legendre functions
c                        of the first kind of order m and degree m
c               ipdnorm: array of exponents of the first derivatives
c                        of associated Legendre functions of the first
c                        kind of order m and degree m
c               pnorm  : array of characteristics of the associated
c                        Legendre functions of the first kind of order m
c                        and degree m
c               ipnorm : array of exponents of the associated Legendre
c                        functions of the first kind of order m and
c                        degree m
c
c     input/output:
c               alpha  : array of coefficients in the recursion
c                        formula for the associated Legendre functions
c               beta   : array of coefficients in the recursion
c                        formula for the associated Legendre functions
c               gamma  : array of coefficients in the recursion
c                        formula for the associated Legendre functions
c               coefa  : array of coefficients in the expression
c                        relating the derivative ratios pdr to the
c                        function ratios pr
c               coefb  : array of coefficients in the expression
c                        relating the derivative ratios pdr to the
c                        function ratios pr
c               coefc  : array of coefficients in the expression
c                        relating the derivative ratios pdr to the
c                        function ratios pr
c               coefd  : array of coefficients in the expression
c                        relating the derivative ratios pdr to the
c                        function ratios pr
c               coefe  : array of coefficients in the expression
c                        relating the derivative ratios pdr to the
c                        function ratios pr
c
c        use param
c
c  real(knd) scalars and arrays
        real(knd) adec,ajterm,am2p1,anden1,anden2,an2tnp1,bargs,coef,
     1            den,rm,rm2,temp1,temp2,temp3,ten,term,teste,testeo,
     2            ta,tb,tc,t1,t2,t3
        real(knd) alpha(maxp),barg(maxt),beta(maxp),coefa(maxp),
     1            coefb(maxp),coefc(maxp),coefd(maxp),coefe(maxp),
     2            gamma(maxp),pdnorm(maxt),pdr(maxt,maxp),pdr1(maxp),
     3            pr(maxt,maxp),pnorm(maxt)
c
c  integer array
        integer ipdnorm(maxt),ipnorm(maxt)
c
        ten=10.0e0_knd
        adec=ten**(-ndec+2)
        nex=range(barg(1))-1
        nfac=nex/3
        if(2*(nfac/2).ne.nfac) nfac=nfac-1
        teste=ten**nfac
        testeo=1.0e0_knd/teste
        rm=m
        rm2=m*m
        am2p1=m+m+1
        m2=m+m
        m2m1=m2-1
        mm1=m-1
        mm2=m-2
        msqp1=2*m*m+1
        coef=1.0e0_knd
        if(iopd.eq.4) coef=-1.0e0_knd
c
c  calculate the coefficients alpha(j), beta(j), and gamma(j) for
c  the three term recursion relating the Legendre function ratios
c
c              m                m
c   pr(k,j) = p    (barg(k)) / p  (barg(k))
c              m+j-1            m+j-3
c
c  and calculate the coefficients coefa(j), coefb(j), coefc(j),
c  coefd(j), and coefe(j) in the expression used to calculate
c  ratios of Legendre function derivatives
c
c               m                 m
c   pdr(k,j) = p'    (barg(k)) / p'  (barg(k))
c               m+j-1             m+j-3
c
        if(limcsav.ge.lim) go to 30
          do 10 j=limcsav+3,lim+2
      	  n=m+j-3
      	  n2=n+n
      	  n2p3=n2+3
      	  n2p1=n2+1
      	  n2m1=n2-1
      	  nmmm2=n-mm2
      	  nmmm1=n-mm1
      	  npmm1=n+mm1
      	  npm=n+m
      	  npmm1=n+mm1
          npmp1=n+m+1
          npmp2=npmp1+1
      	  an2tnp1=2*n*real(n+1,knd)
      	  anden1=nmmm2*real(nmmm1,knd)
      	  anden2=n2m1*anden1
      	  alpha(j)=real(n2p3,knd)*n2p1/anden1
      	  beta(j)=real(n2p1,knd)*(real(msqp1,knd)-an2tnp1)/anden2
      	  gamma(j)=-real(n2p3,knd)*real(npm,knd)*real(npmm1,knd)/
     1              anden2
          coefa(j)=-real(npmp2,knd)/real(nmmm1,knd)
          coefb(j)=real(n2p3,knd)*(n+2)/anden1
          coefc(j)=-real(npmp1,knd)*npmp2/anden1
          coefd(j)=real(npmp1,knd)/real(nmmm2,knd)
          coefe(j)=-real(n+1,knd)*n2p3/anden1
10        continue
        gamma(3)=0.0e0_knd
        gamma(4)=0.0e0_knd
        term=1.0e0_knd
        iterm=0
        if(m.lt.2) go to 30
          do jm=2,m
      	  term=(jm+jm-1)*term
      	    if(term.gt.teste) then
      	    term=term*testeo
      	    iterm=iterm+nfac
            end if
          end do
        jterm=int(log10(term))
        term=term*(ten**(-jterm))
        iterm=iterm+jterm
30      continue
c
c   calculate the ratios of successive Legendre functions of the same
c   parity using the three term recursion relationship
c
c   pr(k,j) = coef*alpha(j)*barg(k)*barg(k) + beta(j) +
c             gamma(j)/pr(k,j-2)
c
c   where coef = -1 when computing functions for r2leg, = +1 otherwise
c
          do 140 k=1,narg
      	  pnorm(k)=term
          ipnorm(k)=iterm
          pdnorm(k)=term
          ipdnorm(k)=iterm
c
c   define the first two ratios equal to unity and (2m+1)*barg(k)
          pr(k,1)=1.0e0_knd
          pr(k,2)=am2p1*barg(k)
          jdelta=1
          if(abs(barg(k)).lt.adec) jdelta=2
          bargs=barg(k)*barg(k)
            do 40 j=3,lim+2,jdelta
            pr(k,j)=coef*alpha(j)*bargs+beta(j)+(gamma(j)/pr(k,j-2))
40          continue
c
c   calculate the corresponding ratios of first derviatives using the
c   relationship (except for eta equal to zero or unity, where special
c   expressions are used)
c
c              (coefa(j)+coef*coefb(j)*barg(k)*barg(k))*pr(k,j)+coefc(j)
c   pdr(k,j) = ----------------------------------------------------
c                  pr(k,j)+coef*coefd(j)+coef*coefe(j)*barg(k)*barg(k)
c
c   where coef = -1 when computing functions for r2leg, = +1 otherwise
c
c   for abs(eta) <= 0.1, calculate the ratios instead using a recursion
c   relation involving successive ratios (except for eta equal to 0,
c   where special expressions are used)
c
          if(iopd.eq.0.or.iopd.eq.2) go to 120
          pdr(k,1)=1.0e0_knd
          pdr(k,2)=1.0e0_knd
          if(abs(barg(k)).ge.adec) go to 50
          pdr(k,2)=am2p1
            do j=4,lim+2,2
            pdr(k,j)=-real(m2m1+j,knd)/(j-2)
            end do
          go to 140
50        if(abs(abs(barg(k))-coef).ge.adec) go to 70
          if(m.eq.0) go to 60
          if(m.ne.2) go to 130
          pdr(k,1)=-2.0e0_knd*barg(k)
          go to 80
60        temp1=1.0e0_knd
          temp2=3.0e0_knd
          pdr(k,2)=1.0e0_knd
          pdr(k,3)=3.0e0_knd*barg(k)
            do j=4,lim+2
            temp3=temp2+real(j-1,knd)
            pdr(k,j)=temp3/temp1
            temp1=temp2
            temp2=temp3
            end do
          go to 140
70        if(m.ne.0) go to 80
          pdr(k,1)=1.0e0_knd
          pdr(k,2)=1.0e0_knd
          pdr(k,3)=3.0e0_knd*barg(k)
          jlow=4
          go to 90
80        pdr(k,2)=am2p1*(+coef*(rm+1.0e0_knd)*bargs-1.0e0_knd)/
     1             (rm*barg(k))
          jlow=3
90        continue
          if(abs(barg(k)).le.0.1e0_knd) go to 110
            do 100 j=jlow,lim+2
            den=(pr(k,j)+coefd(j)+coef*coefe(j)*bargs)
            if(den.eq.0.0e0_knd) den=1.0e-50_knd
            pdr(k,j)=((coefa(j)+coef*coefb(j)*bargs)*pr(k,j)+
     1               coefc(j))/den
            if(iopd.eq.3.and.abs(pdr(k,j)).lt.1.0e-5_knd) go to 110
100        continue
         go to 120
110      continue
         if(m.ne.0) pdr1(1)=pdr(k,2)
         if(m.eq.0) pdr1(2)=pdr(k,3)
           do j=jlow-1,lim+1
           n=j+m-1
           t3=coef
           if(coef.eq.-1.0e0_knd.and.2*(j/2).eq.j) t3=1.0e0_knd
           t1=coef*bargs-1.0e0_knd
           t2=(n*n*t1+rm2)
           ta=j*t2
           tb=(n+n+1)*t3*barg(k)*(t2+n*t1)
           tc=-(n+m)*(t2+(n+n+1)*t1)
           pdr1(j)=(tb+tc/pdr1(j-1))/ta
           end do
           do j=jlow,lim+2
           pdr(k,j)=pdr1(j-2)*pdr1(j-1)
           end do
120       continue
          if(m.eq.0.or.iopd.eq.2.or.iopd.eq.3.or.iopd.eq.4) go to 140
          if(abs(abs(barg(k))-1.0e0_knd).lt.adec) go to 130
          ajterm=rm*log10(1.0e0_knd-bargs)/2.0e0_knd
          jterm=int(ajterm)
          ipnorm(k)=ipnorm(k)+jterm
          pnorm(k)=pnorm(k)*(ten**(ajterm-real(jterm,knd)))
          if(iopd.eq.0) go to 140
          ajterm=log10(rm*abs(barg(k)))+(rm-2.0e0_knd)*
     1           log10(1.0e0_knd-bargs)/2.0e0_knd
          jterm=int(ajterm)
          ipdnorm(k)=ipdnorm(k)+jterm
          pdnorm(k)=-pdnorm(k)*(ten**(ajterm-real(jterm,knd)))
          if(barg(k).lt.0.0e0_knd) pdnorm(k)=-pdnorm(k)
          go to 140
130       pnorm(k)=0.0e0_knd
          ipnorm(k)=0
          if(m.ne.2) pdnorm(k)=0.0e0_knd
          if(m.ne.2) ipdnorm(k)=0
140       continue
        return
        end
c
c
        subroutine qleg (m,lnum,limq,maxmp,maxq,x,ndec,nex,iflagl1,qdr,
     1                   qdqr,qdm1,iqdm1,qdl,iqdl,qr,qm1,iqm1,ql,iql,
     2                   termpq,itermpq,qr1,qdr1,qm0,qdm0)
c
c  purpose:     To calculate ratios of successive associated Legendre
c               functions of the second kind for given c,x, and m.
c               to calculate corresponding ratios of their first
c               derivatives. To calculate the characteristics and
c               exponents of both the Legendre functions of the second
c               kind and their first derivatives.
c
c  parameters:
c
c     input:    m      : m
c               lnum   : number of l values desired (=lmax+1);
c                        also equal to the dimension of the arrays
c                        ql, iql, qdl, and iqdl
c               limq   : the number of associated Legendre function
c                        ratios calculated for given m,lnum,c,ndec,
c                        and x
c               maxmp  : dimension of qdqr array
c               maxq   : dimension of qr and qdr arrays
c               x      : radial coordinate x
c               ndec   : number of decimal digits in real(knd)
c                        arithmetic
c               nex    : maximum exponend in real(knd) arithmetic
c               iflagl1: equal to 1 if Legendre function ratios
c                        used in subroutine r2leg1 will be computed,
c                        equal to zero otherwise
c     output:   qdr    : array of ratios of derivatives of successive
c                        associated Legendre functions of the second
c                        kind
c               qdqr   : array of ratios of derivatives of associated
c                        Legendre functions of the second kind to the
c                        corresponding Legendre function for degrees
c                        from -m to m-1
c               qdm1   : characteristic of the first derivative of
c                        the associated Legendre function of the second
c                        kind with order m and degree m-1, scaled by
c                                       -m/2
c                        (2m-1)!!(x*x+1)
c               iqdm1  : exponent corresponding to qdm1
c               qdl    : array of characteristics of the first
c                        derivatives of the associated Legendre
c                        functions of the second kind with order m
c                        and degrees from m to m+lnum-1, scaled by
c                                       -m/2
c                        (2m-1)!!(x*x+1)
c               iqdl   : array of exponents corresponding to qdl
c               qr     : array of ratios of successive associated
c                        Legendre functions of the second kind
c               qm1    : characteristic of the associated Legendre
c                        function of the second kind with order m
c                                                                 -m/2
c                        and degree m-1, scaled by (2m-1)!!(x*x+1)
c               iqm1   : exponent corresponding to qm1
c               ql     : array of characteristics of the associated
c                        Legendre function of the second kind with
c                        order m and degrees from m to m+lnum-1
c                                                 -m/2
c                        scaled by (2m-1)!!(x*x+1)
c               iql    : array of exponents corresponding to ql
c               termpq : characteristic of the relative size of the
c                        maximum terms in the positive degree q series
c                        and the p series used to calculate r2 and r2d
c                        in subroutine r2leg
c               itermpq: exponent corresponding to termpq
c               qr1    : array of ratios of successive associated
c                        Legendre functions of the second kind
c                        reindexed for use in the Baber and Hasse
c                        expansion
c               qdr1   : array of ratios of derivatives of successive
c                        associated Legendre functions of the second
c                        kind reindexed for use in the Baber and Hasse
c                        expansion
c               qm0    : characteristic of the associated Legendre
c                        function of the second kind with order m
c                                                                 -m/2
c                        and degree 0, scaled by (2m-1)!!(x*x+1)
c               qdm0   : characteristic of the first derivative of
c                        the associated Legendre function of the second
c                        kind with order m and degree 0, scaled by
c                                       -m/2
c                        (2m-1)!!(x*x+1)
c
c        use param
c
c  real(knd) scalars and arrays
        real(knd) ajm,dec,qdm1,qlow,qlow0,qlow1,qmid,qmid0,qmid1,qm1,
     1            qmm,qupp,qupp0,qupp1,q00,q0m,q1m,q11,qmmm1,rin,rin1,
     2            rin2,rin3,rin4,rm,rmsq,ten,term,termpq,teste,testeo,
     3            tjm,tm,tmr,x,xang,xfac,xc,x1d,xsqr
        real(knd) qdl(lnum),qdr(maxq),ql(lnum),qr(maxq),qdqr(maxmp)
        real(knd) qr1(maxq),qdr1(maxq),qdm0,qm0
c
c  integer arrays
        integer iqdl(lnum),iql(lnum)
c
c  Note that the factor i involved in these functions is suppressed
c  and taken into account in calculations in the subroutine r2leg
c
        ten=10.0e0_knd
        dec=ten**(-ndec)
        nfac=nex/3
        if(2*(nfac/2).ne.nfac) nfac=nfac-1
        teste=ten**(nfac)
        testeo=1.0e0_knd/teste
        rm=m
        rmsq=rm*rm
        tm=rm+rm
        tmr=tm/(tm+1.0e0_knd)
        x1d=x*x+1.0e0_knd
        xsqr=sqrt(x1d)
        xang=asin(x/xsqr)
        q00=-atan(1.0e0_knd/x)
        if(m.eq.0) qm0=q00
        iflag=0
        if(x.lt.0.001e0_knd) iflag=1
c
c                m
c  For x < 0.1: Q   is calculated by forward recursion in m starting
c                m
c                   0       1
c  with values for Q   and Q  obtained from closed form expressions,
c                   0       1
c                           -m/2
c  scaled by (2m-1)!!(x*x+1).    For x >= 0.1: It is calculated using
c                                                n     n-1      n-2
c  forward recursion on the expression relating Q  to Q    and Q  .
c                                                m     m        m
c                          0
c  The starting value for Q  is obtained by reverse recursion on
c                          m
c                           0     0        0
c  the expression relating Q  to Q    and Q    and normalizing the
c                           n     n+1      n+2
c                                    0                  1
c  result using the known value for Q  . Similarly for Q  .
c                                    0                  m
c
        iqlterm=0
        if(m.ne.0) go to 10
        ql(1)=q00
        qr(1)=x+1.0e0_knd/q00
        go to 30
10      q11=-x1d*q00-x
        if(m.ne.1) go to 20
        ql(1)=q11
        qr(3)=3.0e0_knd*x-2.0e0_knd/q11
        go to 30
20      continue
        if(x.ge.0.1e0_knd) go to 25
        qlow=q00
        qupp=q11
        if(iflag.eq.1) qmmm1=3.0e0_knd*x*q11-2.0e0_knd
          do jm=1,m-1
          ajm=jm
          tjm=real(jm+jm,knd)/real(jm+jm+1,knd)
          ql(1)=(-x1d-tjm)*qupp-tjm*x1d*qlow
            if(iflag.eq.1) then
            qmmm1=-((ajm+ajm+2)/(ajm+ajm+1))*x1d*qmmm1+x*ql(1)
            end if
          qlow=qupp
          qupp=ql(1)
          end do
        if(iflag.eq.1) qr(m+m+1)=qmmm1/ql(1)
        go to 30
25      limm=-int(ndec/(log10(xsqr-x)))
        qupp0=0.0e0_knd
        qmid0=1.0e0_knd
        qupp1=0.0e0_knd
        qmid1=1.0e0_knd
          do jn=limm+m,m,-1
          rin1=jn+1
          rin2=jn+2
          rin3=jn+jn+3
          if(2*(jn/2).ne.jn) rin3=-rin3
          qlow0=(rin3*x*qmid0-rin2*qupp0)/rin1
          qlow1=(rin3*x*qmid1-rin1*qupp1)/rin2
          qupp0=qmid0
          qupp1=qmid1
          qmid0=qlow0
          qmid1=qlow1
          end do
        q0m=qlow0
        q1m=qlow1
        iqlterm=0
          do jn=m-1,1,-1
          rin1=jn+1
          rin2=jn+2
          rin3=jn+jn+3
          if(2*(jn/2).ne.jn) rin3=-rin3
          qlow0=(rin3*x*qmid0-rin2*qupp0)/rin1
          qlow1=(rin3*x*qmid1-rin1*qupp1)/rin2
            if(abs(qlow0).gt.teste) then
            qlow0=qlow0*testeo
            qlow1=qlow1*testeo
            qmid0=qmid0*testeo
            qmid1=qmid1*testeo
            iqlterm=iqlterm+nfac
            end if
          qupp0=qmid0
          qupp1=qmid1
          qmid0=qlow0
          qmid1=qlow1
          end do
        rin3=qlow0
        qlow0=3.0e0_knd*x*qmid0-2.0e0_knd*qupp0
        q1m=(q11/qlow1)*q1m
        q0m=(q00/qlow0)*q0m
        qlow=q0m
        qmid=q1m
          do j=1,m-1
          rin1=j+j
          rin2=j+j+1
          rin3=(m+j)*(m-j+1)
          rin4=(j+j+1)*(j+j-1)
          qupp=-(rin1*x/rin2)*qmid+(rin3*x1d/rin4)*qlow
          qlow=qmid
          qmid=qupp
          end do
          if(2*(m/4).ne.(m/2)) qupp=-qupp
          ql(1)=qupp
30      iql(1)=int(log10(abs(ql(1))))
        ql(1)=ql(1)*(ten**(-iql(1)))
        iql(1)=iql(1)-iqlterm
c
c                        m    m
c  the ratios qr(k+m) = Q  / Q    , k=m+limq to k=m+1 are calculated
c                        k    k-1
c
c  for x > 0.001 using backward recursion from qr(mxqr+2m) =
c
c   m        m
c  Q      / Q       = x - sqrt(x*x+1), where the last quantity
c   mxqr+m   mxqr+m-1
c
c  is the asymptotic limit of the ratio as mxqr approaches infinity.
c  Otherwise the ratios are calculated using forward recursion from
c  the ratio qr(m+m+1)
c
        if(iflag.eq.1) go to 40
        mxqr=limq-int(ndec/(log10(xsqr-x)))
        if(mxqr.lt.2*limq) mxqr=2*limq
        qupp=x-xsqr
          do jn=mxqr+m,limq+m+1,-1
          rin=jn
          qlow=-(rin+rm-1.0e0_knd)/(x*(rin+rin-1.0e0_knd)
     1         -(rin-rm)*qupp)
          qupp=qlow
          end do
        qr(limq+m+m)=qupp
          do jn=limq+m,m+2,-1
          rin=jn
          qr(jn+m-1)=-(rin+rm-1.0e0_knd)/(x*(rin+rin-1.0e0_knd)
     1               -(rin-rm)*qr(jn+m))
          end do
        go to 50
40      continue
          do jn=m+2,limq+m
          rin=jn
          qr(jn+m)=((rin+rin-1.0e0_knd)*x+(rin+rm-1.0e0_knd)/
     1             qr(jn+m-1))/(rin-rm)
          end do
50      continue
c
c                              m     m
c  calculate ratios qr(k+m) = Q   / Q   ,k=m-1 to k=-m+1 using
c                              k-1   k
c
c                                       m      m
c  backward recursion from qr(m+m-1) = Q    / Q     = x
c                                       m-2    m-1
c
        if(m.eq.0) go to 120
        qr(m+m-1)=x
        if(m.eq.1) go to 70
          do 60 jn=m-1,2-m,-1
          rin=jn
          qr(jn+m-1)=(x*(rin+rin-1.0e0_knd)
     1               +((rin-rm)/qr(jn+m)))/(rin+rm-1.0e0_knd)
          if(qr(jn+m-1).eq.0.0e0_knd) qr(jn+m-1)=dec
60        continue
70      continue
c
c                  m
c  calculation of Q    , m > 0 by forward division of qr ratios
c                  m-1
c
c                 m
c  starting with Q  calculated from its closed form expression.
c                 0
c
          if(2*(m/2).eq.m) then
          qm1=sin(rm*xang)
          qm0=qm1
          qdm0=rm*cos(rm*xang)/x1d
          else
          qm1=cos(rm*xang)
          qm0=qm1
          qdm0=-rm*sin(rm*xang)/x1d
          end if
        xfac=0.5e0_knd*rm*log10(x1d)
        ixfac=int(xfac)
        qm1=qm1*(ten**(xfac-ixfac))
        term=1.0e0_knd
        iterm=0
        if(m.lt.2) go to 90
          do jm=2,m
          ajm=jm
          term=term*(ajm-1.0e0_knd)/(ajm+ajm-1.0e0_knd)
            if(term.lt.testeo) then
            term=term*teste
            iterm=iterm-nfac
            end if
          end do
90      qm1=qm1*term
        jterm=int(log10(abs(qm1)))
        qm1=qm1*ten**(-jterm)
        iqm1=ixfac+iterm+jterm
          if(2*(m/2).eq.m) then
          if(2*(m/4).ne.m/2) qm1=-qm1
          else
          if(2*((m-1)/4).ne.(m-1)/2) qm1=-qm1
          end if
        if(m.lt.2) go to 110
          do jm=1,m-1
          qm1=qm1/qr(jm+m)
            if(abs(qm1).gt.teste) then
            qm1=qm1*testeo
            iqm1=iqm1+nfac
            end if
          end do
110     continue
        iterm=int(log10(abs(qm1)))
        qm1=qm1*ten**(-iterm)
        iqm1=iqm1+iterm
120     continue
c
c  calculation of ratios of the first derivatives of q with respect
c  to x for degrees >= m using the relationship:
c
c                  m    m      [kx]qr(k+m)+(k+m)
c     qdr(k+m) = Q'  / Q'   =  ----------------- , k=m+1 to k=m+lim
c                  k    k-1    [(k-m)]qr(k+m)-kx
c
c
          do jm=m+1,m+limq
          ajm=jm
          qdr(jm+m)=(ajm*x*qr(jm+m)+(ajm+rm))/((ajm-rm)*qr(jm+m)-
     1               ajm*x)
          end do
c
c                   m                       m
c  calculation of q'    from the value for q
c                   m-1                     m-1
c                       m
c  and calculation of q'   , k=0 to lnum-1, from the value
c                       m+k
c        m
c  for q'
c        m
c
        if(m.gt.0) go to 140
        qdl(1)=1.0e0_knd/x1d
        iqdl(1)=0
        qdm0=qdl(1)
        go to 150
140     qdm1=-rm*x*qm1/x1d
        iterm=int(log10(abs(qdm1)))
        qdm1=qdm1*(ten**(-iterm))
        iqdm1=iqm1+iterm
        qdl(1)=rm*(x*ql(1)-2.0e0_knd*qm1*(ten**(iqm1-iql(1))))/x1d
        iqdl(1)=iql(1)
150     continue
        m2m1=m+m-1
          do jl=2,lnum
          ql(jl)=ql(jl-1)*qr(m2m1+jl)
          iql(jl)=iql(jl-1)
            if(abs(ql(jl)).gt.teste) then
            ql(jl)=ql(jl)*testeo
            iql(jl)=iql(jl)+nfac
            end if
            if(abs(ql(jl)).lt.testeo) then
            ql(jl)=ql(jl)*teste
            iql(jl)=iql(jl)-nfac
            end if
          qdl(jl)=qdl(jl-1)*qdr(m2m1+jl)
          iqdl(jl)=iqdl(jl-1)
            if(abs(qdl(jl)).gt.teste) then
            qdl(jl)=qdl(jl)*testeo
            iqdl(jl)=iqdl(jl)+nfac
            end if
            if(abs(qdl(jl)).lt.testeo) then
            qdl(jl)=qdl(jl)*teste
            iqdl(jl)=iqdl(jl)-nfac
            end if
          end do
          do jl=1,lnum
          iterm=int(log10(abs(ql(jl))))
          ql(jl)=ql(jl)*ten**(-iterm)
          iql(jl)=iql(jl)+iterm
          iterm=int(log10(abs(qdl(jl))))
          qdl(jl)=qdl(jl)*ten**(-iterm)
          iqdl(jl)=iqdl(jl)+iterm
          end do
        termpq=rm*log10(xsqr)
        itermpq=int(termpq)
        termpq=ten**(termpq-itermpq)
c
c  Calculation of ratios of the first derivatives of q with respect
c  to x for degrees from 0 to m-1 via backward recursion using a
c  relation developed from the traditional recursion relations.
c  Here
c                 m      m
c     qdr(k+m) = Q'   / Q'  , k=m-1 to k=1
c                 k-1    k
c
c  The backward recursion is started with the value for qdr(m+m-1) =
c  (x*x*(m-1)-1)/(m*x).
c
        if(m.eq.0) go to 180
        qdr(m+m-1)=(x*x*(rm-1.0e0_knd)-1.0e0_knd)/(rm*x)
        if(qdr(m+m-1).eq.0.0e0_knd) qdr(m+m-1)=dec
        if(m.lt.3) go to 180
          do jn=m-1,2,-1
          rin=jn
          term=rin*rin*x1d-rmsq
            if(term.eq.0.0e0_knd) then
            xc=x*(1.0e0_knd+dec)
            term=rin*rin*(xc*xc+1.0e0_knd)-rmsq
            end if
          qdr(jn+m-1)=(x*(jn+jn-1)*(jn*(jn-1)*x1d-rmsq)+
     1                ((jn-m)*((jn-1)*(jn-1)*x1d-rmsq))/
     2                qdr(jn+m))/((jn+m-1)*term)
          if(qdr(jn+m-1).eq.0.0e0_knd) qdr(jn+m-1)=dec
          end do
180     continue
c
c  Calculation of ratios of the first derivative of q with respect
c  to x to the corresponding value for q for degrees from -m to m-1
c  Here
c                  m      m
c     qdqr(k+m) = Q'   / Q   , k=-m+1 to k=m
c                  k-1    k-1
c
c
        if(m.eq.0) go to 190
        qdqr(1)=((rm-1.0e0_knd)*x+(rm+rm-1.0e0_knd)/qr(1))/x1d
          do jn=2,m+m
          rin=jn
          qdqr(jn)=(x*(-rm+rin-1.0e0_knd)-(rin-1.0e0_knd)*qr(jn-1))/
     1              x1d
          end do
190     continue
        if(iflagl1.eq.0) go to 200
c
c  Modification of The ratios qr and qdr to obtain those used in the
c  Baber and Hasse expansion.
c                                     m    m
c  The resulting ratios are qr1(k) = Q  / Q    , k=1 to limq
c                                     k    k-1
c                 m    m
c  and qdr1(k) = Q  / Q    , k=1 to limq
c                 k    k-1
c
          do j=1,limq-m
          qr1(m+j)=qr(m+m+j)
          qdr1(m+j)=qdr(m+m+j)
          end do
          if(m.gt.1) then
            do j=1,m-1
            qr1(j)=-1.0e0_knd/qr(m+j)
            qdr1(j)=-1.0e0_knd/qdr(m+j)
            end do
          end if
          if(m.ne.0) then
          qr1(m)=-(ql(1)/qm1)*(ten**(iql(1)-iqm1))
          qdr1(m)=-(qdl(1)/qdm1)*(ten**(iqdl(1)-iqdm1))
          end if
          if(2*(m/2).eq.m) then
          if(2*(m/4).ne.m/2) qm0=-qm0
          if(2*(m/4).ne.m/2) qdm0=-qdm0
          else
          if(2*((m-1)/4).ne.(m-1)/2) qm0=-qm0
          if(2*((m-1)/4).ne.(m-1)/2) qdm0=-qdm0
          end if
200     continue
        return
        end
c
c
       subroutine pint(cc,m,lnum,x,limint,maxint,maxlp,maxmp,lim1max,
     1                 ndec,nex,ngqs,ngau,wg,xg,step1,nstep1,step2,
     2                 nstep2,step3,nstep3,limi,rpint1,rpint2,pint1,
     3                 pint2,pint3,pint4,norme,pnorm,ipnorm,coefme,
     4                 coefmo)
c
c  purpose:     To calculate integrals of the product of associated
c               Legendre functions and kernels containing spherical
c               Neumann functions and a window function. Four
c               different kernel functions are involved leading to
c               integrals of four different types. The integrals are
c               calculated using Gaussian quadrature.
c
c  parameters:
c
c     input:    cc     : complex c
c               m      : m
c               lnum   : number of l values desired
c               x      : x
c               limint : number of integrals of each of the four types
c                        required
c               maxint : dimension of the integral arrays
c               maxlp  : dimension of characteristic and exponent
c                        arrays of integrals
c               maxmp  : dimension of the spherical Neumann function
c                        array
c               lim1max: dimension of the spherical Bessel function
c                        array used in computing spherical Neumann
c                        functions
c               ndec   : number of decimal digits in real(knd)
c                        arithmetic
c               nex    : maximum exponent in real(knd) arithmetic
c               ngqs   : total number of steps
c               ngau   : order of Gaussian quadrature for substeps
c               wg     : ngau Gaussian quadrature weighting factors
c               xg     : corresponding Gaussian quadrature arguments
c               step1  : size of step 1
c               nstep1 : number of substeps in step1
c               step2  : size of step 2
c               nstep2 : number of substeps in step2
c               step3  : size of step 3
c               nstep3 : number of substeps in step3
c
c     output:   limi   : highest order of Legendre function for which
c                        integrals are calculated, since higher order
c                        integrals could overflow. series in subroutine
c                        r2int will be limited to order limi
c               rpint1 : array of ratios of successive integrals of
c                        the same parity of the first type (l-m even)
c                        or of the third type (l-m odd)
c               rpint2 : array of ratios of successive integrals of
c                        the same parity of the second type (l-m even)
c                        or of the fourth type (l-m odd)
c               pint1  : array of scaled values for the integrals of
c                        the first type
c               pint2  : array of scaled values for the integrals of
c                        the second type
c               pint3  : array of scaled values for the integrals of
c                        the third type
c               pint4  : array of scaled values for the integrals of
c                        the fourth type
c               norme  : scaling exponent for the spherical Neumann
c                        functions of order m
c               pnorm  : array of characteristics for the scaling
c                        factors used for the associated Legendre
c                        functions
c               ipnorm : array of exponents for the scaling factors
c                        used for the associated Legendre functions
c               coefme : coefficient used to multiply r2 to get one of
c                        the two contributions to r2d when l-m is even
c               coefmo : coefficient used to multiply r2 to get one of
c                        the two contributions to r2d when l-m is odd
c
c        use param
c
c  real(knd) scalars and arrays
        real(knd) amo2,an,argb,arn,bn,coef,coefme,coefmo,coefo,dec,
     1            emo2,etai,etaism1,etcoef1,etcoef2,factor,factor2,
     2            f2exp,rm,rn,sargb,step1,step2,step3,substep1,
     3            substep2,substep3,ten,term,tf2,term1,term2,test,
     4            teste,testeo,test1,test2,twom,twomi,x,x2,xsp1,zi,zl,
     5            zu
        real(knd) alpha(maxint),beta(maxint),p(maxint),pnorm(maxlp),
     1            step(ngqs),wg(ngau),xg(ngau)
c
c  complex(knd) scalars and arrays
        complex(knd) arg,cc,coef1,coef2,coef4,darg,darg2,sbes1,sbes2,
     1               sneuna,sneunb,sneunc,sneu1,sneu2,sneu3,ynorm1,
     2               ynorm2,ynorm3
        complex(knd) pint1(maxint),pint2(maxint),pint3(maxint),
     1               pint4(maxint),rpint1(maxint),rpint2(maxint),
     2               sbesf(lim1max)
c
c  integer array
        integer ipnorm(maxlp)
c
        ten=10.0e0_knd
        emo2=(0.5e0_knd)*m
        nfac=nex/3
        if(2*(nfac/2).ne.nfac) nfac=nfac-1
        teste=ten**(nfac)
        testeo=1.0e0_knd/teste
        test1=ten**(nex-3)
        ifac=nex-30
        factor=ten**(-ifac)
        test=1.0e0_knd/factor
        x2=x*x
        xsp1=x2+1.0e0_knd
        dec=ten**(-ndec-1)
        rm=m
        amo2=0.5e0_knd*rm
        lim=limint
        coefme=rm*x/xsp1
        coefmo=(rm*x2+xsp1)/(x*xsp1)
          substep1=step1/nstep1
            do j=1,nstep1
            step(j)=substep1
            end do
          substep2=step2/nstep2
            do j=1,nstep2
            step(nstep1+j)=substep2
            end do
          substep3=step3/nstep3
            do j=1,nstep3
            step(nstep1+nstep2+j)=substep3
            end do
c
c  calculation of scaling factors for the associated Legendre functions
	pnorm(1)=1.0e0_knd
	pnorm(2)=1.0e0_knd
	ipnorm(1)=0
	ipnorm(2)=0
        if(m.eq.0) go to 50
          do 40 n=1,m
          an=n+n
          bn=n+n+1
          pnorm(1)=pnorm(1)*an
          pnorm(2)=pnorm(2)*bn
          iterm1=int(log10(pnorm(1)))
          iterm2=int(log10(pnorm(2)))
          pnorm(1)=pnorm(1)*ten**(-iterm1)
          pnorm(2)=pnorm(2)*ten**(-iterm2)
          ipnorm(1)=ipnorm(1)+iterm1
          ipnorm(2)=ipnorm(2)+iterm2
40	  continue
50	twom=m+m
        pnorm(3)=pnorm(1)*(twom+2)/2
        iterm3=int(log10(pnorm(3)))
        pnorm(3)=pnorm(3)*ten**(-iterm3)
        ipnorm(3)=iterm3+ipnorm(1)
        if(lnum.lt.4) go to 70
          do 60 il=4,lnum,2
          pnorm(il)=pnorm(il-2)*(twom+il-1)/(il-1)
          pnorm(il+1)=pnorm(il-1)*(twom+il)/(il)
          iterm1=log10(pnorm(il))
          iterm2=log10(pnorm(il+1))
          ipnorm(il)=ipnorm(il-2)+iterm1
          ipnorm(il+1)=ipnorm(il-1)+iterm2
          pnorm(il)=pnorm(il)*ten**(-iterm1)
          pnorm(il+1)=pnorm(il+1)*ten**(-iterm2)
60	  continue
70	continue
c
c  calculation of the coefficients in the recursion relation used
c  for the scaled associated Legendre functions
	alpha(1)=(twom+1.0e0_knd)*pnorm(1)/pnorm(2)
        alpha(1)=alpha(1)*ten**(ipnorm(1)-ipnorm(2))
        beta(1)=0.0e0_knd
        alpha(2)=(twom+3.0e0_knd)*pnorm(2)/(pnorm(3)*2.0e0_knd)
        alpha(2)=alpha(2)*ten**(ipnorm(2)-ipnorm(3))
        beta(2)=-(twom+1.0e0_knd)/(twom+2.0e0_knd)
          do 80 il=3,lim+2
          alpha(il)=alpha(il-2)*(twom+il-1)*(twom+il+il-1)*
     1    (il-2)/((il-1)*(twom+il)*(twom+il+il-5))
	  beta(il)=-(twom+il-1)/(twom+il)
80	  continue
c
          do 90 il=1,lim+1,2
          pint1(il)=(0.0e0_knd,0.0e0_knd)
          pint2(il)=(0.0e0_knd,0.0e0_knd)
          pint3(il+1)=(0.0e0_knd,0.0e0_knd)
          pint4(il+1)=(0.0e0_knd,0.0e0_knd)
90	  continue
c
c  calculation of the scaling exponents for the spherical Neumann
c  functions required for the four types of integrals
        twomi=1.0e0_knd
        if(m.eq.0) go to 110
          do 100 n=1,m
	  twomi=twomi*real(n+n-1,knd)/real(n+n,knd)
100       continue
110	continue
        arg=cc*sqrt(xsp1)
        darg=(1.0e0_knd,0.0e0_knd)/arg
        ynorm1=-cos(arg)*darg
        ynorm2=ynorm1*darg-sin(arg)*darg
        normy=0
          do 120 n=3,m+3
          rn=n+n-3
          arn=n+n-1
          ynorm3=-ynorm1+darg*rn*ynorm2
            if(abs(ynorm3).gt.teste) then
            ynorm3=ynorm3*testeo
            ynorm2=ynorm2*testeo
            normy=normy+nfac
            end if
          ynorm1=ynorm2
          ynorm2=ynorm3
120       continue
          iterm=int(log10(abs(ynorm3)))
          normy=normy+iterm
          norme=normy
c
c  Gaussian quadrature integration loops. first dividing integrand
c  into ngqs steps
            do 180 k=1,ngqs
              if(k.eq.1) then
              zl=0.0e0_knd
              zu=step(1)
              end if
              if(k.gt.1) then
              zl=zu
              zu=zl+step(k)
              end if
          etcoef1=(zl+zu)/2.0e0_knd
          etcoef2=step(k)/2.0e0_knd
          coef=step(k)
c
c  Gaussian quadrature integration over each step
            do 170 i=1,ngau
 	    zi=etcoef1+xg(i)*etcoef2
            etai=1.0e0_knd-zi
            etaism1=zi*(2.0e0_knd-zi)
            argb=x2+etaism1
            term2=xsp1*etaism1*etaism1/argb
            term1=sqrt(term2)
            sargb=sqrt(argb)
            coefo=1.0e0_knd/sargb
            arg=cc*sargb
            darg=(1.0e0_knd,0.0e0_knd)/arg
            sneu1=-cos(arg)*darg
            sneu2=term1*darg*(sneu1-sin(arg))
            norm=ifac
            lown=3
            limb=2*int(abs(real(arg))+abs(aimag(arg)))
            lim1=2*limb+20
            if(abs(aimag(arg)).lt.2.0e0_knd.or.m.lt.3) go to 130
            sbesf(lim1)=arg/real(lim1+lim1+1,knd)
              do n=lim1-1,1,-1
              rn=real(n+n+1,knd)
              sbesf(n)=1.0e0_knd/(rn*darg-sbesf(n+1))
              end do
            darg2=darg*darg
            sneuna=sneu1
            sneunb=sneu2/term1
            sbes1=sbesf(1)*sin(arg)*darg
            if(limb.gt.m+2) limb=m+2
              do n=2,limb
              sbes2=sbesf(n)*sbes1
              sneunc=(sbes2*sneunb-darg2)/sbes1
              sbes1=sbes2
                if(n.ne.limb) then
                sneuna=sneunb
                sneunb=sneunc
                end if
              end do
              nf2=int(real(limb-2,knd)*abs(log10(term1)))/nfac+1
              f2exp=real(limb-2,knd)/real(nf2,knd)
              factor2=1.0e0_knd
              tf2=term1**(f2exp)
                do ll=1,nf2
                factor2=factor2*tf2
                  if(factor2.gt.teste) then
                  factor2=factor2*testeo
                  norm=norm+nfac
                  end if
                  if(factor2.lt.testeo) then
                  factor2=factor2*teste
                  norm=norm-nfac
                  end if
                end do
            sneu1=sneuna*factor2
            sneu2=term1*sneunb*factor2
            sneu3=term2*sneunc*factor2
              if(limb.eq.m+2) then
              go to 150
              else

              sneu2=sneu3
              lown=limb+2
              end if
130         continue
              do n=lown,m+3
              rn=real(n+n-3,knd)
              sneu3=-term2*sneu1+term1*darg*rn*sneu2
              if(n.eq.m+3) go to 150
                if(abs(sneu3).gt.teste) then
                sneu3=sneu3*testeo
                sneu2=sneu2*testeo
                norm=norm+nfac
                end if
                if(abs(sneu3).lt.testeo) then
                sneu3=sneu3*teste
                sneu2=sneu2*teste
                norm=norm-nfac
                end if
              sneu1=sneu2
              sneu2=sneu3
140           end do
150         continue
            iterm=int(log10(abs(sneu3)))
            term=ten**(-iterm)
            sneu3=sneu3*term
            sneu2=sneu2*term
            sneu1=sneu1*term
            norm=norm+iterm
            iexp=-norme+norm
            if(iexp.le.-nex+10-ifac) go to 160
              if(iexp.gt.-nex+10) then
              term=ten**iexp
              sneu3=sneu3*term
              sneu2=sneu2*term
              sneu1=sneu1*term
              coef1=coef*sneu1*wg(i)
              coef2=(coef/term1)*coefo*sneu2*wg(i)
              coef4=(coef/term2)*coefo*coefo*etai*sneu3*wg(i)
              iflag=0
              else
              term=ten**(iexp+ifac)
              sneu3=sneu3*term
              sneu2=sneu2*term
              sneu1=sneu1*term
              iflag=1
              end if
            p(1)=twomi*factor
            p(2)=alpha(1)*etai*p(1)
              if(iflag.eq.0) then
              pint1(1)=pint1(1)+coef1*p(1)
              pint2(1)=pint2(1)+coef2*p(1)
              pint4(2)=pint4(2)+coef4*p(2)
              end if
              do il=2,limi,2
              p(il+1)=alpha(il)*etai*p(il)+beta(il)*p(il-1)
              p(il+2)=alpha(il+1)*etai*p(il+1)+beta(il+1)*p(il)
                if(iflag.eq.0) then
                pint1(il+1)=pint1(il+1)+coef1*p(il+1)
                pint2(il+1)=pint2(il+1)+coef2*p(il+1)
                pint4(il+2)=pint4(il+2)+coef4*p(il+2)
                  if(abs(pint1(il+1)).gt.test1.or.abs(pint2(il+1)).gt.
     1             test1.or.abs(pint4(il+2)).gt.test1) then
                  limi=il
                  go to 160
                  end if
                end if
                if(abs(p(il+2)).gt.test) then
                p(il+1)=p(il+1)*factor
                p(il+2)=p(il+2)*factor
                  if(iflag.eq.0) then
                  coef1=coef1*test
                  coef2=coef2*test
                  coef4=coef4*test
                  else
                  coef1=coef*sneu1*wg(i)
                  coef2=(coef/term1)*coefo*sneu2*wg(i)
                  coef4=(coef/term2)*coefo*coefo*etai*sneu3*wg(i)
                  iflag=0
                  end if
                end if
              end do
160         continue
170         continue
180	  continue
190	continue
          do 200 il=1,limi-1,2
          pint3(il+1)=(pint2(il+2)-beta(il+1)*pint2(il))
     1                /alpha(il+1)
200	  continue
210     continue
c        write(40,220) ngau,substep1,substep2,substep3,limi
c220     format(15x,'Gauss quad. order =',i5,'; step sizes = ',f12.10,
c     1         ', ',f12.10,', 'f12.10,'.',/,15x,'integrals for ',i5,
c     2         ' lowest order Legendre functions will be used for r2.')
c
c  calculation of ratios of integrals for ease in compution of r2 and
c  r2d in subroutine r2int
        rpint1(1)=(0.0e0_knd,0.0e0_knd)
        rpint1(2)=(0.0e0_knd,0.0e0_knd)
        rpint2(1)=(0.0e0_knd,0.0e0_knd)
        rpint2(2)=(0.0e0_knd,0.0e0_knd)
          do 240 il=3,limi,2
          rpint1(il)=(pint1(il)*(twom+il-1))/(pint1(il-2)*(il-1))
          rpint2(il)=(pint2(il)*(twom+il-1))/(pint2(il-2)*(il-1))
          rpint1(il+1)=(pint3(il+1)*(twom+il))/(pint3(il-1)*(il))
          rpint2(il+1)=(pint4(il+1)*(twom+il))/(pint4(il-1)*(il))
240	  continue
        limi=limi-1
        return
        end
c
c
        subroutine sphbes (cc,x,limj,maxj,maxlp,sbesf,sbesdf,sbesn,
     1                     ibese,sbesdr)
c
c  purpose:     To calculate ratios of successive first kind spherical
c               Bessel functions of the same parity for given c and x.
c               to calculate corresponding ratios of their first
c               derivatives. To calculate the characteristics and
c               exponents of both the Bessel functions of the first
c               kind and their first derivatives.
c
c  parameters:
c
c     input:    cc     : complex c
c               x      : argument of spherical Bessel functions
c               limj   : the number of spherical Bessel function
c                        ratios calculated for given lnum,c,ndec,
c                        and maximum m desired
c               maxj   : dimension of sbesf vector
c               maxlp  : the number of scale factors that are
c                        calculated
c
c     output:   sbesf  : ratios of successive first kind spherical
c                        Bessel functions of the same parity
c               sbesdf : ratios of first derivatives of successive
c                        first kind spherical Bessel functions of the
c                        same parity
c               sbesn  : characteristics for the spherical
c                        Bessel functions
c               ibese  : exponents for the spherical
c                        Bessel functions
c               sbesdr : ratios of first derivatives of spherical Bessel
c                        functions to the corresponding spherical
c                        spherical functions
c
c        use param
c
c  real(knd) scalars
        real(knd) adj,ci,cm,cr,rn,ten,x
c
c  complex(knd) scalars and arrays
        complex(knd) cc,cx,stemp0,stemp1
        complex(knd) sbesdf(maxj),sbesdr(maxlp),sbesf(maxj),
     1               sbesn(maxlp)
c
c  integer array
        integer ibese(maxlp)
c
        ten=10.0e0_knd
        cx=cc*x
        ci=aimag(cc)
        cr=real(cc)
        cm=abs(cc)
        lim1=2*cm*x+20
c
c  compute first kind Bessel function ratios
c        sbesf(k)= j(n=k,c*x)/j(n=k-1,c*x)
c        sbesn(k)= (j(n=k-1),c*x))*10.0e0_knd**(-ibese(k))
c
        if(int(abs(cx)).lt.limj.or.aimag(cc).ne.0.0e0_knd) go to 20
c
c  for c*x >= limj and aimag(cc) = 0, use forward recursion to
c  get fcn. ratios:
c       j(n+1,c*x)/j(n,c*x)=(2*n+1)/(c*x)-1/(j(n,c*x)/j(n-1,c*x))
c
        stemp0=sin(cx)/cx
        sbesf(1)=(stemp0/cx-cos(cx)/cx)/stemp0
          do 10 n=1,limj-1
          rn=n+n+1
          sbesf(n+1)=(rn/cx)-(1.0e0_knd/sbesf(n))
10        continue
        go to 60
20      continue
c
c  for c*x < lim or aimag(cc) unequal to zero, use backward recursion
c  to get fcn. ratios:
c       j(n,c*x)/j(n-1,c*x) = 1/( (2*n+1)/(c*x) - j(n+1,c*x)/j(n,c*x) )
c
        stemp0=0.0e0_knd
        if(lim1.le.limj) go to 40
          do 30 n=lim1,limj,-1
          rn=n+n+1
          stemp1=1.0e0_knd/(rn/cx-stemp0)
          stemp0=stemp1
30        continue
40      sbesf(limj)=stemp0
          do 50 n=limj-1,1,-1
          rn=n+n+1
          sbesf(n)=1.0e0_knd/(rn/cx-sbesf(n+1))
50        continue
60      continue
c
c  for all c*x, calculate the amplitude and sign scale
c  factors by forward operation on the Bessel function
c  ratios.
        nex=range(rn)-1
          if(0.43e0_knd*x*ci.gt.nex/3) then
          adj=x*ci/log(10.0e0_knd)
          iadj=int(adj)
          sbesn(1)=(0.0e0_knd,1.0e0_knd)*
     1               exp((0.0e0_knd,-1.0e0_knd)*cr*x)/(2.0e0_knd*cx)
          sbesn(1)=sbesn(1)*(10.0e0_knd**(adj-iadj))
          iterm=int(log10(abs(sbesn(1))))
          sbesn(1)=sbesn(1)*(10.0e0_knd**(-iterm))
          ibese(1)=iadj+iterm
          sbesn(2)=sbesn(1)/cx+(0.0e0_knd,1.0e0_knd)*sbesn(1)
          iterm=log10(abs(sbesn(2)))
          sbesn(2)=sbesn(2)*(10.0e0_knd**(-iterm))
          ibese(2)=ibese(1)+iterm
          else
          stemp0=sin(cx)/cx
          stemp1=stemp0/cx-cos(cx)/cx
          ibese(1)=log10(abs(stemp0))
          sbesn(1)=stemp0*ten**(-ibese(1))
          if(abs(sin(cx)).lt.0.5e0_knd.and.abs(cx).gt.1.0e0_knd)
     1       go to 70
          sbesn(2)=sbesn(1)*sbesf(1)
          ibese(2)=log10(abs(sbesn(2)))
          sbesn(2)=sbesn(2)*ten**(-ibese(2))
          ibese(2)=ibese(2)+ibese(1)
          go to 80
70        ibese(2)=log10(abs(stemp1))
          sbesn(2)=stemp1*ten**(-ibese(2))
          sbesf(1)=stemp1/stemp0
80        continue
          end if
           do 90 n=3,maxlp
          sbesn(n)=sbesn(n-1)*sbesf(n-1)
          ibese(n)=log10(abs(sbesn(n)))
          sbesn(n)=sbesn(n)*ten**(-ibese(n))
          ibese(n)=ibese(n)+ibese(n-1)
90        continue
c
c  calculate the ratios of the first derivatives of successive
c  Bessel functions using corresponding function ratios
          do 100 n=1,limj
          rn=n-1
          sbesdf(n)=(cx-(rn+2.0e0_knd)*sbesf(n))/(rn-cx*sbesf(n))
100       continue
c
c  calculate the ratios of the first derivative to the corresponding
c  spherical Bessel function
          do 110 n=1,maxlp
          rn=n-1
          sbesdr(n)=(rn/cx)-sbesf(n)
110       continue
c
c  calculate the ratios of successive functions and derivatives
c  of the same parity
          do 120 n=limj,2,-1
          sbesf(n)=sbesf(n-1)*sbesf(n)
          sbesdf(n)=sbesdf(n-1)*sbesdf(n)
120       continue
        return
        end
c
c
        subroutine sphneu (cc,x,limn,maxn,maxlp,limbes,sneuf,sneun,
     1                     ineue,sneudf,sneudr)
c
c  purpose:     to calculate ratios of spherical Neumann functions
c               and ratios of their first derivatives for given c and x.
c               to calculate the Neumann function characteristics
c               and exponents. to calculate ratios of the first
c               derivatives of the corresponding Neumann functions.
c
c  parameters:
c
c     input:    cc     : complex c
c               x      : argument of spherical Neumann functions
c               limn   : the number of spherical Neumann function
c                        ratios calculated for given lnum,c,ndec,
c                        and maximum m desired
c               maxn   : dimension of sneuf and sneudf arrays
c               maxlp  : the number of values of scale factors
c                        that are calculated
c               limbes : dimension of spherical Bessel function ratios
c                        calculated for use in calculating Neumann
c                        functions
c
c     output:   sneuf  : ratios of successive spherical Neumann
c                        functions of the same parity
c               sneun  : characteristic for the spherical
c                        Neumann functions
c               ineue  : exponent for the spherical
c                        Neumann functions
c               sneudf : ratios of first derivatives of successive
c                        spherical Neumann functions of the same parity
c               sneudr : ratios of first derivatives of spherical
c                        Neumann functions to the corresponding
c                        function
c
c        use param
c
c  real(knd) scalars
        real(knd) adj,ci,cr,rn,rnn,test,x
c
c  complex scalars and arrays
        complex(knd) cc,cx,cx2,sbes1,sbes2,stemp0,stemp1
        complex(knd) sneudf(maxn),sneudr(maxlp),sneuf(maxn),
     1               sneun(maxn),sbesf(limbes)
c
c  integer arrays
        dimension ineue(maxn)
c
c  compute first kind ratios of Neumann functions and ratios
c  of their first derivatives
c
c        sneuf(k)=y(n=k,c*x)/y(n=k-2,c*x)
c        sneun(k)=(y(n=k-1),c*x)*10.e0_knd**(-ineue(k))
c        sneudf(k)=y'(n=k,c*x)/y'(n=k-2,c*x)
c        sneudr(k)=(y'(n=k-1),c*x)/y(n=k-1),c*x))
c
c  calculate j ratios below turning point by backward recursion
c
c       j(n,c*x)/j(n-1,c*x) = 1/( (2*n+1)/(c*x) - j(n+1,c*x)/j(n,c*x) )
c
        cx=cc*x
        cx2=cx*cx
        cr=real(cc)
        ci=aimag(cc)
        test=1.0e0_knd/abs(cx)
        itest=-int(log10(abs(cx)))
        limb=2*int(abs(real(cx))+abs(aimag(cx)))
        lim1=2*limb+20
        if(limb.lt.2) limb=2
        sbesf(lim1)=cx/(lim1+lim1+1)
          do 10 n=lim1-1,1,-1
          rn=real(n+n+1)
          sbesf(n)=1.0e0_knd/(rn/cx-sbesf(n+1))
10      continue
c
c  use relation with j's to compute y's from order zero
c  to the turning point. compute derivative ratios
c  at same time.
c
        nex=range(rn)-1
          if(0.43e0_knd*x*ci.le.nex/3) go to 30
          ndec=precision(rn)
          adj=x*ci/log(10.0e0_knd)
          iadj=int(adj)
          sbes1=(0.0e0_knd,1.0e0_knd)*
     1            exp((0.0e0_knd,-1.0e0_knd)*cr*x)/(2.0e0_knd*cx)
          sbes1=sbes1*(10.0e0_knd**(adj-iadj))
          iterm=int(log10(abs(sbes1)))
          ibes1=iadj+iterm
          sbes1=sbes1*(10.0e0_knd**(-iterm))
          sneun(1)=(0.0e0_knd,1.0e0_knd)*sbes1
          ineue(1)=ibes1
          sbes2=sbes1/cx+(0.0e0_knd,1.0e0_knd)*sbes1
          sneuf(1)=sbes2/sbes1
          iterm=int(log10(abs(sbes2)))
          ibes2=ibes1+iterm
          sbes2=sbes2*(10.0e0_knd**(-iterm))
          sneun(2)=(0.0e0_knd,1.0e0_knd)*sbes2
          ineue(2)=ibes2
        sneudf(1)=-(cx-2.0e0_knd*sneuf(1))/(cx*sneuf(1))
        if(limb.gt.limn) limb=limn
        j=1
        sbes1=sbes2
        ibes1=ibes2
          do n=2,limb
          j=j+1
          rn=real(n)
          rnn=real(n-1)
          sbes2=sbesf(n)*sbes1
          iterm=int(log10(abs(sbes2)))
          sbes2=sbes2*(10.0e0**(-iterm))
          ibes2=ibes1+iterm
          sneun(n+1)=sbesf(n)*sneun(n)
          if(ibes1+ineue(n).lt.ndec+10)
     1    sneun(n+1)=sneun(n+1)-(10.0e0_knd**(-ibes1-ineue(n)))/
     2               (cx2*sbes1)
          sneuf(n)=sneun(n+1)/sneun(n)
          iterm=int(log10(abs(sneun(n+1))))
          sneun(n+1)=sneun(n+1)*(10.0e0_knd**(-iterm))
          ineue(n+1)=ineue(n)+iterm
          sneudf(n)=(cx-(rn+1.0e0_knd)*sneuf(n))/(rnn-cx*sneuf(n))
          if(ibes1+ibes2.lt.itest) go to 20
          sbes1=sbes2
          ibes1=ibes2
          end do
20      limb=max(j,2)
        go to 70
30      continue
        stemp0=-cos(cx)/cx
        stemp1=(stemp0-sin(cx))/cx
        sneuf(1)=stemp1/stemp0
        sneun(1)=stemp0
        sneun(2)=stemp1
        sbes1=sbesf(1)*sin(cx)/cx
        sneudf(1)=-(cx-2.0e0_knd*sneuf(1))/(cx*sneuf(1))
        if(limb.gt.limn) limb=limn
        j=1
          do 40 n=2,limb
          j=j+1
          rn=real(n)
          rnn=real(n-1)
          sbes2=sbesf(n)*sbes1
          sneun(n+1)=sbesf(n)*sneun(n)-1.0e0_knd/(cx2*sbes1)
          sneuf(n)=sneun(n+1)/sneun(n)
          sneudf(n)=(cx-(rn+1.0e0_knd)*sneuf(n))/(rnn-cx*sneuf(n))
          if(abs(sbes1+sbes2).lt.test) go to 50
          sbes1=sbes2
40        continue
50      limb=max(j,2)
c
c  calculate characteristics and exponents for the Neumann functions
c  up to and including the turning point
        ineue(1)=int(log10(abs(stemp0)))
        sneun(1)=stemp0*10.0e0_knd**(-ineue(1))
        ineue(2)=int(log10(abs(stemp1)))
        sneun(2)=stemp1*10.0e0_knd**(-ineue(2))
          do 60 n=3,limb
          ineue(n)=int(log10(abs(sneun(n))))
          sneun(n)=sneun(n)*10.0e0_knd**(-ineue(n))
60        continue
70      continue
c
c  use forward recursion from breakpoint to compute function ratios
c
c       y(n+1,c*x)/y(n,c*x)=(2*n+1)/(c*x)-1/(y(n,c*x)/y(n-1,c*x))
c
c  compute derivative ratios at same time using function ratios.
        if(limb.eq.limn) go to 90
          do 80 n=limb+1,limn
          rn=qfloat(n-1)
          rnn=qfloat(n+n-1)
          sneuf(n)=rnn/cx-1.0e0_knd/sneuf(n-1)
          sneudf(n)=(cx-(rn+2.0e0_knd)*sneuf(n))/(rn-cx*sneuf(n))
80        continue
90      continue
          sneuf(limn+1)=(0.0e0_knd,0.0e0_knd)
          sneuf(limn+2)=(0.0e0_knd,0.0e0_knd)
          sneudf(limn+1)=(0.0e0_knd,0.0e0_knd)
          sneudf(limn+2)=(0.0e0_knd,0.0e0_knd)
c
c  calculate the characteristics and exponents for Neumann
c  functions beyond the turning point by forward operation
c  on the Neumann function ratios:
        if(limb+1.gt.maxlp) go to 110
          do 100 n=limb+1,maxlp
          sneun(n)=sneun(n-1)*sneuf(n-1)
          ineue(n)=int(log10(abs(sneun(n))))
          sneun(n)=sneun(n)*10.0e0_knd**(-ineue(n))
          ineue(n)=ineue(n)+ineue(n-1)
100       continue
110     continue
c
c  calculate the ratios of the first derivatives to the corresponding
c  spherical Neumann functions
          do 120 n=1,maxlp
          rn=real(n-1)
          sneudr(n)=(rn/cx)-sneuf(n)
120       continue
c
c  calculate the ratios of successive functions and derivatives
c  of the same parity
          do 130 n=limn,2,-1
          sneuf(n)=sneuf(n-1)*sneuf(n)
          sneudf(n)=sneudf(n-1)*sneudf(n)
130       continue
        return
        end
c
c
	  end module vb_oblate