      module vb_prolate
	  use regime
	  use elemfunc

	  private

	  public cprofcn

	  contains

	  subroutine cprofcn(cc,m,lnum,ioprad,x1,r1c,ir1e,r1dc,ir1de,r2c,&
                   ir2e,r2dc,ir2de,iopang,iopnorm,narg,arg,&
                   s1,is1e,s1d,is1de, enr,maxd)
!
!  subroutine version of the fortran program cprofcn101
!
!      September 2015
!
!  developed by arnie lee van buren
!  www.mathieuandspheroidalwavefunctions.com
!
!  purpose:     To calculate the first and second kind prolate
!               radial functions and their first derivatives
!               with respect to x for a range of values of l
!               for specified values of m, x and complex c.
!               To calculate the first kind oblate angular
!               functions and their first derivatives with
!               respect to eta for a range of values of l
!               and eta for specified values of m and complex
!               c. Here c = qreal(c) + iqimag(c), where qreal(c)
!               is positive and qimag(c) can be positive or
!               negative. [For convenience in the following
!               descriptive comments about coblfcn, it is
!               assumed that qimag(c) is positive]
!
!
!  Csubpro is written in fortran 90. It uses 128-bit arithmetic
!  (real*16 and complex*32) to maximize the accuracy. It was
!  developed around profcn which obtains spheroidal function
!  values for real values of the size parameter c (=kd/2, where
!  k is the wavenumber and d is the interfocal distance of the
!  elliptical cross-sectionse of the prolate spheroid. Csubpro
!  provides function values for c complex [= qreal(c) + iqimag(c)],
!  where the imaginary part qimag(c) often accounts for losses
!  in wave propagation. Accurate results can be obtained over large
!  parameter ranges with qimag(c) as large as and often larger than
!  qreal(c). The resulting accuracy is comparable to that obtained
!  using profcn [with c equal to qreal(c)] when qimag(c) is no
!  greater than about 20.
!
!  Even though csubpro provides accurate results over wide parameter
!  ranges, there are regions where the accuracy of the resulting
!  radial functions of the second kind and their first derivatives
!  may fall below 10 digits even when using real*16 arithmetic. use
!  of real*8 and complex*16 arithmetic would widen these regions. For
!  profcn where c is real, these regions tend to occur where x is
!  between about 1.001 and 1.1, c is greater than about 1000, m is
!  greater than about 50, and l-m is intermediate. For c complex the
!  same behavior is expected. Other ranges of similarly reduced accuracy
!  may be found when c is complex, expecially when the imaginary part of
!  c has a large magnitude.
!
!  The minimum number of accurate digits desired for the radial
!  spheroidal functions of the second kind is an integer minacc.
!  Both the radial functions of the first kind and the angular
!  functions are computed to ndec digits. Minacc is set below
!  to a value of 16 (equivalent to real*8 precision in most
!  cases). However, if desired, accuracies of 25 digits or more can
!  often be obtained. Ten digits of accuracy are adequate for many
!  applications. In this case setting minacc equal to 10 may be
!  appropriate for regions where profcn may not be able to provide
!  accuracies as great as 16 decimal digits. Setting minacc larger
!  than needed will likely result in increased computation time. This
!  increase can be considerable for the difficult regions discussed
!  above. Setting minacc less than 5 will lead to problems in running
!  profcn due to program logic.
!
!  Angular functions values for some values of eta have reduced accuracy
!  when qreal(c) and/or the magnitude of qimag(c) are nonsmall. The
!  calculated values tend to be less accurate the larger qreal(c), the
!  smaller the value of l - m (for values less than the 2*qreal(c)/pi),
!  and the closer eta is to unity (i.e., the closer theta is to zero
!  degrees). They also tend to be less accurate the larger the magnitude
!  of qimag(c). This decreased accuracy occurs for values of l-m near
!  the breakpoint 2*qreal(c)/pi and for values of eta near zero. The
!  loss of accuracy (in decimal digits) of the angular functions is due
!  to subtraction error and is accompanied by a proportional decrease in
!  the magnitude of the angular function relative to its corresponding
!  associated Legendre function. The decrease in magnitude almost always
!  results in a corresponding reduction in the magnitude of their
!  contribution to the solution of physical problems involving oblate
!  spheroidal functions. Thus the lower accuracy in some of the angular
!  functions almost always has negligible impact on calculated solutions.
!
!  The eigenvalues for c real are traditionally ordered in increasing
!  numerical value to match the behavior of the corresponding Legendre
!  functions. When qimag(c) is small, the eigenvalues can still be
!  ordered in increasing real part. When qimag(c) is greater than about
!  3, however, the ordering of the eigenvalues is more complicated
!  Some of the eigenvalues are similar to low order oblate eigenvalues.
!  In the oblate case when c is real and non-small, neighboring low
!  order oblate eigenvalues beginning with l = m and l = m + 1 are
!  nearly identical. The larger the value of qimag(c), the greater the
!  number of paired eigenvalues. These eigenvalues do not fit easily
!  into the csubpro eigenvalue sequence. The existence of oblate-like
!  eigenvalues is apparently related to the fact that the recursion
!  relation for the oblate angular function expansion coefficients
!  can be obtained from the corresponding prolate recursion relation
!  by replacing c with -ic. Accurate estimates of these oblate-like
!  eigenvalues are given by the standard asymptotic approximations
!  for oblate eigenvalues with c replaced by -i times the prolate value
!  for c, i.e., by qimag(c)-iqreal(c).[If qimag(c) is chosen negative,
!  c is replaced by -qimag(c)+iqreal(c).] The oblate-like eigenvalues
!  are identified by their close numerical agreement with these
!  estimates. I include all of the eigenvalues where both an even and
!  an odd order value are within 1% of an asymptotic estimate. I now
!  order the selected eigenvalues by placing the two paired values with
!  the largest imaginary part in the center of the sequence. I then
!  arrange the others so that their imaginary values decreases to both
!  sides of the center with the two smallest imaginary values on the
!  two ends of the sequence. I make sure that even and odd eigenvalues
!  are interlaced. I then order the eigenvalues that were not oblate-
!  like in increasing real part. When ordered this way, their imaginary
!  part increases with increasing order until a maximum value is
!  obtained. This maximum value is very close to the imaginary part of
!  the two smallest oblate-like eigenvalues. I then insert the oblate-
!  like eigenvalues to one side or the other of this maximum value,
!  making sure the resulting sequence have interlaced even and odd
!  eigenvalues. This choice for ordering the eigenvalues is somewhat
!  arbitrary, but it appears to be a reasonable one.
!
!  Csubpro is designed around the number of decimal digits ndec and
!  the maximum exponent nex available in 128-bit arithmetic on the
!  user's computer. They have been set to 33 and 4930 respectively.
!  If different values are appropriate, changes can be made in
!  the two statements specifiying ndec and nex that appear below
!  following the first parameter declaration statements.
!
!  Csubpro is written to avoid overflow and underflow and can provide
!  radial function values far outside the dynamic range of real*16
!  arithmetic. It does this by expressing each radial function in
!  scientific notation and giving it in two parts. The first part,
!  called the characteristic, is a complex*32 number whose magnitude
!  ranges from unity to just less than 10.0q0. It contains 33 decimal
!  digit values for both the real and the imaginary parts of the
!  function value. The second part is the integer exponent (power 10).
!  For example, the radial function of the first kind is given by r1c
!  and ir1e. The corresponding complex*32 radial function value is then
!  equal to r1c*(10.0q0**ir1e). The radial function values can fall
!  outside the dynamic range of complex*32 arithmetic at large orders,
!  especially when c has a small magnitude. Here the functions of the
!  first kind underflow and the corresponding functions of the second
!  kind overflow.
!
!
!    Input and output parameters from the subroutine call statement are
!    defined below:
!
!          cc     : desired complex value of the size parameter
!                   (= kd/2, where k = complex wavenumber and d =
!                   interfocal length) (complex*32)
!
!          m      : desired value for m (integer)
!
!          lnum   : number of values desired for l (integer)
!
!          ioprad : (integer)
!                 : =0 if radial functions are not computed
!                 : =2 if radial functions of both kinds and
!                      their first derivatives are computed
!
!          x1     : value of the radial coordinate x - 1.0q0
!                  (a nominal value of 10.0q0 can be entered
!                  for x1 if ioprad = 0) (real*16)
!
!          r1c   :  complex*32 vectors of length lnum containing the
!          r1dc     characteristics for the radial functions of the
!                   first kind r1 and their first derivatives
!
!          ir1e   : integer vectors of length lnum containing the
!          ir1de    exponents corresponding to r1c and r1dc
!
!          r2c    : complex*32 vectors of length lnum containing the
!          r2dc     characteristics for the radial functions of the
!                   second kind r2 and their first derivatives
!
!          ir2e   : integer vectors of length lnum containing the
!          ir2de    exponents corresponding to r2c and r2dc
!
!          iopang : (integer)
!                 : =0 if angular functions are not computed
!                 : =1 if angular functions of the first kind
!                      are computed
!                 : =2 if angular functions of the first kind and
!                      their first derivatives are computed
!
!          iopnorm: (integer)
!                 : =0 if not scaled. The angular functions have
!                      the same norm as the corresponding associated
!                      legendre function [i.e., we use the Meixner
!                      -Schafke normalization scheme.] This norm
!                      becomes very large as m becomes very large.
!                 : =1 if angular functions of the first kind
!                      (and their first derivatives if computed)
!                      are scaled by the square root of the
!                      normalization of the corresponding
!                      associated Legendre function. The resulting
!                      scaled angular functions have unity norm.
!
!          narg   : number of values of the angular coordinate eta for
!                   which angular functions are calculated (integer)
!                   must be set equal to 1 if iopang=0
!
!          arg:     vector containing the values of eta for which
!                   angular functions are desired (real*16)
!
!          s1,s1d : complex*32 two-dimensional arrays s1(lnum,narg)
!                   and s1d(lnum,narg) that contain narg calculated
!                   angular functions and their first derivatives
!                   for each of the lnum values of l.
!                   For example, s1(10,1) is the angular function for
!                   l = m +10 -1 and the first value of eta given by
!                   arg(1)
!
!  Output data is written to four files: fort.20, fort.30, fort.40,
!  and fort.50
!
!   fort.20
!
!     This file contains values for all radial functions that have
!     been calculated.
!     The first line in the file contains the values for x, c (input as
!     the parameter cc), and m, formatted as follows (see statement 120
!     in subroutine main):
!
!                x      : e38.30
!                c      : e38.30, e38.30
!                m      : i5
!
!     Each subsequent set of 3 lines contains values for successive
!     orders l, beginning with l = m). The order l and the eigenvalue
!     are given in the first line. The next line contains the complex
!     value for the radial function of the first kind r1 followed by the
!     complex value for its first derivative r1d. Each complex value is
!     given by the characteristic of the real part, followed by the
!     characteristic of the imaginary part, followed by the integer
!     exponent that applies to both part. The third line contains the
!     corresponding values for the radial function of the second kind
!     r2 and its first derivative r2d. At the end of the third line is
!     the integer accuracy, equal to the estimated number of accurate
!     decimal digits in the radial functions. This is usually measured
!     using the wronskian. The values for r2 or r2d are set equal to
!     zero when the accuracy naccr is equal to zero.
!
!   fort.30
!
!     This file fort.30 contains values for all angular functions
!     (of the first kind) that have been calculated. Its first line
!     contains the values for c and m. The second nd line contains the
!     value for the first l (=m).
!     This is followed by a series of narg lines. Each line contains
!     a desired value of angle (ioparg = 0) or angular coordinate eta
!     (ioparg =1) followed by the corresponding angular functions Again
!     the functions are given by the real characteristic, followed by
!     the imaginary characteristic, followed by the common integer
!     If first derivatives of the angular functions are desired, the
!     next line give the corresponding first derivative values. [Only
!     16 decimal digits are written for each characteristic to fort.30.
!     If more digits are desired, the formats 1090 and 1100 can be
!     modified.]
!     The estimated accuracy of s1 (and that of s1d if it is requested
!     are given as a two digit integer following the angular function
!     or its first derivative. Each integer is an estimate of the
!     number of decimal digits of accuracy in the angular function
!     (and its first derivative when iopang = 2). It is a conservative
!     estimate based on the calculated subtraction error in the series
!     calculations of the functions. When the accuracy estimate
!     is equal to 0, the corresponding function values are set equal
!     to zero. (i2; see statements 1090 and 1100.
!
!   fort.40 and fort.50
!
!     These files are diagnostic files that contain information
!     about specific techniques used and numbers of terms required
!     for the radial function and angular function calculations,
!     respectively. They are annotated and should be self
!     explanatory. Calculated Values are given in full real*16 and
!     complex*32 precision.
!
!   The user may desire values for the d coefficients that appear in
!   the expression for the angular functions as well as in many of the
!   expressions used to calculate the radial functions. Ratios of these
!   coefficients are calculated in the subroutine dnorm. See the
!   comments given at the beginning of dnorm.
!
!  scalars and arrays
!
			complex(knd) :: cc, r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum)
			complex(knd) :: s1(lnum, narg), s1d(lnum, narg)
			complex(knd), allocatable, dimension(:,:) :: enr
			real*16 sgn
		complex*32 d01
		integer id01, m, lnum, ioprad, iopang, iopnorm, narg,&
	maxd
        real(knd) arg(narg),c,x1,xt
        character*10 time,date,zone
        integer dt(8)
        integer ir1e(lnum),ir1de(lnum),ir2e(lnum),ir2de(lnum)
		integer is1e(lnum, narg), is1de(lnum, narg)
!
!  ndec and nex should be set to conform with the user's computer.
!
!     ndec: number of decimal digits available in r*16 arithmetic.
!     nex:  maximum exponent available in r*16 arithmetic.
!
!	write(*, *) 'beginning of vB'
		is1de = 0
		is1de = 0

        ndec=33
        nex=4930
!
!  minimum desired accuracy. set equal to 16; can be changed if desired
        minacc=16
!
        ioparg=1
        mmin=m
        minc=0
        mnum=1
!
        call date_and_time(date,time,zone,dt)
!
!  open output files
        open(20, file='fort.20')
        open(30, file='fort.30')
        open(40, file='fort.40')
        open(50, file='fort.50')
!
!  set array dimensions
        c=cqabs(cc)
        maxc=2*iqint((x1+1.0q0)*(qabs(qreal(cc))+qabs(qimag(cc))))+25
        maxm=mmin+minc*(mnum-1)
        maxint=lnum+3*ndec+iqint(c)+5
        maxj=maxint+maxm
        maxp=maxint
        maxn=maxp+maxm
        maxn=max(maxn,maxc)
        maxpdr=4*ndec+5
        neta=993
        ngau=200
        if(ioprad /= 2) go to 10
        lnump=max(lnum+maxm,500)
        if(x1.ge.0.001q0) maxn=2*(lnump*(-18.5q0-20.q0*qlog10(x1))+&
                         5*ndec+4*maxm+c+01000)+maxm+5
        if(x1.gt.0.08q0) maxn=2*(lnump*(0.5q0-3.0q0*qlog10(x1))+&
                        5*ndec+4*maxm+c+01000)+maxm+5
        if(x1.gt.1.0q0) maxn=2*(lnump*0.5q0+5*ndec+4*maxm+c+00500)+&
                       maxm+5
        maxn=max(maxn,maxc)
        maxp=max(maxn,maxp)
        if(x1.lt.1.0q-3) ngau=200-50*iqint(qlog10(x1)-1.0q-30)
        if(x1.lt.1.0q-10) ngau=250-50*iqint(qlog10(x1)-1.0q-30)
        if(x1.lt.1.0q-11) ngau=1200
        if(x1.lt.1.0q-12) ngau=2400
        if(x1.le.0.5q0) maxpdr=maxpdr+iqint(2.q0*c+100.0q0*x1)+400
10      maxq=maxint+maxm+maxm
        maxd=maxp/2+1
        maxdr=maxpdr/2+1
        maxp=max(maxp,maxpdr)
        maxlp=lnum+maxm+5
        maxmp=maxm+5
        maxt=1
        jnenmax=10
        if(iopang.ne.0) maxt=narg
!

!		write(*,*) 'calculated dimension'
!		write(*, *) 'maxd = ', maxd
		allocate(enr(lnum, 0:maxd))
		enr = 0q0
!		write(*, *) 'maxd = ', maxd
        call main (mmin,minc,mnum,lnum,c,cc,ioprad,iopang,iopnorm,&
              minacc,x1,ngau,ioparg,arg,narg,neta,maxd,&
              maxdr,maxint,maxj,maxlp,maxm,maxmp,maxn,maxp,&
              maxpdr,maxq,maxt,jnenmax,ndec,nex,r1c,ir1e,r1dc,&
              ir1de,r2c,ir2e,r2dc,ir2de,s1,s1d,dt,enr(:,1:maxd))
!		write(*,*) 'we are out of main'
!
        end

!
!
!
        subroutine main (mmin,minc,mnum,lnum,c,cc,ioprad,iopang,iopnorm,&
                   minacc,x1,ngau,ioparg,arg,narg,neta,maxd,&
                   maxdr,maxint,maxj,maxlp,maxm,maxmp,maxn,maxp,&
                   maxpdr,maxq,maxt,jnenmax,ndec,nex,qr1,ir1,&
                   qr1d,ir1d,qr2,ir2,qr2d,ir2d,s1,s1d,dt,enr)
!
!  purpose:     to coordinate the calculation of both the prolate
!               spheroidal radial and angular functions and their
!               first derivatives using various algorithms.
!
!  scalars
        real(kind=16) aj1,aj2,apcoef,apcoefn,api,c,coefme,coefmo,coefn,&
          dec,deta,etaval,eta1,factor,pcoef,pcoefe,pcoefet,&
          pcoefn,pcoefo,pcoef1,pdcoefe,pdcoefet,pdcoefo,pi,qdml,&
          qml,rm,rm2,sgn,termpq,x,xb,xbninp,x1,blterm,testw,&
          testw1,testw2
        complex(kind=16) cc,cmfnorm,cfnorm,cr2lc,cr2dlc,c2,c4,d01,dneg,&
             eigval,r1c,r1cb,r1dc,r1dcb,r2c,r2dc,r2ec,r2dec,r2ic,&
             r2dic,r2lc,r2dlc,r2nc,r2dnc,wronc,wront
!
!  arrays with dimension lnum
        integer iqdl(lnum),iql(lnum),ifajo(lnum)
        real*16 qdl(lnum),ql(lnum)
        complex*32 fajo(lnum),eigt(lnum),eigst(lnum)
!
!  complex*32 arrays with dimension lnum
        complex*32 qr1(lnum),qr1d(lnum),qr2(lnum),qr2d(lnum)
!
!  integer arrays with dimension lnum
        integer ir1(lnum),ir1d(lnum),ir2(lnum),ir2d(lnum)
!
!  complex*32 arrays with dimensions lnum and narg
        complex*32 s1(lnum,narg),s1d(lnum,narg)
!
!  arrays with dimension maxd
!  passed
        complex*32 enr(lnum,maxd)
!  created
		complex*32, allocatable, dimension(:) :: bliste, gliste,&
             blisto, glisto
!
!  arrays with dimension maxdr
        complex*32, allocatable, dimension(:) :: drhor
!
!  arrays with dimension maxint
        complex*32, allocatable, dimension(:) :: pint1, pint2, pint3,&
             pint4, rpint1, rpint2
!
!  arrays with dimension maxj
        complex*32, allocatable, dimension(:) :: sbesf, sbesf2, sbesdf,&
             sbesdf2
!
!  arrays with dimension maxlp
        integer, allocatable, dimension(:) :: ibese, ibese2, ipnormint
        real*16, allocatable, dimension(:) :: pnormint
        complex*32, allocatable, dimension(:) :: sbesdr, sbesdr2,&
             sbesn, sbesn2
        complex*32, allocatable, dimension(:) :: sneudre
		complex*32, allocatable, dimension(:,:) :: sneudrsv
        complex*32, allocatable, dimension(:) :: sneudr
!
!  arrays with dimension maxmp
        complex*32, allocatable, dimension(:) :: enrneg
        integer, allocatable, dimension(:) :: norme
!
!  arrays with dimension maxn
        integer, allocatable, dimension(:) :: ineue, ineuee
		integer, allocatable, dimension(:,:) :: ineuesv

		complex*32, allocatable, dimension(:) :: sneufe, sneudfe,&
          sneune
		complex*32, allocatable, dimension(:,:) :: sneufsv, sneudfsv,&
			sneunsv
        complex*32, allocatable, dimension(:) :: sneun, sneuf, sneudf
!
!  arrays with dimension given by maxp
        real*16, allocatable, dimension(:) :: alpha, beta, coefa, coefb,&
          coefc, coefd, coefe, gamma, pdratt, pratb, pratt,&
          prat1
		real*16, allocatable, dimension(:,:) :: pdr, pdrat, pr, prat,&
           pratbsv, prattsv, pdratsv
!
!  arrays with dimension maxpdr
        real*16, allocatable, dimension(:) :: prx, pdrx
!
!  arrays with dimension maxq
        real*16, allocatable, dimension(:) :: qr, qdr
!
!  arrays with dimension maxt
        real*16 arg(maxt)
		real*16, allocatable, dimension(:) :: barg, etainp, pdnorm,&
          pdnorma, pnorm, pnorma, pdtempe,&
          pdtempo, ptempe, ptempo, xin, xlninp
        complex*32, allocatable, dimension(:) :: s1c, s1dc
        integer, allocatable, dimension(:) :: ipdnorm, ipdnorma, ipnorm,&
            ipnorma, ipdtempe, ipdtempo,&
            iptempe, iptempo, is1e, is1de,&
            naccs, naccsd
!
!  arrays with dimension neta
        real*16, allocatable, dimension(:) :: eta, xbn, xln
!
!  arrays with dimension ngau
        real*16, allocatable, dimension(:) :: wr, xr
!
!  miscellaneous integer arrays
        integer nees(100),naccsav(100)
		integer, allocatable, dimension(:) :: neeb, limpsv,&
            limnsv, jelimsv
!
        character*10 time,date,zone
!
        integer dt(8),jetasav(1000)


		dt = 0
		jetasav = 0
		nees = 0
		naccsav = 0
        iqdl = 0
		iql = 0
		ifajo = 0
        qdl = 0
		ql = 0
		fajo = 0
		eigt = 0
		eigst = 0

!  ALLOCATION
		allocate(bliste(maxd))
		allocate(gliste(maxd))
		allocate(blisto(maxd))
		allocate(glisto(maxd))

		allocate(drhor(maxdr))

		allocate(pint1(maxint))
		allocate(pint2(maxint))
		allocate(pint3(maxint))
		allocate(pint4(maxint))
		allocate(rpint1(maxint))
		allocate(rpint2(maxint))

		allocate(sbesf(maxj))
		allocate(sbesf2(maxj))
		allocate(sbesdf(maxj))
		allocate(sbesdf2(maxj))

		allocate(ibese(maxlp))
		allocate(ibese2(maxlp))
		allocate(ipnormint(maxlp))

		allocate(pnormint(maxlp))

		allocate(sbesdr(maxlp))
		allocate(sbesdr2(maxlp))
		allocate(sbesn(maxlp))
		allocate(sbesn2(maxlp))

		allocate(sneudre(maxlp))
		allocate(sneudrsv(jnenmax,maxlp))
        allocate(sneudr(maxlp))

		allocate(enrneg(maxmp))
        allocate(norme(maxmp))

        allocate(ineue(maxn))
		allocate(ineuee(maxn))
		allocate(ineuesv(jnenmax,maxn))

		allocate(sneufe(maxn))
		allocate(sneudfe(maxn))
		allocate(sneune(maxn))

		allocate(sneufsv(jnenmax,maxn))
		allocate(sneudfsv(jnenmax,maxn))
		allocate(sneunsv(jnenmax,maxn))

        allocate(sneun(maxn))
		allocate(sneuf(maxn))
		allocate(sneudf(maxn))

        allocate(alpha(maxp))
		allocate(beta(maxp))
		allocate(coefa(maxp))
		allocate(coefb(maxp))
		allocate(coefc(maxp))
		allocate(coefd(maxp))
		allocate(coefe(maxp))
		allocate(gamma(maxp))
		allocate(pdratt(maxp))
		allocate(pratb(maxp))
		allocate(pratt(maxp))
		allocate(prat1(maxp))

		allocate(pdr(maxt,maxp))
		allocate(pdrat(maxt,maxp))
		allocate(pr(maxt,maxp))
		allocate(prat(maxt,maxp))
		allocate(pratbsv(jnenmax,maxp))
		allocate(prattsv(jnenmax,maxp))
		allocate(pdratsv(jnenmax,maxp))

		allocate(prx(maxpdr))
		allocate(pdrx(maxpdr))

		allocate(qr(maxq))
		allocate(qdr(maxq))

		allocate(barg(maxt))
		allocate(etainp(maxt))
		allocate(pdnorm(maxt))
		allocate(pdnorma(maxt))
		allocate(pnorm(maxt))
		allocate(pnorma(maxt))
		allocate(pdtempe(maxt))
		allocate(pdtempo(maxt))
		allocate(ptempe(maxt))
		allocate(ptempo(maxt))
		allocate(xin(maxt))
		allocate(xlninp(maxt))

        allocate(s1c(maxt))
		allocate(s1dc(maxt))

        allocate(ipdnorm(maxt))
		allocate(ipdnorma(maxt))
		allocate(ipnorm(maxt))
		allocate(ipnorma(maxt))
		allocate(ipdtempe(maxt))
		allocate(ipdtempo(maxt))
		allocate(iptempe(maxt))
		allocate(iptempo(maxt))
		allocate(is1e(maxt))
		allocate(is1de(maxt))
		allocate(naccs(maxt))
        allocate(naccsd(maxt))

		allocate(eta(neta))
		allocate(xbn(neta))
		allocate(xln(neta))

		allocate(wr(ngau))
		allocate(xr(ngau))

		allocate(neeb(jnenmax))
		allocate(limpsv(jnenmax))
		allocate(limnsv(jnenmax))
		allocate(jelimsv(jnenmax))

		bliste = 0
		gliste = 0
		blisto = 0
		glisto = 0

		drhor = 0

		pint1 = 0
		pint2 = 0
		pint3 = 0
		pint4 = 0
		rpint1 = 0
		rpint2 = 0

		sbesf = 0
		sbesf2 = 0
		sbesdf = 0
		sbesdf2 = 0

		ibese = 0
		ibese2 = 0
		ipnormint = 0

		pnormint = 0

		sbesdr = 0
		sbesdr2 = 0
		sbesn = 0
		sbesn2 = 0

		sneudre = 0
		sneudrsv = 0
        sneudr = 0

		enrneg = 0
        norme = 0

        ineue = 0
		ineuee = 0
		ineuesv = 0

		sneufe = 0
		sneudfe = 0
		sneune = 0

		sneufsv = 0
		sneudfsv = 0
		sneunsv = 0

        sneun = 0
		sneuf = 0
		sneudf = 0

        alpha = 0
		beta = 0
		coefa = 0
		coefb = 0
		coefc = 0
		coefd = 0
		coefe = 0
		gamma = 0
		pdratt = 0
		pratb = 0
		pratt = 0
		prat1 = 0

		pdr = 0
		pdrat = 0
		pr = 0
		prat = 0
		pratbsv = 0
		prattsv = 0
		pdratsv = 0

		prx = 0
		pdrx = 0

		qr = 0
		qdr = 0

		barg = 0
		etainp = 0
		pdnorm = 0
		pdnorma = 0
		pnorm = 0
		pnorma = 0
		pdtempe = 0
		pdtempo = 0
		ptempe = 0
		ptempo = 0
		xin = 0
		xlninp = 0

        s1c = 0
		s1dc = 0

        ipdnorm = 0
		ipdnorma = 0
		ipnorm = 0
		ipnorma = 0
		ipdtempe = 0
		ipdtempo = 0
		iptempe = 0
		iptempo = 0
		is1e = 0
		is1de = 0
		naccs = 0
		naccsd = 0

		eta = 0
		xbn = 0
		xln = 0

		wr = 0
		xr = 0

		neeb = 0
		limpsv = 0
		limnsv = 0
		jelimsv = 0

!       write(*,*) 'allocated everything'
!		write(*,*) 'jnenmax =', jnenmax
!		write(*,*) 'ngau =', ngau
!		write(*,*) 'neta =', neta
!		write(*,*) 'maxt =', maxt
!		write(*,*) 'maxpdr =', maxpdr
!		write(*,*) 'maxp =', maxp
!		write(*,*) 'maxq =', maxq
!		write(*,*) 'maxn =', maxn
!		write(*,*) 'maxlp =', maxlp
!		write(*,*) 'maxmp =', maxmp
!		write(*,*) 'maxj =', maxj
!		write(*,*) 'maxint =', maxint
!		write(*,*) 'maxdr =', maxdr
!		write(*,*) 'maxd =', maxd

10      dec=10.q0**(-ndec-1)
        if(ioprad.ne.0) x=x1+1.q0
        jtest=ndec-minacc-2
        pi=3.1415926535897932384626433832795028841971q0
        api=pi/180.q0
        c2=cc*cc
        c4=c2*c2
        nbp=iqint(2.0q0*(qabs(qreal(cc))+qabs(qimag(cc)))/3.14q0)
        imax=1.2*max(25,nbp/2)+5
        lical=imax+imax
!
!  begin loops
          if(iopang.eq.0) go to 20
            do jarg=1,narg
            barg(jarg)=arg(jarg)
            end do
20        continue
          igau=0
          if(ioprad.ne.0) write(40,30) x,cc
30        format(1x,'x = ',e38.30,/,1x,'c = ',e38.30,e38.30)
          time1=60*dt(6)+dt(7)+0.001*dt(8)
          if(ioprad.eq.2) wront=1.0q0/(cc*x1*(x1+2.0q0))
40        continue
            maxe=max(100,nbp+nbp)+10
            ibflag1=0
!
!         mnum is equal to 1 and mmin is equal to the input value of m
!
            do 900 mi=1,mnum
!			write(*,*) 'mi = ', mi
            m=mmin+minc*(mi-1)
            m2=m+m
            if(iopang.ne.0) write(50,50) cc,m
50          format(1x,'c = ',e38.30,e38.30,'; m = ',i5)
            if(ioprad.ne.0) write(40,60) m
60          format(1x,'m = ',i5)
            if(iopang.ne.0) write(30,70) cc,m
70          format(1x,e38.30,e38.30,i5)
            rm=qfloat(m)
            rm2=qfloat(m+m)
            icounter=0
            iopleg=0
            iopneu=0
            iopeta=0
            iopint=1
            jintm=0
            iopd=3
            limcsav=0
            if(ioprad.ne.2) go to 80
            if(x1.le.0.4q0.and.c.le.10.0q0) iopleg=1
            if(x1.gt.0.4q0.and.c.le.10.0q0) iopneu=1
            if(x1.le.0.4q0.and.c.le.15.0q0.and.minacc.le.16) iopleg=1
            if(x1.gt.0.4q0.and.c.le.20.0q0.and.minacc.le.16) iopneu=1
            if(iopleg.eq.1.or.iopneu.eq.1) iopint=0
            ioppsum=1
            iopqnsum=1
            if(m.eq.0) iopqnsum=0
            nee=1
            jnen=0
            incnee=64
            nacceta=0
            msearch=0
80          continue
            if(iopang.eq.0) go to 90
            limps1=lnum+3*ndec+iqint(c)
            if((limps1+3).gt.maxp) limps1=maxp-3
            iopd=0
            if(iopang.eq.2) iopd=1
            call pleg(m,limps1,maxp,limcsav,iopd,ndec,barg,narg,maxt,&
                pr,pdr,pdnorm,ipdnorm,pnorm,ipnorm,alpha,beta,&
                gamma,coefa,coefb,coefc,coefd,coefe)
            limcsav=limps1
            iopd=3
90          if(ioprad.eq.0.or.mi.ne.1) go to 100
            limj=lnum+3*ndec+iqint(c)+maxm
            xb=qsqrt(x1*(x1+2.q0))
            call sphbes(cc,xb,limj,maxj,maxlp,sbesf,sbesdf,sbesn,ibese,&
                  sbesdr)
            if(qimag(cc).ne.0.0q0) call sphbes(cc,x,limj,maxj,maxlp,&
                    sbesf2,sbesdf2,sbesn2,ibese2,sbesdr2)
100         iflag=0
            iflagp=0
            ibflag2=0
            legflag=0
            jflagleg=0
            naccleg=0
            legstart=max(nbp,iqint(.05q0*rm*qreal(cc)))
            nflag=0
            lowacc=ndec
            lowtest=minacc
            nacctest=minacc
            naccintp=0
            nacclegp=0
            naccneup=0
            naccr=minacc
            ietacount=0
            naccsub=0
            incnflag=0
            iplflag=0
            factor=1.0q0
            iflagpc=1
            iflagbesb=0
            iopbesb=0
110         continue
            if(ioprad.ne.0) write(20,120) x,cc,m
120         format(1x,e38.30,e38.30,e38.30,i5)
!           write(45,120) x,cc,m
            time1=60*dt(6)+dt(7)+0.001*dt(8)
            if(qabs(qimag(cc)).gt.50.0q0) iopeta=1
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc cycle by l
              do 850 li=1,lnum
!			  write(*,*) 'li = ', li
              l=m+(li-1)
              if(iopang.ne.0) write(30,140) l
140           format(1x,i6)
              if(iopang.ne.0) write(50,150) l
150           format(1x,'l = ',i6)
              ix=l-m-2*((l-m)/2)
              iopnee=0
              naccrsav=naccr
              naccr=-1
              nsav=ndec
              limdrad=3*ndec+iqint(c)
              if(ioprad.ne.0.and.li.ne.1) limdrad=jbes+jbes+20+&
                                            iqint(qsqrt(c))
              if(iopint.ne.0.and.li.ne.1.and.jintm.gt.jbes)&
            limdrad=jintm+jintm+20+iqint(qsqrt(c))
              limdbesb=3*ndec+iqint(c)+li
              if(iflagbesb.eq.1.and.iopbesb.eq.0.and.&
            limdbesb.gt.limdrad) limdrad=limdbesb
              limdang=3*ndec+iqint(c)
              if(iopang.ne.0.and.li.ne.1) limdang=jang+jang+20+&
                                            iqint(qsqrt(c))
              if(iopang.eq.0) limd=limdrad
              if(ioprad.eq.0) limd=limdang
              if(iopang.ne.0.and.ioprad.ne.0) limd=max(limdang,limdrad)
              if(li.eq.1) limmf=limdang
              if(li.gt.1) limmf=jmf+jmf+20+c/25
              limd=max(limd,limmf)
              if(ioprad.ne.2) go to 155
              if(iopleg.eq.1) limdleg=l-m+3*ndec+iqint(c)
              if(iopleg.eq.2) limdleg=jleg+jleg+20+iqint(qsqrt(c))
              if(iopleg.ne.0) limd=max(limd,limdleg)
              limdneu=limd
              lplus=max(l,500)
              if(x1.ge.0.001q0) limdneu=2*((lplus)*(-18.5q0-&
                                  20.q0*qlog10(x1))&
                                  +5*ndec+4*m+c+01000)
              if(x1.gt.0.08q0) limdneu=2*((lplus)*(0.5q0-&
                                  3.0q0*qlog10(x1))+&
                                  5*ndec+4*m+c+01000)
              if(x1.gt.1.0q0) limdneu=2*((lplus)*0.5q0+5*ndec+4*m+c+&
                                00500)
              if(iopneu.eq.2.and.naccneu.gt.0)&
               limdneu=jneu+jneu+20+iqint(qsqrt(c))
              if(iopneu.ne.0) limd=max(limd,limdneu)
              limdeta=limd
              if(x1.ge.0.001q0) limdeta=2*((lplus)*(-18.5q0-20.q0*&
                                qlog10(x1))+5*ndec+4*m+c+01000)
              if(x1.gt.0.08q0) limdeta=2*((lplus)*(0.5q0-3.0q0*&
                                qlog10(x1))+5*ndec+4*m+c+01000)
              if(x1.gt.1.0q0) limdeta=2*((lplus)*0.5q0+5*ndec+4*m+c+&
                                00500)
              if(iopeta.eq.3.and.naccrsav.gt.minacc)&
                        limdeta=jeta+jeta+500+c/10
              if(iopeta.eq.3.and.naccrsav.le.minacc)&
                        limdeta=jeta+jeta+500+c
              if(iopeta.ne.0) limd=max(limd,limdeta)
155           continue
              if(limd.gt.maxp) limd=maxp
!
!        obtain estimates for the eigenvalues to be used as starting
!        values for the bouwkamp procedure
              if(l.eq.m) call geteig(m,cc,nbp,lnum,maxe,eigst)
!
!        use bouwkamp procedure to obtain accurate eigenvalues
              if(l.eq.m) ienre=(3*ndec+iqint(c))/2
              if(l.eq.m) jlowe=1
              if(l.eq.m) limdle=2
              if(l.eq.m+1) ienro=(3*ndec+iqint(c))/2
              if(l.eq.m+1) jlowo=1
              if(l.eq.m+1) limdlo=3
!
!  compute the coeficients in the bouwkamp method
              if(ix.eq.1) go to 160
!
!  beta coefficients (bliste) for l-m even
              if(limdle.gt.limd) go to 163
              j=jlowe
                do 158 i=limdle,limd,2
                i2=i+i
                blterm=qfloat(i)*qfloat(i-1)*qfloat(m2+i)*&
                 qfloat(m2+i-1)/(qfloat(m2+i2-1)*&
                 qfloat(m2+i2-1)*qfloat(m2+i2-3)*&
                 qfloat(m2+i2+1))
                bliste(j)=qcmplx(blterm,0.0q0)*c4
                j=j+1
158             continue
!
!  gamma coeficients (gliste) for l-m even
              j=jlowe
                do 159 i=limdle-1,limd+1,2
                i2=i+i
                gliste(j)=qfloat(m+i-1)*qfloat(m+i)+0.5q0*c2*(1.0q0-&
                    qfloat(m2*m2-1)/(qfloat(m2+i2-3)*&
                    qfloat(m2+i2+1)))
                j=j+1
159             continue
              go to 163
160           continue
!
!  beta coefficients (blist0) for l-m odd
              if(limdlo.gt.limd) go to 163
              j=jlowo
                do 161 i=limdlo,limd,2
                i2=i+i
                blterm=qfloat(i)*qfloat(i-1)*qfloat(m2+i)*&
                 qfloat(m2+i-1)/(qfloat(m2+i2-1)*&
                 qfloat(m2+i2-1)*qfloat(m2+i2-3)*&
                 qfloat(m2+i2+1))
                blisto(j)=qcmplx(blterm,0.0q0)*c4
                j=j+1
161             continue
!
!  gamma coeficient (glist0) for l-m odd
              j=jlowo
                do 162 i=limdlo-1,limd+1,2
                i2=i+i
                glisto(j)=qfloat(m+i-1)*qfloat(m+i)+0.5q0*c2*(1.0q0-&
                    qfloat(m2*m2-1)/(qfloat(m2+i2-3)*&
                    qfloat(m2+i2+1)))
              j=j+1
162           continue
163           continue
              if(li.le.lical) eigval=eigst(li)
              if(li.gt.lical) eigval=4.0q0*eigt(li-1)-&
              6.0q0*eigt(li-2)+4.0q0*eigt(li-3)-eigt(li-4)
              itestm=ndec
              idigc=ndec

              if(ix.eq.0) call conver (l,m,cc,limd,ndec,maxd,bliste,&
                                 gliste,lical,ioprad,ienre,eigval,&
                                 enr(li,:),idigc,itestm)
              if(ix.eq.1) call conver (l,m,cc,limd,ndec,maxd,blisto,&
                                 glisto,lical,ioprad,ienro,eigval,&
                                 enr(li,:),idigc,itestm)

              eigt(li)=eigval
              minacc=min(minacc,idigc)
              if(ioprad.ne.0) write(20,165) l,eigval
              if(ioprad.ne.0) write(40,170) eigval
165           format(1x,'l =',i5,2x,'eigenvalue =',e24.16,e24.16)
170           format(10x,'eigenvalue =',e40.31,e40.31)
              if(ix.eq.1) go to 176
              if(limdle.gt.limd) go to 177
              limdle=limd+2
              if(2*(limd/2).ne.limd) limdle=limd+1
              jlowe=limd/2+1
              go to 177
176           if(limdlo.gt.limd) go to 177
              limdlo=limd+1
              if(2*(limd/2).ne.limd) limdlo=limd+2
              jlowo=(limd-1)/2+1
177           call dnorm (l,m,cc,ndec,limd,maxd,enr(li,:),sgn,d01,id01,&
                     cmfnorm,cfnorm,jmf,nsubmf,jfla,nsubf)

             jmf=max(jmf,jfla)
              if(ioprad.eq.0) go to 720
              if(l.eq.m.and.nsubmf.gt.jtest) iopint=1
              if(l.eq.m.and.nsubmf.gt.jtest.and.iopneu.ne.0) iopneu=0
              if(l.eq.m.and.nsubmf.gt.jtest.and.iopleg.ne.0) iopleg=0
!
!  determine prolate radial functions of the first kind
              if(li.eq.1) limr1=3*ndec+iqint(c)
              if(li.ne.1) limr1=jbesa+jbesa+20+c/25
              call r1besa(l,m,cc,x1,limr1,ndec,maxd,enr(li,:),maxj,&
             maxlp,nex,iflag,sbesf,sbesdf,sbesn,ibese,sbesdr,&
             d01,id01,r1c,ir1e,r1dc,ir1de,cfnorm,jbesa,factor,&
			   nsubr1a)

              naccr1a=min(idigc,itestm)-max(nsubr1a,nsubf)
              write(40,180) r1c,ir1e,r1dc,ir1de
180           format(10x,'r1 = ', f35.31,f35.31,i6,/,10x,&
               'r1d = ',f35.31,f35.31,i6)
              ichoicer1=1
              iflagbesb=0
              if(naccr1a.lt.minacc+5) iflagbesb=1
              naccr1=naccr1a
              if(naccr1a.gt.minacc) then
              iopbesb=0
              jbesb=0
              go to 186
              end if
!
              if(iflagpc.eq.0) go to 185
              prat1(1)=1.0q0
              prat1(2)=rm2+1.0q0
                do jp=3,limj-maxm
                aj1=qfloat(jp-1)
                aj2=qfloat(jp-2)
                prat1(jp)=(rm2+aj1)*(rm2+aj2)/(aj1*aj2)
                end do
              pcoefn=x1*(x1+2.0q0)/(x*x)
              apcoefn=(rm/2.q0)*qlog10(pcoefn)
              ipcoefn=iqint(apcoefn)
              pcoefn=10.0q0**(apcoefn-ipcoefn)
              iflagpc=0
185           continue
!
              if(li.eq.1.or.iopbesb.eq.0) limr1=3*ndec+iqint(c)
              if(li.ne.1.and.iopbesb.eq.1) limr1=jbesb+jbesb+20+c/25
              iopbesb=1
              if(qimag(cc).ne.0.0q0) call r1besb(l,m,cc,x1,limr1,ndec,&
             maxd,maxlp,maxj,maxp,enr(li,:),sbesf2,sbesn2,ibese2,&
             sbesdf2,sbesdr2,prat1,pcoefn,ipcoefn,cmfnorm,r1cb,&
             ir1eb,r1dcb,ir1deb,jbesb,nsubr1b)

              write(40,180) r1cb,ir1eb,r1dcb,ir1deb
              naccr1b=min(idigc,itestm)-max(nsubr1b,nsubmf)
              if(naccr1a.ge.naccr1b) go to 186
              naccr1=naccr1b
              r1c=r1cb
              ir1e=ir1eb
              r1dc=r1dcb
              ir1de=ir1deb
              ichoicer1=2
186           jbes=max(jbesa,jbesb)
              if(ioprad.eq.2) ir2est=iqint(qlog10(cqabs(wront)))-ir1de+1
!
!  determine prolate radial functions of the second kind
!
!  calculation using integration technique
              if(ioprad.ne.2) go to 680
              if(iopint.eq.0) go to 230
              if(iopint.eq.2) go to 190
              limint=lnum+3*ndec+iqint(c)
              if(igau.eq.0) call gauss(ndec,ngau,xr,wr)
              igau=1
              ngqs=10
              if(c.gt.2000.0q0) ngqs=ngqs*(c/2000.0q0)*&
                                     (c/2000.0q0)
              call pint(cc,m,lnum,x1,limint,maxint,maxlp,maxmp,ndec,&
                  wr,xr,ngau,ngqs,rpint1,rpint2,pint1,&
                  pint2,pint3,pint4,norme,pnormint,ipnormint,&
                  coefme,coefmo)
190           continue
              if(iopint.eq.1) limint=3*ndec+iqint(c)
              if(iopint.eq.2) limint=jintm+jintm+20+iqint(qsqrt(c))
              call r2int(l,m,cc,limint,ndec,maxd,enr(li,:),d01,id01,&
                   maxint,maxmp,nex,maxlp,rpint1,rpint2,&
                   pint1,pint2,pint3,pint4,norme,pnormint,&
                   ipnormint,coefme,coefmo,ir2est,r2ic,ir2ie,&
                   r2dic,ir2die,jint,coefn,icoefn)

              iopint=2
              if(jint.gt.jintm) jintm=jint
              if(iopint.eq.3.and.jint.lt.jintm) jint=jintm
              wronc=r1c*r2dic*10.0q0**(ir1e+ir2die)-r2ic*r1dc*&
              10.0q0**(ir2ie+ir1de)
              naccint=-iqint(qlog10(cqabs((wronc-wront)/wront)+dec))
              if(naccint.lt.0) naccint=0
              testw1=qabs(qreal((r1c*r2dic)*10.0q0**(ir1e+ir2die)))
              testw2=qabs(qreal((r2ic*r1dc)*10.0q0**(ir2ie+ir1de)))
              testw=testw1
              if(testw2.gt.testw1) testw=testw2
              naccsubri=-iqint(qlog10(qabs(qreal(wronc)/testw)+dec))
              testw1=qabs(qimag((r1c*r2dic)*10.0q0**(ir1e+ir2die)))
              testw2=qabs(qimag((r2ic*r1dc)*10.0q0**(ir2ie+ir1de)))
              testw=testw1
              if(testw2.gt.testw1) testw=testw2
              naccsubii=-iqint(qlog10(qabs(qimag(wronc)/testw)+dec))
              naccsubi=min(naccsubri,naccsubii)
              if(naccsubi.lt.0) naccsubi=0
              if(ichoicer1.eq.1) naccint=min(naccint+naccsubi,idigc-&
                                       nsubf)
              if(ichoicer1.eq.2) naccint=min(naccint+naccsubi,idigc-&
                                       nsubmf)
              naccout=naccint
              if(naccintp-naccint.gt.8) naccint=naccintp
              if(naccout.lt.naccr) go to 200
              naccr=naccout
              r2c=r2ic
              ir2e=ir2ie
              r2dc=r2dic
              ir2de=ir2die
200           continue
              istartr2=1
              if(naccint.gt.minacc.and.naccintp.gt.minacc) then
              iopleg=0
              iopneu=0
              iopeta=0
              istartr2=0
              end if
              if(naccout.ge.minacc) then
              iopneu=0
              iopeta=0
              end if
              if(naccint.eq.0) iopint=0
              naccintp=naccout
210           if(naccsubi.eq.0) write(40,215) naccout,r2ic,ir2ie,&
                                        r2dic,ir2die
215           format(15x,'accuracy =',i3,&
               ' decimal digits.'/,10x,'r2 = ', f35.31,f35.31,&
               2x,i6,/,10x,'r2d = ',f35.31,f35.31,2x,i6)
              if(naccsubi.gt.0) write(40,220) naccout,naccsubi,r2ic,&
                                        ir2ie,r2dic,ir2die
220           format(15x,'accuracy =',i3,&
               ' decimal digits, adjusted for',i3,' subtraction',&
               ' digits in forming Wronskian',/,10x,'r2 = ',&
               f35.31,f35.31,2x,i6,/,10x,'r2d = ',f35.31,f35.31,&
               2x,i6)
230           continue
!
!  calculation using Legendre expansion and joining factor
              naccleg=0
              if(iopleg.eq.0) go to 360
              if(iopleg.eq.2.or.jflagleg.eq.1) go to 310
              jflagleg=1
              limdr=c+ndec+50.q0*x1+200
              if(limdr.gt.maxdr) limdr=maxdr
              if(ioppsum.eq.0) go to 250
              xin(1)=x
              limpleg=limdr+limdr
              call pleg(m,limpleg,maxp,limcsav,iopd,ndec,xin,1,maxt,&
                  prat,pdrat,pdnorma,ipdnorma,pnorma,ipnorma,&
                  alpha,beta,gamma,coefa,coefb,coefc,coefd,coefe)
              limcsav=max(limcsav,limpleg)
                do jj=1,limpleg
                prx(jj)=prat(1,jj)
                pdrx(jj)=pdrat(1,jj)
                end do
250           limq=lnum+3*ndec+iqint(c)
              call qleg(m,lnum,limq,maxq,x1,ndec,qdr,qdml,iqdml,qdl,&
                  iqdl,qr,qml,iqml,ql,iql,termpq,itermpq)
              fajo(1)=cc/(rm2-1.0q0)
              ifajo(1)=0
              if(m.eq.0) go to 280
                do im=1,m
                fajo(1)=fajo(1)*qfloat(im+im)/cc
                if(cqabs(fajo(1)).lt.1.q+10) go to 260
                fajo(1)=fajo(1)*(1.q-10)
                ifajo(1)=ifajo(1)+10
260             continue
                if(cqabs(fajo(1)).gt.1.q-10) go to 270
                fajo(1)=fajo(1)*(1.q+10)
                ifajo(1)=ifajo(1)-10
270             continue
                end do
280           continue
              fajo(2)=-cc*fajo(1)/(rm2-3.0q0)
              ifajo(2)=ifajo(1)
                do jl=3,lnum-1,2
                fajo(jl)=fajo(jl-2)*(qfloat(jl+m+m-1)/qfloat(jl-2))
                ifajo(jl)=ifajo(jl-2)
                if(cqabs(fajo(jl)).lt.1.0q10) go to 290
                fajo(jl)=fajo(jl)*1.0q-10
                ifajo(jl)=ifajo(jl)+10
290             fajo(jl+1)=fajo(jl-1)*(qfloat(jl+m+m-1)/qfloat(jl))
                ifajo(jl+1)=ifajo(jl-1)
                if(cqabs(fajo(jl+1)).lt.1.0q10) go to 300
                fajo(jl+1)=fajo(jl+1)*1.0q-10
                ifajo(jl+1)=ifajo(jl+1)+10
300             end do
              if(2*(lnum/2).eq.lnum.or.lnum.eq.2) go to 310
              fajo(lnum)=fajo(lnum-2)*qfloat(lnum+m+m-1)/qfloat(lnum-2)
              ifajo(lnum)=ifajo(lnum-2)
310           continue
              limleg=l-m+3*ndec+iqint(c)
              limdr=c+ndec+50.0q0*x1+200
              if(iopleg.eq.2) limleg=jleg+jleg+20+iqint(qsqrt(c))
              if(iopleg.eq.2) limdr=jlegp+10+iqint(0.5q0*qsqrt(c))
              if(limdr.gt.maxdr) limdr=maxdr
              call dalt(l,m,cc,limdr,maxdr,maxmp,ioppsum,eigval,&
                  enrneg,drhor,dneg,idneg)
              call r2leg(l,m,cc,x1,lnum,limleg,limdr,iflagp,&
                   ndec,maxd,maxmp,maxpdr,maxdr,maxq,enr(li,:),&
                   enrneg,drhor,d01,id01,dneg,idneg,cfnorm,&
                   nsubf,cmfnorm,nsubmf,prx,pdrx,qdr,qdml,iqdml,&
                   qdl,iqdl,qr,qml,iqml,ql,iql,fajo,ifajo,termpq,&
                   itermpq,ioppsum,iopqnsum,sgn,r1dc,ir1de,r2lc,&
                   ir2le,r2dlc,ir2dle,jleg,jlegp)

              wronc=r1c*r2dlc*10.0q0**(ir1e+ir2dle)-r2lc*r1dc*&
               10.0q0**(ir2le+ir1de)
              naccleg=-iqint(qlog10(cqabs((wronc-wront)/wront)+dec))
              if(naccleg.lt.0) naccleg=0
              testw1=qabs(qreal((r1c*r2dlc)*10.0q0**(ir1e+ir2dle)))
              testw2=qabs(qreal((r2lc*r1dc)*10.0q0**(ir2le+ir1de)))
              testw=testw1
              if(testw2.gt.testw1) testw=testw2
              naccsubrl=-iqint(qlog10(qabs(qreal(wronc)/testw))+dec)
              testw1=qabs(qimag((r1c*r2dlc)*10.0q0**(ir1e+ir2dle)))
              testw2=qabs(qimag((r2lc*r1dc)*10.0q0**(ir2le+ir1de)))
              testw=testw1
              if(testw2.gt.testw1) testw=testw2
              naccsubil=-iqint(qlog10(qabs(qimag(wronc)/testw)+dec))
              naccsubl=min(naccsubrl,naccsubil)
              if(naccsubl.lt.0) naccsubl=0
              if(naccsubl.lt.2) go to 315
              if(ichoicer1.eq.1) naccsubl=min(naccsubl,idigc-nsubf)
              if(ichoicer1.eq.2) naccsubl=min(naccsubl,idigc-nsubmf)
              if(naccsubl.lt.0) naccsub1=0
              if(ichoicer1.eq.1) naccleg=min(naccleg+naccsubl,idigc-&
                                       nsubf)
              if(ichoicer1.eq.2) naccleg=min(naccleg+naccsubl,idigc-&
                                       nsubmf)
              nagg1=-iqint(qlog10(qabs((cqabs(r2lc*(10.0q0**ir2le))-&
              cqabs(r1c*(10.0q0**ir1e)))/cqabs(r1c*&
              (10.0q0**ir1e)))+dec))
              nagg2=-iqint(qlog10(qabs((cqabs(r2dlc*(10.0q0**ir2dle))-&
              cqabs(r1dc*(10.0q0**ir1de)))/cqabs(r1dc*&
              (10.0q0**ir1de)))+dec))
              nagg=max(nagg1,nagg2)
              if(nagg.lt.0) nagg=0
              if(ichoicer1.eq.1.and.nagg.gt.idigc-nsubf)&
            nagg=idigc-nsubf
              if(ichoicer1.eq.2.and.nagg.gt.idigc-nsubmf)&
            nagg=idigc-nsubmf
              if(naccleg.lt.nagg) naccleg=nagg
315           naccout=naccleg
              if(nacclegp-naccleg.gt.8) naccleg=nacclegp
              nacclegp=naccleg
              if(naccout.le.naccr) go to 320
              naccr=naccout
              r2c=r2lc
              ir2e=ir2le
              r2dc=r2dlc
              ir2de=ir2dle
320           continue
              if(naccout.ge.naccrsav) then
              iopleg=2
              else
              if(iopleg.eq.1.and.l.ne.m) iopleg=0
              legstart=l+naccrsav-naccout
              end if
              if(naccout.gt.minacc) then
              iopleg=2
              iopneu=0
              iopeta=0
              iopint=0
              end if
340           continue
              if(naccsubl.lt.2) write(40,350) naccout,r2lc,ir2le,&
                                        r2dlc,ir2dle
350           format(15x,'accuracy =',i3,&
               ' decimal digits.'/,10x,'r2 = ', f35.31,f35.31,&
               2x,i6,/,10x,'r2d = ',f35.31,f35.31,2x,i6)
              if(naccsubl.ge.2) write(40,355) naccout,naccsubl,r2lc,&
                                        ir2le,r2dlc,ir2dle
355           format(15x,'accuracy =',i3,&
               ' decimal digits, adjusted for',i3,' subtraction',&
               ' digits in forming Wronskian',/,10x,'r2 = ',&
               f35.31,f35.31,2x,i6,/,10x,'r2d = ',f35.31,f35.31,&
               2x,i6)
360           continue
!
!  calculation using conventional Neumann expansion (eta=1)
              if(iopneu.eq.0) go to 420
              if(iopneu.eq.2) go to 380
              if(ibflag1.eq.1) go to 370
              ibflag1=1
              lnump=max(lnum+maxm,500)
              limn=2*(lnump*(-18.5q0-20.q0*qlog10(x1))+&
              5*ndec+4*m+c+01000)+maxm
              if(x1.gt.0.08q0) limn=2*(lnump*(0.5q0-3.0q0*qlog10(x1))+&
                             5*ndec+4*m+c+01000)+maxm
              if(x1.gt.1.0q0) limn=2*(lnump*0.5q0+5*ndec+4*m+c+00500)+&
                             maxm
              if(limn.gt.maxn) limn=maxn
              limbesf=2*iqint(qreal(cc*x)+qabs(qimag(cc*x)))+25
              call sphneu(cc,x,limn,maxn,maxlp,limbesf,sneuf,sneun,&
                    ineue,sneudf,sneudr)
370           if(ibflag2.eq.1) go to 380
              ibflag2=1
              lp=max(lnum+m,500)
              limp=2*(lp*(-18.5q0-20.q0*qlog10(x1))+&
             5*ndec+4*m+c+01000)
              if(x1.gt.0.08q0) limp=2*(lp*(0.5q0-3.0q0*qlog10(x1))+&
                             5*ndec+4*m+c+01000)
              if(x1.gt.1.0q0) limp=2*(lp*0.5q0+5*ndec+4*m+c+00500)
              if(limp.gt.maxp) limp=maxp
              prat1(1)=1.0q0
              prat1(2)=rm2+1.0q0
                do jp=3,limp
                aj1=qfloat(jp-1)
                aj2=qfloat(jp-2)
                prat1(jp)=(rm2+aj1)*(rm2+aj2)/(aj1*aj2)
                end do
              pcoefn=x1*(x1+2.0q0)/(x*x)
              apcoefn=(rm/2.q0)*qlog10(pcoefn)
              ipcoefn=iqint(apcoefn)
              pcoefn=10.0q0**(apcoefn-ipcoefn)
380           continue
              lplus=max(l,500)
              limneu=2*((lplus)*(-18.5q0-20.0q0*qlog10(x1))+&
               5*ndec+4*m+c+00500)
              if(x1.ge.0.1q0) limneu=2*((lplus)*(0.5q0-qlog10(x1))+&
                               5*ndec+4*m+c+00500)
              if(x1.gt.1.0q0) limneu=2*((lplus)*0.5q0+5*ndec+4*m+c+&
                               00500)
              if(iopneu.eq.2.and.naccneu.gt.0) limneu=jneu+jneu+20+&
                                                iqint(qsqrt(c))
              if(limneu.gt.limp-2) limneu=limp-2
              call r2neu(l,m,cc,x1,limneu,ndec,maxd,maxlp,maxn,&
                   maxp,minacc,enr(li,:),sneuf,sneun,ineue,sneudf,&
                   sneudr,prat1,pcoefn,ipcoefn,cmfnorm,r1dc,ir1de,&
                   r2nc,ir2ne,r2dnc,ir2dne,jneu,naccsub)

              wronc=r1c*r2dnc*10.q0**(ir1e+ir2dne)-r2nc*r1dc*&
              10.q0**(ir2ne+ir1de)
              naccneu=-iqint(qlog10(cqabs((wronc-wront)/wront)+dec))
              if(naccneu.lt.0) naccneu=0
              testw1=qabs(qreal((r1c*r2dnc)*10.0q0**(ir1e+ir2dne)))
              testw2=qabs(qreal((r2nc*r1dc)*10.0q0**(ir2ne+ir1de)))
              testw=testw1
              if(testw2.gt.testw1) testw=testw2
              naccsubrn=-iqint(qlog10(qabs(qreal(wronc)/testw))+dec)
              testw1=qabs(qimag((r1c*r2dnc)*10.0q0**(ir1e+ir2dne)))
              testw2=qabs(qimag((r2nc*r1dc)*10.0q0**(ir2ne+ir1de)))
              testw=testw1
              if(testw2.gt.testw1) testw=testw2
              naccsubin=-iqint(qlog10(qabs(qimag(wronc)/testw)+dec))
              naccsubn=min(naccsubrn,naccsubin)
              if(naccsubn.lt.0) naccsubn=0
              if(naccsubn.lt.2) go to 385
              if(ichoicer1.eq.2) naccsubn=min(naccsubn,idigc-nsubmf)
              if(naccsubn.lt.0) naccsubn=0
              if(ichoicer1.eq.1) naccneu=min(naccneu+naccsubn,idigc-&
                                       nsubf)
              if(ichoicer1.eq.2) naccneu=min(naccneu+naccsubn,&
                                       idigc-nsubmf)
              if(naccneu.lt.0) naccneu=0
385           naccout=naccneu
              if(naccneup-naccneu.gt.8) naccneu=naccneup
              naccneup=naccneu
              if(naccout.le.naccr) go to 390
              naccr=naccout
              r2c=r2nc
              ir2e=ir2ne
              r2dc=r2dnc
              ir2de=ir2dne
390           continue
              if(naccneu.ge.minacc) then
              iopneu=2
              iopeta=0
              end if
              if(naccneu.gt.minacc) then
              nflag=1
              if(li.gt.nbp) iopint=0
              end if
              if(iopeta.eq.0.and.naccr.lt.minacc) then
              nflag=0
              iopeta=1
              nee=max(nee,993-incnee*(minacc-naccneu+1))
              end if
              if(nflag.eq.1) go to 400
              if(naccneu.lt.naccr) iopneu=0
400           IF(naccsubn.lt.2) write(40,410) naccout,r2nc,ir2ne,&
                                        r2dnc,ir2dne
410           format(15x,'Wronskian accuracy =',i3,&
               ' decimal digits.'/,10x,'r2 = ', f35.31,f35.31,2x,&
               i5,/,10x,'r2d = ',f35.31,f35.31,2x,i5)
              if(naccsubn.ge.2) write(40,415) naccout,naccsubn,r2nc,&
                                        ir2ne,r2dnc,ir2dne
415           format(15x,'accuracy =',i3,&
               ' decimal digits, adjusted for',i3,' subtraction',&
               ' digits in forming Wronskian',/,10x,'r2 = ',&
               f35.31,f35.31,2x,i6,/,10x,'r2d = ',f35.31,f35.31,&
               2x,i6)
420           continue
!
!  calculation using the variable eta expansion
              if(iopeta.eq.0) go to 670
                do 430 inn=1,100
                nees(inn)=0
430             naccsav(inn)=0
              inen=0
              neemark=nee
              naccetamax=0
              neemax=nee
              naccnmax=0
              nacctemp=0
              netatry=1
              if(iopeta.gt.1) go to 440
              kounter=0
                do jnet=1,neta
                eta(jnet)=qcos((neta+1-jnet)*.5q0*pi/(neta+1))
                xbn(jnet)=qsqrt(x1*(x1+2.0q0)+eta(jnet)*eta(jnet))
                xln(jnet)=eta(jnet)*(x1+1.0q0)/xbn(jnet)
                end do
              iopeta=2
440           if(iopeta.eq.3) go to 540
              etaval=eta(nee)
450           xbninp=xbn(nee)
              netainp=1
              etainp(1)=eta(nee)
              xlninp(1)=xln(nee)
              lplus=max(l,500)
              limn=2*((lplus)*(-18.5q0-20.0q0*qlog10(x1))+&
             5*ndec+10*incnee+4*m+c+01000)+m
              if(x1.gt.0.08q0) limn=2*((lplus)*(0.5q0-&
           3.0q0*qlog10(x1))+10*ndec+10*incnee+4*m+c+01000)+m
              if(x1.gt.1.0q0) limn=2*((lplus)*0.5q0+5*ndec+10*incnee+&
                             4*m+c+00500)+m
              if(limn.gt.maxn) limn=maxn-2
              limp=limn-m
              if(jnen.eq.0) go to 510
              jnenlim=jnen
              if(jnen.gt.jnenmax) jnenlim=jnenmax
              limplim=limp
              limnlim=limn
                do 500 jn=1,jnenlim
                if(nee.ne.neeb(jn)) go to 500
                if(limplim.gt.limpsv(jn)) limplim=limpsv(jn)
                if(limnlim.gt.limnsv(jn)) limnlim=limnsv(jn)
                  do 460 je=1,limplim
                  pratb(je)=pratbsv(jn,je)
                  pratt(je)=prattsv(jn,je)
460               pdratt(je)=pdratsv(jn,je)
                  do 470 je=1,limnlim
                  sneufe(je)=sneufsv(jn,je)
470               sneudfe(je)=sneudfsv(jn,je)
                  jelim=maxlp
                  if(maxlp.gt.limn+1) jelim=limn+1
                  if(jelim.gt.jelimsv(jn)) jelim=jelimsv(jn)
                  do 480 je=1,jelim
                  sneune(je)=sneunsv(jn,je)
                  sneudre(je)=sneudrsv(jn,je)
480               ineuee(je)=ineuesv(jn,je)
                write(40,490) etaval
490             format(8x,'r2eta: reused expansion functions for eta ='&
                 ,f13.9,'.')
                go to 530
500             continue
510           continue
              jnen=jnen+1
              jnencur=jnen-(jnenmax*int((jnen-1)/jnenmax))
              neeb(jnencur)=nee
              limbesf=2*iqint(qabs(qreal(cc*xbninp))+&
                qabs(qimag(cc*xbninp)))+25
              call sphneu(cc,xbninp,limn,maxn,maxlp,limbesf,sneufe,&
                    sneune,ineuee,sneudfe,sneudre)
                do je=1,limn
                sneufsv(jnencur,je)=sneufe(je)
                sneudfsv(jnencur,je)=sneudfe(je)
                limnsv(jnencur)=limn
                end do
                jelim=maxlp
                if(maxlp.gt.limn+1) jelim=limn+1
                do 520 je=1,jelim
                sneunsv(jnencur,je)=sneune(je)
                sneudrsv(jnencur,je)=sneudre(je)
520             ineuesv(jnencur,je)=ineuee(je)
              jelimsv(jnencur)=jelim
              call pleg(m,limp,maxp,limcsav,iopd,ndec,xlninp,netainp,&
                  maxt,prat,pdrat,pdnorma,ipdnorma,pnorma,ipnorma,&
                  alpha,beta,gamma,coefa,coefb,coefc,coefd,coefe)
              limcsav=max(limcsav,limp)
                do je=1,limp
                pratt(je)=prat(1,je)
                pdratt(je)=pdrat(1,je)
                prattsv(jnencur,je)=pratt(je)
                pdratsv(jnencur,je)=pdratt(je)
                limpsv(jnencur)=limp
                end do
              limpd=2*(lnum+iqint(c)+ndec)
              if(limpd.gt.limp) limpd=limp
              iopd=2
              call pleg(m,limpd,maxp,limcsav,iopd,ndec,etainp,netainp,&
                  maxt,prat,pdrat,pdnorma,ipdnorma,pnorma,ipnorma,&
                  alpha,beta,gamma,coefa,coefb,coefc,coefd,coefe)
              iopd=3
                do je=1,limpd
                pratb(je)=prat(1,je)
                pratbsv(jnencur,je)=pratb(je)
                end do
              pratb(limpd+1)=0.0q0
              pratb(limpd+2)=0.0q0
              pratbsv(jnencur,limpd+1)=0.0q0
              pratbsv(jnencur,limpd+2)=0.0q0
530           continue
              pcoefe=((x1*(x1+2.0q0))/(x1*(x1+2.0q0)+eta(nee)**2))
              apcoef=(rm/2.0q0)*qlog10(pcoefe)
              ipcoefe=iqint(apcoef)
              pcoefe=10.0q0**(apcoef-ipcoefe)
              pcoefo=pcoefe*pratt(2)/pratb(2)
              ipcoefo=ipcoefe
              pdcoefe=pcoefe
              if(m.ne.0) pdcoefe=pcoefe*rm*xln(nee)/(xln(nee)*xln(nee)-&
                            1.0q0)
              ipdcoefe=ipcoefe
              pdcoefo=pdcoefe*pdratt(2)/pratb(2)
              ipdcoefo=ipdcoefe
              if(li.lt.3) go to 540
                do jl=3,li+ix,2
                pcoefe=pcoefe*pratt(jl)/pratb(jl)
                iterm=qlog10(qabs(pcoefe))
                pcoefe=pcoefe*10.0q0**(-iterm)
                ipcoefe=ipcoefe+iterm
                pdcoefe=pdcoefe*pdratt(jl)/pratb(jl)
                iterm=qlog10(qabs(pdcoefe))
                pdcoefe=pdcoefe*10.0q0**(-iterm)
                ipdcoefe=ipdcoefe+iterm
                end do
              continue
              if(li.lt.4) go to 540
                do jl=4,li+1-ix,2
                pcoefo=pcoefo*pratt(jl)/pratb(jl)
                iterm=qlog10(qabs(pcoefo))
                pcoefo=pcoefo*10.0q0**(-iterm)
                ipcoefo=ipcoefo+iterm
                pdcoefo=pdcoefo*pdratt(jl)/pratb(jl)
                iterm=qlog10(qabs(pdcoefo))
                pdcoefo=pdcoefo*10.0q0**(-iterm)
                ipdcoefo=ipdcoefo+iterm
                end do
540           continue
              if(ix.eq.0) go to 550
              pcoefet=pcoefo
              ipcoefet=ipcoefo
              pcoefo=pcoefo*pratt(li+2)/pratb(li+2)
              iterm=iqint(qlog10(qabs(pcoefo)))
              pcoefo=pcoefo*10.0q0**(-iterm)
              ipcoefo=ipcoefo+iterm
              pdcoefet=pdcoefo
              ipdcoefet=ipdcoefo
              pdcoefo=pdcoefo*pdratt(li+2)/pratb(li+2)
              iterm=iqint(qlog10(qabs(pdcoefo)))
              pdcoefo=pdcoefo*10.0q0**(-iterm)
              ipdcoefo=ipdcoefo+iterm
              go to 560
550           pcoefet=pcoefe
              ipcoefet=ipcoefe
              pcoefe=pcoefe*pratt(li+2)/pratb(li+2)
              iterm=iqint(qlog10(qabs(pcoefe)))
              pcoefe=pcoefe*10.0q0**(-iterm)
              ipcoefe=ipcoefe+iterm
              pdcoefet=pdcoefe
              ipdcoefet=ipdcoefe
              pdcoefe=pdcoefe*pdratt(li+2)/pratb(li+2)
              iterm=iqint(qlog10(qabs(pdcoefe)))
              pdcoefe=pdcoefe*10.0q0**(-iterm)
              ipdcoefe=ipdcoefe+iterm
560           continue
              lplus=max(l,500)
              limeta=2*((lplus)*(-18.5q0-20.0q0*qlog10(x1))+&
               5*ndec+4*m+c+01000)
              if(x1.gt.0.08q0) limeta=2*((lplus)*(0.50q0-&
              3.0q0*qlog10(x1))+5*ndec+4*m+c+01000)
              if(x1.gt.1.0q0) limeta=2*((lplus)*0.5q0+5*ndec+4*m+c+&
                               00500)
              if(iopeta.eq.3.and.naccrsav.gt.minacc)&
                        limeta=jeta+jeta+500+c/10
              if(iopeta.eq.3.and.naccrsav.le.minacc)&
                        limeta=jeta+jeta+500+c
              if(iopeta.eq.2) limeta=max(limeta,jeta+jeta+500+iqint(c))
              if(limeta.gt.limp-2) limeta=limp-2
              if(limeta.gt.limd) limeta=limd
              call r2eta(l,m,cc,x1,etaval,nee,limeta,ndec,maxd,maxlp,&
                   maxn,maxp,minacc,lowtest,enr(li,:),sneufe,&
                   sneune,ineuee,sneudfe,sneudre,pdratt,pratb,&
                   pratt,pcoefet,ipcoefet,pdcoefet,ipdcoefet,&
                   nsubf,nsubmf,idigc,ichoicer1,naccr1,r1c,ir1e,&
                   r1dc,ir1de,naccsub,naccr,r2ec,ir2ee,r2dec,&
                   ir2dee,nacceta,nacciop,jeta,iopnee,neemark,&
                   naccnmax,naccsube)

              netatry=netatry+1
              naccetas=nacceta
580           if(nacciop.eq.0) write(40,590)
590           format(15x,'r2eta accuracy is calculated using the',&
               ' wronskian.')
              if(nacciop.eq.1) write(40,600)
600           format(15x,'r2eta accuracy is set equal to estimated',&
               ' numerator accuracy.')
              if(nacciop.eq.2) write(40,605)
605           format(15x,'r2eta accuracy set equal to the smaller'&
                 ' of the estimted numerator and denominator',&
                 ' accuracies.')
              if(naccetas.lt.3.and.nacciop.eq.1) naccetas=0
              iopeta=3
              if(naccetas.gt.nacctemp) nacctemp=naccetas
              if(naccetas.lt.naccr) go to 610
              naccr=min(naccetas,naccr1)
              r2c=r2ec
              ir2e=ir2ee
              r2dc=r2dec
              ir2de=ir2dee
610           continue
620           if(nacciop.ne.0.or.naccsube.lt.2) write(40,625)&
           naccetas,etaval,nee,r2ec,ir2ee,r2dec,ir2dee
625           format(15x,'accuracy = ',i3,' decimal digits; eta',&
               ' = ',f12.9,'; nee = ',i4,/,10x,'r2 = ', f35.31,&
               f35.31,i5,/,10x,'r2d = ',f35.31,f35.31,i5)
              if(nacciop.eq.0.and.naccsube.ge.2) write(40,630)&
           naccetas,naccsube,etaval,nee,r2ec,ir2ee,r2dec,ir2dee
630           format(15x,'accuracy =',i3,' decimal digits,'&
                ' adjusted for',i3,' subtraction digits in'&
                ' forming Wronskian;'/,15x,' eta = ',f12.9,';'&
                ' nee = ',i4,/,10x,'r2 = ',f35.31,f35.31,2x,i6,&
                /,10x,'r2d = ',f35.31,f35.31,2x,i6)
              if(naccetas.ge.minacc) ietacount=ietacount+1
              if(naccetas.ge.naccetamax) neemax=nee
              if(naccetas.gt.naccetamax) naccetamax=naccetas
              if(ietacount.ge.5) incnflag=1
              IF(naccetas.ge.minacc.and.iplflag.eq.0) nee=nee-incnee
              if(nee.lt.1) nee=1
              if(naccetas.ge.minacc) iplflag=1
              if(naccetas.ge.minacc) go to 660

              iopeta=2
              if(iplflag.eq.1.and.incnflag.eq.1.and.netatry.eq.2)&
               iopnee=0
              ietacount=0
              if(iopnee.eq.0) go to 650
              if(iopnee.eq.2) nee=neemax-incnee
              if(iopnee.eq.1) nee=neemark-incnee
              if(nee.lt.1) nee=1
              if(iopnee.eq.2.and.naccetas.lt.lowtest-1) go to 640
              incnee=8
              if(x1.ge.0.05q0) incnee=16
              if(x1.ge.0.1q0) incnee=32
640           if(nacctemp.ge.lowtest-1) msearch=1
              if(msearch.eq.0) iopeta=0
              if(msearch.eq.0) lowtest=lowtest-1
              go to 670
650           if(nee.eq.neta) go to 660
              nee=nee+incnee
              if(nee.gt.neta) nee=neta
              if(msearch.ne.0) kounter=0
              go to 440
660           continue

              iopeta=3
              if(naccetas.lt.minacc.and.nee.eq.neta)&
               nee=nee-incnee
              if(naccetas.lt.minacc) iopeta=2
              if(nee.ne.neta) msearch=1
              if(naccetas.ge.minacc) kounter=kounter+1
              if(kounter.ge.(2*incnee).and.msearch.ne.0)&
               incnee=2*incnee
              if(incnee.gt.64) incnee=64
              if(x1.lt.0.1q0.and.incnee.gt.32) incnee=32
              if(iopint.ne.0.and.naccetas.lt.lowacc) iopeta=0
              if(iopeta.eq.0) nacctest=naccetas
              if(naccetas.lt.minacc) iplflag=0
670           if(naccr.gt.0) go to 680
              naccr=0
              r2c=1.q0
              ir2e=-ndec
              r2dc=1.q0
              ir2de=-ndec
680           if(ioprad.eq.2.and.nacciop.eq.0) write(20,690) r1c,ir1e,&
               r1dc,ir1de,r2c,ir2e,r2dc,ir2de,naccr
690           format(10x,2(f18.15,f18.15,i5,2x),/,10x,&
               2(f18.15,f18.15,i5,2x),i2,' w')
              if(ioprad.eq.2.and.nacciop.ne.0) write(20,700) r1c,ir1e,&
               r1dc,ir1de,r2c,ir2e,r2dc,ir2de,naccr
700           format(10x,2(f18.15,f18.15,i5,2x),/,10x,&
               2(f18.15,f18.15,i5,2x),i2,' e')
              if(ioprad.eq.1) write(20,710) r1c,ir1e,r1dc,ir1de
710           format(10x,2(f18.15,f18.15,i5,2x))

              qr1(li)=r1c
              ir1(li)=ir1e
              qr1d(li)=r1dc
              ir1d(li)=ir1de
              qr2(li)=r2c
              ir2(li)=ir2e
              qr2d(li)=r2dc
              ir2d(li)=ir2de
              if(lowacc.gt.naccr) lowacc=naccr
              if(istartr2.eq.1) then
              if(ndec-nsubmf.ge.minacc.and.ndec-nsubf.ge.minacc.and.&
           x1.le.0.4q0.and.iopleg.eq.0.and.l.ge.legstart)&
           iopleg=1
              if(ndec-nsubmf.ge.minacc.and.x1.ge.(0.001q0).and.&
           iopneu.eq.0.and.iopleg.ne.2) iopneu=1
              if(iopeta.eq.0.and.x1.ge.0.001q0.and.iopneu.eq.0.and.&
           iopleg.ne.2) iopeta=1
              end if
              if(nee.eq.neta.and.iopeta.ne.0) iopneu=1
              if(nee.eq.neta) iopeta=0
              naccsub=max(naccsubi,naccsubl,naccsubn,naccsube)
!             if(naccr.lt.minacc) write(45,*) l,naccr
720           if(iopang.eq.0) go to 850
!
!  determine first kind prolate angular function
              if(l.eq.m) lims1=3*ndec+iqint(c)
              if(l.ne.m) lims1=jang+jang+20+iqint(qsqrt(c))
              if(lims1.gt.maxp) lims1=maxp

              call s1leg(l,m,cc,iopang,iopnorm,barg,narg,lims1,ndec,&
                   maxt,maxd,maxp,enr(li,:),sgn,pr,pdr,pdnorm,&
                   ipdnorm,pnorm,ipnorm,pdtempe,ipdtempe,&
                   pdtempo,ipdtempo,ptempe,iptempe,ptempo,&
                   iptempo,s1c,is1e,s1dc,is1de,naccs,naccsd,jang)

                do 810 jarg=1,narg
                if(ioparg.eq.0.and.iopang.eq.1) write(50,730)&
                 arg(jarg),naccs(jarg)
                if(ioparg.eq.0.and.iopang.eq.2) write(50,735)&
                 arg(jarg),naccs(jarg),naccsd(jarg)
730             format(1x,'theta = ',e40.31,'  accuracy = ',i3,&
                 ' digits.')
735             format(1x,'theta = ',e40.31,'  accuracies = ',i3,&
                 ' and',i3,' digits.')
                if(ioparg.eq.1.and.iopang.eq.1) write(50,740)&
                 barg(jarg),naccs(jarg)
                if(ioparg.eq.1.and.iopang.eq.2) write(50,745)&
                 barg(jarg),naccs(jarg),naccsd(jarg)
740             format(1x,'eta = ',e40.31,'  accuracy = ',i3,&
                 ' digits.')
745             format(1x,'eta = ',e40.31,'  accuracies = ',i3,&
                 ' and',i3,' digits.')
                if(ioparg.eq.1.and.iopang.eq.1) write(30,750)&
                barg(jarg),s1c(jarg),is1e(jarg),naccs(jarg)
                if(ioparg.eq.1.and.iopang.eq.2) write(30,760)&
                barg(jarg),s1c(jarg),is1e(jarg),naccs(jarg),&
                s1dc(jarg),is1de(jarg),naccsd(jarg)
                if(ioparg.eq.0.and.iopang.eq.1) write(30,750)&
                arg(jarg),s1c(jarg),is1e(jarg),naccs(jarg)
                if(ioparg.eq.0.and.iopang.eq.2) write(30,760)&
                arg(jarg),s1c(jarg),is1e(jarg),naccs(jarg),&
                s1dc(jarg),is1de(jarg),naccsd(jarg)
                if(iopang.eq.1) write(50,770) s1c(jarg),is1e(jarg)
                if(iopang.eq.2) write(50,780) s1c(jarg),is1e(jarg),&
                s1dc(jarg),is1de(jarg)
750             format(1x,e30.15,2x,f18.15,f18.15,2x,i5,2x,i2)
760             format(1x,e30.15,2x,f18.15,f18.15,2x,i5,2x,i2/,&
                 33x,f18.15,f18.15,2x,i5,2x,i2)
770             format(10x,'s1 = ',f40.31,f40.31,2x,i5)
780             format(10x,'s1 = ',f40.31,f40.31,2x,i5,/,&
                 10x,'s1d = ',f40.31,f40.31,2x,i5)
                s1(li,jarg)=s1c(jarg)*(10.0q0**is1e(jarg))
                if(iopang.eq.2) s1d(li,jarg)=s1dc(jarg)*&
                                       (10.0q0**is1de(jarg))
!			write(*,*) 'here'
810             continue
!			write(*,*) 'here'
850           continue
!			write(*,*) 'here'
900       continue
          call date_and_time(date,time,zone,dt)
          time2=60*dt(6)+dt(7)+.001*dt(8)
          caltime=time2-time1
          if(caltime.lt.0.0) caltime=caltime+3600.0
          if(ioprad.ne.0) write(40,910) caltime
          if(ioprad.eq.0) write(50,910) caltime
910       format(3x,'elapsed time for this run was ',f8.1,' seconds.')

!  DEALLOCATION
!		write(*,*) 'deallocation'
		deallocate(bliste)
!				write(*,*) 'deallocation 0'

		deallocate(gliste)
		deallocate(blisto)
		deallocate(glisto)

		deallocate(drhor)

		deallocate(pint1)
		deallocate(pint2)
		deallocate(pint3)
		deallocate(pint4)
		deallocate(rpint1)
		deallocate(rpint2)
!		write(*,*) 'deallocation 1'

		deallocate(sbesf)
		deallocate(sbesf2)
		deallocate(sbesdf)
		deallocate(sbesdf2)

		deallocate(ibese)
		deallocate(ibese2)
		deallocate(ipnormint)

		deallocate(pnormint)

		deallocate(sbesdr)
		deallocate(sbesdr2)
		deallocate(sbesn)
		deallocate(sbesn2)

!		write(*,*) 'deallocation 2'

		deallocate(sneudre)
		deallocate(sneudrsv)
        deallocate(sneudr)

		deallocate(enrneg)
        deallocate(norme)

        deallocate(ineue)
		deallocate(ineuee)
		deallocate(ineuesv)

		deallocate(sneufe)
		deallocate(sneudfe)
		deallocate(sneune)
!		write(*,*) 'deallocation 3'

		deallocate(sneufsv)
		deallocate(sneudfsv)
		deallocate(sneunsv)

        deallocate(sneun)
		deallocate(sneuf)
		deallocate(sneudf)

        deallocate(alpha)
		deallocate(beta)
		deallocate(coefa)
		deallocate(coefb)
		deallocate(coefc)
		deallocate(coefd)
		deallocate(coefe)
		deallocate(gamma)
		deallocate(pdratt)
		deallocate(pratb)
		deallocate(pratt)
		deallocate(prat1)
!				write(*,*) 'deallocation 4'

		deallocate(pdr)
		deallocate(pdrat)
		deallocate(pr)
		deallocate(prat)
		deallocate(pratbsv)
		deallocate(prattsv)
		deallocate(pdratsv)

		deallocate(prx)
		deallocate(pdrx)

		deallocate(qr)
		deallocate(qdr)

		deallocate(barg)
		deallocate(etainp)
		deallocate(pdnorm)
		deallocate(pdnorma)
		deallocate(pnorm)
		deallocate(pnorma)
		deallocate(pdtempe)
		deallocate(pdtempo)
		deallocate(ptempe)
		deallocate(ptempo)
		deallocate(xin)
		deallocate(xlninp)
!				write(*,*) 'deallocation 5'

        deallocate(s1c)
		deallocate(s1dc)

        deallocate(ipdnorm)
		deallocate(ipdnorma)
		deallocate(ipnorm)
		deallocate(ipnorma)
		deallocate(ipdtempe)
		deallocate(ipdtempo)
		deallocate(iptempe)
		deallocate(iptempo)
		deallocate(is1e)
		deallocate(is1de)
		deallocate(naccs)
		deallocate(naccsd)
!				write(*,*) 'deallocation 6'

		deallocate(eta)
		deallocate(xbn)
		deallocate(xln)

		deallocate(wr)
		deallocate(xr)

		deallocate(neeb)
		deallocate(limpsv)
		deallocate(limnsv)
		deallocate(jelimsv)
!				write(*,*) 'deallocation 7'

        return
        end
!
!
!
        subroutine s1leg (l,m,c,iopang,iopnorm,barg,narg,lims1,ndec,&
                    maxt,maxd,maxp,enr,sgn,pr,pdr,pdnorm,ipdnorm,&
                    pnorm,ipnorm,pdtempe,ipdtempe,pdtempo,&
                    ipdtempo,ptempe,iptempe,ptempo,iptempo,s1c,&
                    is1e,s1dc,is1de,naccs,naccsd,jang)
!
!  purpose:     to calculate the prolate angular functions of the first
!               kind and their first derivatives with respect to eta.
!
!  parameters:
!
!     input :   l       : l
!               m       : m
!               c       : c
!               iopang  : index = 1 when angular functions of the
!                         first kind are calculated; = 2 when the
!                         first derivatives with respect to eta are
!                         also calculated
!               iopnorm : = 1 when the angular functions (and
!                         first derivatives) are scaled by the
!                         square root of the normalization of the
!                         corresponding Legendre function, giving them
!                         unity norm; iopnorm = 0 otherwise
!               barg    : array of eta values for which angular
!                         functions are desired
!               narg    : number of eta values
!               lims1   : approximately twice the maximum number
!                         of terms available to be taken in the
!                         sums
!               ndec    : number of decimal digits available in
!                         r*16 arithmetic
!               maxt    : dimension of barg, pdnorm, ipdnorm, pnorm,
!                         ipnorm, pdtempe, ipdtempe, pdtempo, ipdtempo,
!                         ptempe, iptempe, ptempo, iptempo, s1c, is1e,
!                         s1dc, is1de, and naccs arrays. first dimension
!                         of the doubly dimensioned arrays pr and pdr
!               maxd    : dimension of enr array
!               maxp    : second dimension of pr and pdr arrays
!               enr     : array of d coefficient ratios
!               sgn     : sign of d constant multiplying the
!                         associated Legendre of order m and
!                         degree l in the series for s1 and s1d
!               pr      : array of ratios of successive first kind
!                         associated Legendre functions of the same
!                         parity
!               pdr     : array of ratios of successive derivatives of
!                         first kind associated Legendre functions of
!                         the same parity
!               pdnorm  : array of characteristics of the first
!                         derivatives of associated Legendre functions
!                         of the first kind of order m and degree m
!               ipdnorm : array of exponents corresponding to pdnorm
!               pnorm   : array of characteristics of the associated
!                         Legendre functions of the first kind of order
!                         m and degree m
!               ipnorm  : array of exponents corresponding to pnorm
!               pdtempe : storage array of characteristics of the ratio
!                         of the first derivative of the associated
!                         Legendre function of order m and degree l - 2
!                         or l - 1, depending on whether l - m is even
!                         or odd, to the first derivative of the
!                         function of order m and degree m
!               ipdtempe: array of exponents corresponding to pdtempe
!               pdtempo : storage array of characteristics of the ratio
!                         of the first derivative of the associated
!                         Legendre function of order m and degree l - 2
!                         or l - 1, depending on whether l - m is odd
!                         or even, to the first derivtive of the
!                         function of order m and degree m
!               idptempo: array of exponents corresponding to pdtempo
!               ptempe  : storage array of characteristics of the ratio
!                         of the associated Legendre function of order
!                         m and degree l - 2 or l - 1, depending on
!                         whether l - m is even or odd, to the function
!                         of order m and degree m
!               iptempe : array of exponents corresponding to ptempe
!               ptempo  : storage array of characteristics of the ratio
!                         of the associated Legendre function of order
!                         m and degree l - 2 or l - 1, depending on
!                         whether l - m is odd or even, to the
!                         function of order m and degree m
!               iptempo : array of exponents corresponding to ptempo
!
!
!     output:   s1c    : array of characteristics of prolate
!                        angular functions of the first kind
!               is1e   : array of exponents of prolate angular
!                        functions of the first kind
!               s1dc   : array of characteristics of derivative with
!                        respect to eta of prolate angular functions
!                        of the first kind
!               is1de  : array of exponents of derivatives with respect
!                        to eta of prolate angular functions of first
!                        kind
!               naccs  : array of integer estimates of the number of
!                        accurate decimal digits in the values obtained
!                        for s1
!               naccsd : array of integer estimates of the number of
!                        accurate decimal digits in the values obtained
!                        for s1d
!               jang   : maximum value of the index j in the forward
!                        sum for r1 and r1d, i.e., the highest enr(j)
!                        used
!
!  real*16 scalars and arrays
        real*16 adec,aj,aj2,dcon,dec,factor,fterm,rm2,rm2m1,rm2m3,&
          rm2p1,sgn
        complex*32 c,coef,dnew,dnewd,dnum,dold,doldd,s1,s1d
        real*16 barg(maxt),pdr(maxt,maxp),pdnorm(maxt),&
          pnorm(maxt),pr(maxt,maxp),pdtemp(maxt),ptemp(maxt),&
          pdtempe(maxt),ptempe(maxt),pdtempo(maxt),ptempo(maxt)
        complex*32 enr(maxd),s1c(maxt),s1dc(maxt)
!
!  integer arrays
        dimension ipdnorm(maxt),ipnorm(maxt),ipdtemp(maxt),iptemp(maxt),&
            ipdtempe(maxt),iptempe(maxt),ipdtempo(maxt),&
            iptempo(maxt),is1de(maxt),is1e(maxt),naccs(maxt),&
            naccsd(maxt)
!
        dec=10.q0**(-ndec-1)
        dcon=dec
        adec=1000.q0*dec
        rm2=qfloat(m+m)
        rm2m1=qfloat(m+m-1)
        rm2p1=qfloat(m+m+1)
        rm2m3=qfloat(m+m-3)
        if(l.gt.(m+1)) go to 30
          do 20 k=1,narg
          if(pnorm(k).eq.0.q0) go to 20
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
          if(pnorm(k).eq.0.q0) go to 100
          if(ix.ne.0) go to 60
          ptempe(k)=ptempe(k)*pr(k,l-m+1)
          if(qabs(ptempe(k)).lt.1.q+10) go to 40
          ptempe(k)=ptempe(k)*(1.q-10)
          iptempe(k)=iptempe(k)+10
40        ptemp(k)=ptempe(k)
          iptemp(k)=iptempe(k)
          if(qabs(barg(k)).lt.adec) go to 100
          pdtempe(k)=pdtempe(k)*pdr(k,l-m+1)
          if(qabs(pdtempe(k)).lt.1.q+10) go to 50
          pdtempe(k)=pdtempe(k)*(1.q-10)
          ipdtempe(k)=ipdtempe(k)+10
50        pdtemp(k)=pdtempe(k)
          ipdtemp(k)=ipdtempe(k)
          go to 100
60        if(qabs(barg(k)).lt.adec) go to 80
          ptempo(k)=ptempo(k)*pr(k,l-m+1)
          if(qabs(ptempo(k)).lt.1.q+10) go to 70
          ptempo(k)=ptempo(k)*(1.q-10)
          iptempo(k)=iptempo(k)+10
70        ptemp(k)=ptempo(k)
          iptemp(k)=iptempo(k)
80        pdtempo(k)=pdtempo(k)*pdr(k,l-m+1)
          if(qabs(pdtempo(k)).lt.1.q+10) go to 90
          pdtempo(k)=pdtempo(k)*(1.q-10)
          ipdtempo(k)=ipdtempo(k)+10
90        pdtemp(k)=pdtempo(k)
          ipdtemp(k)=ipdtempo(k)
100        continue
110     continue
        lim=lims1/2-ix
!
!  compute the normalizing factor
        dnum=(1.0q0,0.0q0)
        coef=(1.0q0,0.0q0)
        jlow=l-m+2
        jterm=lm2
          do 120 j=jlow,lims1,2
          aj=j
          aj2=aj+aj
          jterm=jterm+1
          coef=coef*(aj+rm2)*enr(jterm)*(aj+rm2m1)*enr(jterm)&
         *(aj2+rm2m3)/(aj*(aj-1.q0)*(aj2+rm2p1))
          dnum=dnum+coef
          if(cqabs(coef/dnum).lt.dcon) go to 130
120       continue
130     jlow=l-m
        jn=jterm
        if(jlow.lt.2) go to 150
        coef=(1.0q0,0.0q0)
        jterm=lm2
        j=jlow
          do 140 jj=2,jlow,2
          aj=j
          aj2=aj+aj
          coef=coef*aj*(aj-1.q0)*(aj2+rm2p1)/((aj2+rm2m3)*enr(jterm)&
         *enr(jterm)*(aj+rm2)*(aj+rm2m1))
          jterm=jterm-1
          j=j-2
          dnum=dnum+coef
          if(cqabs(coef/dnum).lt.dcon) go to 150
140       continue
150     dnum=(1.q0,0.0q0)/cqsqrt(dnum)
        write(50,155) jn,lim
155     format(5x,'Meixner-Schafke normalization series converged in '&
           ,i6,' terms; ',i6,' terms available.')
        jlow=lm2+1
        jang=0
!
!  compute the associated Legendre function normalization factor
        factor=1.0q0
        ifactor=0
        if(iopnorm.eq.0) go to 210
        if(m.eq.0) go to 170
          do 160 j=1,m
          aj=j
          factor=factor*(aj+aj)*(aj+aj-1.q0)
          if(factor.lt.1.q100) go to 160
          factor=factor*1.q-100
          ifactor=ifactor+100
160       continue
170     if(l.eq.m) go to 190
          do 180 j=1,l-m
          aj=j
          factor=factor*(rm2+aj)/(aj)
          if(factor.lt.1.q100) go to 180
          factor=factor*1.q-100
          ifactor=ifactor+100
180       continue
190     factor=factor*2.q0/(l+l+1.q0)
        factor=qsqrt(factor)
        ifactor=ifactor/2
        iterm=iqint(qlog10(factor))
        factor=factor*(10.q0**(-iterm))
        ifactor=ifactor+iterm
        write(50,200)
200     format(5x,'s1 is normalized to have unit norm.')
210     continue
        if(iopnorm.eq.0) write(50,215)
215     format(5x,'s1 has the same normalization as the',&
         ' corresponding Legendre function.')


!
!  compute the angular function s1
          do 380 k=1,narg
          if(pnorm(k).eq.0.q0) go to 220
          if((ix.eq.1).and.(qabs(barg(k)).lt.adec)) go to 220
          if(((qabs(qabs(barg(k))-1.q0)).lt.adec).and.(m.ne.0)) go to 220
          go to 230
220       s1c(k)=(0.0q0,0.0q0)
          is1e(k)=0
          naccs(k)=ndec
          go to 300
230       dold=(1.0q0,0.0q0)
          s1=dold
          fterm=1.0q0
            do 240 j=jlow,lim
            dnew=dold*enr(j)*pr(k,j+j+ixx2)
            s1=s1+dnew
            if(cqabs(dnew).gt.fterm) fterm=cqabs(dnew)
            if(cqabs(dnew/s1).lt.dcon) go to 250
            dold=dnew
240         continue
250       if(j.gt.jang) jang=j
          write(50,260) barg(k),j
260       format(8x,'s1 calculation for eta = ',f13.8,' converged in ',&
           i6,' terms')
          if(lm2.lt.1) go to 280
          dold=(1.q0,0.0q0)
          j=lm2
            do 270 jj=1,lm2
            dnew=dold/(pr(k,j+j+ixx2)*enr(j))
            s1=s1+dnew
            if(cqabs(dnew).gt.fterm) fterm=cqabs(dnew)
            if(cqabs(dnew/s1).lt.dcon) go to 280
            dold=dnew
            j=j-1
270         continue
280       s1c(k)=s1*dnum*(ptemp(k)*pnorm(k)*sgn/factor)
          if(s1c(k).ne.(0.0q0,0.0q0)) iterm=iqint(qlog10(cqabs(s1c(k))))
          if(s1c(k).eq.(0.0q0,0.0q0)) iterm=0
          s1c(k)=s1c(k)*(10.0q0**(-iterm))
          is1e(k)=iptemp(k)+ipnorm(k)+iterm-ifactor
          if(cqabs(s1c(k)).ge.1.q0) go to 290
          s1c(k)=s1c(k)*10.0q0
          is1e(k)=is1e(k)-1
290       if(s1.eq.(0.0q0,0.0q0)) naccs(k)=0
          if(s1.ne.(0.0q0,0.0q0)) naccs(k)=ndec-2-&
        qlog10(cqabs(fterm/s1))
          if(naccs(k).lt.0) naccs(k)=0
          if(naccs(k).gt.0) go to 300
          naccs(k)=0
          s1c(k)=(0.0q0,0.0q0)
          is1e(k)=0
          s1dc(k)=(0.0q0,0.0q0)
          is1de(k)=0
          go to 380
!
!       compute the first derivative of the angular function when
!       iopang equals 2
300       if(iopang.ne.2) go to 380
          if(pnorm(k).eq.0.q0) go to 310
          if((ix.eq.0).and.(qabs(barg(k)).lt.adec)) go to 310
          if(((qabs(qabs(barg(k))-1.q0)).lt.adec).and.(m.ne.0)&
        .and.(m.ne.2)) go to 310
          go to 320
310       s1dc(k)=(0.0q0,0.0q0)
          is1de(k)=0
          naccsd(k)=ndec
          go to 370
320       doldd=(1.0q0,0.0q0)
          s1d=doldd
          fterm=1.0q0
          if(l.eq.0) s1d=(0.0q0,0.0q0)
            do 330 j=jlow,lim
            dnewd=doldd*enr(j)*pdr(k,j+j+ixx2)
            s1d=s1d+dnewd
            if(cqabs(dnewd).gt.fterm) fterm=cqabs(dnewd)
            if(cqabs(dnewd/s1d).lt.dcon) go to 340
            doldd=dnewd
330         continue
340       if(lm2.lt.1) go to 360
          doldd=(1.0q0,0.0q0)
          j=lm2
          ja=lm2
          if(m.eq.0.and.ix.eq.0) ja=lm2-1
          if(ja.eq.0) go to 360
            do 350 jj=1,ja
            dnewd=doldd/(pdr(k,j+j+ixx2)*enr(j))
            s1d=s1d+dnewd
            if(cqabs(dnewd).gt.fterm) fterm=cqabs(dnewd)
            if(cqabs(dnewd/s1d).lt.dcon) go to 360
            doldd=dnewd
            j=j-1
350         continue
360       s1dc(k)=s1d*dnum*(pdtemp(k)*pdnorm(k)*sgn/factor)
          if(s1d.eq.(0.0q0,0.0q0)) naccsd(k)=0
          if(s1d.ne.(0.0q0,0.0q0)) naccsd(k)=ndec-2-&
        qlog10(cqabs(fterm/s1d))
          if(naccsd(k).lt.0) naccsd(k)=0
          if(s1dc(k).ne.(0.0q0,0.0q0))&
         iterm=iqint(qlog10(cqabs(s1dc(k))))
          if(s1dc(k).eq.(0.0q0,0.0q0)) iterm=0
          s1dc(k)=s1dc(k)*10.0q0**(-iterm)
          is1de(k)=ipdtemp(k)+ipdnorm(k)+iterm-ifactor
          if(cqabs(s1dc(k)).ge.1.0q0) go to 370
          s1dc(k)=s1dc(k)*10.0q0
          is1de(k)=is1de(k)-1
370       continue
380       continue
        return
        end
!
!
!
        subroutine r1besa (l,m,c,x1,limr1,ndec,maxd,enr,maxj,maxlp,&
                     nex,iflag,sbesf,sbesdf,sbesn,ibese,sbesdr,&
                     d01,id01,r1c,ir1e,r1dc,ir1de,dfnorm,jbesa,&
                     factor,nsubr1a)
!
!  purpose:     to calculate the prolate radial function of the
!               first kind and its first derivative with respect
!               to x, using an expansion of spherical Bessel
!               functions of the first kind with argument
!               c*sqrt(x*x-1).
!
!  parameters:
!
!     input:    l      : l
!               m      : m
!               c      : complex c
!               x1     : x-1
!               limr1  : approximately twice the maximum number of
!                        terms available to be taken in the series
!               ndec   : number of decimal digits available in
!                        r*16 arithmetic
!               maxd   : dimension of enr array
!               enr    : d coefficient ratios
!               maxj   : dimension of sbesf and sbesdf arrays
!               maxlp  : maximum  l value desired; dimension
!                        of the sbesdr, sbesn, and ibese arrays
!               nex    : maximum exponent available in real*16
!                        arithmetic
!               iflag  : integer = 1 if forward series not needed;
!                        =0 if the forward series is computed
!               sbesf  : array of ratios of successive first kind
!                        spherical Bessel functions of the same parity
!               sbesdf : array of ratios of successive derivatives of
!                        first kind spherical Bessel functions of the
!                        same parity
!               sbesn  : array of characteristics for Bessel functions
!               ibese  : array of exponents corresponding to sbesn
!               sbesdr : value of ratio of first derivative of
!                        spherical Bessel function to the corresponding
!                        Bessel function
!               d01    : characteristic of the ratio of the first d
!                        coefficient, either d0 or d1, to the d
!                        coefficient for n = l - m
!               id01   : exponent for d01
!
!     output  : r1c    : characteristic of prolate radial function
!                        of the first kind
!               ir1e   : exponent of prolate radial function of the
!                        first kind
!               r1dc   : characteristic of derivative with respect
!                        to x of prolate radial function of the first
!                        kind
!               ir1de  : exponent of derivative with respect to x of
!                        prolate radial function of the first kind
!               dfnorm : Flammer normalization sum for the d constants;
!                        computed below as the denominator series
!               jbesa  : maximum value of the index j in the forward
!                        sum for r1 and r1d, i.e., the highest enr(j)
!                        used
!               factor : coefficient used in computing one of the two
!                        contributions to r1d for l odd and m equal to
!                        0. The value used for the previous odd l is
!                        input to the subroutine. The value for the
!                        current odd value of l is computed and output
!                        to be availabe for the next odd value of l.
!               nsubr1a: maximum number of digits of subtraction error
!                        in calculating r1 and r1d using this subroutine
!
!  real*16 and complex*32 scalars and arrays
        real*16 con,dec,factor,r1rp,r1ip,r1drp,r1dip,rj,teste,testeo,x1
        complex*32 c,dfnorm,dnew,dnewd,dold,doldd,d01,r1bot,r1c,r1d,&
             r1dc,r1dstore,r1sig,r1temp,r1top,r1dtop,r1topd,&
             r1dtopd,term,termd
        complex*32 enr(maxd),sbesdf(maxj),sbesdr(maxlp),sbesf(maxj),&
             sbesn(maxlp)
!
!  integer array
        dimension ibese(maxlp)
!
!  convergence ratio dec is set according to the requested accuracy
        dec=10.q0**(-ndec-1)
        lm2=(l-m)/2
!
!  ix=0 for l-m even, ix=1 for l-m odd
        ix=l-m-2*lm2
        lim=limr1/2-ix
        con=(x1+1.q0)/(qsqrt(x1)*qsqrt((x1+2.q0)))
        nfac=nex-ndec
        teste=10.q0**nfac
        testeo=1.q0/teste
        ir1tope=0
        mml=m+m-1+ix
        iflagd=0
        if(x1.lt.0.1q0.and.ix.eq.1.and.m.eq.0) iflagd=1
        if(iflagd.eq.1.and.l.ne.1) factor=factor*qfloat(l)/qfloat(l-1)
!
!
!  compute radial function of the first kind r1 and its first
!  derivative r1d
!
!  forward summation of numerator series
        dold=(1.0q0,0.0q0)
        doldd=(1.0q0,0.0q0)
        r1rp=1.0q0
        r1drp=1.0q0
        r1ip=0.0q0
        r1dip=0.0q0
        r1top=dold
        r1dtop=doldd
        if(iflag.eq.1) go to 50
          do 20 j=lm2+1,lim
          jj=j+j+ix
          rj=qfloat(jj)
          dnew=dold*enr(j)*sbesf(jj+m)*qfloat((jj+mml))/qfloat((jj-ix))
          dnewd=doldd*enr(j)*sbesdf(jj+m)*qfloat((jj+mml))/&
          qfloat((jj-ix))
          r1top=r1top+dnew
          r1dtop=r1dtop+dnewd
          if(qreal(dnew).gt.0.0q0) r1rp=r1rp+qreal(dnew)
          if(qreal(dnewd).gt.0.0q0) r1drp=r1drp+qreal(dnewd)
          if(qimag(dnew).gt.0.0q0) r1ip=r1ip+qimag(dnew)
          if(qimag(dnewd).gt.0.0q0) r1dip=r1dip+qimag(dnewd)
          if((cqabs(dnew/r1top)+cqabs(dnewd/r1dtop)).lt.dec) go to 30
          dold=dnew
          doldd=dnewd
20        continue
30      continue
        jbesa=min(j,lim)
        if(iflagd.eq.0.or.l.ne.1) go to 40
        r1topd=r1top-1.0q0
        r1dtopd=r1dtop-1.0q0
40      continue
50	continue
!
!  backward summation of numerator series
        if (lm2.lt.1) go to 80
        dold=(1.0q0,0.0q0)
        doldd=(1.0q0,0.0q0)
          do 70 j=lm2,1,-1
          jj=j+j+ix
          rj=qfloat(jj)
          dnew=dold*(jj-ix)/(sbesf(jj+m)*qfloat((jj+mml))*enr(j))
          dnewd=doldd*(jj-ix)/(sbesdf(jj+m)*qfloat((jj+mml))*enr(j))
          r1top=r1top+dnew
          r1dtop=r1dtop+dnewd
          if(qreal(dnew).gt.0.0q0) r1rp=r1rp+qreal(dnew)
          if(qreal(dnewd).gt.0.0q0) r1drp=r1drp+qreal(dnewd)
          if(qimag(dnew).gt.0.0q0) r1ip=r1ip+qimag(dnew)
          if(qimag(dnewd).gt.0.0q0) r1dip=r1dip+qimag(dnewd)
          if((cqabs(dnew/r1top)+cqabs(dnewd/r1dtop)).lt.dec) go to 80
          if(cqabs(r1top).lt.teste) go to 60
          r1top=r1top*testeo
          r1rp=r1rp*testeo
          r1ip=r1ip*testeo
          dnew=dnew*testeo
          ir1tope=ir1tope+nfac
          r1dtop=r1dtop*testeo
          r1drp=r1drp*testeo
          r1dip=r1dip*testeo
          dnewd=dnewd*testeo
60        dold=dnew
          doldd=dnewd
70        continue
          if(jj.ne.3) iflagd=0
          if(iflagd.eq.0) go to 80
          r1topd=r1top-dnew
          r1dtopd=r1dtop-dnewd
          go to 90
80        continue
          iflagd=0
90        nsubr=0
          if(r1rp.ne.0.0q0) nsubr=-qlog10(qabs(qreal(r1top)/r1rp)+dec)
          if(nsubr.lt.0) nsubr=0
          nsubi=0
          if(r1ip.ne.0.0q0) nsubi=-&
             qlog10(qabs(qimag(r1top)/r1ip)+dec)
          if(nsubi.lt.0) nsubi=0
          nsubdr=0
          if(r1drp.ne.0.0q0)&
          nsubdr=-qlog10(qabs(qreal(r1dtop)/r1drp)+dec)
          if(nsubdr.lt.0) nsubdr=0
          nsubdi=0
          if(qimag(r1dtop).ne.0.0q0) nsubdi=-&
                        qlog10(qabs(qimag(r1dtop)/r1dip)+dec)
          if(nsubdi.lt.0) nsubdi=0
          r1bot=dfnorm
!
!  compute r1 and r1d
        r1temp=r1top*sbesn(l+1)/r1bot
        if(ix.eq.1) r1temp=r1temp*con
        iterm=qlog10(cqabs(r1temp))
        ir1e=ir1tope+ibese(l+1)+iterm
        r1c=r1temp*(10.0q0**(-iterm))
        if(cqabs(r1c).ge.1.0q0) go to 100
        r1c=r1c*10.q0
        ir1e=ir1e-1
100     if(iflagd.eq.1) r1temp=r1temp*r1topd/r1top
        if(iflagd.eq.1) r1dtop=r1dtopd
        r1d=r1dtop*sbesn(l+1)*c*con*sbesdr(l+1)/r1bot
        r1dstore=r1d*con
        if(ix.eq.1) r1d=r1d*con-r1temp/(x1*(x1+1.q0)*(x1+2.q0))
        ndsubr=0
        if(ix.eq.1) ndsubr=qlog10(qabs(qreal(r1dstore)/qreal(r1d))+dec)
        if(ndsubr.lt.0) ndsubr=0
        ndsubi=0
        if(ix.eq.1.and.qimag(r1d).ne.0.0q0)&
       ndsubi=qlog10(qabs(qimag(r1dstore)/qimag(r1d))+dec)
        if(ndsubi.lt.0) ndsubi=0
        if(iflagd.eq.0) go to 110
        term=x1*(x1+2.0q0)*sbesdr(2)*sbesn(2)
        termd=term-sbesn(3)*(10.0q0**(ibese(3)-ibese(2)))
        ndsub1r=qlog10(qabs(qreal(term)/qreal(termd))+dec)
        if(ndsub1r.lt.0) ndsub1r=0
        ndsub1i=0
        if(qimag(termd).ne.0.0q0)&
          ndsub1i=qlog10(qabs(qimag(term)/qimag(termd))+dec)
        if(ndsub1i.lt.0) ndsub1i=0
        termd=termd*(c*d01/(x1*(x1+2.0q0)*r1bot*factor))*&
       (10.0q0**(id01+ibese(2)-ibese(l+1)-ir1tope))
        r1d=r1d+termd
        ndsub2r=qlog10(qabs(qreal(termd)/qreal(r1d))+dec)
        if(ndsub2r.lt.0) ndsub2r=0
        ndsub2i=0
        if(qimag(r1d).ne.0.0q0)&
          ndsub2i=qlog10(qabs(qimag(termd)/qimag(r1d))+dec)
        if(ndsub2i.lt.0) ndsub2i=0
        ndsubr=ndsubr+ndsub1r+ndsub2r
        if(ndsubr.lt.0) ndsubr=0
        if(ndsubr.gt.ndec) ndsub=ndec
        ndsubi=ndsubi+ndsub1i+ndsub2i
        if(ndsubi.lt.0) ndsubi=0
        if(ndsubi.gt.ndec) ndsubi=ndec
110     continue
        nsubdr=nsubdr+ndsubr
        nsubdi=nsubdi+ndsubi
        nsub=max(nsubr,nsubi)
        nsubd=max(ndsubr,ndsubi)
        if(iflag.eq.0) write(40,120) jbesa,lim,nsub,nsubd
120     format(8x,'r1besa: numerator conv. in ',i6,' terms with',i6,&
         ' available;',/,15x,i3,' digits of subtraction error',&
         ' in r1,',i3,' digits in r1d')
        if(iflag.eq.1) write(40,130) jbesa,lim,nsub,nsubd
130     format(8x,'r1besa: numerator conv. in ',i6,' terms with',i6,&
         ' available;',i3,' digits of'/,15x,' subtraction error',&
         ' in r1,',i3,' digits in r1d. forward sum not required.')
        if(ix.eq.1.and.ndsubr.ne.0.or.ndsubi.ne.0) write(40,140) ndsubr,&
                                                          ndsubi
140     format(16x,'r1d subt. error includes',i3,' digits',&
         ' from adding its two components.')
        nsubr1a=max(nsubi,nsubr,nsubdi,nsubdr)
        if(nsubr1a.gt.ndec) nsubr1a=ndec
        iterm=qlog10(cqabs(r1d))
        ir1de=ir1tope+ibese(l+1)+iterm
        r1dc=r1d*(10.0q0**(-iterm))
        if(cqabs(r1dc).ge.1.0q0) go to 150
        r1dc=r1dc*10.q0
        ir1de=ir1de-1
150     continue
        mfac=ir1e-ibese(l+1)
        if(mfac.gt.(ndec+5)) iflag=1
        if(mfac.le.(ndec+5)) iflag=0
        return
        end
!
!
!
        subroutine r1besb (l,m,c,x1,limbes,ndec,maxd,maxlp,maxj,maxp,&
                     enr,sbesf,sbesn,ibese,sbesdf,sbesdr,prat1,&
                     pcoefn,ipcoefn,dmfnorm,r1c,ir1e,r1dc,ir1de,&
                     jbesb,nsubr1b)
!
!  purpose:     to calculate the prolate radial function of the
!               first kind and its first derivative with respect
!               to x, using the traditional expansions in terms of
!               spherical Bessel functions.
!
!  parameters:
!
!     input:    l      : l
!               m      : m
!               c      : c
!               x1     : x-1
!               limbes : maximum number of terms to be taken in the
!                        series summations for r1 and r1d
!               ndec   : number of decimal digits available in
!                        r*16 arithmetic
!               maxd   : dimension of enr array
!               maxlp  : maximum  l value desired; dimension
!                        of the sbesn and ineue arrays
!               maxj   : dimension of sbesf and sbesdr arrays
!               maxp   : dimension of prat1 array
!               enr    : array of ratios of successive d coefficients
!               sbesf  : array of ratios of successive spherical Bessel
!                        functions of the same parity
!               sbesn  : array of characteristics for Bessel functions
!               ibese  : array of exponents for Bessel functions
!               sbesdf : array of ratios of successive first derivatives
!                        of spherical Bessel functions of same parity
!               sbesdr : array of ratios of first derivatives of Bessel
!                        functions to the corresponding Bessel functions
!               prat1  : array of ratios of successive coefficients in
!                        r1 and r1d sum
!               pcoefn : characteristic of coefficient for term in both
!                        r1 and r1d sums that contains Bessel function
!                        of order l + m
!               ipcoefn: exponent (to the base 10) corresponding to
!                        pcoefn
!               dmfnorm: Morse-Feshbach normalization factor of the
!                        d coefficients. equal to the reciprocal of the
!                        value of the d constant d(n = l - m) using
!                        this normalization for the angular functions
!
!     output:   r1c    : characteristic of prolate radial function
!                        of the first kind
!               ir1e   : exponent of prolate radial function of the
!                        first kind
!               r1dc   : characteristic of derivative with respect
!                        to x of prolate radial function of the first
!                        kind
!               ir1de  : exponent of derivative with respect to x of
!                        prolate radial function of the first kind
!               jbesb  : index of term where convergence is
!                        achieved for r1 or for r1d, whichever term is
!                        larger
!               nsubr1b: maximum number of digits of subtraction error
!                        in calculating r1 and r1d using this subroutine
!
!  real*16 scalars and arrays
        real*16 dcon,dec,dnewi,dnewr,dnewdi,dnewdr,pcoefn,rm,rm2,&
          r1dcoef,r1dip,r1ip,r1dci,r1dcr,r1drp,r1rp,x1
        complex*32 c,dmfnorm,dnew,dnewd,dold,doldd,r1c,r1dc,&
             r1dc1,r1dc2,r1dtemp,r1temp
        real*16 prat1(maxp)
        complex*32 enr(maxd),sbesdr(maxlp),sbesn(maxlp),&
             sbesf(maxj),sbesdf(maxj)
!
!  integer arrays
        dimension ibese(maxlp)
!
        rm=qfloat(m)
        dec=10.0q0**(-ndec)
        dcon=dec
        r1dcoef=rm/((x1+1.q0)*(x1+2.q0)*x1)
        rm2=rm*2.q0
        lm2=(l-m)/2
!
!  ix = 0 for l-m even; ix = 1 for l-m odd
        ix=l-m-2*lm2
        lim=limbes/2-ix
!
!  compute radial function of the first kind
!
!  backward series
        r1temp=(1.0q0,0.0q0)
        r1dtemp=(1.0q0,0.0q0)
        r1rp=1.0q0
        r1ip=0.0q0
        r1drp=1.0q0
        r1dip=0.0q0
        if (lm2.lt.1) go to 20
        dold=(1.0q0,0.0q0)
        doldd=(1.0q0,0.0q0)
          do 10 j=lm2,1,-1
          jj=j+j+ix
          dnew=-dold/(sbesf(jj+m)*prat1(jj+1)*enr(j))
          dnewd=-doldd/(sbesdf(jj+m)*prat1(jj+1)*enr(j))
          r1temp=r1temp+dnew
          r1dtemp=r1dtemp+dnewd
          dnewr=qreal(dnew)
          if(dnewr.gt.0.0q0) r1rp=r1rp+dnewr
          dnewi=qimag(dnew)
          if(dnewi.gt.0.0q0) r1ip=r1ip+dnewi
          dnewdr=qreal(dnewd)
          if(dnewdr.gt.0.0q0) r1drp=r1drp+dnewdr
          dnewdi=qimag(dnewd)
          if(dnewdi.gt.0.0q0) r1dip=r1dip+dnewdi
          dold=dnew
          doldd=dnewd
10      continue
20      continue
!
!  forward series
        dold=(1.0q0,0.0q0)
        doldd=(1.0q0,0.0q0)
          do 30 j=lm2+1,lim
          jj=j+j+ix
          dnew=-dold*enr(j)*sbesf(jj+m)*prat1(jj+1)
          dnewd=-doldd*enr(j)*sbesdf(jj+m)*prat1(jj+1)
          r1temp=r1temp+dnew
          r1dtemp=r1dtemp+dnewd
          dnewr=qreal(dnew)
          if(dnewr.gt.0.0q0) r1rp=r1rp+dnewr
          dnewi=qimag(dnew)
          if(dnewi.gt.0.0q0) r1ip=r1ip+dnewi
          dnewdr=qreal(dnewd)
          if(dnewdr.gt.0.0q0) r1drp=r1drp+dnewdr
          dnewdi=qimag(dnewd)
          if(dnewdi.gt.0.0q0) r1dip=r1dip+dnewdi
          if((cqabs(dnew/r1temp)+cqabs(dnewd/r1dtemp)).lt.dcon) go to 40
          dold=dnew
          doldd=dnewd
30        continue
40      continue
        jbesb=min(j,lim)
!
!  combining results to form the radial function characteristics
!  r1c and r1dc and corresponding exponents ir1e and ir1de
        r1c=r1temp*sbesn(l+1)*pcoefn/dmfnorm
        iterm=iqint(qlog10(cqabs(r1c)))
        ir1e=ibese(l+1)+ipcoefn+iterm
        r1c=r1c*10.0q0**(-iterm)
        if(cqabs(r1c).ge.1.0q0) go to 50
        r1c=r1c*10.0q0
        ir1e=ir1e-1
50	continue
        isubr=-iqint(qlog10(qabs(qreal(r1temp)/r1rp)+dec))
        if(isubr.lt.0) isubr=0
        isubi=0
        if(r1ip.ne.0.0q0) isubi=-iqint(qlog10(qabs(qimag(r1temp)&
                             /r1ip)+dec))
        if(isubi.lt.0) isubi=0
        isubdr=-iqint(qlog10(qabs(qreal(r1dtemp)/r1drp)+dec))
        if(isubdr.lt.0) isubdr=0
        isubdi=0
        if(r1dip.ne.0.0q0) isubdi=-iqint(qlog10(qabs(qimag(r1dtemp)&
                               /r1dip)+dec))
        if(isubdi.lt.0) isubdi=0
        r1dc1=r1dcoef*r1c
        r1dc2=(c*r1dtemp*sbesn(l+1)*sbesdr(l+1)*pcoefn/&
       dmfnorm)*10.0q0**(ibese(l+1)+ipcoefn-ir1e)
        r1dc=r1dc1+r1dc2
        r1dcr=qreal(r1dc)
        r1dci=qimag(r1dc)
        isubdr1=-iqint(qlog10(qabs(r1dcr/qreal(r1dc1))+dec))
        isubdi1=0
        if(qimag(r1dc1).ne.0.0q0) isubdi1=-iqint(qlog10(qabs(r1dci&
                               /qimag(r1dc1))+dec))
        if(isubdi1.lt.0) isubdi1=0
        if(isubdr1.lt.0) isubdr1=0
        isubdr=isubdr+isubdr1
        isubdi=isubdi+isubdi1
        isub=max(isubr,isubi)
        isubd=max(isubdr,isubdi)
        nsubr1b=max(isubr,isubi,isubdr,isubdi)
        if(nsubr1b.gt.ndec) nsubr1b=ndec
        write(40,60) jbesb,lim,isub,isubd
60      format(8x,'r1besb: numerator conv. in ',i6,' terms with',i6,&
         ' available.',i3,' digits of sub. error',/,15x,&
         ' in r1 numerator. ',i3,' digits of sub. error in',&
         ' r1d numerator.')
        iterm=iqint(qlog10(cqabs(r1dc)))
        ir1de=ir1e+iterm
 	r1dc=r1dc*10.0q0**(-iterm)
        if(cqabs(r1dc).ge.1.0q0) go to 70
        r1dc=r1dc*10.0q0
        ir1de=ir1de-1
70	continue
        return
        end
!

!
!
        subroutine r2int (l,m,c,limint,ndec,maxd,enr,d01,id01,&
                    maxint,maxmp,nex,maxlp,rpint1,rpint2,&
                    pint1,pint2,pint3,pint4,norme,pnorm,ipnorm,&
                    coefme,coefmo,ir2est,r2c,ir2e,r2dc,ir2de,jint,&
                    coefn,icoefn)
!
!
!  purpose:     to calculate values of the radial function of the
!               second kind and its first derivative using an integral
!               representation of the radial functions in terms of the
!               angular function of the first kind together with a
!               Neumann function kernal. the angular function is
!               expanded in a series of associated Legendre functions.
!               gaussian quadrature is used (in subroutine pint) to
!               evaluate the resulting integrals involving associated
!               Legendre functions times the Neumann function kernel.
!               this subroutine r2int performs the summation of the
!               integrals times d constants to obtain r2 and r2d.
!
!  parameters:
!
!     input:    l      : l
!               m      : m
!               c      : c
!               limint : approximately twice the maximum number of
!                        terms available to be taken in the series
!               ndec   : number of decimal digits available in
!                        r*16 arithmetic
!               maxd   : dimension of enr array
!               enr    : d coefficient ratios
!               d01    : characteristic of the first d constant,
!                        either d0 or d1, depending on whether l-m
!                        is even or odd
!               id01   : exponent (base 10) of the first d constant
!               maxint : dimension of pint and rpint arrays
!               maxmp  : dimension of norme array
!               nex    : maximum exponent available in real*16
!                        arithmetic
!               maxlp  : maximum  l value desired; dimension
!                        of the pnorm and ipnorm arrays
!               rpint1 : arrays of ratios of successive integrals of
!                        either the first or the third kind, depending
!                        on whether l-m is even or odd
!               rpint2 : array of ratios of successive integrals of
!                        either the second or the fourth kind,
!                        depending on whether l-m is even or odd
!               pint1  : array of scaled values for the integrals of
!                        the first kind
!               pint2  : array of scaled values for the integrals of
!                        the second kind
!               pint3  : array of scaled values for the integrals of
!                        the third kind
!               pint4  : array of scaled values for the integrals of
!                        the fourth kind
!               norme  : array of exponents used to scale the Neumann
!                        function of order m involved in the integrals
!               pnorm  : array of characteristics of the scaling factors
!                        used for the associated Legendre functions in
!                        the integrals to avoid overflow
!               ipnorm : array of exponents (base 10) corresponding to
!                        pnorm
!               coefme : coefficient used to multiply the r2 sum to
!                        obtain its contribution to r2d when l-m is even
!               coefmo : coefficient used to multiply the r2 sum to
!                        obtain its contribution to r2d when l-m is odd
!               ir2est : estimate of the exponent (base 10) of the
!                        radial function of the second kind (estimated
!                        using calculted values for r1d and the
!                        theoretical wronskian
!
!     output:   r2c    : characteristic of prolate radial function
!                        of the first kind
!               ir2e   : exponent of prolate radial function of the
!                        first kind
!               r2dc   : characteristic of derivative with respect
!                        to x of prolate radial function of the first
!                        kind
!               ir2de  : exponent of derivative with respect to x of
!                        prolate radial function of the first kind
!               jint   : maximum value of the index j in the forward
!                        sum for r2 and r2d, i.e., the highest enr(j)
!                        used
!               coefn  : characteristic of a coefficient calculated
!                        for l = m and used for higher values of l
!               icoefn : exponent for coefn
!
!  real*16 scalars and arrays
        real*16 aj,arr,coefa,coefme,coefmo,coefn,dec,ri,&
          rm,rm2,r2dposi,r2dposr,r2posi,r2posr
        complex*32 c,coefl,dnew,dnewd,dold,doldd,d01,r2c,r2dc,r2dtemp,&
             r2temp
        complex*32 enr(maxd),pint1(maxint),pint2(maxint),&
          pint3(maxint),pint4(maxint),rpint1(maxint),&
          rpint2(maxint)
        real*16 pnorm(maxlp)
!
!  integer arrays
        dimension norme(maxmp),ipnorm(maxlp)
!
        rm=qfloat(m)
        rm2=rm+rm
        lm2=(l-m)/2
        ix=l-m-2*lm2
        ixx=ix-1
        ixx2=ixx+2
        lim=limint/2-ix
!
!  compute the leading coefficient
        if(l.gt.m) go to 20
        icoefn=norme(m+1)
        coefn=0.5q0
        if(m.eq.0) go to 20
          do 10 i=1,m
          ri=qfloat(i)
	  coefn=coefn/(ri+ri)
          iterm=iqint(qlog10(qabs(coefn)))
          coefn=coefn*10.0q0**(-iterm)
10        icoefn=icoefn+iterm
20      continue
        if(ix.eq.0) coefa=(rm2+1.q0)*coefn
        if(ix.eq.1) coefa=(rm2+3.q0)*coefn
        if((ix.eq.0).and.(2*(lm2/2).ne.lm2)) coefa=-coefa
        if((ix.eq.1).and.(2*((l-m-1)/4).ne.(l-m-1)/2)) coefa=-coefa
	coefl=coefa/d01
        icoefl=-id01+icoefn
        dec=10.q0**(-ndec-1)
        jlow=lm2+1
!
!  compute the integrals of s
!
!  forward summation for series for r2 and r2d
        dold=(1.0q0,0.0q0)
        doldd=(1.0q0,0.0q0)
        r2dtemp=doldd
        r2temp=dold
        r2posr=1.0q0
        r2posi=0.0q0
        r2dposr=1.0q0
        r2dposi=0.0q0
          do 30 j=jlow,lim
          dnew=dold*enr(j)*rpint1(j+j+ixx2)
          dnewd=doldd*enr(j)*rpint2(j+j+ixx2)
          r2temp=r2temp+dnew
          r2dtemp=r2dtemp+dnewd
          if(qreal(dnew).gt.0.q0) r2posr=r2posr+qreal(dnew)
          if(qreal(dnewd).gt.0.q0) r2dposr=r2dposr+qreal(dnewd)
          if(qimag(dnew).gt.0.q0) r2posi=r2posi+qimag(dnew)
          if(qimag(dnewd).gt.0.q0) r2dposi=r2dposi+qimag(dnewd)
          if((cqabs(dnew/r2temp)+cqabs(dnewd/r2dtemp)).lt.dec) go to 40
          dold=dnew
          doldd=dnewd
30        continue
!
!  backward summation for series for r2 and r2d
40      jint=j
        if(jint.gt.lim) jint=lim
        if(lm2.eq.0) go to 60
        dold=(1.0q0,0.0q0)
        doldd=(1.0q0,0.0q0)
        j=lm2
          do 50 jj=1,lm2
          dnew=dold/(rpint1(j+j+ixx2)*enr(j))
          dnewd=doldd/(rpint2(j+j+ixx2)*enr(j))
          r2temp=r2temp+dnew
          r2dtemp=r2dtemp+dnewd
          if(qreal(dnew).gt.0.q0) r2posr=r2posr+qreal(dnew)
          if(qreal(dnewd).gt.0.q0) r2dposr=r2dposr+qreal(dnewd)
          if(qimag(dnew).gt.0.q0) r2posi=r2posi+qimag(dnew)
          if(qimag(dnewd).gt.0.q0) r2dposi=r2dposi+qimag(dnewd)
          if((cqabs(dnew/r2temp)+cqabs(dnewd/r2dtemp)).lt.dec) go to 60
          dold=dnew
          doldd=dnewd
          j=j-1
50        continue
60      continue
        isubr=-iqint(qlog10(qabs(qreal(r2temp)/r2posr)+dec))
        isubi=0
        if(r2posi.ne.0.0q0) isubi=-iqint(qlog10(qabs(qimag(r2temp)&
                             /r2posi)+dec))
        isubdr=-iqint(qlog10(qabs(qreal(r2dtemp)/r2dposr)+dec))
        isubdi=0
        if(r2dposi.ne.0.0q0) isubdi=-iqint(qlog10(qabs(qimag(r2dtemp)&
                               /r2dposi)+dec))
	r2temp=r2temp*coefl*pnorm(l-m+1)
        if(ix.eq.0) r2temp=r2temp*pint1(l-m+1)
        if(ix.eq.1) r2temp=r2temp*pint3(l-m+1)
        iterm=iqint(qlog10(cqabs(r2temp)))
        ir2e=iterm+ipnorm(l-m+1)+icoefl
        r2c=r2temp*10.0q0**(-iterm)
        if(cqabs(r2c).ge.1.0q0) go to 70
        r2c=r2c*10.0q0
        ir2e=ir2e-1
70	r2dtemp=-r2dtemp*coefl*pnorm(l-m+1)*c
        if(ix.eq.0) r2dtemp=r2dtemp*pint2(l-m+1)+r2temp*coefme
        if(ix.eq.1) r2dtemp=r2dtemp*pint4(l-m+1)+r2temp*coefmo
        if(ix.eq.0) jsubdr=-iqint(qlog10(qabs(qreal(r2dtemp)/&
                    qreal(r2temp*coefme))+dec))
        jsubdi=0
        if(ix.eq.0.and.qimag(r2temp).ne.0.0q0) jsubdi=-&
                    iqint(qlog10(qabs(qimag(r2dtemp)/&
                    qimag(r2temp*coefme))+dec))
        if(ix.eq.1) jsubdr=-iqint(qlog10(qabs(qreal(r2dtemp)/&
                    qreal(r2temp*coefmo))+dec))
        if(ix.eq.1.and.qimag(r2temp).ne.0.0q0) jsubdi=-&
                    iqint(qlog10(qabs(qimag(r2dtemp)/&
                    qimag(r2temp*coefmo))+dec))
        if(jsubdr.lt.0) jsubdr=0
        if(jsubdi.lt.0) jsubdi=0
        isubdr=isubdr+jsubdr
        isubdi=isubdi+jsubdi
        isub=max(isubr,isubi)
        isubd=max(isubdr,isubdi)
        write(40,80) jint,lim,isub,isubd
80      format(8x,'r2int: converged in ',i6,' terms; 'i6,&
         ' available; ',i3,' digits of subtraction error',/,&
         15x,'in r2, ',i3,' digits of subtraction error in r2d.')
        jterm=iqint(qlog10(cqabs(r2dtemp)))
        ir2de=jterm+ipnorm(l-m+1)+icoefl
        r2dc=r2dtemp*10.0q0**(-jterm)
        if(cqabs(r2dc).ge.1.0q0) go to 90
        r2dc=r2dc*10.0q0
        ir2de=ir2de-1
90      continue
        return
        end
!
!
!
        subroutine r2leg (l,m,c,x1,lnum,limleg,limdr,iflagp,&
                    ndec,maxd,maxmp,maxpdr,maxdr,maxq,enr,&
                    enrneg,drhor,d01,id01,dneg,idneg,dfnorm,&
                    nsubf,dmfnorm,nsubmf,prx,pdrx,qdr,qdml,iqdml,&
                    qdl,iqdl,qr,qml,iqml,ql,iql,fajo,ifajo,termpq,&
                    itermpq,ioppsum,iopqnsum,sgn,r1dc,ir1de,r2c,&
                    ir2e,r2dc,ir2de,jleg,jlegp)
!
!  purpose:     to evaluate the prolate radial function of the
!               second kind and its first derivative with respect
!               to x using the traditional expansion in associated
!               Legendre functions.
!
!  parameters:
!
!     input :   l       : l
!               m       : m
!               c       : c
!               x1      : x-1
!               lnum    : number of l values desired
!               limleg  : approximately twice the maximum number
!                         of terms available to be taken in qsum,
!                         (sum involving q's time d constants)
!               limdr   : maximum number of terms available to be
!                         taken in psum (sum involving p's time
!                         d rho constants)
!               iflagp  : integer flag set = 1 if psum series converges
!                         fully; set = 0 otherwise
!               ndec    : number of decimal digits available in
!                         r*16 arithmetic
!               maxd    : dimension of enr array
!               maxmp   : dimension of enrneg array
!               maxpdr  : dimension of prx and pdrx arrays
!               maxdr   : dimension of drhor array
!               maxq    : dimension of qr and qdr arrays
!               enr     : array of d coefficient ratios
!               enrneg  : array of d coefficient ratios with
!                         negative subscripts
!               drhor   : array of d rho coefficient ratios
!               d01     : characteristic of the ratio of the first d
!                         constant with nonnegative subscript, either
!                         d0 or d1 depending on whether l-m is even or
!                         odd, to the d constant with subscript l - m
!               id01    : exponent (base 10) corresponding to d01
!               dneg    : characteristic of the ratio of the first d
!                         coefficient with negative subscript (either
!                         -1 or -2, depending on whether l - m is odd
!                         or even) to the d constant with subscript
!                         l - m
!               idneg   : exponent corresponding to dneg
!               dfnorm  : Flammer normalization sum divided by the d
!                         constant with subscript l - m
!               dmfnorm : Morse-Feshbach normalization factor of the
!                         d coefficients. equal to the reciprocal of the
!                         value of the d constant d(n = l - m) using
!                         this normalization for the angular functions
!               prx     : ratios of successive Legendre functions of
!                         the first kind of the same parity
!               pdrx    : ratios of successive first derivatives of
!                         Legendre functions of the first kind of the
!                         same parity
!               qdr     : ratios of first derivatives of successive
!                         Legendre functions of the second kind
!               qdml    : characteristic of the first derivative of
!                         the associated Legendre function of the second
!                         kind with order m and degree m-1
!               iqdml   : exponent corresponding to qdml
!               qdl     : characteristic of the first derivative of
!                         the associated Legendre function of the second
!                         kind with order m and degree m
!               iqdl    : exponent corresponding to qdl
!               qr      : array of ratios of successive associated
!                         Legendre functions of the second kind
!               qml     : characteristic of the associated Legendre
!                         function of the second kind with order m
!                         and degree m-1
!               iqml    : exponent corresponding to qml
!               ql      : characteristic of the associated Legendre
!                         function of the second kind with order m and
!                         degree m
!               iql     : exponent corresponding to ql
!               fajo    : characteristic of the joining factor of the
!                         second kind
!               ifajo   : exponent corresponding to fajo
!               termpq  : characteristic of the relative size of the
!                         maximum terms in the positive degree q series
!                         and the p series used to calculate r2 and r2d
!               itermpq : exponent corresponding to termpq
!               ioppsum : integer flag = 0 if psum need not be computed
!                         since its contribution to r2 and r2d is
!                         negligible; = 1 if psum is computed
!               iopqnsum: integer flag = 0 if qnsum need not be computed
!                         since its contribution to r2 and r2d is
!                         negligible; = 1 if qnsum is computed
!               sgn     : sign of the d coefficients with subscript
!                         l - m
!               r1dc    : charcteristic of corresponding first
!                         derivative of the radial function of the first
!                         kind
!               ir1de   : exponent of corresponding first derivative of
!                         the radial function of the first kind
!
!     output:   r2c     : characteristic of prolate
!                         radial function of the second kind
!               ir2e    : exponent of prolate radial function of the
!                         second kind
!               r2dc    : characteristic of derivative with
!                         respect to x of prolate radial function
!                         of the second kind
!               ir2de   : exponent of derivative with respect to x of
!                         prolate radial function of second kind
!               jleg    : maximum number of terms taken in qsum
!               jlegp   : maximum number of terms taken in psum
!
!  real*16 scalars and arrays
        real*16 dconp,dconq,dconqn,dec,qdml,qml,rm,sgn,sumr,sumi,&
          termpq,test,testd,testm,testdm,tm,x1
        complex*32 c,dfnorm,dmfnorm,dneg,dnegjf,dnew,dnewd,dold,doldd,&
          d01,psum,pdsum,pstest,qndsum,qdsum,qnsum,qsum,&
          r1dc,r2c,r2dc,spsum,spdsum
        real*16 prx(maxpdr),pdrx(maxpdr),qdl(lnum),qdr(maxq),ql(lnum),&
          qr(maxq)
        real*16 r2qposr,r2dqposr,r2qposi,r2dqposi,r2qnposr,r2dqnposr,&
          r2qnposi,r2dqnposi,r2pposr,r2dpposr,r2pposi,r2dpposi,&
          sr2pposr,sr2dpposr,sr2pposi,sr2dpposi
        complex*32 drhor(maxdr),enr(maxd),enrneg(maxmp),fajo(lnum)
!
!  integer arrays
        dimension ifajo(lnum),iqdl(lnum),iql(lnum)
!
        dec=10.q0**(-ndec-1)
        lm2=(l-m)/2
        ix=l-m-2*lm2
        imxp=m+m+ix
        ixx=1-ix
        lim1=limleg/2-ix
        lim2=limdr-1
        rm=qfloat(m)
        tm=rm+rm
        dconq=dec
        dconqn=dec
        dconp=dec
        dnegjf=dneg*d01
        if(m.eq.0) dnegjf=d01
        idnegjf=idneg+id01
        if(m.eq.0) idnegjf=id01
        fajo(l-m+1)=fajo(l-m+1)*dmfnorm*dnegjf/dfnorm
        iterm=iqint(qlog10(cqabs(fajo(l-m+1))))
        fajo(l-m+1)=fajo(l-m+1)*10.0q0**(-iterm)
        ifajo(l-m+1)=ifajo(l-m+1)+idnegjf+iterm
        nsubqr=0
        nsubqi=0
        nsubqdr=0
        nsubqdi=0
        nsubqnr=0
        nsubqni=0
        nsubqndr=0
        nsubqndi=0
        nsubpr=0
        nsubpi=0
        nsubpdr=0
        nsubpdi=0
!
!  begin calculation of series for r2
!
!  calculate d*q sum over positive n using pyramid summation
!
!  backward summation
        qsum=(1.q0,0.0q0)
        qdsum=(1.q0,0.0q0)
        r2qposr=1.0q0
        r2qposi=0.0q0
        r2dqposr=1.0q0
        r2dqposi=0.0q0
        if(lm2.lt.1) go to 20
        dold=(1.0q0,0.0q0)
        doldd=(1.0q0,0.0q0)
        j=lm2
          do 10 jj=1,lm2
          dnew=dold/(qr(j+j+imxp)*qr(j+j+imxp-1)*enr(j))
          qsum=qsum+dnew
          dnewd=doldd/(qdr(j+j+imxp)*qdr(j+j+imxp-1)*enr(j))
          qdsum=qdsum+dnewd
          if(qreal(dnew).gt.0.q0) r2qposr=r2qposr+qreal(dnew)
          if(qreal(dnewd).gt.0.q0) r2dqposr=r2dqposr+qreal(dnewd)
          if(qimag(dnew).gt.0.q0) r2qposi=r2qposi+qimag(dnew)
          if(qimag(dnewd).gt.0.q0) r2dqposi=r2dqposi+qimag(dnewd)
          if((cqabs(dnew/qsum)+cqabs(dnewd/qdsum)).lt.dconq) go to 20
          dold=dnew
          doldd=dnewd
          j=j-1
10        continue
20      continue
!
!  forward summation
        jlow=lm2+1
        dold=(1.0q0,0.0q0)
        doldd=(1.0q0,0.0q0)
          do 30 j=jlow,lim1
          dnew=dold*enr(j)*qr(j+j+imxp)*qr(j+j+imxp-1)
          qsum=qsum+dnew
          dnewd=doldd*enr(j)*qdr(j+j+imxp)*qdr(j+j+imxp-1)
          qdsum=qdsum+dnewd
          if(qreal(dnew).gt.0.q0) r2qposr=r2qposr+qreal(dnew)
          if(qreal(dnewd).gt.0.q0) r2dqposr=r2dqposr+qreal(dnewd)
          if(qimag(dnew).gt.0.q0) r2qposi=r2qposi+qimag(dnew)
          if(qimag(dnewd).gt.0.q0) r2dqposi=r2dqposi+qimag(dnewd)
          if((cqabs(dnew/qsum)+cqabs(dnewd/qdsum)).lt.dconq) go to 40
          dold=dnew
          doldd=dnewd
30        continue
40      continue
        jleg=j
        if(jleg.gt.lim1) jleg=lim1
        itestqsum=-iqint(qlog10(cqabs(dnew/qsum)+dec))
!
        if(qreal(qsum).ne.0.0q0) subqr=iqint(qlog10(qabs(r2qposr/qreal(qsum)+dec)))
        if(nsubqr.lt.0) nsubqr=0
        if(qimag(qsum).ne.0.0q0)&
subqi=iqint(qlog10(qabs(r2qposi/qimag(qsum)+dec)))
        if(nsubqi.lt.0) nsubqi=0
        nsubq=max(nsubqi,nsubqr)
        if(qreal(qdsum).ne.0.0q0)&
subqdr=iqint(qlog10(qabs(r2dqposr/qreal(qdsum)+dec)))
        if(nsubqdr.lt.0) nsubqdr=0
        if(qimag(qdsum).ne.0.0q0)&
subqdi=iqint(qlog10(qabs(r2dqposi/qimag(qdsum)+dec)))
        if(nsubqdi.lt.0) nsubqdi=0
        nsubqd=max(nsubqdi,nsubqdr)
!
        qsum=qsum*ql(l-m+1)/(fajo(l-m+1)*termpq)
        iterm=iqint(qlog10(cqabs(qsum)))
        qsum=qsum*(10.0q0**(-iterm))
        iqsum=iql(l-m+1)-ifajo(l-m+1)-itermpq+iterm
        qdsumsav=qdsum
        qdsum=qdsum*qdl(l-m+1)/(fajo(l-m+1)*termpq)
        iterm=iqint(qlog10(cqabs(qdsum)))
        qdsum=qdsum*(10.0q0**(-iterm))
        iqdsum=iqdl(l-m+1)-ifajo(l-m+1)-itermpq+iterm
!
!  calculate d*q sum over negative n
        qnsum=(0.0q0,0.0q0)
        qndsum=(0.0q0,0.0q0)
        r2qnposr=0.0q0
        r2dqnposr=0.0q0
        r2qnposi=0.0q0
        r2dqnposi=0.0q0
        iqnsum=0
        iqndsum=0
        nsubqn=0
        nsubqnd=0
        j2=0
        nmterm=0
        if(iopqnsum.eq.0.or.m.eq.0) go to 90
        nmterm=m
        qnsum=enrneg(m)
        qndsum=enrneg(m)
        j2=1
        if(ix.eq.1) go to 50
        qnsum=qnsum*qr(m+m-1)
        qndsum=qndsum*qdr(m+m-1)
50      continue
        if(qreal(qnsum).gt.0.0q0) r2qnposr=qreal(qnsum)
        if(qreal(qndsum).gt.0.0q0) r2dqnposr=qreal(qndsum)
        if(qimag(qnsum).gt.0.0q0) r2qnposi=qimag(qnsum)
        if(qimag(qndsum).gt.0.0q0) r2dqnposi=qimag(qndsum)
        if(m.eq.1) go to 80
        dold=qnsum
        doldd=qndsum
          do 60 j=2,m
          dnew=dold*enrneg(m-j+1)*qr(imxp-j-j+1)*qr(imxp-j-j+2)
          qnsum=qnsum+dnew
          dnewd=doldd*enrneg(m-j+1)*qdr(imxp-j-j+1)*qdr(imxp-j-j+2)
          qndsum=qndsum+dnewd
          if(qreal(dnew).gt.0.q0) r2qnposr=r2qnposr+qreal(dnew)
          if(qreal(dnewd).gt.0.q0) r2dqnposr=r2dqnposr+qreal(dnewd)
          if(qimag(dnew).gt.0.q0) r2qnposi=r2qnposi+qimag(dnew)
          if(qimag(dnewd).gt.0.q0) r2dqnposi=r2dqnposi+qimag(dnewd)
          dold=dnew
60        doldd=dnewd
70      j2=j
        if(j2.gt.m) j2=m
80      continue
        if(qreal(qnsum).ne.0.0q0)&
subqnr=iqint(qlog10(qabs(r2qnposr/qreal(qnsum)+dec)))
        if(nsubqnr.lt.0) nsubqnr=0
        if(qimag(qnsum).ne.0.0q0)&
subqni=iqint(qlog10(qabs(r2qnposi/qimag(qnsum)+dec)))
        if(nsubqni.lt.0) nsubqni=0
        nsubqn=max(nsubqni,nsubqnr)
        if(qreal(qndsum).ne.0.0q0)&
subqndr=iqint(qlog10(qabs(r2dqnposr/qreal(qndsum)+dec)))
        if(nsubqndr.lt.0) nsubqndr=0
        if(qimag(qndsum).ne.0.0q0)&
subqndi=iqint(qlog10(qabs(r2dqnposi/qimag(qndsum)+dec)))
        if(nsubqndi.lt.0) nsubqndi=0
        nsubqnd=max(nsubqndi,nsubqndr)
!
        qnsum=qnsum*qml*d01/(fajo(l-m+1)*termpq)
        iterm=iqint(qlog10(cqabs(qnsum)))
        qnsum=qnsum*(10.0q0**(-iterm))
        iqnsum=iqml+id01-ifajo(l-m+1)-itermpq+iterm
        qnsum=qnsum*(10.0q0**(iqnsum-iqsum))
        qndsumsav=qndsum
        qndsum=qndsum*qdml*d01/(fajo(l-m+1)*termpq)
        iterm=iqint(qlog10(cqabs(qndsum)))
        qndsum=qndsum*(10.0q0**(-iterm))
        iqndsum=iqdml+id01-ifajo(l-m+1)-itermpq+iterm
        qndsum=qndsum*(10.0q0**(iqndsum-iqdsum))
90      continue
!
!       calculate d(rho|n)*p summation
        psum=(0.0q0,0.0q0)
        pdsum=(0.0q0,0.0q0)
        ipsum=0
        ipdsum=0
        r2pposr=0.0q0
        r2dpposr=0.0q0
        r2pposi=0.0q0
        r2dpposi=0.0q0
        nsubp=0
        nsubpd=0
        jlegp=0
        itestpsum=0
        if(ioppsum.eq.0) go to 160
        psum=prx(ixx+1)*drhor(1)
        pdsum=pdrx(ixx+1)*drhor(1)
        dold=psum
        doldd=pdsum
        if(m.ne.0.or.ix.ne.1) go to 100
        pdsum=(0.0q0,0.0q0)
        doldd=drhor(1)
100     continue
        spsum=psum
        spdsum=pdsum
        testm=1.0q0
        testdm=1.0q0
        iflagp=0
        jlegpf=1
        jlegpd=1
        if(qreal(psum).gt.0.0q0) r2pposr=qreal(psum)
        if(qreal(pdsum).gt.0.0q0) r2dpposr=qreal(pdsum)
        if(qimag(psum).gt.0.0q0) r2pposi=qimag(psum)
        if(qimag(pdsum).gt.0.0q0) r2dpposi=qimag(pdsum)
        sr2pposr=r2pposr
        sr2pposi=r2pposi
        sr2dpposr=r2dpposr
        sr2dpposi=r2dpposi
          do 130 j=2,lim2
          dnew=dold*drhor(j)*prx(j+j-ix)
          psum=psum+dnew
          dnewd=doldd*drhor(j)*pdrx(j+j-ix)
          pdsum=pdsum+dnewd
          test=cqabs(dnew/psum)
          testd=cqabs(dnewd/pdsum)
          if(qreal(dnew).gt.0.q0) r2pposr=r2pposr+qreal(dnew)
          if(qreal(dnewd).gt.0.q0) r2dpposr=r2dpposr+qreal(dnewd)
          if(qimag(dnew).gt.0.q0) r2pposi=r2pposi+qimag(dnew)
          if(qimag(dnewd).gt.0.q0) r2dpposi=r2dpposi+qimag(dnewd)
          if(test.gt.testm.or.test.eq.0.q0) go to 110
          testm=test
          spsum=psum
          sr2pposr=r2pposr
          sr2pposi=r2pposi
          jlegpf=j
110       if(testd.gt.testdm.or.testd.eq.0.q0) go to 120
          testdm=testd
          spdsum=pdsum
          sr2dpposr=r2dpposr
          sr2dpposi=r2dpposi
          jlegpd=j
120       continue
          if((test+testd).lt.dconp) go to 140
          dold=dnew
          doldd=dnewd
130       continue
        go to 150
140     continue
        iflagp=1
150     psum=spsum
        pdsum=spdsum
        r2pposr=sr2pposr
        r2pposi=sr2pposi
        r2dpposr=sr2dpposr
        r2dpposi=sr2dpposi
        jlegp=max(jlegpf,jlegpd)
        test=max(testm,testdm)
        itestpsum=-iqint(qlog10(test+dec))
        if(itestpsum.gt.ndec) itestpsum=ndec
        if(itestpsum.lt.0) itestpsum=0
        if(qreal(psum).ne.0.0q0)&
subpr=iqint(qlog10(qabs(r2pposr/qreal(psum)+dec)))
        if(nsubpr.lt.0) nsubpr=0
        if(qimag(psum).ne.0.0q0)&
subpi=iqint(qlog10(qabs(r2pposi/qimag(psum)+dec)))
        if(nsubpi.lt.0) nsubpi=0
        nsubp=max(nsubpi,nsubpr)
        if(qreal(pdsum).ne.0.0q0)&
subpdr=iqint(qlog10(qabs(r2dpposr/qreal(pdsum)+dec)))
        if(nsubpdr.lt.0) nsubpdr=0
        if(qimag(pdsum).ne.0.0q0)&
subpdi=iqint(qlog10(qabs(r2dpposi/qimag(pdsum)+dec)))
        if(nsubpdi.lt.0) nsubpdi=0
        nsubpd=max(nsubpdi,nsubpdr)
!
        psum=psum*dnegjf*termpq/fajo(l-m+1)
        iterm=0
        if(psum.ne.(0.0q0,0.0q0)) iterm=iqint(qlog10(cqabs(psum)))
        psum=psum*(10.0q0**(-iterm))
        ipsum=idnegjf+itermpq-ifajo(l-m+1)+iterm
        psum=psum*(10.0q0**(ipsum-iqsum))
        pdsum=pdsum*dnegjf*termpq/fajo(l-m+1)
        if(m.ne.0) then
        pdsum=pdsum*rm*(x1+1.q0)/(x1*(x1+2.q0))
        end if
        iterm=0
        if(pdsum.ne.0.q0) iterm=iqint(qlog10(cqabs(pdsum)))
        pdsum=pdsum*(10.0q0**(-iterm))
        ipdsum=idnegjf+itermpq-ifajo(l-m+1)+iterm
        pdsum=pdsum*(10.q0**(ipdsum-iqdsum))
!
160     r2c=qsum+qnsum+psum
        nsubpcor=0
        if(ioppsum.eq.1) nsubpcor=-iqint(qlog10(cqabs((psum/r2c))))
        nsubp=nsubp-nsubpcor
        if(nsubp.lt.0) nsubp=0
        nsubqcor=-iqint(qlog10(cqabs((qsum/r2c))))
        nsubq=nsubq-nsubqcor
        if(nsubq.lt.0) nsubq=0
        nsubqncor=0
        if(iopqnsum.eq.1) nsubqncor=-iqint(qlog10(cqabs((qnsum/r2c))))
        nsubqn=nsubqn-nsubqncor
        if(nsubqn.lt.0) nsubqn=0
        nsub=max(nsubq,nsubp,nsubqn)
        nsub=max(nsub,nsubf,nsubmf)
        if(nsub.gt.ndec) nsub=ndec
!
        nqns=0
        if(qnsum.ne.(0.0q0,0.0q0))&
      nqns=iqint(qlog10(cqabs(qnsum/r2c)))
        if(nqns.lt.(-ndec-1)) iopqnsum=0
        iterm=iqint(qlog10(cqabs(r2c)))
        r2c=r2c*(10.0q0**(-iterm))
        ir2e=iqsum+iterm
        if(cqabs(r2c).ge.1.q0) go to 170
        r2c=r2c*10.0q0
        ir2e=ir2e-1
170     continue
!
        r2dc=qdsum+qndsum+pdsum
        nsubpdcor=0
        if(ioppsum.eq.1) nsubpdcor=-iqint(qlog10(cqabs((pdsum/r2dc))))
        nsubpd=nsubpd-nsubpdcor
        if(nsubpd.lt.0) nsubpd=0
        nsubqdcor=-iqint(qlog10(cqabs((qdsum/r2dc))))
        nsubqd=nsubqd-nsubqdcor
        if(nsubqd.lt.0) nsubqd=0
        nsubqndcor=0
        if(iopqnsum.eq.1) nsubqndcor=&
                             -iqint(qlog10(cqabs((qndsum/r2dc))))
        nsubqnd=nsubqnd-nsubqndcor
        if(nsubqnd.lt.0) nsubqnd=0
        nsubd=max(nsubqd,nsubpd,nsubqnd)
        nsubd=max(nsubd,nsubf,nsubmf)
        if(nsubd.gt.ndec) nsubd=ndec
!
        iterm=iqint(qlog10(cqabs(r2dc)))
        r2dc=r2dc*(10.0q0**(-iterm))
        ir2de=iqdsum+iterm
        if(cqabs(r2dc).ge.1.0q0) go to 180
        r2dc=r2dc*10.0q0
        ir2de=ir2de-1
180     continue
        if(ioppsum.eq.1) write(40,190) itestqsum,jleg,itestpsum,&
                                 jlegp,lim1,lim2,nsub,nsubd
190     format(8x,'r2leg: qsum converged to',i3,' digits in ',i6,&
         ' terms; psum to',i3,' digits in',i6,' terms.',&
         /,15x,i6,' and ',i6,' terms available.',i3,&
         ' digits of sub. error in r2, ',i3,' digits in r2d.')
        if(ioppsum.eq.0) write(40,200) itestqsum,jleg,lim1,nsub,nsubd
200     format(8x,'r2leg: qsum converged to',i3,' digits in ',i6,&
         ' terms; ',i6,' terms available.',/,14x,i3,' digits of',&
         ' sub. error in r2;',i3,' digits of sub. error in r2d.'&
          ,/,15x,'psum is negligible and not calculated.')
        nps=0
        if(psum.ne.(0.0q0,0.0q0))&
                nps=iqint(qlog10(cqabs(psum/r2c)))
        if(nps.lt.(-ndec-1)) ioppsum=0
        return
        end
!
!
!
        subroutine r2neu (l,m,c,x1,limneu,ndec,maxd,maxlp,maxn,&
                    maxp,minacc,enr,sneuf,sneun,ineue,sneudf,&
                    sneudr,prat1,pcoefn,ipcoefn,dmfnorm,&
                    r1dc,ir1de,r2c,ir2e,r2dc,ir2de,jneu,naccsub)
!
!  purpose:     to calculate the prolate radial function of the
!               second kind and its first derivative with respect
!               to x, using the traditional expansions in terms of
!               spherical Neumann functions.
!
!  parameters:
!
!     input:    l      : l
!               m      : m
!               c      : complex c
!               x1     : x-1
!               limneu : maximum number of terms to be taken in the
!                        series summations for r2 and r2d
!               ndec   : number of decimal digits available in
!                        r*16 arithmetic
!               maxd   : dimension of enr array
!               maxlp  : maximum  l value desired; dimension
!                        of the sneun, sneudn, ineue, and ineude arrays
!               maxn   : dimension of sneuf and sneudf arrays
!               maxp   : dimension of prat1 array
!               minacc : number of decimal digits of desired accuracy
!                        of the resulting radial functions
!               enr    : array of ratios of successive d coefficients
!               sneuf  : array of ratios of successive spherical Neumann
!                        functions of the same parity
!               sneun  : array of characteristics for Neumann functions
!               ineue  : array of exponents for Neumann functions
!               sneudf : array of ratios of successive first derivatives
!                        of spherical Neumann functions of same parity
!               sneudr : array of ratios of first derivatives of Neumann
!                        functions to the corresponding functions
!               prat1  : array of ratios of successive coefficients in
!                        r2 and r2d sum
!               pcoefn : characteristic of coefficient for term in both
!                        r2 and r2d sums that contains Neumann function
!                        of order l + m
!               ipcoefn: exponent (to the base 10) corresponding to
!                        pcoefn
!               dmfnorm: Morse-Feshbach normalization factor of the
!                        d coefficients. equal to the reciprocal of the
!                        value of the d constant d(n = l - m) using
!                        this normalization for the angular functions
!               r1dc   : charcteristic of corresponding first
!                        derivative of the radial function of the first
!                        kind
!               ir1de  : exponent of corresponding first derivative of
!                        the radial function of the first kind
!
!     output:   r2c    : characteristic of prolate radial function
!                        of the second kind
!               ir2e   : exponent of prolate radial function of the
!                        second kind
!               r2dc   : characteristic of derivative with respect
!                        to x of prolate radial function of the second
!                        kind
!               ir2de  : exponent of derivative with respect to x of
!                        prolate radial function of the second kind
!               jneu   : index of term where best convergence is
!                        achieved for r2 or for r2d, whichever term is
!                        larger
!               naccsub: maximum subtraction error in forming Wronskian
!                        for previous value of l
!
!
!  real*16 and complex*32 scalars and arrays
        real*16 dconb,dconf,dconi,dec,pcoefn,rm,rm2,r2dcoef,r2est,&
          r2test,sumdpi,sumdpr,sumpi,sumpr,sumdpit,sumdprt,&
          sumpit,sumprt,test,testd,testdm,testm,x1
        complex*32 c,dmfnorm,dnew,dnewd,dold,doldd,r1dc,r2c,r2dc,&
             r2dtemp,r2temp,sr2temp,sr2dtemp,sumcoef
        real*16 prat1(maxp)
        complex*32 enr(maxd),sneudr(maxlp),sneun(maxn),&
             sneuf(maxn),sneudf(maxn)
!
!  integer arrays
        dimension ineue(maxn)
!
        rm=qfloat(m)
        dconf=10.0q0**(-minacc-4)
        dconi=10.0q0**(ndec+2)
        dec=10.0q0**(-ndec)
        sumcoef=(10.0q0**(-ir1de-ineue(l+1)-ipcoefn+naccsub))/&
          (c*x1*(x1+2.0q0)*r1dc*sneun(l+1)*pcoefn)
        r2est=cqabs(sumcoef*dmfnorm)
        dconb=r2est/dconi
        r2test=r2est*dconi
        r2dcoef=rm/((x1+1.q0)*(x1+2.q0)*x1)
        rm2=rm*2.q0
        lm2=(l-m)/2
!
!  ix = 0 for l-m even; ix = 1 for l-m odd
        ix=l-m-2*lm2
        lim=limneu/2-ix
!
!  compute radial function of the second kind
!
!  backward series
        r2temp=(1.0q0,0.0q0)
        sumpr=1.0q0
        sumpi=0.0q0
        r2dtemp=(1.0q0,0.0q0)
        sumdpr=1.0q0
        sumdpi=0.0q0
        if (lm2.lt.1) go to 20
        dold=(1.0q0,0.0q0)
        doldd=(1.0q0,0.0q0)
          do 10 j=lm2,1,-1
          jj=j+j+ix
          dnew=-dold/(sneuf(jj+m)*prat1(jj+1)*enr(j))
          dnewd=-doldd/(sneudf(jj+m)*prat1(jj+1)*enr(j))
          r2temp=r2temp+dnew
          if(qreal(dnew).gt.0.0q0) sumpr=sumpr+qreal(dnew)
          if(qimag(dnew).gt.0.0q0) sumpi=sumpi+qimag(dnew)
          r2dtemp=r2dtemp+dnewd
          if(qreal(dnewd).gt.0.0q0) sumdpr=sumdpr+qreal(dnewd)
          if(qimag(dnewd).gt.0.0q0) sumdpi=sumdpi+qimag(dnewd)
          if(cqabs(dnew/r2temp)+cqabs(dnewd/r2dtemp).lt.dconb) go to 20
          dold=dnew
          doldd=dnewd
10      continue
20      continue
!
!  forward series
        dold=(1.0q0,0.0q0)
        doldd=(1.0q0,0.0q0)
        testm=1.0q0
        testdm=1.0q0
        sr2temp=r2temp
        sumprt=sumpr
        sumpit=sumpi
        sr2dtemp=r2dtemp
        sumdprt=sumdpr
        sumdpit=sumdpi
        js=lim
        jds=lim
          do 70 j=lm2+1,lim
          jj=j+j+ix
          dnew=-dold*enr(j)*sneuf(jj+m)*prat1(jj+1)
          dnewd=-doldd*enr(j)*sneudf(jj+m)*prat1(jj+1)
          r2temp=r2temp+dnew
          if(qreal(dnew).gt.0.0q0) sumpr=sumpr+qreal(dnew)
          if(qimag(dnew).gt.0.0q0) sumpi=sumpi+qimag(dnew)
          r2dtemp=r2dtemp+dnewd
          if(qreal(dnewd).gt.0.0q0) sumdpr=sumdpr+qreal(dnewd)
          if(qimag(dnewd).gt.0.0q0) sumdpi=sumdpi+qimag(dnewd)
          test=cqabs(dnew/r2temp)
          testd=cqabs(dnewd/r2dtemp)
!         if(cqabs(r2temp).gt.r2test) go to 80
          if(test.lt.testm) go to 30
          go to 40
30        testm=test
          sr2temp=r2temp
          sumprt=sumpr
          sumpit=sumpi
          js=j
40        continue
          if(testd.lt.testdm) go to 50
          go to 60
50        testdm=testd
          sr2dtemp=r2dtemp
          sumdprt=sumdpr
          sumdpit=sumdpi
          jds=j
60        continue
          if ((test+testd).lt.dconf) go to 80
          dold=dnew
          doldd=dnewd
70        continue
80      r2temp=sr2temp
        r2dtemp=sr2dtemp
        sumpr=sumprt
        sumpi=sumpit
        sumdpr=sumdprt
        sumdpi=sumdpit
90      continue
        jneu=max(js,jds)
        if(j.ge.lim) jneu=jneu+qint(cqabs(c))
        jtestm=-iqint(qlog10(testm+dec))
        if(jtestm.lt.0) jtestm=0
        if(jtestm.gt.ndec) jtestm=ndec
        jtestdm=-iqint(qlog10(testdm+dec))
        if(jtestdm.lt.0) jtestdm=0
        if(jtestdm.gt.ndec) jtestdm=ndec
        naccns1r=-iqint(qlog10(qabs(qreal(r2temp)/sumpr)+dec))
        naccns1i=0
        if(sumpi.ne.0.0q0) naccns1i=-iqint(qlog10(qabs(qimag(r2temp)&
                            /sumpi)+dec))
        if(naccns1r.lt.0) naccns1r=0
        if(naccns1r.gt.ndec) naccns1r=ndec
        if(naccns1i.lt.0) naccns1i=0
        if(naccns1i.gt.ndec) naccns1i=ndec
        naccns2r=-iqint(qlog10(qabs(qreal(r2dtemp)/sumdpr)+dec))
        naccns2i=0
        if(sumdpi.ne.0.0q0) naccns2i=-iqint(qlog10(qabs(qimag(r2dtemp)&
                            /sumdpi)+dec))
        if(naccns2r.lt.0) naccns2r=0
        if(naccns2r.gt.ndec) naccns2r=ndec
        if(naccns2i.lt.0) naccns2i=0
        if(naccns2i.gt.ndec) naccns2i=ndec
        naccns1=max(naccns1r,naccns1i)
        naccns2=max(naccns2r,naccns2i)
        naccns=max(naccns1,naccns2)
        naccn1=min(jtestm,ndec-naccns1)
        naccn2=min(jtestdm,ndec-naccns2)
        naccn=min(naccn1,naccn2)
        if(naccn.lt.0) naccn=0
        write(40,100) j,lim,js,jtestm,naccns1,jds,jtestdm,naccns2,&
                naccn1,naccn2
100     format(8x,'r2neu: numerator converged in ',i6,' terms; ',i6,&
         ' terms available. best r2 at ',i6,' terms with',&
         /15x,'convergence to',i3,' digits;',i3,' digits',&
         ' subtr. error. best r2d at ',i6,' terms with',/,15x,&
         'convergence to ',i2,' digits;',i3,' digits subtr.',&
         ' error. estimated numerator accuracy is ',i2,' digits.'&
         ' for r2 and',i2,' digits fot r2d')
!
!  combining results to form the radial function characteristics
!  r2c and r2dc and corresponding exponents ir2e and ir2de
        r2c=r2temp*sneun(l+1)*pcoefn/dmfnorm
        iterm=iqint(qlog10(cqabs(r2c)))
        ir2e=ineue(l+1)+ipcoefn+iterm
        r2c=r2c*10.0q0**(-iterm)
        if(cqabs(r2c).ge.1.0q0) go to 110
        r2c=r2c*10.0q0
        ir2e=ir2e-1
110	continue
        r2dc=r2dcoef*r2c+(c*r2dtemp*sneun(l+1)*sneudr(l+1)*pcoefn/&
       dmfnorm)*10.0q0**(ineue(l+1)+ipcoefn-ir2e)
        iterm=iqint(qlog10(cqabs(r2dc)))
        ir2de=ir2e+iterm
 	r2dc=r2dc*10.0q0**(-iterm)
        if(cqabs(r2dc).ge.1.0q0) go to 120
        r2dc=r2dc*10.0q0
        ir2de=ir2de-1
120	continue
        return
        end
!
!
!
        subroutine r2eta (l,m,c,x1,eta,nee,limeta,ndec,maxd,maxlp,&
                    maxn,maxp,minacc,lowtest,enr,sneuf,sneun,&
                    ineue,sneudf,sneudr,pdratt,pratb,pratt,&
                    pcoefn,ipcoefn,pdcoefn,ipdcoefn,nsubf,nsubmf,&
                    idigc,ichoicer1,naccr1,r1c,ir1e,r1dc,ir1de,&
                    naccsub,naccr,r2c,ir2e,r2dc,ir2de,nacceta,&
                    nacciop,jeta,iopnee,neemark,naccnmax,naccsube)
!
!  purpose:     to calculate the prolate radial function of the
!               second kind and its first derivative with respect
!               to x, using an expansion of spherical Neumann
!               functions.
!
!  parameters:
!
!     input:   l        : l
!              m        : m
!              c        : complex c
!              x1       : x-1
!              eta      : value for eta used in calculation
!              nee      : index in the array of eta values in the main
!                         program that corresponds to the value of eta
!                         used in r2eta calculations
!              limeta   : maximum number of terms available in the sums
!                         for r2 and r2d
!              ndec     : number of decimal digits available in
!                         r*16 arithmetic
!              maxd     : dimension of enr array
!              maxlp    : maximum  l value desired; dimension
!                         of the sneudr arrays
!              maxn     : dimension of sneuf, sneudf, sneun and ineue
!                         arrays
!              maxp     : dimension of pdratt, pratb, and pratt arrays
!              minacc   : minimum number of accurate decimal digits
!                         that are requested
!              lowtest  : minimum accuracy of calculated radial
!                         functions of the second kind for previous
!                         values of l
!              enr      : array of ratios of successive d coefficients
!              sneuf    : array of ratios of successive spherical
!                         Neumann functions of the same parity
!              sneun    : array of characteristics for Neumann functions
!              ineue    : array of exponents corresponding to sneun
!              sneudf   : array of ratios of successive first
!                         derivatives of spherical Neumann functions of
!                         the same parity
!              sneudr   : array of ratios of first derivatives of the
!                         spherical Neumann functions to the
!                         corresponding functions
!              pdratt   : array of ratios of successive first
!                         derivatives of the associated Legendre
!                         functions of the first kind of the same parity
!                         (used in numerator series)
!              pratb    : array of ratios of successive associated
!                         Legendre functions of the first kind of the
!                         same parity (used in denominator series)
!              pratt    : array of ratios of successive associated
!                         Legendre functions of the first kind of the
!                         same parity (used in numerator series)
!              pcoefn   : characteristic of the ratio of the numerator
!                         and denominator associated Legendre functions
!                         of the first kind of order m and degree l
!              ipcoefn  : exponent corresponding to pcoefn
!              pdcoefn  : characteristic of the ratio of the first
!                         derivative of the associated Legendre function
!                         of the first kind in the numerator and the
!                         associated Legendre function of the first kind
!                         in the denominator, both of order m and
!                         degree l
!              ipdcoefn : exponent corresponding to pdcoefn
!              nsubf    : number of decimal digits of subtraction error
!                         in computing the Flammer normalization factor
!              nsubmf   : number of decimal digits of subtraction error
!                         in computing the Morse-Feshbach normalization
!                         factor
!              idigc    : number of decimal digits of convergence of the
!                         eigenvalue using the Bouwkamp procedure
!              ichoicer1: integer equal to 1 if the radial functions
!                         of the first kind are computed using the
!                         expansion with eta set equal to zero; equal
!                         to 2 if the traditional expansion with eta
!                         set equal to unity is used
!              naccr1   : estimate of number of decimal digits of
!                         accuracy of the radial function of the first
!                         kind (and its first derivative) calculated
!                         using the subtraction errors in their
!                         calculation
!              r1c      : first kind radial function (calculated in
!                         either r1besa or r1besb)
!              irie     : exponent corresponding to r1c
!              r1dc     : characteristic of the first derivative with
!                         respect to x of the radial function of the
!                         first kind (calculated in r1besa or r1besb)
!              ir1de    : exponent corresponding to r1dc
!              naccsub  : maximum subtraction error in computing the
!                         Wronskian using the calculated radial function
!                         values from all of the methods used for the
!                         previous value of l
!              naccr    : maximum accuracy obtained for the radial
!                         functions of the second kind for the current
!                         value of l using other methods as well as
!                         previous tries using r2eta
!
!     output:  r2c      : characteristic of the prolate radial function
!                         of the second kind
!              ir2e     : exponent of the prolate radial function of the
!                         second kind
!              r2dc     : characteristic of the first derivative with
!                         respect to x of the prolate radial function
!                         of the second kind
!              ir2de    : exponent corresponding to r2dc
!              nacceta  : estimated number of accurate decimal digits in
!                         r2 and r2d. computed from the wronskian if
!                         nacciop = 0; estimated from the series if
!                         nacciop = 1
!              nacciop  : integer flag = 1 if the denominator in the
!                         expressions for r2 and r2d is computed using
!                         the theoretical wronskian and known values for
!                         r1 and r1d. nacciop = 0 otherwise
!              jeta     : maximum number of terms taken in the numerator
!                         sums for r2 and r2d
!              iopnee   : integer flag = 0 if none of the values used
!                         previously for eta in the variable eta method
!                         has led to estimated accuracy in the numerator
!                         of at least 3 digits (where iopnee is set to
!                         1) or has led to a subtraction error larger
!                         than either ndec-lowtest digits or ndec-naccr
!                         digits in either the r2 or the r2d numerator
!                         series (where iopnee is set to 2)
!              neemark  : index for the eta value in the array storing
!                         eta values in the main program that
!                         corresponds to the last eta value used for
!                         which iopnee = 0
!              naccnmax : maximum numerator accuracy (in decimal digits)
!                         obtained for the current value of l
!              naccsube : corresponding Wronskian subtraction error for
!                         the variable eta method
!
!
!  real*16 and complex*32 scalars and arrays
        real*16 ca,dcon,dconi,dec,eta,etas,pcoefn,pdcoefn,&
          rm,rm2,r2dcoef1,r2est,r2test,sumdpi,sumdnpr,sumdnpi,&
          sumdpr,sumnpi,sumnpr,test,testd,testdm,testm,testw,&
          testw1,testw2,txi,txr,txdi,txdr,xet,xets,x1
        complex*32 c,denom,dnew,dnewd,dnewd1,dnewd2,dnewsum,dnewdsum,&
             dold,doldd,doldd1,doldd2,r1c,r1dc,r2c,r2dc,r2dcoef2,&
             r2dtemp,r2temp,reld12,sr2dtemp,sr2temp,sumcoef,&
             wronc,wront
        real*16 pratb(maxp),pratt(maxp),pdratt(maxp)
        complex*32 enr(maxd),sneudr(maxlp),sneun(maxn),sneuf(maxn),&
             sneudf(maxn)
!
!  integer arrays
        dimension ineue(maxn)
!
        ca=cqabs(c)
        dec=10.0q0**(-ndec-2)
        rm=qfloat(m)
        dcon=10.0q0**(-minacc-4)
        dconi=10.0q0**(ndec+2)
        etas=eta*eta
        xet=qsqrt(x1*(x1+2.0q0)+etas)
        xets=xet*xet
        sumcoef=(10.0q0**(-ir1de-ineue(l+1)-ipcoefn+naccsub))/&
          (c*x1*(x1+2.0q0)*r1dc*sneun(l+1)*pcoefn)
        r2dcoef1=-eta*(1.q0-etas)/(xets*xet)
        r2dcoef2=c*(x1+1.q0)/xet
        reld12=(r2dcoef2/r2dcoef1)*sneudr(l+1)*(pcoefn/&
         pdcoefn)*10.0q0**(ipcoefn-ipdcoefn)
        rm2=rm*2.0q0
        lm2=(l-m)/2
!
!  ix = 0 for l-m even; ix = 1 for l-m odd
        ix=l-m-2*lm2
        lim=limeta/2-ix
!
!  compute radial function of the second kind and its first derivative
!
!  backward series for denominator
        denom=(1.0q0,0.0q0)
        sumdpr=1.0q0
        sumdpi=0.0q0
        if (lm2.lt.1) go to 20
        dold=(1.0q0,0.0q0)
          do 10 j=lm2,1,-1
          jj=j+j+ix
          dnew=dold/(pratb(jj+1)*enr(j))
          denom=denom+dnew
          if(qreal(dnew).gt.0.0q0) sumdpr=sumdpr+qreal(dnew)
          if(qimag(dnew).gt.0.0q0) sumdpi=sumdpi+qimag(dnew)
          if(cqabs(dnew/denom).lt.dec) go to 20
          dold=dnew
10        continue
20      continue
!
!  forward series for denominator
        dold=(1.0q0,0.0q0)
          do 30 j=lm2+1,lim
          jj=j+j+ix
          dnew=dold*enr(j)*pratb(jj+1)
          denom=denom+dnew
          if(qreal(dnew).gt.0.0q0) sumdpr=sumdpr+qreal(dnew)
          if(qimag(dnew).gt.0.0q0) sumdpi=sumdpi+qimag(dnew)
          if(cqabs(dnew/denom).lt.dec) go to 40
          dold=dnew
30        continue
40      continue
        jden=j
        nsubdenr=-iqint(qlog10(qabs(qreal(denom)/sumdpr)+dec))
        nsubdeni=0
        if(sumdpi.ne.0.0q0) nsubdeni=-iqint(qlog10(qabs(qimag(denom)/&
                               sumdpi)+dec))
        if(nsubdenr.gt.ndec) nsubdenr=ndec
        if(nsubdeni.gt.ndec) nsubdeni=ndec
        naccdr=ndec-max(2,iqint(qlog10(ca)))-nsubdenr
        naccdi=ndec-max(2,iqint(qlog10(ca)))-nsubdeni
        nsubden=max(nsubdenr,nsubdeni)
        if(naccdr.lt.0) naccdr=0
        if(naccdi.lt.0) naccdi=0
        naccd=min(naccdr,naccdi)
        r2est=cqabs(sumcoef*denom)
        r2test=r2est*dconi
!
!  backward series for numerator
        dold=(1.0q0,0.0q0)
        doldd1=(1.0q0,0.0q0)
        doldd2=reld12
        r2temp=(1.0q0,0.0q0)
        sumnpr=1.0q0
        sumnpi=0.0q0
        sumdnpr=0.0q0
        sumdnpi=0.0q0
        r2dtemp=doldd2
        if(l.ne.0) r2dtemp=r2dtemp+doldd1
        if(qreal(r2dtemp).gt.0.q0) sumdnpr=qreal(r2dtemp)
        if(qimag(r2dtemp).gt.0.q0) sumdnpi=qimag(r2dtemp)
        if(lm2.lt.1) go to 60
          do 50 j=lm2,1,-1
          jj=j+j+ix
          dnew=-dold/(sneuf(jj+m)*pratt(jj+1)*enr(j))
          dnewd1=-doldd1/(sneuf(jj+m)*pdratt(jj+1)*enr(j))
          dnewd2=-doldd2/(sneudf(jj+m)*pratt(jj+1)*enr(j))
          r2temp=r2temp+dnew
          if(qreal(dnew).gt.0.0q0) sumnpr=sumnpr+qreal(dnew)
          if(qimag(dnew).gt.0.0q0) sumnpi=sumnpi+qimag(dnew)
          dnewd=dnewd1+dnewd2
          r2dtemp=r2dtemp+dnewd
          if(qreal(dnewd).gt.0.0q0) sumdnpr=sumdnpr+qreal(dnewd)
          if(qimag(dnewd).gt.0.0q0) sumdnpi=sumdnpi+qimag(dnewd)
          if(cqabs(dnew/r2temp)+cqabs(dnewd/r2dtemp).lt.dec) go to 60
          dold=dnew
          doldd1=dnewd1
          doldd2=dnewd2
50        continue
60      continue
        if(m.eq.0.and.jj.eq.2) then
        r2dtemp=r2dtemp-dnewd1
        if(qreal(dnewd1).gt.0.0q0) sumdnpr=sumdnpr-qreal(dnewd1)
        if(qimag(dnewd1).gt.0.0q0) sumdnpi=sumdnpi-qimag(dnewd1)
        end if
!
!  forward series for numerator
        dold=(1.0q0,0.0q0)
        doldd1=(1.0q0,0.0q0)
        doldd2=reld12
        test=1.0q0
        testd=1.0q0
        testm=1.0q0
        testdm=1.0q0
        js=lm2
        jds=lm2
        sr2temp=r2temp
        sr2dtemp=r2dtemp
        txr=sumnpr
        txi=sumnpi
        txdr=sumdnpr
        txdi=sumdnpi
        dnewsum=(0.0q0,0.0q0)
        dnewdsum=(0.0q0,0.0q0)
        doldd=reld12
        if(l.ne.0) doldd=reld12+(1.0q0,0.0q0)
        kount=0
        kountd=0
          do 110 j=lm2+1,lim-1
          jj=j+j+ix
          kount=kount+1
          kountd=kountd+1
          dnew=-dold*enr(j)*sneuf(jj+m)*pratt(jj+1)
          dnewd1=-doldd1*enr(j)*sneuf(jj+m)*pdratt(jj+1)
          dnewd2=-doldd2*enr(j)*sneudf(jj+m)*pratt(jj+1)
          if(((qreal(dnew)/qreal(dold)).le.0.0q0).or.(kount.eq.100))&
          go to 70
          dnewsum=dnewsum+dnew
          go to 80
70        r2temp=r2temp+dnewsum
          kount=0
          if(cqabs(r2temp).gt.r2test) go to 120
          if(qreal(dnewsum).gt.0.0q0) sumnpr=sumnpr+qreal(dnewsum)
          if(qimag(dnewsum).gt.0.0q0) sumnpi=sumnpi+qimag(dnewsum)
          if(cqabs(dnewsum).ne.0.0q0) test=cqabs(dnewsum/r2temp)
          dnewsum=dnew
          if(test.ge.testm) go to 80
          testm=test
          sr2temp=r2temp
          js=j
          txr=sumnpr
          txi=sumnpi
80        dnewd=dnewd1+dnewd2
          if(((qreal(dnewd)/qreal(doldd)).le.0.0q0).or.(kountd.eq.100))&
          go to 90
          dnewdsum=dnewdsum+dnewd
          go to 100
90        r2dtemp=r2dtemp+dnewdsum
          kountd=0
          if(qreal(dnewdsum).gt.0.0q0) sumdnpr=sumdnpr+qreal(dnewdsum)
          if(qimag(dnewdsum).gt.0.0q0) sumdnpi=sumdnpi+qimag(dnewdsum)
          if(cqabs(dnewdsum).ne.0.0q0) testd=cqabs(dnewdsum/r2dtemp)
          dnewdsum=dnewd
          if(testd.ge.testdm) go to 100
          testdm=testd
          sr2dtemp=r2dtemp
          jds=j
          txdr=sumdnpr
          txdi=sumdnpi
100       if ((test+testd).lt.dcon) go to 120
          if((dnew.eq.(0.0q0,0.0q0)).and.(dnewd.eq.(0.0q0,0.0q0)))&
           go to 120
          dold=dnew
          doldd1=dnewd1
          doldd2=dnewd2
          doldd=dnewd1+dnewd2
110       continue
120     r2temp=sr2temp
        r2dtemp=sr2dtemp
130     jmax=j
        jeta=js
        if(jds.gt.js) jeta=jds
        if(jden.gt.jeta) jeta=jden
        testdm=cqabs(r2dtemp)*testdm
        continue
        jtestm=-iqint(qlog10(qabs(testm)))
        if(jtestm.lt.0) jtestm=0
        if(jtestm.gt.ndec) jtestm=ndec
        jtestdm=-iqint(qlog10(cqabs(testdm/r2dtemp)))
        if(jtestdm.lt.0) jtestdm=0
        if(jtestdm.gt.ndec) jtestdm=ndec
        naccns1r=-iqint(qlog10(qabs(qreal(r2temp)/txr)+dec))
        naccns1i=0
        if(txi.ne.0.0q0) naccns1i=-iqint(qlog10(qabs(qimag(r2temp)&
                            /txi)+dec))
        if(naccns1r.lt.0) naccns1r=0
        if(naccns1r.gt.ndec) naccns1r=ndec
        if(naccns1i.lt.0) naccns1i=0
        if(naccns1i.gt.ndec) naccns1i=ndec
        naccns2r=-iqint(qlog10(qabs(qreal(r2dtemp)/txdr)+dec))
        naccns2i=0
        if(txdi.ne.0.0q0) naccns2i=-iqint(qlog10(qabs(qimag(r2dtemp)&
                            /txdi)+dec))
        if(naccns2r.lt.0) naccns2r=0
        if(naccns2r.gt.ndec) naccns2r=ndec
        if(naccns2i.lt.0) naccns2i=0
        if(naccns2i.gt.ndec) naccns2i=ndec
        naccns1=max(naccns1r,naccns1i)
        naccns2=max(naccns2r,naccns2i)
        naccns=max(naccns1,naccns2)
        naccn1=min(jtestm,ndec-naccns1)
        naccn2=min(jtestdm,ndec-naccns2)
        naccn=min(naccn1,naccn2)
        if(naccn.lt.0) naccn=0
        naccnmaxp=naccnmax
        if(naccn.gt.naccnmax) naccnmax=naccn
        if(naccnmax.gt.naccn.and.naccnmax.ge.3) iopnee=1
        if(ndec-naccns.lt.lowtest) iopnee=2
        if(ndec-naccns.lt.naccr) iopnee=2
        if(iopnee.eq.0) neemark=nee
!
!
!  combining results to form the radial function characteristics
!  r2c and r2dc and corresponding exponents ir2e and ir2de
        r2c=r2temp*sneun(l+1)*pcoefn/denom
        iterm=0
        if(r2c.ne.(0.0q0,0.0q0)) iterm=iqint(qlog10(cqabs(r2c)))
        ir2e=ineue(l+1)+ipcoefn+iterm
        r2c=r2c*10.0q0**(-iterm)
        r2dc=(r2dcoef1*r2dtemp*sneun(l+1)*pdcoefn/denom)*&
       10.0q0**(ineue(l+1)+ipdcoefn-ir2e)
        iterm=0
        if(r2dc.ne.(0.0q0,0.0q0)) iterm=iqint(qlog10(cqabs(r2dc)))
        ir2de=ir2e+iterm
 	r2dc=r2dc*10.0q0**(-iterm)
        write(40,140) jmax,jden,lim,js,jtestm,naccns1,jds,jtestdm,&
               naccns2,naccn,naccd
140     format(8x,'r2eta: numerator, denominator converged in ',&
         i6,' ,',i6,' terms; ',i6,' terms available.',/,&
         15x,'best r2 at ',i6,' terms with convergence to',i3,&
         ' digits;',i3,' digits subtr. error.',/,15x,&
         'best r2d at ',i6,' terms with convergence to ',i3,&
         ' digits;',i3,' digits subtr. error.',/,15x,&
         'estimated numerator and denominator accuracy is ',i4,&
         ' and',i4,' digits.')
        wronc=r1c*r2dc*(10.0q0**(ir1e+ir2de))-r2c*r1dc*&
              (10.q0**(ir2e+ir1de))
        wront=1.0q0/(c*x1*(x1+2.q0))
        nacceta=-iqint(qlog10(cqabs((wronc-wront)/wront)+dec))
        if(nacceta.lt.0) nacceta=0
        testw1=qabs(qreal((r1c*r2dc)*10.0q0**(ir1e+ir2de)))
        testw2=qabs(qreal((r2c*r1dc)*10.0q0**(ir2e+ir1de)))
        testw=testw1
        if(testw2.gt.testw1) testw=testw2
        naccsubre=-iqint(qlog10(qabs(qreal(wronc)/testw))+dec)
        testw1=qabs(qimag((r1c*r2dc)*10.0q0**(ir1e+ir2de)))
        testw2=qabs(qimag((r2c*r1dc)*10.0q0**(ir2e+ir1de)))
        testw=testw1
        if(testw2.gt.testw1) testw=testw2
        naccsubie=-iqint(qlog10(qabs(qimag(wronc)/testw)+dec))
        naccsube=min(naccsubre,naccsubie)
        if(ichoicer1.eq.1) naccsube=min(naccsube,idigc-&
                                 nsubf)
        if(ichoicer1.eq.2) naccsube=min(naccsube,idigc-&
                                 nsubmf)
        naccsube=min(naccsube,idigc-nsubden)
        if(naccsube.lt.0) naccsube=0
        if(naccsube.lt.2) go to 145
        if(ichoicer1.eq.1) nacceta=min(nacceta+naccsube,idigc-&
                                 nsubf,idigc-nsubden)
        if(ichoicer1.eq.2) nacceta=min(nacceta+naccsube,idigc-&
                                 nsubmf,idigc-nsubden)
        if(nacceta.lt.0) nacceta=0
145     continue
        nacciop=0
        if(naccn.ge.naccr1.and.naccd.ge.naccr1) then
          nacceta=max(min(naccn,naccd),nacceta)
          nacciop=2
        end if
        if((naccn-naccsube).gt.naccd.and.nacceta.lt.minacc.and.&
               naccd.lt.minacc) nacciop=1
        if(nacciop.ne.1) go to 160
        write (4,150)
150     format(15x,'denominator calculated using wronskian.')
        nacceta=naccn-naccsube
        r2c=r2c*wront/wronc
        iterm=0
        if(r2c.ne.(0.0q0,0.0q0)) iterm=iqint(qlog10(cqabs(r2c)))
        ir2e=ir2e+iterm
        r2c=r2c*10.0q0**(-iterm)
        r2dc=r2dc*wront/wronc
        iterm=0
        if(r2dc.ne.(0.0q0,0.0q0)) iterm=iqint(qlog10(cqabs(r2dc)))
        ir2de=ir2de+iterm
        r2dc=r2dc*10.0q0**(-iterm)
160     if(cqabs(r2c).ge.1.0q0) go to 170
        r2c=r2c*10.0q0
        ir2e=ir2e-1
170     continue
        if(cqabs(r2dc).ge.1.0q0) go to 180
        r2dc=r2dc*10.0q0
        ir2de=ir2de-1
180     continue
        if(nacciop.eq.0.and.naccn.gt.nacceta) naccnmax=max(naccnmaxp,&
                                                      nacceta)
        return
        end
!
!
!
        subroutine geteig (m,c,nbp,lnum,maxe,eigst)
!
!  purpose:     to determine estimates of the eigenvalues for
!               orders up to and beyond the break point for a given m
!               and complex c. These estimates will be starting values
!               for the bouwkamp procedure.
!
!  parameters:
!
!     input:    m     : m
!               c     : complex c
!               nbp   : breakpoint value
!               lnum  : number of values of l for which functions
!                       are desired
!               maxe  : dimension of matrix elements
!     output:   eigst : array of eigenvalue estimates
!
!  real*16 scalars and complex*32 scalars and arrays
        complex*32 c,c2,c4,d(maxe),e(maxe),f(maxe),eigst(lnum),g(maxe)
        real*16 cqr,cqi,c2m,xl,cm,em

!
        m2=m+m
        em=m
        c2=c*c
        c4=c2*c2
        cqi=qabs(qimag(c))
        cqr=qabs(qreal(c))
        cm=cqabs(c)
        c2m=cqabs(c2)
!  obtain starting eigenvalues eig(i) for the Bouwkamp procedure
!
            lime=max(100,nbp+nbp)
            imax=1.2*max(25,nbp/2)+5
              do 10 i=1,lime
              xl=qfloat(m+i+i-2)
              d(i)=xl*(xl+1.0q0)+c2*(2.0q0*xl*(xl+1.0q0)-2.0q0*em*em-&
             1.0q0)/((2.0q0*xl-1.0q0)*(2.0q0*xl+3.0q0))
              d(i)=d(i)/c2m
10            continue
            nm1=lime-1
              do 20 i=1,nm1
              xl=qfloat(m+i+i-2)
              e(i)=c2*(1.0q0/(2.0q0*xl+3.0q0))*qsqrt(((xl+2.0q0+em)*&
             (xl+1.0q0+em)*(xl+2.0q0-em)*(xl+(1.0q0)-em))/&
             ((2.0q0*xl+5.0q0)*(2.0q0*xl+1.0q0)))
              e(i)=e(i)/c2m
20            continue
            call tridiag(d,e,lime,maxe)
              do 30 i=1,lime
              f(i)=d(i)*c2m
30            continue
              do 40 i=1,lime
              xl=qfloat(m+i+i-1)
              d(i)=xl*(xl+1.0q0)+c2*(2.0q0*xl*(xl+1.0q0)-2.0q0*em*em-&
             1.0q0)/((2.0q0*xl-1.0q0)*(2.0q0*xl+3.0q0))
              d(i)=d(i)/c2m
40            continue
            nm1=lime-1
              do 50 i=1,nm1
              xl=qfloat(m+i+i-1)
              e(i)=c2*(1.0q0/(2.0q0*xl+3.0q0))*qsqrt(((xl+2.0q0+em)*&
             (xl+1.0q0+em)*(xl+2.0q0-em)*(xl+(1.0q0)-em))/&
             ((2.0q0*xl+5.0q0)*(2.0q0*xl+1.0q0)))
              e(i)=e(i)/c2m
50            continue
            call tridiag(d,e,lime,maxe)
              do 60 i=1,lime
              g(i)=d(i)*c2m
60            continue
            call eigorder(c,m,lime,maxe,nbp,f,g)
            lnumo2=lnum/2
            iupp=min(imax,lnumo2)
            if(iupp.eq.0) go to 80
              do 70 i=1,iupp
              eigst(i+i-1)=f(i)
              eigst(i+i)=g(i)
70            continue
            lnumo2=lnum/2
80          if(2*lnumo2.ne.lnum.and.imax.gt.lnumo2)&
           eigst(lnum)=f(lnumo2+1)
            return
            end
!
!
!
        subroutine tridiag(d,e,n,np)
!
!  Purpose:     To calculate the eigenvalues of a symmetric tridiagonal
!               matrix by the implicit ql method and order them in
!               ascending value.
!
! Subroutine tridiag is a translation by Brian Gladman of the algol
! procedure imtql2, num. math. 12, 377-383(1968) by martin and
! wilkinson,as modified  in num. math. 15, 450(1970) by dubrulle,
! handbook for auto. comp.,vol.ii-linear algebra, 241-248(1971). It
! has been converted to complex arithmetic for use in cprofcn.
! The option and code for finding the eigenvectors has been removed.
!
!  parameters:
!
!     input:    d  : diagonal elements of the tridiagonal matrix
!               e  : the subdiagonal elements of the matrix
!               n  : order of the matrix
!               np : dimension of the vectors d and e
!
!     output:   d  : eigenvalues in ascending value.
!
        complex*32 b,c,f,g,p,r,s
        complex*32 d(np),e(np)
        real*16 dp
        integer ev_no
!
        do ev_no = 1, n
!  iterate to find the next eigenvalue
          do it = 0, 30
!  look for a small sub-diagonal element
             do i = ev_no, n - 1
             dp = cqabs(d(i)) + cqabs(d(i + 1))
             if(dp + cqabs(e(i)).eq.dp) exit
             end do
!
!  end the iteration if we have an eigenvalue
!          if(i.eq.ev_no) exit
!  form an implicit shift
          g = (d(ev_no + 1) - d(ev_no)) / (qcmplx(2.0q0,0.0q0)*e(ev_no))
          r = diag(g, qcmplx(1.0q0,0.0q0))
          if(qreal(g)*qreal(r).ge.0.0q0)&
       g = d(i) - d(ev_no) + e(ev_no) / (g + r)
          if(qreal(g)*qreal(r).lt.0.0q0)&
       g = d(i) - d(ev_no) + e(ev_no) / (g - r)
          s = qcmplx(1.0q0,0.0q0)
          c = qcmplx(1.0q0,0.0q0)
          p = qcmplx(0.0q0,0.0q0)
!
!  perform a plane rotation followed by Givens
!  rotations to restore tridiagonal form
            do j = i, ev_no + 1, -1
            f = s * e(j - 1)
            b = c * e(j - 1)
            r = diag(f, g)
            e(j) = r
!  recover from underflow
            if(qreal(r).eq.0.0q0.and.qimag(r).eq.0.0q0) exit
            s = f / r
            c = g / r
            g = d(j) - p
            r = (d(j - 1) - g) * s + qcmplx(2.0q0,0.0q0) * c * b
            p = s * r
            d(j) = g + p
            g = c * r - b
!
            end do
          d(j) = d(j) - p
          if(j.eq.ev_no) e(j) = g
          e(i) = qcmplx(0.0q0,0.0q0)
          end do
        end do
        return
        end
!
!
!
        function diag(x, y)
        complex*32 diag,x,y
        if(cqabs(x).ge.cqabs(y)) then
          diag = x * cqsqrt(qcmplx(1.0q0,0.0q0) + (y / x) ** 2)
        else
          diag = y * cqsqrt(qcmplx(1.0q0,0.0q0) + (x / y) ** 2)
        end if
        return
        end
!
!
!
        subroutine eigorder(c,m,n,np,nbp,f,g)
!
!
!  purpose:     To order the eigenvalues estimates so that the
!               resulting radial function magnitudes have a similar
!               behavior with increasing order as the corresponding
!               spherical Bessel functions
!  parameters:
!
!     input:    c     : complex c
!               m     : m
!               n     : number of even eigenvalue estimates obtained
!                       from the tridiagonal matrix; also the number
!                       of odd eigenvalue estimates
!               np    : dimension of both the array of even estimates
!                       and the array of odd estimates
!               nbp   : breakpoint value
!
!     output:   f     : array of ordered even eigenvalue estimates
!               g     : array of ordered odd eigenvalue estimates
!
!  integer, real*16 and complex*32 scalars and arrays
		implicit none
        complex*32 c,f(np),g(np),fm(np),gm(np),fs(np),gs(np),p,testp
        integer nflage(n),nflago(n), m, n, np, nbp, iflag, il, ilow,&
	i, imax, isav, j, jlowf, jlowg, k, kk, limp, limpr, nmax,&
  limpc
        real*16 test1,test2,cr,ci,cm2,cr2,ci2,testpr,testpi
!
		fm = 0
		gm = 0
		fs = 0
		gs = 0
		nflage = 0
		nflago = 0
        limp = 0
        imax=1.2*max(25,nbp/2)+5
        cr=qreal(c)
        ci=qimag(c)
        cr2=cr*cr
        ci2=ci*ci
        cm2=cr2+ci2
!
!  order even eigenvalues in ascending real part
         do i = 1,n-1
         k = i
         p=f(i)
           do j = i + 1, n
           if(qreal(f(j)).lt.qreal(p)) then
           k = j
           p = f(j)
           end if
           end do
         if(k.ne.i) then
         f(k)=f(i)
         f(i)=p
         end if
         end do
!
!  order odd eigenvalues in ascending real part
         do i = 1,n-1
         k = i
         p=g(i)
           do j = i + 1, n
           if(qreal(g(j)).lt.qreal(p)) then
           k = j
           p = g(j)
           end if
           end do
         if(k.ne.i) then
         g(k)=g(i)
         g(i)=p
         end if
         end do
         do k=1,imax
         end do
!
!  locate oblate-like eigenvalues
         nmax=iqint(qabs(ci))/2
         if(nmax.eq.0) go to 50
         limp=0
           do j=1,nmax
           testpr=cr2-ci2+qabs(ci)*(2*(m+2*j-1))-((m+j)*(2*j-1)-j+1)*&
            (1.0q0+(qabs(ci)/cm2)*(m+2*j-1)/4.0q0)
           testpi=2*cr*(qabs(ci)-m-2*j+1)-((m+j)*(2*j-1)-j+1)*&
            (cr/cm2)*(m+2*j-1)/4.0q0
           testpi=sign(testpi,ci)
           testp=qcmplx(testpr,testpi)
             do i=1,imax
             if(qabs((qreal(f(i))-testpr)/qreal(f(i))).gt.0.01q0)&
          go to 10
             if(qabs((qimag(f(i))-testpi)/qimag(f(i))).gt.0.01q0)&
          go to 10
             limp=limp+1
             fs(limp)=f(i)
             isav=i
               do kk=i,imax-limp
               f(kk)=f(kk+1)
               end do
             go to 20
10           end do
           go to 50
20           il=max(isav-1,1)
             do i=il,isav
             if(qabs((qreal(g(i))-testpr)/qreal(g(i))).gt.0.01q0)&
           go to 30
             if(qabs((qimag(g(i))-testpi)/qimag(g(i))).gt.0.01q0)&
           go to 30
             gs(limp)=g(i)
               do kk=i,imax-limp
               g(kk)=g(kk+1)
               end do
             go to 40
30           continue
             end do
           gs(limp)=g(isav)
           do kk=isav,imax-limp
           g(kk)=g(kk+1)
           end do
           go to 50
40         continue
           end do
50       continue
!
!  order oblate-like even eigenvalues for insertion in eigenvalue
!  sequence
         if(limp.eq.0) go to 160
         limpc=limp/2+1
         limpr=(limp+1)/2
           do i = 1, limpr
			if (limpc+i-1 <= np .and. i+i-1 <= np) then
           fm(limpc+i-1)=fs(i+i-1)
		   endif
           end do
         if(limp.eq.1) go to 60
           do i = 1, limp-limpr
		   if (limpc-i <= np .and. i+i <= np) then
           fm(limpc-i)=fs(i+i)
		   endif
           end do
60       continue
!
!  order oblate-like odd eigenvalues for insertion in eigenvalue
!  sequence
         limpc=(limp+1)/2
         limpr=limp/2
           do i = 1, limp-limpr
		   if (limpc-i+1 <= np .and. i+i-1 <= np ) then
           gm(limpc-i+1)=gs(i+i-1)
		   endif
           end do
           if(limp.eq.1) go to 70
           do i = 1, limpr
		   if (limpc+i <= np .and. i+i <= np) then
           gm(limpc+i)=gs(i+i)
		   endif
           end do
70       continue
!
!  insert paired eigenvalues in proper location in sequence
         iflag=1
         if(limp.ne.(2*(limp/2))) iflag=-iflag
         ilow=1
80       continue
           do i = ilow,imax
           if(nflage(i).ne.1) go to 90
             do j = i,imax-1
             f(j)=f(j+1)
             nflage(j)=nflage(j+1)
             end do
           ilow=i
           go to 80
90         end do
         ilow=1
100      continue
           do i = ilow,imax
           if(nflago(i).ne.1) go to 110
             do j = i,imax-1
             g(j)=g(j+1)
             nflago(j)=nflago(j+1)
             end do
           ilow=i
           go to 100
110        end do
!
           do j=1,imax-limp-1
           if(qabs(qimag(g(j)))-qabs(qimag(f(j))).lt.0.0q0) exit
           if(qabs(qimag(f(j+1)))-qabs(qimag(g(j))).lt.0.0q0) go to 130
           end do
!		   write(*,*) 'j = ', j
         test1=4.0q0*qimag(g(j-1))-6.0q0*qimag(f(j-1))+&
         4.0q0*qimag(g(j-2))-qimag(f(j-2))
         test2=4.0q0*qimag(g(j))-6.0q0*qimag(f(j+1))+&
         4.0q0*qimag(g(j+1))-qimag(f(j+2))
         if(qabs(qimag(f(j))-test1).gt.qabs(qimag(f(j))-test2))go to 120
         jlowf=j+1
         jlowg=j
         iflag=-iflag
         go to 150
120      jlowf=j
         jlowg=j
         go to 150
130      test1=4.0q0*qimag(f(j))-6.0q0*qimag(g(j-1))+&
         4.0q0*qimag(f(j-1))-qimag(g(j-2))
         test2=4.0q0*qimag(f(j+1))-6.0q0*qimag(g(j+1))+&
         4.0q0*qimag(f(j+2))-qimag(g(j+2))
         if(qabs(qimag(g(j))-test1).gt.qabs(qimag(g(j))-test2))go to 140
         jlowf=j+1
         jlowg=j+1
         go to 150
140      jlowf=j+1
         jlowg=j
         iflag=-iflag
150      continue
!			write(*,*) 'imax=', imax, 'limp=', limp, 'jlowf=', jlowf
           do i =1,imax-jlowf-limp+1
		   if (imax-i+1 <= np .and. imax-limp-i+1 <= np) then
           f(imax-i+1)=f(imax-limp-i+1)
		   endif
           end do
           do i =1,imax-jlowg-limp+1
		   if (imax-i+1 > 0 .and. imax-i+1 <= np .and.&
		   imax-limp-i+1 > 0 .and. imax-limp-i+1 <= np) then
			g(imax-i+1)=g(imax-limp-i+1)
		   endif
           end do
         k=jlowf-1
           do i = 1,limp
           k=k+1
           if(iflag.eq.1) f(k)=fm(i)
           if(iflag.eq.-1) f(k)=fm(limp-i+1)
           end do
         k=jlowg-1
           do i = 1,limp
           k=k+1
           if(iflag.eq.-1) g(k)=gm(limp-i+1)
           if(iflag.eq.1) g(k)=gm(i)
           end do
160      continue
         return
         end
!
!
!
        subroutine conver (l,m,c,limd,ndec,maxd,blist,glist,&
                     lical,ioprad,ienr,eigval,enr,idigc,itestm)
!
!  purpose:     to determine the eigenvalue using the boukwamp method.
!  parameters:
!
!     input:    l     : l
!               m     : m
!               c     : complex c
!               limd  : number of enr values computed
!               ndec  : number of decimal digits available in r*16
!                       arithmetic
!               maxd  : dimension of enr,blist,glist arrays
!               blist : coefficients used in computing d coefficients
!               glist : coefficients used in computing d coefficients
!               eigval: estimated value of the eigenvalue
!               lical : number of eigenvalue estimates that are
!                       believed accurate enough so that the Bouwkamp
!                       procedure is not required to refine the estimate
!               ioprad: integer equal to 0 if radial functions are not
!                       desired; equal to 1 if radial functions of the
!                       first kind only are required; equal to 2 if
!                       radial functions of both the first and second
!                       kinds are desired
!               ienr  : number of d coeficient ratios required for
!                       convergence of the sum involved in the
!                       Bouwkamp eigenvalue procedure for the previous
!                       value of l
!
!     output:   eigval: converged eigenvalue
!               enr   : array of limd values
!               idigc : number of decimal digits of convergence
!                       for the Boouwkamp procedure
!               itestm: maximum number of decimal digits of agreement
!                       for one of the d coefficient ratios calculated
!                       forward and for the same d coefficient ratio
!                       calculated backward using recursion
!               ienr  : number of d coeficient ratios required for
!                       convergence of the sum involved in the
!                       Bouwkamp eigenvalue procedure for the current
!                       value of l
!
!  real*16 scalars and complex*32 scalars and arrays
        real*16 dec,eigdec,eigtest
        complex*32 c,cora,corb,de,dl,eigval,enrc
        complex*32 blist(maxd),eigst,enr(maxd),glist(maxd),enrf(maxd)
!
		enrf = 0
        eigst=eigval
        dec=10.0q0**(-ndec-1)
        eigdec=dec*100.q0
        ncsave=0
        li=l-m+1
        lm2=(l-m)/2
        limdb=2*ienr+20
        if(l.eq.m.or.l.eq.m+1) limdb=2*ienr
        if(limdb.gt.limd) limdb=limd

!		write(*,*) 'conver limd = ', limd, 'limbd =', limdb
!
!  begin bouwkamp procedure
        ix=l-m-2*lm2
        ifc=1
        lim2=limdb/2-ix
        iglim=lim2+1
        irio=lm2+1
        iw1=lm2+2
        li=l-m+1
10      enr(1)=eigval-glist(1)
!		write(*,*) 'lim2=', lim2
        if(irio.lt.2) go to 30
!
!  evaluate the continued fraction
          do 20 i=1,irio-1
          enr(i+1)=-blist(i)/enr(i)-glist(i+1)+eigval
20        continue
30      enr(lim2)=-blist(lim2)/(glist(iglim)-eigval)
        iw15=lim2-1
        ip=iw1+iw15
!
!  evaluate the continued fraction
          do 40 i=iw1,iw15
          ipi=ip-i
          enr(ipi)=-blist(ipi)/(glist(ipi+1)-eigval+enr(ipi+1))
40        continue
50      enrc=-blist(irio)/(glist(irio+1)-eigval+enr(irio+1))
        de=enrc*enrc/blist(irio)
        corb=de
!
!  compute denominator in the bouwkamp
          do 60 i=iw1,lim2
          de=enr(i)*enr(i)/blist(i)*de
          corb=corb+de
          if(cqabs(de/corb).lt.dec) go to 70
60        continue
70      ienr=i
        de=qcmplx(1.0q0,0.0q0)
        cora=de
!
!  compute the denominator in the bouwkamp
          do 80 i=1,lm2
          de=blist(irio-i)/(enr(irio-i)*enr(irio-i))*de
          cora=cora+de
          if(cqabs(de/cora).lt.dec) go to 90
80        continue
!
!  compute the correction to the eigenvalue
90      dl=(enrc-enr(irio))/(cora+corb)
        eigval=dl+eigval
!
!  eigenvalue accurate enough?
        eigtest=cqabs(dl/eigval)
        if(eigtest.lt.eigdec) go to 100
        nc=-iqint(qlog10(eigtest+dec))
        if(nc.le.ncsave+2) go to 100
        ncsave=nc
        ifc=ifc+1
        if(ifc.lt.10) go to 10
100     continue
        idigc=-iqint(qlog10(eigtest+dec))
        if(idigc.gt.ndec) idigc=ndec
        ncon=-iqint(qlog10(cqabs((eigst-eigval)/eigval)+dec))
        if(ncon.gt.ndec) ncon=ndec
        jflagc=0
        if(idigc.lt.ndec-2.and.ncon.ge.idigc-2.and.li.le.lical) then
          eigval=eigst
          jflagc=1
        end if
        idcoef=1
!
!  calculate the d constant ratios (enr)
110     enr(1)=eigval-glist(1)
        if(lm2.lt.2) go to 125
          do 120 i=1,lm2-1
          enr(i+1)=-blist(i)/enr(i)-glist(i+1)+eigval
120       continue
125     continue
        if(lm2.ne.0) enrc=-blist(lm2)/enr(lm2)-glist(lm2+1)+eigval
        if(lm2.eq.0) enrc=enr(1)
        lim2=limd/2-ix
        enr(lim2)=-blist(lim2)/(glist(lim2+1)-eigval)
        ilow=lm2+1
        iupp=lim2-1
        ip=ilow+iupp
!		write(*,*) '2 lim2 =', lim2
          do 130 i=ilow,iupp
          ipi=ip-i
          enr(ipi)=-blist(ipi)/(glist(ipi+1)-eigval+enr(ipi+1))
130       continue
        itestm=-iqint(qlog10(cqabs((enrc-enr(lm2+1))/enr(lm2+1))+&
           dec))
        if(itestm.gt.ndec) itestm=ndec
        if(jflagc.eq.0) write(40,140) l,eigst,idigc,ifc
140     format(1x,'l = ',i4,' eigen. estimate',e40.31,e40.31,&
         /,10x,'conv. to',i3,' dig. at ifc =',i3)
        if(jflagc.eq.1) write(40,150) l,eigst,idigc,ifc
150     format(1x,'l = ',i4,' eigen. estimate',e40.31,e40.31,&
         /,10x,'conv. to',i3,' dig. at ifc =',i3,'. Estimate is',&
         ' likely more accurate than converged value and will',&
         ' be used.')
        if(idcoef.eq.0.and.ioprad.ne.0) write(40,200) lim2,itestm,lm2+1
        if(idcoef.eq.0.and.ioprad.eq.0) write(50,200) lim2,itestm,lm2+1
        if(idcoef.eq.0) go to 220
        if(lm2.eq.0) go to 170
          do 160 i=lm2,1,-1
160       enr(i)=-blist(i)/(glist(i+1)-eigval+enr(i+1))
170     continue
        enrf(1)=eigval-glist(1)
        itestm=1
!		write(*,*) '3 lim2 = ', lim2
          do 180 n=2,lim2
          enrf(n)=-blist(n-1)/enrf(n-1)-glist(n)+eigval
          itest=-iqint(qlog10(cqabs((enrf(n)-enr(n))/enrf(n))+&
           dec))
          if(itest.gt.ndec) itest=ndec
          if(itest.le.itestm) go to 180
          itestm=itest
!		  write(*,*) 'assigning nlim, n = ', n
          nlim=n
          if(itestm.ge.ndec-1) go to 190
180       continue
190     continue
        idigc=itestm
        if(ioprad.ne.0) write(40,200) lim2,itestm,nlim
        if(ioprad.eq.0) write(50,200) lim2,itestm,nlim
200     format(15x,i6,' "d" coefs.: forward and backward',&
         ' recursion match to ',i3,' digits at n = ',i3)
!			write(*,*) 'nlim = ', nlim
          do 210 n=1,min(nlim, size(enr))
          enr(n)=enrf(n)
210       continue
220     continue
        return
        end
!
!
!
        subroutine dnorm (l,m,c,ndec,limd,maxd,enr,sgn,d01,id01,&
                     dmfnorm,dfnorm,jmf,nsubmf,jfla,nsubf)
!
!  purpose:     to compute unscaled d coefficient ratios using scaled
!               ratios, to compute the Morse-Feshbach normalization
!               factor for the d coefficients, to compute the Flammer
!               normalization factor for the d coefficients, to compute
!               the characteristic and exponent of the first d constant,
!               either d0 or d1, and to compute the sign of the d
!               constant for n=l-m.
!
!
!  parameters:
!
!     input:    l       : l
!               m       : m
!               c       : complex c
!               ndec    : number of decimal digits available in
!                         r*16 arithmetic
!               limd    : approximately twice the maximum number
!                         of terms available to be taken in the sum
!               maxd    : dimension of enr array
!               enr     : array of scaled d coefficient ratios
!
!     output:   enr     : array of unscaled d coefficient ratios
!                         enr(i) = ratio of the d coefficient with
!                         subscript 2*i+ix to the d coefficient with
!                         subscript 2*(i-1)+ix. Here ix =0 when l-m is
!                         even and ix=1 when l-m is odd.
!                         If the user needs the d coefficent ratios,
!                         they are available below right before
!                         statement 20.
!               sgn     : sign of the d coefficient for n=l-m
!               d01     : characteristic of ratio of first d
!                         constant (either d0 or d1, depending on
!                         whether l-m is even or odd) to the d
!                         constant of order equal to l-m
!               id01    : exponent corresponding to d01
!               dmfnorm : Morse-Feshbach normalization factor of the
!                         d coefficients. equal to the reciprocal of
!                         the value of the d constant d(n = l - m) using
!                         this normalization for the angular functions
!               dfnorm  : Flammer normalization factor of the
!                         d coefficients. equal to the reciprocal of
!                         the value of the d constant d(n = l - m) using
!                         this normalization for the angular functions
!               jmf     : maximum index of enr required for convergence
!                         of dmfnorm
!               nsubmf  : number of decimal digits of subtraction error
!                         incurred in calculating dmfnorm
!               jfla    : maximum index of enr required for convergence
!                         of dfnorm
!               nsubf   : number of decimal digits of subtraction error
!                         incurred in calculating dfnorm
!
!  real*16 scalars and complex*32 scalars and array
        real*16 aj,arr,dec,dfnormi,dfnormr,dmfnormi,dmfnormr,ea,enrr,&
          diterm,rm2,sgn,sumip,sumrp,termi,termr
        complex*32 c,csq,decinv,dfnorm,dmfnorm,d01,term,en
        complex*32 enr(maxd)
!
        rm2=qfloat(m+m)
        dec=10.0q0**(-ndec-1)
        decinv=qcmplx(10.0q0**(ndec+1),0.0q0)
        csq=c*c
        lm2=(l-m)/2
        ix=l-m-2*lm2
        lim2=limd/2-ix
        mml=m+m-1+ix
        sgn=1.0q0
!       en=(1.0q0,0.0q0)
          do 20 i=1,lim2
          arr=qfloat(ix+i+i)
          ea=arr+arr+rm2
          if(i.gt.lm2) go to 10
          enrr=enr(i)
          if(enrr.lt.0.0q0)  sgn=sgn*(-1.0q0)
10        enr(i)=(ea-1.0q0)*(ea+1.0q0)*enr(i)/((arr+rm2)*&
             (arr+rm2-1.0q0)*csq)
20        continue
!          en=(1.0q0,1.0q0)
!          do 15 i=1,lim2
!          en=en*enr(i)
!          if(cqabs(en).lt.dec) go to 16
!15        continue
!16        continue
!
!  compute the Morse-Feshbach normalization factor
!    forward summation of series
        term=(1.0q0,0.0q0)
        dmfnorm=term
        sumrp=1.0q0
        sumip=0.0q0
        jlow=l-m+2
        jterm=lm2
          do 30 j=jlow,limd,2
          aj=qfloat(j)
          jterm=jterm+1
          term=term*(aj+rm2)*enr(jterm)*(aj+rm2-1.0q0)/&
         (aj*(aj-1.0q0))
          termr=term
          if(termr.gt.0.0q0) sumrp=sumrp+termr
          termi=qimag(term)
          if(termi.gt.0.0q0) sumip=sumip+termi
          dmfnorm=dmfnorm+term
          if(cqabs(term).lt.dec) go to 40
30        continue
40      jlow=l-m
        jmf=jterm
!    backward summation of series
        if(jlow.lt.2) go to 60
        term=(1.0q0,0.0q0)
        jterm=lm2
          do 50 j=jlow,2,-2
          aj=qfloat(j)
          term=term*aj*(aj-1.0q0)/((aj+rm2&
         -1.0q0)*(aj+rm2)*enr(jterm))
          termr=term
          if(termr.gt.0.0q0) sumrp=sumrp+termr
          termi=qimag(term)
          if(termi.gt.0.0q0) sumip=sumip+termi
          jterm=jterm-1
          dmfnorm=dmfnorm+term
          if(cqabs(term).lt.dec) go to 60
50        continue
60      continue
        dmfnormr=dmfnorm
        jsubr=0
        if(dmfnormr*sumrp.ne.0.0q0)&
        jsubr=iqint(qlog10(qabs(sumrp/dmfnormr)+dec))
        if(jsubr.lt.0) jsubr=0
        jsubi=0
        dmfnormi=qimag(dmfnorm)
        if(dmfnormi*sumip.ne.0.0q0)&
subi=iqint(qlog10(qabs(sumip/dmfnormi)+dec))
        if(jsubi.lt.0) jsubi=0
        nsubmf=max(jsubr,jsubi)
        if(nsubmf.gt.ndec) nsubmf=ndec
        write(40,70) jmf,nsubmf
70      format(28x'Morse-Feshbach norm. converged in ',&
         i6,' terms with ',i3,' digits subt. error.')
!
!  compute the Flammer normalization factor
!    forward summation of series
        term=(1.0q0,0.0q0)
        dfnorm=term
        sumrp=1.0q0
        sumip=0.0q0
          do 80 j=lm2+1,lim2
          jj=j+j+ix
          term=-term*enr(j)*qfloat((jj+mml))/&
         (qfloat(jj-ix))
          termr=term
          if(termr.gt.0.0q0) sumrp=sumrp+termr
          termi=qimag(term)
          if(termi.gt.0.0q0) sumip=sumip+termi
          dfnorm=dfnorm+term
          if(cqabs(term/dfnorm).lt.dec) go to 90
80        continue
90      continue
        jfla=min(j,lim2)
!
!    backward summation of series
        if(lm2.lt.1) go to 110
        term=(1.0q0,0.0q0)
          do 100 j=lm2,1,-1
          jj=j+j+ix
          term=-term*(jj-ix)/(qfloat(jj+mml)&
         *enr(j))
          termr=term
          if(termr.gt.0.0q0) sumrp=sumrp+termr
          termi=qimag(term)
          if(termi.gt.0.0q0) sumip=sumip+termi
          dfnorm=dfnorm+term
          if(cqabs(term/dfnorm).lt.dec) go to 110
100       continue
110     continue
        dfnormr=dfnorm
        jsubr=0
        if(dfnormr*sumrp.ne.0.0q0)&
        jsubr=iqint(qlog10(qabs(sumrp/dfnormr)+dec))
        if(jsubr.lt.0) jsubr=0
        jsubi=0
        dfnormi=qimag(dfnorm)
        if(dfnormi*sumip.ne.0.0q0)&
subi=iqint(qlog10(qabs(sumip/dfnormi)+dec))
        if(jsubi.lt.0) jsubi=0
        nsubf=max(jsubr,jsubi)
        if(nsubf.gt.ndec) nsubf=ndec
        write(40,120) jfla,nsubf
120     format(28x,'Flammer norm. converged in ',i6,' terms; ',&
         'with ',i2,' digits subt. error.')
!
!  compute the d0(c|ml) or d1(c|ml)
        id01=0
   	d01=(1.0q0,0.0q0)
        if(lm2.eq.0) go to 140
          do 130 kjl=1,lm2
          kkjl=lm2-kjl+1
    	  d01=d01/enr(kkjl)
          if(cqabs(d01).gt.dec) go to 130
          d01=d01*decinv
          id01=id01-ndec-1
130       continue
        iterm=iqint(qlog10(cqabs(d01)))
        d01=d01*(10.0q0**(-iterm))
!		write(*,*) 'this must be final d01 =', d01
        id01=id01+iterm
140     continue
        return
        end
!
!
!
        subroutine dalt (l,m,c,limdr,maxdr,maxmp,ioppsum,eigval,&
                   enrneg,drhor,dneg,idneg)
!
!  purpose:     to calculate d ratios with negative subscripts
!               and d-rho ratios.
!  parameters:
!
!     input:    l       : l
!               m       : m
!               c       : complex c
!               limdr   : number of ratios of successive d-rho
!                         coefficients calculated
!               maxdr   : dimension of drhor array
!               maxmp   : dimension of enrneg array
!               ioppsum : integer index = 0 if no d rho coefficients
!                         are calculated (psum not needed for r2leg)
!               eigval  : eigenvalue
!
!     output:   enrneg  : array of d coefficient ratios with
!                         negative subscripts
!               drhor   : array of d rho coefficient ratios
!               dneg    : chracteristic of the first d coefficient
!                         with negative index
!               idneg   : exponent (base 10) of the first d coefficient
!                         with negative index
!
!  real*16 scalars and complex*32 scalars and arrays
        real*16 t,uterm,wterm
        complex*32 c,dneg,eigval,vterm
        complex*32 enrneg(maxmp),drhor(maxdr)
!
!  if l-m is even, ix=0; if l-m is odd, ix=1
        ix=(l-m)-2*((l-m)/2)
        t=qfloat(m+m-ix-ix)
!
!  calculate ratios of d constants with negative subscripts
        if(m.eq.0) go to 30
          do 10 i=1,m+1
          enrneg(i)=(0.0q0,0.0q0)
10        continue
        n=2-2*m+ix
        call un(n,m,uterm)
        call vnm2(n,m,c,eigval,vterm)
!
!       enrneg(k) = { d(-2m+2k-2)/d(-2m+2k), l-m even    }
!                   { d(-2m+1+2k-2)/d(-2m+1+2k), l-m odd }
!
!       calculations continue up to and including
!       enrneg(k=m) = { d(-2)/d(0), l-m even }
!                     { d(-1)/d(1), l-m odd  }
!
        enrneg(1)=-uterm/vterm
        dneg=enrneg(1)
        idneg=0
          do 20 i=2,2*m,2
          ii=i-2*m+ix
          if (ii.ge.0) go to 30
          n=ii+2
          call un(n,m,uterm)
          call vnm2(n,m,c,eigval,vterm)
          call wnm4(n,m,wterm)
          j=i/2
          enrneg(j+1)=-uterm/(wterm*enrneg(j)+vterm)
          dneg=dneg*enrneg(j+1)
          if(cqabs(dneg).gt.1.0q-10) go to 20
          dneg=dneg*1.0q10
          idneg=idneg-10
20        continue
        iterm=iqint(qlog10(cqabs(dneg)))
        dneg=dneg*(10.0q0**(-iterm))
        idneg=idneg+iterm
30      if(ioppsum.eq.0) go to 60
!
!  calculate ratios of d rho constants
!
!       drhor(k-m) = { d(rho|2k)/d(rho|2k-2), l-m even  }
!                    { d(rho|2k-1)/d(rho|2k-3), l-m odd }
!
          do 40 i=1,limdr
          drhor(i)=(0.0q0,0.0q0)
40        continue
!
          do 50 i=2*limdr,6,-2
          n=4-i+ix-m-m
          ii=(i-2)/2
          call un(n,m,uterm)
          call vnm2(n,m,c,eigval,vterm)
          call wnm4(n,m,wterm)
          drhor(ii)=-uterm/(wterm*drhor(ii+1)+vterm)
50        continue
        n=-2*m+ix
        call vnm2(n,m,c,eigval,vterm)
        call wnm4(n,m,wterm)
        ii=ii-1
!
!       the final value of ii is 1;
!       drhor(1) has a special value:
!       drhor(1) = { d(rho|2m+2)/d(-2m), l-m even  }
!                  { d(rho|2m+1)/d(-2m+1), l-m odd }
!
        drhor(ii)=1.0q0/((t-1.0q0)*(t+1.0q0)*(wterm*drhor(ii+1)+&
              vterm))
        if(ix.eq.1) drhor(ii)=-drhor(ii)
60      continue
        return
        end
!
!
!
        subroutine un (n,m,uterm)
!
!  purpose:     to evaluate the expression u  for given n and m.
!                                           n
!
!  parameters:
!
!     input:    n      : subscript of u
!               m      : m
!
!     output:   uterm  : u
!                         n
!
!  real*16 scalars
        real*16 r,rn,uterm
!
        rn=qfloat(n)
        r=qfloat(n+m+m)
        uterm=r*(r-1.q0)/((rn+r+1.q0)*(rn+r-1.q0))
        return
        end
!
!
!
        subroutine vnm2 (n,m,c,eigval,vterm)
!
!  purpose:     to evaluate the expression v    for given n, m, c, and
!               eigval.                     n-2
!
!  parameters:
!
!     input:    n      : subscript of v
!               m      : m
!               c      : c
!               eigval : eigenvalue
!
!     output:   vterm  : v
!                         n-2
!
!  real*16 and complex*32 scalars
        real*16 r,rm
        complex*32 c,eigval,vterm
!
        rm=qfloat(m)
        r=qfloat(n+m-2)
        vterm=(2.q0*r*(r+1.q0)-2.q0*rm*rm-1.q0)/((2.q0*r+3.q0)*&
         (2.q0*r-1.q0))+(r*(r+1.q0)-eigval)/(c*c)
        return
        end
!
!
!
        subroutine wnm4 (n,m,wterm)
!
!  purpose:     to evaluate the expression w    for given n and m.
!                                           n-4
!  parameters:
!
!     input:    n      : subscript of w
!               m      : m
!
!     output:   wtermm : w
!                         n-4
!
!  real*16 scalars
        real*16 r,rm,wterm
!
        rm=qfloat(m)
        r=qfloat(n-4)
        wterm=(r+2.q0)*(r+1.q0)/((r+r+rm+rm+3.q0)*(r+r+rm+rm+1.q0))
        return
        end
!
!
!
	subroutine gauss (ndec,n,x,w)
!
!  purpose:     To evaluate the coordinates and weighting factors
!               for an nth order Gaussian quadrature
!
!  parameters:
!
!     input:    ndec: number of decimal digits in real*16 arithmetic;
!                     usually equal to 33
!               n   : order of quadrature
!
!     output:   x   : coordinate values for quadrature
!               w   : weighting factors
!
!  real*16 scalars and arrays
        real*16 delta,der,pi,s,t,test,u,v,z
        real*16 x(n),w(n)
!
        test=10.0q0**(-ndec)
        imax=(n+1)/2
        pi=3.1415926535897932384626433832795028841971q0
          do 40 i=1,imax
	  z=qcos(pi*(i-0.25q0)/(n+0.5q0))
            do 20 j=1,30
            u=0.0q0
	    v=1.0q0
	      do 10 k=1,n
	      t=u
              u=v
	      v=((k+k-1)*z*u-(k-1)*t)/k
10   	      continue
            s=z*z-1.0q0
	    der=n*(z*v-u)/s
	    delta=-v/der-0.5q0*v*v*((n*n*s-n*z*z-n)*v+2.0q0*n*z*u)/&
        (der*der*der*s*s)
            z=z+delta
	    if(qabs(delta/z).lt.test) go to 30
20          continue
30        continue
	  x(i)=-z
	  x(n+1-i)=z
	  w(i)=2.0q0/((1.0q0-z*z)*der*der)
	  w(n+1-i)=w(i)
40	  continue
	return
	end
!
!
!
        subroutine pleg (m,lim,maxp,limcsav,iopd,ndec,barg,narg,maxt,pr,&
                   pdr,pdnorm,ipdnorm,pnorm,ipnorm,alpha,beta,&
                   gamma,coefa,coefb,coefc,coefd,coefe)
!
!  purpose:     To calculate ratios of successive associated Legendre
!               functions of the first kind for given arguments barg.
!               to calculate corresponding ratios of their first
!               derivatives. To calculate the characteristics and
!               exponents of both the Legendre functions of the first
!               kind and their first derivatives.
!
!  parameters:
!
!     input:    m      : m
!               lim    : two less than the number of associated Legendre
!                        function ratios calculated for given arguments
!               maxp   : dimension of alpha, beta, gamma, coefa, coefb,
!                        coefc, coefd, and coefe arrays and second
!                        dimension of pr and pdr arrays
!               limcsav: integer equal to the number of coefficients in
!                        each of the arrays alpha, beta, gamma, coefa,
!                        coefb, coefc, coefd, and coefe arrays that
!                        have already been calculated in earlier calls
!                        to pleg for this value of m and will not be
!                        calculated again. [Note that the minimum
!                        array index for the coefficients is 3 and
!                        the maximum array index is limcsav+2]
!               iopd   : integer that is set = 0 if derivatives of
!                        Legendre functions (i.e., their ratios)
!                        are not required when iopang = 1 and the
!                        first derivatives of the angular functions
!                        are not requested.
!                        iopd is set = 1 when iopang = 2 and pleg is
!                        also being used to obtain ratios of first
!                        derivatives of Legendre functions for use in
!                        computing the first derivatives of the angular
!                        functions.
!                        iopd is set = 2 when pleg is being used to
!                        compute ratios of Legendre functions for use in
!                        the calculation of the denominator term used in
!                        calculating the radial functions of the second
!                        kind and their first derivatives in r2eta.
!                        iopd is set = 3 when pleg is being used to
!                        compute ratios of both the Legendre functions
!                        and their first derivatives for use in the
!                        calculation of the numerator terms used
!                        in r2eta to calculate the radial functions of
!                        the second kind and their first deriatives.
!               ndec   : number of decimal digits in real*16 arithmetic
!               barg   : array of narg values of eta for which Legendre
!                        functions are to be calculated
!               narg   : number of specified values of eta in barg array
!               maxt   : dimension of barg array
!
!     output:   pr     : array of ratios of successive first kind
!                        associated Legendre functions of the same
!                        parity
!               pdr    : array of ratios of successive derivatives of
!                        first kind associated Legendre functions of
!                        the same parity
!               pdnorm : array of characteristics of the first
!                        derivatives of associated Legendre functions
!                        of the first kind of order m and degree m
!               ipdnorm: array of exponents of the first derivatives
!                        of associated Legendre functions of the first
!                        kind of order m and degree m
!               pnorm  : array of characteristics of the associated
!                        Legendre functions of the first kind of order m
!                        and degree m
!               ipnorm : array of exponents of the associated Legendre
!                        functions of the first kind of order m and
!                        degree m
!
!     input/output:
!               alpha  : array of coefficients in the recursion
!                        formula for the associated Legendre functions
!               beta   : array of coefficients in the recursion
!                        formula for the associated Legendre functions
!               gamma  : array of coefficients in the recursion
!                        formula for the associated Legendre functions
!               coefa  : array of coefficients in the expression
!                        relating the derivative ratios pdr to the
!                        function ratios pr
!               coefb  : array of coefficients in the expression
!                        relating the derivative ratios pdr to the
!                        function ratios pr
!               coefc  : array of coefficients in the expression
!                        relating the derivative ratios pdr to the
!                        function ratios pr
!               coefd  : array of coefficients in the expression
!                        relating the derivative ratios pdr to the
!                        function ratios pr
!               coefe  : array of coefficients in the expression
!                        relating the derivative ratios pdr to the
!                        function ratios pr
!
!  real*16 scalars and arrays
        real*16 adec,ajterm,am2p1,anden1,anden2,an2tnp1,bargs,den,rm,&
          temp1,temp2,temp3,term
        real*16 alpha(maxp),barg(maxt),beta(maxp),coefa(maxp),&
          coefb(maxp),coefc(maxp),coefd(maxp),coefe(maxp),&
          gamma(maxp),pdnorm(maxt),pdr(maxt,maxp),&
          pr(maxt,maxp),pnorm(maxt)
!
!  integer array
        dimension ipdnorm(maxt),ipnorm(maxt)
!
        adec=10.q0**(-ndec+2)
        rm=qfloat(m)
        am2p1=qfloat(m+m+1)
        m2=m+m
        m2m1=m2-1
        mm1=m-1
        mm2=m-2
        msqp1=2*m*m+1
!
!  calculate the coefficients alpha(j), beta(j), and gamma(j) for
!  the three term recursion relating the Legendre function ratios
!
!              m                m
!   pr(k,j) = p    (barg(k)) / p  (barg(k))
!              m+j-1            m+j-3
!
!  and calculate the coefficients coefa(j), coefb(j), coefc(j),
!  coefd(j), and coefe(j) in the expression used to calculate
!  ratios of Legendre function derivatives
!
!               m                 m
!   pdr(k,j) = p'    (barg(k)) / p'  (barg(k))
!               m+j-1             m+j-3
!
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
      	  an2tnp1=qfloat(2*n)*qfloat(n+1)
      	  anden1=qfloat(nmmm2)*qfloat(nmmm1)
      	  anden2=qfloat(n2m1)*anden1
      	  alpha(j)=qfloat(n2p3)*qfloat(n2p1)/anden1
      	  beta(j)=qfloat(n2p1)*(qfloat(msqp1)-an2tnp1)/anden2
      	  gamma(j)=-qfloat(n2p3)*qfloat(npm)*qfloat(npmm1)/anden2
          coefa(j)=-qfloat(npmp2)/qfloat(nmmm1)
          coefb(j)=qfloat(n2p3*(n+2))/anden1
          coefc(j)=-qfloat(npmp1*npmp2)/anden1
          coefd(j)=qfloat(npmp1)/qfloat(nmmm2)
          coefe(j)=-qfloat((n+1)*n2p3)/anden1
10        continue
        gamma(3)=0.q0
        gamma(4)=0.q0
        term=1.q0
        iterm=0
        if(m.lt.2) go to 30
          do jm=2,m
      	  term=qfloat((jm+jm-1))*term
      	  if(term.lt.(1.q5)) go to 20
      	  term=term*(1.q-5)
      	  iterm=iterm+5
20  	  continue
          end do
30      continue
!
!   calculate the ratios of successive Legendre functions of the same
!   parity using the three term recursion relationship
!
!   pr(k,j) = alpha(j)*barg(k)*barg(k) + beta(j) + gamma(j)/pr(k,j-2)
!
          do 130 k=1,narg
      	  pnorm(k)=term
          ipnorm(k)=iterm
          pdnorm(k)=term
          ipdnorm(k)=iterm
!
!   define the first two ratios equal to unity and (2m+1)*barg(k)
          pr(k,1)=1.q0
          pr(k,2)=am2p1*barg(k)
          jdelta=1
          if(qabs(barg(k)).lt.adec) jdelta=2
          bargs=barg(k)*barg(k)
            do 40 j=3,lim+2,jdelta
            pr(k,j)=alpha(j)*bargs+beta(j)+(gamma(j)/pr(k,j-2))
40          continue
!
!   calculate the corresponding ratios of first derviatives using the
!   relationship (except for eta equal to zero or unity, where special
!   expressions are used)
!
!              (coefa(j)+coefb(j)*barg(k)*barg(k))*pr(k,j)+coefc(j)
!   pdr(k,j) = ----------------------------------------------------
!                  pr(k,j)+coefd(j)+coefe(j)*barg(k)*barg(k)
!
          if(iopd.eq.0.or.iopd.eq.2) go to 110
          pdr(k,1)=1.q0
          pdr(k,2)=1.q0
          if(qabs(barg(k)).ge.adec) go to 50
          pdr(k,2)=am2p1
            do j=4,lim+2,2
            pdr(k,j)=-qfloat(m2m1+j)/qfloat(j-2)
            end do
          go to 130
50        if(qabs(qabs(barg(k))-1.q0).ge.adec) go to 70
          if(m.eq.0) go to 60
          if(m.ne.2) go to 120
          pdr(k,1)=-2.q0*barg(k)
          go to 80
60        temp1=1.q0
          temp2=3.q0
          pdr(k,2)=1.q0
          pdr(k,3)=3.q0*barg(k)
            do j=4,lim+2
            temp3=temp2+qfloat(j-1)
            pdr(k,j)=temp3/temp1
            temp1=temp2
            temp2=temp3
            end do
          go to 130
70        if(m.ne.0) go to 80
          pdr(k,1)=1.q0
          pdr(k,2)=1.q0
          pdr(k,3)=3.q0*barg(k)
          jlow=4
          go to 90
80        pdr(k,2)=am2p1*((rm+1.q0)*bargs-1.q0)/(rm*barg(k))
          jlow=3
90        continue
            do 100 j=jlow,lim+2
            den=(pr(k,j)+coefd(j)+coefe(j)*bargs)
            if(den.eq.0.0q0) den=1.0q-50
            pdr(k,j)=((coefa(j)+coefb(j)*bargs)*pr(k,j)+coefc(j))/den
100        continue
110       continue
          if(m.eq.0.or.iopd.eq.2.or.iopd.eq.3) go to 130
          if(qabs(qabs(barg(k))-1.q0).lt.adec) go to 120
          ajterm=rm*qlog10(1.q0-bargs)/2.q0
          jterm=iqint(ajterm)
          ipnorm(k)=ipnorm(k)+jterm
          pnorm(k)=pnorm(k)*(10.q0**(ajterm-qfloat(jterm)))
          if(iopd.eq.0) go to 130
          ajterm=qlog10(rm*qabs(barg(k)))+(rm-2.q0)*&
           qlog10(1.q0-bargs)/2.q0
          jterm=iqint(ajterm)
          ipdnorm(k)=ipdnorm(k)+jterm
          pdnorm(k)=-pdnorm(k)*(10.q0**(ajterm-qfloat(jterm)))
          if(barg(k).lt.0.q0) pdnorm(k)=-pdnorm(k)
          go to 130
120       pnorm(k)=0.q0
          ipnorm(k)=0
          if(m.ne.2) pdnorm(k)=0.q0
          if(m.ne.2) ipdnorm(k)=0
130       continue
        return
        end
!

!
!
        subroutine qleg (m,lnum,limq,maxq,x1,ndec,qdr,qdml,iqdml,qdl,&
                   iqdl,qr,qml,iqml,ql,iql,termpq,itermpq)
!
!  purpose:     to calculate ratios of successive associated Legendre
!               functions of the second kind for given c,x, and m.
!               to calculate corresponding ratios of their first
!               derivatives. to calculate the characteristics and
!               exponents of both the Legendre functions of the second
!               kind and their first derivatives.
!
!  parameters:
!
!     input:    m      : m
!               lnum   : number of l values desired (=lmax+1);
!                        also equal to the dimension of the arrays
!                        ql, iql, qdl, and iqdl
!               limq   : the number of associated Legendre function
!                        ratios calculated for given m,lnum,c,ndec,
!                        and x1
!               maxq   : dimension of qr and qdr arrays
!               x1     : x - 1.q0
!               ndec   : number of decimal digits in real*16 arithmetic
!
!     output:   qdr    : ratios of derivatives of successive
!                        associated Legendre functions of the second
!                        kind
!               qdml   : characteristic of the first derivative of
!                        the associated Legendre function of the second
!                        kind with order m and degree m-1
!               iqdml  : exponent corresponding to qdml
!               qdl    : characteristic of the first derivative of
!                        the associated Legendre function of the second
!                        kind with order m and degree m
!               iqdl   : exponent corresponding to qdl
!               qr     : ratios of successive associated Legendre
!                        functions of the second kind
!               qml    : characteristic of the associated Legendre
!                        function of the second kind with order m
!                        and degree m-1
!               iqml   : exponent corresponding to qml
!               ql     : characteristic of the associated Legendre
!                        function of the second kind with order m and
!                                                           -m/2
!                        degree m, scaled by (2m-1)!!(x*x-1)
!               iql    : exponent corresponding to ql
!               termpq : characteristic of the relative size of the
!                        maximum terms in the positive degree q series
!                        and the p series used to calculate r2 and r2d
!                        in subroutine r2leg
!               itermpq: exponent corresponding to termpq
!
!  real*16 scalars and arrays
        real*16 ajm,dec,qdml,qlow,qml,qupp,q00,q11,rin,rm,&
          term,termpq,tjm,tm,tmr,x,x1,x1d,xsqr
        real*16 qdl(lnum),qdr(maxq),ql(lnum),qr(maxq)
!
!  integer arrays
        dimension iqdl(lnum),iql(lnum)
!
        dec=10.q0**(-ndec)
        rm=qfloat(m)
        tm=rm+rm
        tmr=tm/(tm+1.q0)
        qflag=0
        x=x1+1.q0
        x1d=(x+1.q0)*x1
        xsqr=qsqrt(x1d)
        mxqrest=limq+ndec*iqint((1.q0-1.q0/qlog10(x-xsqr)))
        if(m.eq.0) mlimq=50000*ndec+limq
        if(m.eq.1) mlimq=12000*ndec+limq
        if(m.eq.2) mlimq=5000*ndec+limq
        if(m.eq.3) mlimq=600*ndec+limq
        if(m.ge.4) mlimq=100*ndec+limq
        if(m.eq.1.and.x1.lt.1.q-9) mlimq=50000*ndec+limq
        mxqr=min(mxqrest,mlimq)
        mxqrpm=mxqr+m
        write(40,5) mxqrpm
5       format(15x,'used backward recursion to calculate ratios of q',&
         ' functions starting at order',i8)
!
!                              m    m
!  calculate ratios qr(k+m) = q  / q    ,k=m+limq to k=m+1 using
!                              k    k-1
!
!                                         m       m               1/2
!  backward recursion from qr(maxm+2m) = q     / q      =x-(x*x-1)
!                                         mxqr+m  mxqr+m-1
!
        qupp=x-xsqr
          do jn=mxqr+m,limq+m+1,-1
          rin=qfloat(jn)
          qlow=(rin+rm-1.q0)/(x*(rin+rin-1.q0)&
         -(rin-rm)*qupp)
          qupp=qlow
          end do
        qr(limq+m+m)=qupp
          do 10 jn=limq+m,m+2,-1
          rin=qfloat(jn)
          qr(jn+m-1)=(rin+rm-1.q0)/(x*(rin+rin-1.q0)&
               -(rin-rm)*qr(jn+m))
10        continue
!
!                              m     m
!  calculate ratios qr(k+m) = q   / q   ,k=m-1 to k=-m+1 using
!                              k-1   k
!
!                                       m      m
!  backward recursion from qr(m-1+m) = q    / q     = x
!                                       m-2    m-1
!
20      if(m.eq.0) go to 100
        qr(m+m-1)=x
        if(m.eq.1) go to 40
          do 30 jn=m-1,2-m,-1
          rin=qfloat(jn)
          qr(jn+m-1)=(x*(rin+rin-1.q0)&
               -((rin-rm)/qr(jn+m)))/(rin+rm-1.q0)
30        continue
40      continue
!
!                  m
!  calculation of q    , m > 0 by forward division of qr ratios
!                  m-1
!
!                 m
!  starting with q  calculated from its closed form expression.
!                 0
!
        qml=rm*qlog10(x+1.q0)-qlog10(2.q0)
        term=rm*qlog10(x1/(x+1.q0))
        if(term.lt.-ndec) go to 50
        qml=qml+qlog10(1.q0-((x1/(x+1.q0))**m))
50      continue
        term=1.q0
        iterm=0
        if(m.lt.2) go to 70
          do jm=2,m
          ajm=qfloat(jm)
          term=term*(ajm-1.q0)/(ajm+ajm-1.q0)
          if(term.gt.dec) go to 60
          term=term/dec
          iterm=iterm-ndec
60        continue
          end do
70      term=qlog10(term)
        qml=qml+term+iterm
        iqml=iqint(qml)
        qml=10.q0**(qml-iqml)
        if(2*(m/2).ne.m) qml=-qml
        if(m.lt.2) go to 90
          do jm=1,m-1
          qml=qml/qr(jm+m)
          if(qabs(qml).gt.dec) go to 80
          qml=qml/dec
          iqml=iqml-ndec
80        continue
          end do
90      continue
        iterm=iqint(qlog10(qabs(qml)))
        qml=qml*10.q0**(-iterm)
        iqml=iqml+iterm
!
!                  m
!  calculation of q   by forward recursion in m starting with values
!                  m
!       0       1
!  for q   and q  obtained from their closed form expressions, scaled
!       0       1
!                    -m/2
!  by (2m-1)!!(x*x-1).
!
100     q00=0.5q0*qlog((x+1.q0)/x1)
        if(m.ne.0) go to 110
        ql(1)=q00
        go to 130
110     q11=x1d*q00-x
        if(m.ne.1) go to 120
        ql(1)=q11
        go to 130
120     qlow=q00
        qupp=q11
          do jm=1,m-1
          tjm=qfloat(jm+jm)/qfloat(jm+jm+1)
          ql(1)=(x1d-tjm)*qupp+tjm*x1d*qlow
          qlow=qupp
          qupp=ql(1)
          end do
130     iql(1)=iqint(qlog10(qabs(ql(1))))
        ql(1)=ql(1)*(10.q0**(-iql(1)))
!
!  calculation of ratios of the first derivatives of q with respect
!  to x, using the relationships:
!
!                  m    m      [kx]qr(k+m)-(k+m)
!     qdr(k+m) = q'  / q'   =  ----------------- , k=m+lim to k=m+1
!                  k    k-1    [(k-m)]qr(k+m)-kx
!
!                  m      m    [(k-m)]qr(k+m)-kx
!     qdr(k+m) = q'   / q'  =  ----------------- , k=m-1 to k=-m+1
!                  k-1    k    [kx]qr(k+m)-(k+m)
!
          do jm=m+1,m+limq
          ajm=qfloat(jm)
          qdr(jm+m)=(ajm*x*qr(jm+m)-(ajm+rm))/((ajm-rm)*qr(jm+m)-ajm*x)
          end do
!
        if(m.eq.0) go to 140
          do jm=1-m,m-1
          ajm=qfloat(jm)
          qdr(jm+m)=(ajm*x*qr(jm+m)-(ajm-rm))/((ajm+rm)*qr(jm+m)-ajm*x)
          end do
140     continue
!
!                   m         m                      m        m
!  calculation of q'    and q'  from the values for q    and q  .
!                   m-1       m                      m-1      m
!
        if(m.gt.0) go to 150
        qdl(1)=-1.q0/x1d
        iqdl(1)=0
        go to 160
150     qdml=-rm*x*qml/x1d
        iterm=iqint(qlog10(qabs(qdml)))
        qdml=qdml*(10.q0**(-iterm))
        iqdml=iqml+iterm
        qdl(1)=rm*(x*ql(1)-2.q0*qml*(10.q0**(iqml-iql(1))))/x1d
        iqdl(1)=iql(1)
160     continue
        m2m1=m+m-1
          do jl=2,lnum
          ql(jl)=ql(jl-1)*qr(m2m1+jl)
          iql(jl)=iql(jl-1)
          if(qabs(ql(jl)).gt.1.q-10) go to 170
          ql(jl)=ql(jl)*1.q10
          iql(jl)=iql(jl)-10
170       qdl(jl)=qdl(jl-1)*qdr(m2m1+jl)
          iqdl(jl)=iqdl(jl-1)
          if(qabs(qdl(jl)).gt.1.q-10) go to 180
          qdl(jl)=qdl(jl)*1.q10
          iqdl(jl)=iqdl(jl)-10
180       end do
        termpq=rm*qlog10(xsqr)
        itermpq=iqint(termpq)
        termpq=10.q0**(termpq-itermpq)
        return
        end
!
!
!
       subroutine pint (c,m,lnum,x1,limint,maxint,maxlp,maxmp,ndec,&
                  wg,xg,ngau,ngqs,rpint1,rpint2,pint1,pint2,&
                  pint3,pint4,norme,pnorm,ipnorm,coefme,coefmo)
!
!  purpose:     to calculate integrals of the product of associated
!               Legendre functions and kernels containing spherical
!               Neumann functions and a window function. four
!               different kernel functions are involved leading to
!               integrals of four different types. the integrals are
!               calculated using gaussian quadrature.
!
!  parameters:
!
!     input:    c      : complex c
!               m      : m
!               lnum   : number of l values desired
!               x1     : x - 1.q0
!               limint : number of integrals of each of the four types
!                        required
!               maxint : dimension of the integral arrays
!               maxlp  : dimension of characteristic and exponent
!                        arrays of integrals
!               maxmp  : dimension of the spherical Neumann function
!                        array
!               ndec   : number of decimal digits in real*16 arithmetic
!               wg     : gaussian quadrature weighting factors
!               xg     : gaussian quadrature arguments
!               ngau   : order of gaussian quadrature
!               ngqs   : number of gaussian quadrature steps the
!                        integrand is divided into
!
!     output:   rpint1 : array of ratios of successive integrals of
!                        the same parity of the first type (l-m even)
!                        or of the third type (l-m odd)
!               rpint2 : array of ratios of successive integrals of
!                        the same parity of the second type (l-m even)
!                        or of the fourth type (l-m odd)
!               pint1  : array of scaled values for the integrals of
!                        the first type
!               pint2  : array of scaled values for the integrals of
!                        the second type
!               pint3  : array of scaled values for the integrals of
!                        the third type
!               pint4  : array of scaled values for the integrals of
!                        the fourth type
!               norme  : array of scaling exponents for the spherical
!                        Neumann functions
!               pnorm  : array of characteristics for the scaling
!                        factors used for the associated Legendre
!                        functions
!               ipnorm : array of exponents for the scaling factors
!                        used for the associated Legendre functions
!               coefme : coefficient for the expression for r2 and r2d
!                        using the integration method (l-m even)
!               coefmo : coefficient for the expression for r2 and r2d
!                        using the integration method (l-m odd)
!
!  real*16 scalars and arrays and complex*32 scalars and arrays
        real*16 ak,amo2,an,argb,arn,bn,coef,coefll,coefme,&
          coefmo,coefo,dec,etai,etaim1,etais,etal,etau,etcoef1,&
          etcoef2,rm,rn,sargb,scal,step,step0,step1,step2,&
          twom,twomi,x1,x1sm1
        complex*32 c,arg,coef1,coef2,coef4,darg,sneu1,sneu2,&
             sneu3,stemp0
        real*16 alpha(maxint),beta(maxint),p(maxint),pnorm(maxlp),&
          wg(ngau),xg(ngau)
        complex*32 pint1(maxint),pint2(maxint),pint3(maxint),&
             pint4(maxint),rpint1(maxint),rpint2(maxint),&
             sneun(maxmp),ynormn(maxmp)
!
!  integer arrays
        dimension norme(maxmp),ipnorm(maxlp)
!
        x1sm1=x1*(x1+2.q0)
        dec=10.q0**(-ndec-1)
        amo2=0.5q0*m
        lim=limint
        step0=1.q0/ngqs
        coefme=(m*(x1+1.q0))/x1sm1
        coefmo=(m*(x1+1.q0)**2+x1sm1)/((x1+1.q0)*x1sm1)
        rm=qfloat(m)
        if(x1.ge.0.1q0) go to 10
        step1=1.q0/(ngqs*(1.q0-3.q0*qlog10(x1)*(1.0q0+3.0q0*&
              qlog10(rm+10.0q0))))
        if(step1.gt.step0) step1=step0
        go to 20
10      step1=step0
20      step2=(1.q0-step1)/(ngqs-1)
        write(40,30) ngau,step1
30      format(10x,'order of gauss quadrature =',i4,'. first step',&
         ' size =',f10.6,'.')
!
!  calculation of scaling factors for the associated Legendre functions
	pnorm(1)=1.0q0
	pnorm(2)=1.0q0
	ipnorm(1)=0
	ipnorm(2)=0
        if(m.eq.0) go to 50
          do 40 n=1,m
          an=qfloat(n+n)
          bn=qfloat(n+n+1)
          pnorm(1)=pnorm(1)*an
          pnorm(2)=pnorm(2)*bn
          iterm1=iqint(qlog10(pnorm(1)))
          iterm2=iqint(qlog10(pnorm(2)))
          pnorm(1)=pnorm(1)*10.0q0**(-iterm1)
          pnorm(2)=pnorm(2)*10.0q0**(-iterm2)
          ipnorm(1)=ipnorm(1)+iterm1
          ipnorm(2)=ipnorm(2)+iterm2
40	  continue
50	twom=qfloat(m+m)
        pnorm(3)=pnorm(1)*(twom+qfloat(2))/qfloat(2)
        iterm3=iqint(qlog10(pnorm(3)))
        pnorm(3)=pnorm(3)*10.q0**(-iterm3)
        ipnorm(3)=iterm3+ipnorm(1)
        if(lnum.lt.4) go to 70
          do 60 il=4,lnum,2
          pnorm(il)=pnorm(il-2)*(twom+il-1)/(il-1)
          pnorm(il+1)=pnorm(il-1)*(twom+il)/(il)
          iterm1=qlog10(pnorm(il))
          iterm2=qlog10(pnorm(il+1))
          ipnorm(il)=ipnorm(il-2)+iterm1
          ipnorm(il+1)=ipnorm(il-1)+iterm2
          pnorm(il)=pnorm(il)*10.0q0**(-iterm1)
          pnorm(il+1)=pnorm(il+1)*10.0q0**(-iterm2)
60	  continue
70	continue
!
!  calculationk of the coefficients in the recursion relation used
!  for the scaled associated Legendre functions
	alpha(1)=(twom+1.q0)*pnorm(1)/pnorm(2)
        alpha(1)=alpha(1)*10.0q0**(ipnorm(1)-ipnorm(2))
        beta(1)=0.0q0
        alpha(2)=(twom+3.0q0)*pnorm(2)/(pnorm(3)*2.0q0)
        alpha(2)=alpha(2)*10.0q0**(ipnorm(2)-ipnorm(3))
        beta(2)=-(twom+1.0q0)/(twom+2.0q0)
          do 80 il=3,lim+2
          alpha(il)=alpha(il-2)*(twom+il-1)*(twom+il+il-1)* &
    (il-2)/((il-1)*(twom+il)*(twom+il+il-5))
	  beta(il)=-(twom+il-1)/(twom+il)
80	  continue
!
          do 90 il=1,lim+2,2
          pint1(il)=(0.0q0,0.0q0)
          pint2(il)=(0.0q0,0.0q0)
          pint3(il+1)=(0.0q0,0.0q0)
          pint4(il+1)=(0.0q0,0.0q0)
90	  continue
!
!  calculation of the scaling exponents for the spherical Bessel
!  functions required for the four types of integrals
        twomi=1.q0
        if(m.eq.0) go to 110
          do 100 n=1,m
100	  twomi=twomi*(n+n-1)/(n+n)
110	continue
        arg=c*qsqrt(x1*(x1+2.0q0))
        darg=1.0q0/arg
        stemp0=-cqcos(arg)*darg
        norme(1)=qlog10(cqabs(stemp0))
        ynormn(1)=stemp0*10.0q0**(-norme(1))
        ynormn(2)=(stemp0-cqsin(arg))*darg
        norme(2)=norme(1)
        ynormn(2)=ynormn(2)*10.0q0**(-norme(1))
!
          do 120 n=3,m+3,2
          rn=qfloat(n+n-3)
          arn=qfloat(n+n-1)
          ynormn(n)=-ynormn(n-2)+darg*rn*ynormn(n-1)
          ynormn(n+1)=-ynormn(n-1)+darg*arn*ynormn(n)
          norme(n+1)=qlog10(cqabs(ynormn(n+1)))
          scal=10.0q0**(-norme(n+1))
          ynormn(n+1)=ynormn(n+1)*scal
          ynormn(n)=ynormn(n)*scal
          norme(n+1)=norme(n+1)+norme(n-1)
          norme(n)=norme(n+1)
120        continue
!
!  gaussian quadrature integration loops. first dividing integrand
!  into ngqs steps
          do 180 k=1,ngqs
          ak=qfloat(k)
          if(k.ne.1) go to 130
          etal=0.q0
          etau=step1
          step=step1
          go to 140
130       etal=step1+(ak-2.q0)*step2
          etau=etal+step2
          step=step2
140       etcoef1=(etau+etal)/2.q0
 	  etcoef2=(etau-etal)/2.q0
!
!  gaussian quadrature integration over each step
	    do 170 i=1,ngau
 	    etai=etcoef1+xg(i)*etcoef2
            etais=etai*etai
            etaim1=1.q0-etais
            argb=x1sm1+etais
            coef=step*((x1sm1*etaim1/argb)**amo2)
            if(coef.lt.dec) go to 190
            sargb=qsqrt(argb)
            coefo=(x1+1.q0)/sargb
            arg=c*sargb
            darg=1.0q0/arg
            stemp0=-cqcos(arg)/arg
            sneun(1)=stemp0
            sneun(1)=sneun(1)*10.0q0**(-norme(1))
            sneun(2)=(stemp0-cqsin(arg))*darg
            sneun(2)=sneun(2)*10.0q0**(-norme(2))
              do 150 n=3,m+3,2
              rn=qfloat(n+n-3)
              arn=qfloat(n+n-1)
              sneun(n)=-sneun(n-2)+darg*rn*sneun(n-1)
              sneun(n+1)=-sneun(n-1)+darg*arn*sneun(n)
              sneun(n+1)=sneun(n+1)*10.0q0**(-norme(n+1)+norme(n-1))
              sneun(n)=sneun(n)*10.0q0**(-norme(n)+norme(n-1))
150           continue
            sneu1=sneun(m+1)
            sneu2=sneun(m+2)*10.0q0**(norme(m+2)-norme(m+1))
            sneu3=sneun(m+3)*10.0q0**(norme(m+3)-norme(m+1))
 	    p(1)=twomi*(etaim1**amo2)
            p(2)=alpha(1)*etai*p(1)
            coef1=coef*sneu1*wg(i)
            coef2=coef*coefo*sneu2*wg(i)
            coef4=coef*coefo*coefo*etai*sneu3*wg(i)
            pint1(1)=pint1(1)+coef1*p(1)
            pint2(1)=pint2(1)+coef2*p(1)
            pint4(2)=pint4(2)+coef4*p(2)
!
              do 160 il=2,lim+1,2
              p(il+1)=alpha(il)*etai*p(il)+beta(il)*p(il-1)
              p(il+2)=alpha(il+1)*etai*p(il+1)+beta(il+1)*p(il)
              pint1(il+1)=pint1(il+1)+coef1*p(il+1)
              pint2(il+1)=pint2(il+1)+coef2*p(il+1)
              pint4(il+2)=pint4(il+2)+coef4*p(il+2)
160	      continue
170         continue
180	  continue
190	continue
          do 200 il=1,lim,2
          pint3(il+1)=(pint2(il+2)-beta(il+1)*pint2(il))&
                /alpha(il+1)
200	  continue
!
!  calculation of ratios of integrals for ease in compution of r2 and
!  r2d in subroutine r2int
        rpint1(1)=(0.0q0,0.0q0)
        rpint1(2)=(0.0q0,0.0q0)
        rpint2(1)=(0.0q0,0.0q0)
        rpint2(2)=(0.0q0,0.0q0)
          do 210 il=3,lim,2
          rpint1(il)=pint1(il)*(twom+il-1)/(pint1(il-2)*(il-1))
          rpint2(il)=pint2(il)*(twom+il-1)/(pint2(il-2)*(il-1))
          rpint1(il+1)=pint3(il+1)*(twom+il)/(pint3(il-1)*(il))
          rpint2(il+1)=pint4(il+1)*(twom+il)/(pint4(il-1)*(il))
210	  continue
        return
        end
!
!
!
        subroutine sphbes (c,x,limj,maxj,maxlp,sbesf,sbesdf,sbesn,&
                     ibese,sbesdr)
!
!  purpose:     to calculate ratios of successive first kind spherical
!               Bessel functions of the same parity for given c and x.
!               to calculate corresponding ratios of their first
!               derivatives. to calculate the characteristics and
!               exponents of both the Bessel functions of the first
!               kind and their first derivatives.
!
!  parameters:
!
!     input:    c      : complex c
!               x      : x
!               limj   : the number of spherical Bessel function
!                        ratios calculated for given lnum,c,ndec,
!                        and maximum m desired
!               maxj   : dimension of sbesf vector
!               maxlp  : the number of scale factors
!                        that are calculated
!
!     output:   sbesf  : ratios of successive first kind spherical
!                        Bessel functions of the same parity
!               sbesdf : ratios of first derivatives of successive
!                        first kind spherical Bessel functions of the
!                        same parity
!               sbesn  : characteristics for the spherical
!                        Bessel functions
!               ibese  : exponents for the spherical
!                        Bessel functions
!               sbesdr : ratios of first derivatives of spherical Bessel
!                        functions to the corresponding spherical
!                        spherical functions
!
!  real*16 scalars and complex*32 scalars and arrays
        real*16 cxr,rn,x
        complex*32 c,cx,f0,f1,r1,stemp0,stemp1
        complex*32 sbesdf(maxj),sbesdr(maxlp),sbesf(maxj),sbesn(maxlp)
!
!  integer array
        dimension ibese(maxlp)
!
        cx=c*x
        cxr=qreal(cx)
        lim1=2*iqint(cqabs(cx))+20
        f0=cqsin(cx)/cx
        f1=(f0-cqcos(cx))/cx
        r1=f1/f0
!
!  compute first kind Bessel function ratios
!        sbesf(k)= j(n=k,c*x)/j(n=k-1,c*x)
!        sbesn(k)= (j(n=k-1),c*x))*10.q0**(-ibese(k))
!
        if (cqabs(cx).lt.limj) go to 20
!
!  for c*x >= lim, use forward recursion to
!  get fcn. ratios:
!       j(n+1,c*x)/j(n,c*x)=(2*n+1)/(c*x)-1/(j(n,c*x)/j(n-1,c*x))
!
        sbesf(1)=r1
          do 10 n=1,limj-1
          rn=qfloat(n+n+1)
          sbesf(n+1)=(rn/cx)-(1.0q0/sbesf(n))
10        continue
        go to 60
20      continue
!
!  for c*x < lim, use backward recursion to
!  get fcn. ratios:
!       j(n,c*x)/j(n-1,c*x) = 1/( (2*n+1)/(c*x) - j(n+1,c*x)/j(n,c*x) )
!
        stemp0=cx/qfloat(2*lim1)
        if(lim1.le.limj) go to 40
          do 30 n=lim1,limj,-1
          rn=qfloat(n+n+1)
          stemp1=1.0q0/(rn/cx-stemp0)
          stemp0=stemp1
30        continue
40      sbesf(limj)=stemp0
          do 50 n=limj-1,1,-1
          rn=qfloat(n+n+1)
          sbesf(n)=1.0q0/(rn/cx-sbesf(n+1))
50        continue
60      continue
!
!  for all c*x, calculate the amplitude and sign scale
!  factors by forward operation on the Bessel function
!  ratios.
        ibese(1)=qlog10(cqabs(f0))
        sbesn(1)=f0*10.0q0**(-ibese(1))
        if(qabs(qsin(cxr)).lt.0.5q0.and.cxr.gt.1.0q0) go to 70
        sbesn(2)=sbesn(1)*sbesf(1)
        ibese(2)=qlog10(cqabs(sbesn(2)))
        sbesn(2)=sbesn(2)*10.0q0**(-ibese(2))
        ibese(2)=ibese(2)+ibese(1)
        go to 80
70      ibese(2)=qlog10(cqabs(f1))
        sbesn(2)=f1*(10.0q0**(-ibese(2)))
80      continue
          do 90 n=3,maxlp
          sbesn(n)=sbesn(n-1)*sbesf(n-1)
          ibese(n)=qlog10(cqabs(sbesn(n)))
          sbesn(n)=sbesn(n)*10.0q0**(-ibese(n))
          ibese(n)=ibese(n)+ibese(n-1)
90        continue
!
!  calculate the ratios of the first derivatives of successive
!  Bessel functions using corresponding function ratios
          do n=1,limj
          rn=qfloat(n-1)
          sbesdf(n)=(cx-(rn+2.0q0)*sbesf(n))/(rn-cx*sbesf(n))
          end do
!
!  calculate the ratios of the first derivative to the corresponding
!  spherical Bessel function
          do n=1,maxlp
          rn=qfloat(n-1)
          sbesdr(n)=(rn/cx)-sbesf(n)
          end do
!
!  calculate the ratios of successive functions and derivatives
!  of the same parity
          do n=limj,2,-1
          sbesf(n)=sbesf(n-1)*sbesf(n)
          sbesdf(n)=sbesdf(n-1)*sbesdf(n)
          end do
        return
        end
!
!
!
        subroutine sphneu (c,x,limn,maxn,maxlp,limbes,sneuf,sneun,&
                     ineue,sneudf,sneudr)
!
!  purpose:     to calculate ratios of spherical Neumann functions
!               and ratios of their first derivatives for given c and x.
!               to calculate the Neumann function characteristics
!               and exponents. to calculate ratios of the first
!               derivatives to corresponding functions.
!
!  parameters:
!
!     input:    c      : complex c
!               x      : x
!               limn   : the number of spherical Neumann function
!                        ratios calculated for given lnum,c,ndec,
!                        and maximum m desired
!               maxn   : dimension of sneuf and sneudf arrays
!               maxlp  : the number of values of scale factors
!                        that are calculated
!               limbes : dimension of spherical Bessel function ratios
!                        calculated for use is calculating Neumann
!                        functions
!
!     output:   sneuf  : ratios of successive spherical Neumann
!                        functions of the same parity
!               sneun  : characteristic for the spherical
!                        Neumann functions
!               ineue  : exponent for the spherical
!                        Neumann functions
!               sneudf : ratios of first derivatives of successive
!                        spherical Neumann functions of the same parity
!               sneudr : ratios of first derivatives of spherical
!                        Neumann functions to the corresponding
!                        function
!
!  real*16 and complex*32 scalars and arrays
        real*16 rn,rnn,test,x
        complex*32 c,cx,cx2,sbes1,sbes2,stemp0,stemp1,sbes
        complex*32 sneudf(maxn),sneudr(maxlp),sneuf(maxn),sneun(maxn)
        complex*32 sbesf(limbes)
!
!  integer arrays
        dimension ineue(maxn)
!
!  compute first kind ratios of Neumann functions and ratios
!  of their first derivatives
!
!        sneuf(k)=y(n=k,c*x)/y(n=k-2,c*x)
!        sneun(k)=(y(n=k-1),c*x)*10.q0**(-ineue(k))
!        sneudf(k)=y'(n=k,c*x)/y'(n=k-2,c*x)
!        sneudr(k)=(y'(n=k-1),c*x)/y(n=k-1),c*x))
!
!  calculate j ratios near break point by backward recursion
!
!       j(n,c*x)/j(n-1,c*x) = 1/( (2*n+1)/(c*x) - j(n+1,c*x)/j(n,c*x) )
!
        cx=c*x
        cx2=cx*cx
        test=1.0q0/cqabs(cx)
        limb=iqint(qabs(qreal(cx))+qabs(qimag(cx)))
        lim1=2*limb+20
        if(limb.lt.2) limb=2
        sbesf(lim1)=cx/(lim1+lim1+1)
          do 10 n=lim1-1,1,-1
          rn=qfloat(n+n+1)
        sbesf(n)=1.0q0/(rn/cx-sbesf(n+1))
10      continue
!
!  use relation with j's to compute y's from order zero
!  to the turning point. compute derivative ratios
!  at same time.

        stemp0=-cqcos(cx)/cx
        stemp1=(stemp0-cqsin(cx))/cx
        sneuf(1)=stemp1/stemp0
        sneun(1)=stemp0
        sneun(2)=stemp1
        sbes1=sbesf(1)*cqsin(cx)/cx
        sneudf(1)=-(cx-2.0q0*sneuf(1))/(cx*sneuf(1))
        if(limb.gt.limn) limb=limn
        j=1
          do 20 n=2,limb
          j=j+1
          rn=qfloat(n)
          rnn=qfloat(n-1)
          sbes2=sbesf(n)*sbes1
          sneun(n+1)=(sbes2*sneun(n)-1.0q0/cx2)/sbes1
          sneuf(n)=sneun(n+1)/sneun(n)
          sneudf(n)=(cx-(rn+1.0q0)*sneuf(n))/(rnn-cx*sneuf(n))
          if(cqabs(sbes1+sbes2).lt.test) go to 30
20        sbes1=sbes2
30      limb=max(j,2)
!
!  calculate characteristics and exponents for the Neumann functions
!  up to and including the turning point
        ineue(1)=iqint(qlog10(cqabs(stemp0)))
        sneun(1)=stemp0*10.0q0**(-ineue(1))
        ineue(2)=iqint(qlog10(cqabs(stemp1)))
        sneun(2)=stemp1*10.0q0**(-ineue(2))
          do 40 n=3,limb
          ineue(n)=iqint(qlog10(cqabs(sneun(n))))
40        sneun(n)=sneun(n)*10.0q0**(-ineue(n))
!
!  use forward recursion from break point to compute function ratios
!
!       y(n+1,c*x)/y(n,c*x)=(2*n+1)/(c*x)-1/(y(n,c*x)/y(n-1,c*x))
!
!  compute derivative ratios at same time using function ratios.
        if(limb.eq.limn) go to 55
          do 50 n=limb+1,limn
          rn=qfloat(n-1)
          rnn=qfloat(n+n-1)
          sneuf(n)=rnn/cx-1.0q0/sneuf(n-1)
50        sneudf(n)=(cx-(rn+2.0q0)*sneuf(n))/(rn-cx*sneuf(n))
55      continue
          sneuf(limn+1)=(0.0q0,0.0q0)
          sneuf(limn+2)=(0.0q0,0.0q0)
          sneudf(limn+1)=(0.0q0,0.0q0)
          sneudf(limn+2)=(0.0q0,0.0q0)
!
!  calculate the characteristics and exponents for Neumann
!  functions beyond the turning point by forward operation
!  on the Neumann function ratios:
        if(limb+1.gt.maxlp) go to 65
          do n=limb+1,maxlp
          sneun(n)=sneun(n-1)*sneuf(n-1)
          ineue(n)=iqint(qlog10(cqabs(sneun(n))))
          sneun(n)=sneun(n)*10.0q0**(-ineue(n))
          ineue(n)=ineue(n)+ineue(n-1)
          end do
65      continue
!
!  calculate the ratios of the first derivatives to the corresponding
!  spherical Neumann functions
          do n=1,maxlp
          rn=qfloat(n-1)
          sneudr(n)=(rn/cx)-sneuf(n)
          end do
!
!     calculate the ratios of successive functions and derivatives

!     of the same parity
          do n=limn,2,-1
          sneuf(n)=sneuf(n-1)*sneuf(n)
          sneudf(n)=sneudf(n-1)*sneudf(n)
          end do
        return
        end
        end module vb_prolate