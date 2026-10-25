c=======================================================================
c     charm_mod.f -- D-meson framework additions to AMPT v2.26t9b
c     (D-meson/src, forked 2026-10-17 from the md5-verified pristine
c     part2/src_stock).  Every switch is read from the environment and
c     defaults to OFF; with all switches unset the binary must reproduce
c     the stock binary byte for byte (gate G0, tests/gate_g0.sh).
c
c     AMPT_HQTRIG=1      HIJING heavy-quark trigger: IHPR2(3)=3, one hard
c                        c-cbar pair forced into a random primary NN
c                        collision of every event (hijing1.383 notes,
c                        7/1994).  pT threshold of the pair from
c                        AMPT_HQPTMIN (GeV, default 0 = inclusive, i.e.
c                        no p0 cut on the triggered pair, as in
c                        Zheng et al. PRC 101 034905).
c     AMPT_HQCLEAN=1     no charm from the untriggered minijets
c                        (MSUB 81/82 off there) and no g->ccbar in
c                        final-state showers (MSTJ(45)=3), so each event
c                        carries exactly the triggered pair.
c     AMPT_HQDIRECT=1    charm quarks enter ZPC with the momenta they had
c                        on the HIJING string BEFORE Lund fragmentation
c                        (Zhang et al. EPJC/arXiv:2403.06099 sec 2.1,
c                        Fan et al. arXiv:2601.12287), formation time
c                        E/mT^2.  The 3-momentum difference to the
c                        string-melting value is given to the light
c                        partons of the same string (energy weighted),
c                        so 3-momentum is conserved string by string.
c     AMPT_SIGHQ=r       sigma(heavy-any)/sigma(light-light) in ZPC
c                        (default 1 = stock).  Implemented in the
c                        closest-approach test d^2 < sigma/pi of all
c                        twelve isco* routines.
c     AMPT_CHARM_OUT=1   write ana/charm.dat (records E,I,L,F,G,H;
c                        format in docs/CHARM_DAT_FORMAT.md).
c     AMPT_LQFRAC=f      fraction of light (u,d) partons sampled into
c                        charm.dat records L/G (default 0.02).
c     AMPT_ECRIT=e       (src_v9) critical energy density in GeV/fm^3 of
c                        the local parton freeze-out of ISOFT=5 (Lin's
c                        "test local freezeout", linana.f LOCAL; stock
c                        value 1.0).  Only read by ISOFT=5.  With ISOFT=5
c                        a parton whose local transverse energy density
c                        falls below 2/3*e keeps the momentum it had then
c                        (later collisions do not change it) and is
c                        hadronised from that freeze-out point.
c                        HQFRZG closes a stock hole: a parton never
c                        frozen before the cascade ends would otherwise
c                        be hadronised with its INITIAL momentum.
c     AMPT_SOUT=1        (src_v11) also write EVERY s/sbar parton into the L
c                        and G records (diagnostic of strange-quark flow).
c     AMPT_HQSHAD=s      (src_v13; off when unset) scales HIJING's impact-
c                        parameter-dependent nuclear shadowing of the
c                        TRIGGERED hard process only (PYSTFU, R -> 1 +
c                        s (R - 1)); minijets and the bulk keep the full
c                        HIJING shadowing.  HIJING gives the gluon
c                        R(x=0.01) ~ 0.46 at the nucleus centre and ~0.6
c                        on average, global nPDF fits (EPPS21, nCTEQ15)
c                        ~0.8-0.9: s ~ 0.4.
c     AMPT_HQLATE=x      (src_v12, ISOFT=5 only; off when unset) charm
c                        decouples later than the bulk: a charm quark
c                        freezes only when the local density falls below
c                        x times the bulk criterion (x=0: never; it then
c                        keeps scattering until the cascade ends), and it
c                        still scatters off FROZEN partons, changing only
c                        its own momentum (no back-reaction on the frozen
c                        bulk, as in Langevin/LBT heavy-quark transport).
c                        A stand-in for the hadronic-phase D diffusion
c                        that ART lacks.  Not energy conserving by the
c                        recoil of those collisions (charm only).
c     AMPT_ZPC_TSTOP=t   (src_v9) hard stop of the whole parton cascade
c                        at lab time t (fm/c), Zhang & Xu PRC 96 044907
c                        sec II A; ported unchanged from the validated
c                        ampt_field/investigating_density_issue patch.
c=======================================================================
      BLOCK DATA HQBLK
      include 'hqctl.inc'
      data ihqtrg,ihqcln,ihqdir,ichout,idstar,ncolhj,ihqiso,ihqfrg,
     &     ihqmtm/0,0,0,0,0,0,0,0,1/
      data hqpr,hqpm,hqmtmx,hqeps,hqpcut/0.5,0.5,1e9,0.05,1e9/
      data hqfbar,hqtilt/-1.,0./
      data hqecr,hqtstp/-1.,-1./
      data ihqsout/0/
      data hqlate/-1./
      data hqshad,ihqtrx/-1.,0/
      data ihqpdf,ihqpst,ihqhrd,nhqcp/0,1,0,2/
      data hqmuf,hqcron,hqxsf,hqsab,hqvio/2.,-1.,1.5,0.,0.,0./
      data ihqdrf/0/
      data ihqshw,hqwsh/0,1./
      data ihqcon,ihqbd/0,0/
      data ihqhad,ihqhds,ihqhlc,nhqhc,hqsgm,hqsgb/0,0,0,0,0,0,7.,10./
      data ihqhrc,ihqhsd/0,24680/
      data hqptmn,sighqr,hqlqf,lqseed/0.,1d0,0.02,12345/
      END

c-----------------------------------------------------------------------
      SUBROUTINE HQINIT
c     read the switches; called from main.f before HIJSET
      include 'hqctl.inc'
      COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      character*32 cv
      SAVE
      call getenv('AMPT_HQTRIG',cv)
      if(cv(1:1).eq.'1') ihqtrg=1
      call getenv('AMPT_HQCLEAN',cv)
      if(cv(1:1).eq.'1') ihqcln=1
      call getenv('AMPT_HQDIRECT',cv)
      if(cv(1:1).eq.'1') ihqdir=1
      call getenv('AMPT_CHARM_OUT',cv)
      if(cv(1:1).eq.'1') ichout=1
      call getenv('AMPT_DSTAR1',cv)
      if(cv(1:1).eq.'1') idstar=1
      call getenv('AMPT_HQISO',cv)
      if(cv(1:1).eq.'1') ihqiso=1
      call getenv('AMPT_HQISO_PCUT',cv)
      if(len_trim(cv).gt.0) read(cv,*) hqpcut
      call getenv('AMPT_TILT',cv)
      if(len_trim(cv).gt.0) read(cv,*) hqtilt
      if(hqtilt.ne.0.) write(6,105) hqtilt
 105  format(' # charm_mod: AMPT_TILT=',f7.3,' fm: light partons shifted',
     &     ' by dx = TILT*eta_s (|eta_s|<=3), charm untilted')
      call getenv('AMPT_ECRIT',cv)
      if(len_trim(cv).gt.0) read(cv,*) hqecr
      if(hqecr.gt.0.) write(6,106) hqecr
 106  format(' # charm_mod: AMPT_ECRIT=',f8.4,' GeV/fm^3 (local parton',
     &     ' freeze-out, used by ISOFT=5 only)')
      call getenv('AMPT_SOUT',cv)
      if(cv(1:1).eq.'1') ihqsout=1
      call getenv('AMPT_HQLATE',cv)
      if(len_trim(cv).gt.0) read(cv,*) hqlate
      if(hqlate.ge.0.) write(6,108) hqlate
 108  format(' # charm_mod: AMPT_HQLATE=',f7.3,': charm freezes below',
     &     ' this x bulk criterion and scatters off frozen partons',
     &     ' (charm momentum only), ISOFT=5')
      call getenv('AMPT_HQSHAD',cv)
      if(len_trim(cv).gt.0) read(cv,*) hqshad
      if(hqshad.ge.0.) write(6,109) hqshad
 109  format(' # charm_mod: AMPT_HQSHAD=',f7.3,': HIJING nuclear',
     &     ' shadowing of the TRIGGERED hard process scaled,',
     &     ' R -> 1 + s (R - 1)')
      call getenv('AMPT_HQPDF',cv)
      if(len_trim(cv).gt.0) read(cv,*) ihqpdf
      call getenv('AMPT_HQPDFSET',cv)
      if(len_trim(cv).gt.0) read(cv,*) ihqpst
      call getenv('AMPT_HQMUF',cv)
      if(len_trim(cv).gt.0) read(cv,*) hqmuf
      call getenv('AMPT_HQXSF',cv)
      if(len_trim(cv).gt.0) read(cv,*) hqxsf
      if(ihqpdf.gt.0) then
c     CTEQ6.1M (LHAPDF set cteq61, member 0); LHAPDF_DATA_PATH must
c     point at the directory holding cteq61/
         call initpdfsetbyname('cteq61')
         call initpdf(0)
         write(6,110) ihqpdf,ihqpst,hqmuf,hqxsf
         if(hqshad.ge.0.) write(6,111)
      endif
 110  format(' # charm_mod: AMPT_HQPDF=',i1,' (1: CTEQ6.1M x EPS09s',
     &     ' NLO set',i3,', 2: CTEQ6.1M only) for the TRIGGERED hard',
     &     ' process; muF^2 =',f5.2,' x PYTHIA Q^2; max safety',f5.2)
 111  format(' # charm_mod: AMPT_HQSHAD ignored (superseded by HQPDF)')
      call getenv('AMPT_HQDREST',cv)
      if(cv(1:1).eq.'1') ihqdrf=1
      if(ihqdrf.eq.1) write(6,113)
 113  format(' # charm_mod: AMPT_HQDREST=1: coalescence distance in the',
     &     ' rest frame of the quarks (Fan et al. Eq. 1)')
      call getenv('AMPT_HQSHW',cv)
      if(cv(1:1).eq.'1') ihqshw=1
      if(ihqshw.eq.1) write(6,114)
 114  format(' # charm_mod: AMPT_HQSHW=1: nuclear modification of the',
     &     ' triggered pair as an event weight (E record field 18)')
      call getenv('AMPT_HQCONST',cv)
      if(cv(1:1).eq.'1') ihqcon=1
      if(ihqcon.eq.1) write(6,115)
 115  format(' # charm_mod: AMPT_HQCONST=1: sum m_Q of the m_inv',
     &     ' criterion with JETSET constituent masses PARF(101-104)',
     &     ' (Fan et al. Eq. 2)')
      call getenv('AMPT_HQBDIST',cv)
      if(cv(1:1).eq.'1') ihqbd=1
      if(ihqbd.eq.1) write(6,116)
      call getenv('AMPT_HQHAD',cv)
      if(cv(1:1).eq.'1') ihqhad=1
      call getenv('AMPT_HQHADM',cv)
      if(len_trim(cv).gt.0) read(cv,*) hqsgm
      call getenv('AMPT_HQHADB',cv)
      if(len_trim(cv).gt.0) read(cv,*) hqsgb
      call getenv('AMPT_HQHADDS',cv)
      if(cv(1:1).eq.'1') ihqhds=1
      call getenv('AMPT_HQHADLC',cv)
      if(cv(1:1).eq.'1') ihqhlc=1
      call getenv('AMPT_HQHADREC',cv)
      if(cv(1:1).eq.'1') ihqhrc=1
      if(ihqhad.eq.1) write(6,117) hqsgm,hqsgb,ihqhds,ihqhlc,ihqhrc
 117  format(' # charm_mod: AMPT_HQHAD=1: isotropic elastic charm-',
     &     'hadron rescattering in ART, sigma =',f6.2,' mb (mesons),',
     &     f6.2,' mb (baryons), He-Fries-Rapp PLB 701 (2011) 445;',
     &     ' Ds',i2,', charm baryons',i2,', partner recoil',i2,
     &     ' (0: probe approximation, bulk untouched)')
 116  format(' # charm_mod: AMPT_HQBDIST=1: baryon d = sqrt(2/3) x',
     &     ' |r_c - (r_q1+r_q2)/2| (equal-mass Jacobi lambda)')
      call getenv('AMPT_HQCRONIN',cv)
      if(len_trim(cv).gt.0) read(cv,*) hqcron
      if(hqcron.ge.0.) write(6,112) hqcron
 112  format(' # charm_mod: AMPT_HQCRONIN delta=',f6.3,': cc-bar pair',
     &     ' kT broadening, Zhang et al. PLB 846 138219 Eqs. 1-3')
      call getenv('AMPT_ZPC_TSTOP',cv)
      if(len_trim(cv).gt.0) read(cv,*) hqtstp
      if(hqtstp.gt.0.) write(6,107) hqtstp
 107  format(' # charm_mod: AMPT_ZPC_TSTOP=',f8.3,' fm/c: parton',
     &     ' cascade stopped at this lab time')
      call getenv('AMPT_HQMTMODE',cv)
      if(len_trim(cv).gt.0) read(cv,*) ihqmtm
      call getenv('AMPT_HQFBAR',cv)
      if(len_trim(cv).gt.0) read(cv,*) hqfbar
      call getenv('AMPT_HQFRAG',cv)
      if(cv(1:1).eq.'1') ihqfrg=1
      call getenv('AMPT_HQPR',cv)
      if(len_trim(cv).gt.0) read(cv,*) hqpr
      call getenv('AMPT_HQPM',cv)
      if(len_trim(cv).gt.0) read(cv,*) hqpm
      call getenv('AMPT_HQMTMAX',cv)
      if(len_trim(cv).gt.0) read(cv,*) hqmtmx
      call getenv('AMPT_HQEPS',cv)
      if(len_trim(cv).gt.0) read(cv,*) hqeps
      call getenv('AMPT_HQPTMIN',cv)
      if(len_trim(cv).gt.0) read(cv,*) hqptmn
      call getenv('AMPT_SIGHQ',cv)
      if(len_trim(cv).gt.0) read(cv,*) sighqr
      call getenv('AMPT_LQFRAC',cv)
      if(len_trim(cv).gt.0) read(cv,*) hqlqf
      if(ihqtrg.eq.1) then
         IHPR2(3)=3
         HIPR1(10)=hqptmn
      endif
      if(ihqcln.eq.1) then
         MSTJ(45)=3
      endif
      if(ichout.eq.1)
     &     OPEN (97, FILE = 'ana/charm.dat', STATUS = 'UNKNOWN')
      if(ihqiso.eq.1) write(6,103) hqpcut
 103  format(' # charm_mod: HQISO=1 only for pair CM momentum below',
     &     e10.3,' GeV')
c     fragmentation into charm baryons: PYTHIA diquark picture,
c     P(qq)/(P(q)+P(qq)) with PARJ(1) = P(qq)/P(q) unless AMPT_HQFBAR set
      if(hqfbar.lt.0.) hqfbar=PARJ(1)/(1.+PARJ(1))
      if(ihqfrg.eq.1) write(6,104) ihqmtm
 104  format(' # charm_mod: HQMTMODE=',i1,' (1: coalesce if mT<mT_max',
     &     ' OR criteria [Fan et al.]; 0: mT<mT_max AND criteria [dm2-6])')
      if(ihqfrg.eq.1) write(6,102) hqpr,hqpm,hqmtmx,hqeps,hqfbar
 102  format(' # charm_mod: HQFRAG=1 p_r=',f6.3,' fm  p_m=',f6.3,
     &     '  mT_max=',e10.3,' GeV  eps_c=',f6.3,' f_baryon=',f6.3)
      if(ihqtrg+ihqcln+ihqdir+ichout+idstar+ihqiso+ihqfrg.gt.0
     &     .or.sighqr.ne.1d0)
     &     write(6,101) ihqtrg,hqptmn,ihqcln,ihqdir,sighqr,ichout,hqlqf,
     &     idstar,ihqiso
 101  format(' # charm_mod: HQTRIG=',i1,' HQPTMIN=',f6.2,' HQCLEAN=',i1,
     &     ' HQDIRECT=',i1,' SIGHQ=',f7.3,' CHARM_OUT=',i1,
     &     ' LQFRAC=',f6.3,' DSTAR1=',i1,' HQISO=',i1)
      RETURN
      END

c-----------------------------------------------------------------------
      double precision FUNCTION HQFAC(i,j)
c     ZPC: factor on cutof2 (=sigma/pi) for the pair (i,j) in /prec2/
      implicit double precision (a-h,o-z)
      parameter (MAXPTN=400001)
      common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &     px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &     xmass(MAXPTN), ityp(MAXPTN)
      include 'hqctl.inc'
      HQFAC=1d0
      if(sighqr.eq.1d0) return
      if(iabs(ityp(i)).eq.4.or.iabs(ityp(j)).eq.4) HQFAC=sighqr
      return
      end

c-----------------------------------------------------------------------
      SUBROUTINE HQCOLL(i,j)
c     ZPC: count 2-body scatterings of charm partons (/prec2/ index)
      parameter (MAXPTN=400001)
      double precision gx,gy,gz,ft,px,py,pz,e,xmass
      common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &     px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &     xmass(MAXPTN), ityp(MAXPTN)
      common /hqcnt/ ncolhq(MAXPTN)
      double precision gxfrz,gyfrz,gzfrz,ftfrz,pxfrz,pyfrz,pzfrz,efrz,
     &     xmfrz,tfrz
      common /frzprc/
     &     gxfrz(MAXPTN), gyfrz(MAXPTN), gzfrz(MAXPTN), ftfrz(MAXPTN),
     &     pxfrz(MAXPTN), pyfrz(MAXPTN), pzfrz(MAXPTN), efrz(MAXPTN),
     &     xmfrz(MAXPTN),
     &     tfrz(302), ifrz(MAXPTN), idfrz(MAXPTN), itlast
      common/anim/nevent,isoft,isflag,izpc
      include 'hqctl.inc'
c     ISOFT=5: a collision involving a frozen parton changes no momentum
c     (zpc.f newmom), so it is not counted; with AMPT_HQLATE an unfrozen
c     charm quark still scatters off a frozen parton (newmom), counted
      if(isoft.eq.5) then
         if(ifrz(i).eq.1.or.ifrz(j).eq.1) then
            if(hqlate.lt.0.) return
            if(.not.((ifrz(i).ne.1.and.iabs(ityp(i)).eq.4).or.
     &           (ifrz(j).ne.1.and.iabs(ityp(j)).eq.4))) return
         endif
      endif
      if(iabs(ityp(i)).eq.4) ncolhq(i)=ncolhq(i)+1
      if(iabs(ityp(j)).eq.4) ncolhq(j)=ncolhq(j)+1
      return
      end

c-----------------------------------------------------------------------
      SUBROUTINE HQFRZG
c     ISOFT=5, end of ZPC (zpstrg, before /frzprc/ replaces /prec2/):
c     a parton that never fell below the critical density before the
c     last collision still carries its INITIAL values in /frzprc/
c     (zpc.f inizpc).  Give it its final values instead, i.e. freeze it
c     at its last interaction, as ISOFT=4 does.  Counts kept for the Z
c     record (nfrz0 = such partons, nfrzt = all partons).
      parameter (MAXPTN=400001)
      double precision gx,gy,gz,ft,px,py,pz,e,xmass
      common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &     px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &     xmass(MAXPTN), ityp(MAXPTN)
      double precision gxfrz,gyfrz,gzfrz,ftfrz,pxfrz,pyfrz,pzfrz,efrz,
     &     xmfrz,tfrz
      common /frzprc/
     &     gxfrz(MAXPTN), gyfrz(MAXPTN), gzfrz(MAXPTN), ftfrz(MAXPTN),
     &     pxfrz(MAXPTN), pyfrz(MAXPTN), pzfrz(MAXPTN), efrz(MAXPTN),
     &     xmfrz(MAXPTN),
     &     tfrz(302), ifrz(MAXPTN), idfrz(MAXPTN), itlast
      common /hqlfz/ nfrz0, nfrzt
      COMMON /PARA1/ MUL
      SAVE
      nfrz0=0
      nfrzt=mul
      do 10 i=1,mul
         if(ifrz(i).eq.1) goto 10
         nfrz0=nfrz0+1
         idfrz(i)=ityp(i)
         gxfrz(i)=gx(i)
         gyfrz(i)=gy(i)
         gzfrz(i)=gz(i)
         ftfrz(i)=ft(i)
         pxfrz(i)=px(i)
         pyfrz(i)=py(i)
         pzfrz(i)=pz(i)
         efrz(i)=e(i)
         xmfrz(i)=xmass(i)
 10   continue
      return
      end

c-----------------------------------------------------------------------
      SUBROUTINE HQSAVE(ISG)
c     HIJING, string-melting branch.  ISG=0: reset for this event and
c     count charm on projectile/target strings (expected 0).  ISG>0:
c     save the charm partons of independent string ISG before HIJFRG.
      PARAMETER (MAXSTR=150001, MAXHQS=200)
      COMMON/HJJET1/NPJ(300),KFPJ(300,500),PJPX(300,500),
     &     PJPY(300,500),PJPZ(300,500),PJPE(300,500),
     &     PJPM(300,500),NTJ(300),KFTJ(300,500),
     &     PJTX(300,500),PJTY(300,500),PJTZ(300,500),
     &     PJTE(300,500),PJTM(300,500)
      COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &     K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &     PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)
      COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      common /hqsav/ nhqs,nhqj1,nhqovf,hqisg(MAXHQS),hqkf(MAXHQS),
     &     hqp(MAXHQS,5)
      integer hqisg,hqkf
      SAVE
      if(ISG.eq.0) then
         nhqs=0
         nhqovf=0
         nhqj1=0
c     charm on projectile/target strings (g->ccbar in showers of jets
c     attached to nucleon strings); string id as istrg0 in hijing:
c     10000+jp (projectile), 20000+jt (target)
         do 10 jp=1,IHNT2(1)
            do 10 l=1,NPJ(jp)
               if(iabs(KFPJ(jp,l)).ne.4) goto 10
               nhqj1=nhqj1+1
               call hqput(10000+jp,KFPJ(jp,l),PJPX(jp,l),PJPY(jp,l),
     &              PJPZ(jp,l),PJPE(jp,l),PJPM(jp,l))
 10      continue
         do 20 jt=1,IHNT2(3)
            do 20 l=1,NTJ(jt)
               if(iabs(KFTJ(jt,l)).ne.4) goto 20
               nhqj1=nhqj1+1
               call hqput(20000+jt,KFTJ(jt,l),PJTX(jt,l),PJTY(jt,l),
     &              PJTZ(jt,l),PJTE(jt,l),PJTM(jt,l))
 20      continue
         return
      endif
      do 30 l=1,NJSG(ISG)
         if(iabs(K2SG(ISG,l)).ne.4) goto 30
         call hqput(ISG,K2SG(ISG,l),PXSG(ISG,l),PYSG(ISG,l),
     &        PZSG(ISG,l),PESG(ISG,l),PMSG(ISG,l))
 30   continue
      return
      end

c-----------------------------------------------------------------------
      SUBROUTINE HQPUT(isg,kf,ax,ay,az,ae,am)
      PARAMETER (MAXHQS=200)
      common /hqsav/ nhqs,nhqj1,nhqovf,hqisg(MAXHQS),hqkf(MAXHQS),
     &     hqp(MAXHQS,5)
      integer hqisg,hqkf
      if(nhqs.ge.MAXHQS) then
         nhqovf=nhqovf+1
         return
      endif
      nhqs=nhqs+1
      hqisg(nhqs)=isg
      hqkf(nhqs)=kf
      hqp(nhqs,1)=ax
      hqp(nhqs,2)=ay
      hqp(nhqs,3)=az
      hqp(nhqs,4)=ae
      hqp(nhqs,5)=am
      return
      end

c-----------------------------------------------------------------------
      SUBROUTINE HQHTOP(npar)
c     end of HTOP: (a) AMPT_HQDIRECT replacement, (b) reset ZPC charm
c     collision counters, (c) charm.dat records E, I, L.
      PARAMETER (MAXSTR=150001, MAXPTN=400001, MAXHQS=200)
      DOUBLE PRECISION  GX0, GY0, GZ0, FT0, PX0, PY0, PZ0, E0, XMASS0
      DOUBLE PRECISION vxp0,vyp0,vzp0,xstrg0,ystrg0,xstrg,ystrg
      DOUBLE PRECISION dlx(MAXHQS),dly(MAXHQS),dlz(MAXHQS),esum,de,
     &     ebef,eaft,w,pp2,bestd,dd
      COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
     &     PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
     &     XMASS0(MAXPTN), ITYP0(MAXPTN)
      COMMON /ilist7/ LSTRG0(MAXPTN), LPART0(MAXPTN)
      common /precpa/vxp0(MAXPTN),vyp0(MAXPTN),vzp0(MAXPTN),
     1     xstrg0(MAXPTN),ystrg0(MAXPTN),
     2     xstrg(MAXPTN),ystrg(MAXPTN),istrg0(MAXPTN),istrg(MAXPTN)
      COMMON/HMAIN1/EATT,JATT,NATT,NT,NP,N0,N01,N10,N11
      COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      COMMON /AREVT/ IAEVT, IARUN, MISS
      include 'hqctl.inc'
      common /hqsav/ nhqs,nhqj1,nhqovf,hqisg(MAXHQS),hqkf(MAXHQS),
     &     hqp(MAXHQS,5)
      integer hqisg,hqkf
      common /hqcnt/ ncolhq(MAXPTN)
      common /hqdia/ nmatch,nfail,dethq
      double precision dethq
      dimension iused(MAXHQS)
      SAVE
c
      do 5 i=1,npar
         ncolhq(i)=0
 5    continue
      nmatch=0
      nfail=0
      dethq=0d0
      nc0=0
      nnan=0
      do 6 i=1,npar
         if(iabs(ITYP0(i)).eq.4) nc0=nc0+1
c     non-finite partons from stock string melting (seen once in 1e6
c     events: a NaN-momentum hadron on a HIJING projectile string);
c     counted in the E record so the analysis can drop the event
         if(.not.(E0(i).eq.E0(i)).or.abs(E0(i)).gt.1d30) nnan=nnan+1
 6    continue
c
      if(ihqdir.eq.1.and.nhqs.gt.0) then
         do 7 k=1,nhqs
            iused(k)=0
            dlx(k)=0d0
            dly(k)=0d0
            dlz(k)=0d0
 7       continue
c     (a1) each melted charm parton takes the momentum of the unused
c     saved parton of the same flavour on the same string that is
c     closest in momentum
         ebef=0d0
         eaft=0d0
         do 20 i=1,npar
            ebef=ebef+E0(i)
            if(iabs(ITYP0(i)).ne.4) goto 20
            isg=istrg0(LSTRG0(i))
            kbest=0
            bestd=1d30
            do 10 k=1,nhqs
               if(iused(k).ne.0.or.hqkf(k).ne.ITYP0(i)
     &              .or.hqisg(k).ne.isg) goto 10
               dd=(PX0(i)-hqp(k,1))**2+(PY0(i)-hqp(k,2))**2
     &              +(PZ0(i)-hqp(k,3))**2
               if(dd.lt.bestd) then
                  bestd=dd
                  kbest=k
               endif
 10         continue
            if(kbest.eq.0) then
               nfail=nfail+1
               goto 20
            endif
            iused(kbest)=1
            nmatch=nmatch+1
            dlx(kbest)=PX0(i)-dble(hqp(kbest,1))
            dly(kbest)=PY0(i)-dble(hqp(kbest,2))
            dlz(kbest)=PZ0(i)-dble(hqp(kbest,3))
            PX0(i)=dble(hqp(kbest,1))
            PY0(i)=dble(hqp(kbest,2))
            PZ0(i)=dble(hqp(kbest,3))
            E0(i)=dsqrt(PX0(i)**2+PY0(i)**2+PZ0(i)**2+XMASS0(i)**2)
c     formation time E/mT^2 of the charm quark itself (hbar c = 0.197)
            FT0(i)=0.197d0*E0(i)/(PX0(i)**2+PY0(i)**2+XMASS0(i)**2)
 20      continue
c     (a2) give each string's momentum difference to its light partons
         do 40 k=1,nhqs
            if(iused(k).eq.0) goto 40
            if(dlx(k).eq.0d0.and.dly(k).eq.0d0.and.dlz(k).eq.0d0)
     &           goto 40
            isg=hqisg(k)
            esum=0d0
            do 30 i=1,npar
               if(iabs(ITYP0(i)).eq.4) goto 30
               if(istrg0(LSTRG0(i)).ne.isg) goto 30
               esum=esum+E0(i)
 30         continue
            if(esum.le.0d0) then
               nfail=nfail+1
               goto 40
            endif
            do 35 i=1,npar
               if(iabs(ITYP0(i)).eq.4) goto 35
               if(istrg0(LSTRG0(i)).ne.isg) goto 35
               w=E0(i)/esum
               PX0(i)=PX0(i)+w*dlx(k)
               PY0(i)=PY0(i)+w*dly(k)
               PZ0(i)=PZ0(i)+w*dlz(k)
               E0(i)=dsqrt(PX0(i)**2+PY0(i)**2+PZ0(i)**2+XMASS0(i)**2)
 35         continue
 40      continue
         do 45 i=1,npar
            eaft=eaft+E0(i)
 45      continue
         dethq=eaft-ebef
         do 46 k=1,nhqs
            if(iused(k).eq.0) nfail=nfail+1
 46      continue
      endif
c
c     src_v14: Cronin broadening of the pair; PYTHIA maximum monitor
      call hqckick(npar)
      call hqvmax
c
c     AMPT_TILT (src_v8, sensitivity study, not a validated physics
c     model): tilt the light-parton source in the x-eta_s plane, Bozek &
c     Wyskiel / Chatterjee & Bozek style.  Partons start at z=0, t=0, so
c     eta_s at formation = momentum rapidity.  Charm (binary-collision
c     profile) is left untilted, as in Chatterjee & Bozek PRL 120 192301.
      if(hqtilt.ne.0.) then
         do 70 i=1,npar
            if(iabs(ITYP0(i)).eq.4) goto 70
            if(E0(i).le.dabs(PZ0(i))) goto 70
            yt=0.5*sngl(dlog((E0(i)+PZ0(i))/(E0(i)-PZ0(i))))
            yt=max(-3.,min(3.,yt))
            GX0(i)=GX0(i)+dble(hqtilt*yt)
 70      continue
      endif
c
      if(ichout.ne.1) return
      if(ihqshw.eq.1) then
         write(97,204) 'E',IAEVT,HINT1(19),NP,NT,NATT,npar,nc0,nhqs,
     &        nhqj1,nhqovf,nmatch,nfail,dethq,HINT1(61),ncolhj,nnan,hqwsh
         write(6,205) IAEVT,khqw,hqwx,hqwq2,hqsab,hqwr,hqwsh
 205     format(' # hqshw ev',i6,' kf',2i4,' x',2e13.5,' Q2',f9.3,
     &        ' sA,sB',2f7.3,' r',2f9.5,' w',f9.6)
      else
         write(97,201) 'E',IAEVT,HINT1(19),NP,NT,NATT,npar,nc0,nhqs,
     &        nhqj1,nhqovf,nmatch,nfail,dethq,HINT1(61),ncolhj,nnan
      endif
      do 60 i=1,npar
         ia=iabs(ITYP0(i))
         if(ia.eq.4) then
            write(97,202) 'I',IAEVT,ITYP0(i),PX0(i),PY0(i),PZ0(i),
     &           XMASS0(i),GX0(i),GY0(i),GZ0(i),FT0(i),istrg0(LSTRG0(i))
         elseif(ia.eq.1.or.ia.eq.2.or.(ia.eq.3.and.ihqsout.eq.1)) then
            if(ia.eq.3.or.hqrand().lt.hqlqf)
     &           write(97,202) 'L',IAEVT,ITYP0(i),PX0(i),PY0(i),PZ0(i),
     &           XMASS0(i),GX0(i),GY0(i),GZ0(i),FT0(i),istrg0(LSTRG0(i))
         endif
 60   continue
 201  format(a1,1x,i7,1x,f8.4,4(1x,i7),6(1x,i4),1x,e12.4,1x,e12.5,
     &     1x,i6,1x,i4)
 204  format(a1,1x,i7,1x,f8.4,4(1x,i7),6(1x,i4),1x,e12.4,1x,e12.5,
     &     1x,i6,1x,i4,1x,f9.6)
 202  format(a1,1x,i7,1x,i6,4(1x,f11.5),4(1x,e12.5),1x,i7)
      return
      end

c-----------------------------------------------------------------------
      real FUNCTION hqrand()
c     private LCG for output subsampling; never touches RANART/NSEED so
c     the physics random sequence is unchanged by AMPT_CHARM_OUT
      include 'hqctl.inc'
      integer*8 ls
      ls=mod(int(lqseed,8)*16807_8,2147483647_8)
      if(ls.le.0) ls=ls+2147483647_8
      lqseed=int(ls)
      hqrand=real(lqseed)/2147483647.
      return
      end

c-----------------------------------------------------------------------
      SUBROUTINE HQOUTF
c     start of PTOH (ZPC finished): charm partons at freeze-out (F) with
c     their ZPC scattering count, and a light-parton sample (G)
      parameter (MAXPTN=400001)
      double precision gx,gy,gz,ft,px,py,pz,e,xmass
      common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &     px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &     xmass(MAXPTN), ityp(MAXPTN)
      common /hqcnt/ ncolhq(MAXPTN)
      COMMON /PARA1/ MUL
      COMMON /AREVT/ IAEVT, IARUN, MISS
      include 'hqctl.inc'
      common /hqlfz/ nfrz0, nfrzt
      common/anim/nevent,isoft,isflag,izpc
      if(ichout.ne.1) return
c     Z: parton-phase end diagnostics (src_v9): ISOFT, partons, partons
c     frozen only by HQFRZG, eps_c (ISOFT=5), lab stop time
      if(isoft.eq.5.or.hqtstp.gt.0.) then
         if(isoft.ne.5) then
            nfrz0=0
            nfrzt=mul
         endif
         write(97,203) 'Z',IAEVT,isoft,nfrzt,nfrz0,hqecr,hqtstp
      endif
 203  format(a1,1x,i7,1x,i2,1x,i7,1x,i7,1x,f9.4,1x,f9.3)
      do 10 i=1,mul
         ia=iabs(ityp(i))
         if(ia.eq.4) then
            write(97,202) 'F',IAEVT,ityp(i),px(i),py(i),pz(i),
     &           xmass(i),gx(i),gy(i),gz(i),ft(i),ncolhq(i)
         elseif(ia.eq.1.or.ia.eq.2.or.(ia.eq.3.and.ihqsout.eq.1)) then
            if(ia.eq.3.or.hqrand().lt.hqlqf)
     &           write(97,202) 'G',IAEVT,ityp(i),px(i),py(i),pz(i),
     &           xmass(i),gx(i),gy(i),gz(i),ft(i),0
         endif
 10   continue
 202  format(a1,1x,i7,1x,i6,4(1x,f11.5),4(1x,e12.5),1x,i7)
      return
      end

c-----------------------------------------------------------------------
      SUBROUTINE HQOUTH
c     end of PTOH: charm hadrons right after coalescence (H) with the
c     charm constituent momentum and the largest constituent separation
c     at the latest constituent freeze-out time (lab frame, fm)
      PARAMETER (MAXSTR=150001)
      DOUBLE PRECISION  PXSGS,PYSGS,PZSGS,PESGS,PMSGS,
     1     GXSGS,GYSGS,GZSGS,FTSGS
      double precision tmx,xq(3),yq(3),zq(3),dmx,dij,pcx,pcy,pcz
      COMMON/HMAIN1/EATT,JATT,NATT,NT,NP,N0,N01,N10,N11
      COMMON/HMAIN2/KATT(MAXSTR,4),PATT(MAXSTR,4)
      COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &     K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &     PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)
      COMMON/SOFT/PXSGS(MAXSTR,3),PYSGS(MAXSTR,3),PZSGS(MAXSTR,3),
     &     PESGS(MAXSTR,3),PMSGS(MAXSTR,3),GXSGS(MAXSTR,3),
     &     GYSGS(MAXSTR,3),GZSGS(MAXSTR,3),FTSGS(MAXSTR,3),
     &     K1SGS(MAXSTR,3),K2SGS(MAXSTR,3),NJSGS(MAXSTR)
      COMMON /ARPRC/ ITYPAR(MAXSTR),
     &     GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     &     PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     &     XMAR(MAXSTR)
      COMMON /AREVT/ IAEVT, IARUN, MISS
      common /hqmod/ ihmode(MAXSTR)
      include 'hqctl.inc'
      SAVE
      if(ichout.ne.1) return
      n=0
      do 50 isg=1,NSG
         if(NJSGS(isg).eq.0) goto 50
         n=n+1
         nq=0
         ic=0
         do 10 l=1,3
            if(K2SGS(isg,l).ne.0) nq=nq+1
            if(iabs(K2SGS(isg,l)).eq.4.and.ic.eq.0) ic=l
 10      continue
         if(ic.eq.0) goto 50
         tmx=FTSGS(isg,1)
         do 20 l=2,nq
            tmx=max(tmx,FTSGS(isg,l))
 20      continue
         do 30 l=1,nq
            xq(l)=GXSGS(isg,l)+PXSGS(isg,l)/PESGS(isg,l)*
     &           (tmx-FTSGS(isg,l))
            yq(l)=GYSGS(isg,l)+PYSGS(isg,l)/PESGS(isg,l)*
     &           (tmx-FTSGS(isg,l))
            zq(l)=GZSGS(isg,l)+PZSGS(isg,l)/PESGS(isg,l)*
     &           (tmx-FTSGS(isg,l))
 30      continue
         dmx=0d0
         do 40 l=1,nq-1
            do 40 m=l+1,nq
               dij=dsqrt((xq(l)-xq(m))**2+(yq(l)-yq(m))**2
     &              +(zq(l)-zq(m))**2)
               dmx=max(dmx,dij)
 40      continue
         pcx=PXSGS(isg,ic)
         pcy=PYSGS(isg,ic)
         pcz=PZSGS(isg,ic)
         write(97,203) 'H',IAEVT,KATT(n,1),PATT(n,1),PATT(n,2),
     &        PATT(n,3),XMAR(n),GXAR(n),GYAR(n),GZAR(n),FTAR(n),
     &        pcx,pcy,pcz,dmx,nq,ihmode(n)
 50   continue
 203  format(a1,1x,i7,1x,i6,4(1x,f11.5),4(1x,e12.5),3(1x,f11.5),
     &     1x,f9.4,1x,i2,1x,i1)
      return
      end

c-----------------------------------------------------------------------
      SUBROUTINE HQNCOL(n)
c     HIJING: store the number of primary NN collisions of the event
      include 'hqctl.inc'
      ncolhj=n
      return
      end

c-----------------------------------------------------------------------
      SUBROUTINE HQDSTR
c     AMPT_DSTAR1 (PTOH, right after the flavour/species loop): set the
c     primordial vector/pseudoscalar ratio of open-charm mesons to 1, as
c     the improved AMPT does (Zheng, Zhang, Shi, Lin, PRC 101 034905,
c     sec. II): the charm-light q-qbar candidates of the event are ordered
c     by excess mass m_inv - m1 - m2; the lower half become D (Ds), the
c     upper half D* (Ds*).  For an odd count the middle one is a D* with
c     probability 1/2 (private LCG, physics random sequence untouched).
      PARAMETER (MAXSTR=150001, MAXC=400)
      DOUBLE PRECISION  PXSGS,PYSGS,PZSGS,PESGS,PMSGS,
     1     GXSGS,GYSGS,GZSGS,FTSGS
      double precision ex(MAXC), xm
      COMMON/HMAIN2/KATT(MAXSTR,4),PATT(MAXSTR,4)
      COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &     K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &     PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)
      COMMON/SOFT/PXSGS(MAXSTR,3),PYSGS(MAXSTR,3),PZSGS(MAXSTR,3),
     &     PESGS(MAXSTR,3),PMSGS(MAXSTR,3),GXSGS(MAXSTR,3),
     &     GYSGS(MAXSTR,3),GZSGS(MAXSTR,3),FTSGS(MAXSTR,3),
     &     K1SGS(MAXSTR,3),K2SGS(MAXSTR,3),NJSGS(MAXSTR)
      COMMON /ARPRC/ ITYPAR(MAXSTR),
     &     GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     &     PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     &     XMAR(MAXSTR)
      include 'hqctl.inc'
      dimension ih(MAXC), iord(MAXC)
      SAVE
      if(idstar.ne.1) return
      n=0
      nh=0
      do 10 isg=1,NSG
         if(NJSGS(isg).eq.0) goto 10
         n=n+1
         if(NJSGS(isg).ne.2) goto 10
         k1=iabs(K2SGS(isg,1))
         k2=iabs(K2SGS(isg,2))
         if((k1.eq.4).eqv.(k2.eq.4)) goto 10
         if(nh.ge.MAXC) goto 10
         xm=dsqrt(max(0d0,(PESGS(isg,1)+PESGS(isg,2))**2
     &        -(PXSGS(isg,1)+PXSGS(isg,2))**2
     &        -(PYSGS(isg,1)+PYSGS(isg,2))**2
     &        -(PZSGS(isg,1)+PZSGS(isg,2))**2))
         nh=nh+1
         ih(nh)=n
         ex(nh)=xm-PMSGS(isg,1)-PMSGS(isg,2)
 10   continue
      if(nh.eq.0) return
c     insertion sort by excess mass
      do 20 i=1,nh
         iord(i)=i
 20   continue
      do 30 i=2,nh
         j=iord(i)
         l=i-1
 25      if(l.ge.1) then
            if(ex(iord(l)).gt.ex(j)) then
               iord(l+1)=iord(l)
               l=l-1
               goto 25
            endif
         endif
         iord(l+1)=j
 30   continue
      nv=nh/2
      if(mod(nh,2).eq.1.and.hqrand().lt.0.5) nv=nv+1
      do 40 i=1,nh
         k=ih(iord(i))
         ka=iabs(KATT(k,1))
         ks=isign(1,KATT(k,1))
         kb=(ka/10)*10
         if(i.gt.nh-nv) then
            kn=ks*(kb+3)
         else
            kn=ks*(kb+1)
         endif
         KATT(k,1)=kn
         ITYPAR(k)=kn
 40   continue
      return
      end

c-----------------------------------------------------------------------
      logical FUNCTION HQISOF(i,j,pp2)
c     ZPC getht: AMPT_HQISO=1 makes every scattering that involves a c or
c     cbar isotropic in the pair CM frame (transport cross section = sigma),
c     light-light scatterings keep the forward-peaked pQCD-like form.  The
c     same uniform deviate is used, so the ZPC random sequence is unchanged.
c     AMPT_HQISO_PCUT (GeV): isotropic only if the pair CM momentum
c     sqrt(pp2) is below the cut (pQCD-like weakening at high energy);
c     default 1e9 = always.
      parameter (MAXPTN=400001)
      double precision gx,gy,gz,ft,px,py,pz,e,xmass,pp2
      common /prec2/gx(MAXPTN),gy(MAXPTN),gz(MAXPTN),ft(MAXPTN),
     &     px(MAXPTN), py(MAXPTN), pz(MAXPTN), e(MAXPTN),
     &     xmass(MAXPTN), ityp(MAXPTN)
      include 'hqctl.inc'
      HQISOF=.false.
      if(ihqiso.ne.1) return
      HQISOF=(iabs(ityp(i)).eq.4.or.iabs(ityp(j)).eq.4)
     &     .and.pp2.lt.dble(hqpcut)**2
      return
      end

c-----------------------------------------------------------------------
      SUBROUTINE HQFRAG
c     end of PTOH (hadron momenta/positions set): coalescence criterion
c     for open-charm MESON and BARYON candidates (exactly one c or cbar)
c     and independent fragmentation of those that fail it (Zhang, Zheng,
c     Shi, Lin arXiv:2403.06099 sec 2.3; Fan, Zhang, Zheng, Shi
c     arXiv:2601.12287):
c       keep coalescence if  dmax < p_r  and  m_inv < sum m_q + p_m (m_H
c       - sum m_q)  and  m_T(H) < mT_max ;
c       otherwise the charm quark fragments: z from Peterson (eps_c), a
c       q2-q2bar pair from vacuum (u:d:s = 1:1:0.3), D = (c q2bar) along
c       the c direction with light-cone fraction z of (E+|p|)_c, spin 1
c       with probability PARJ(13).  The rest of the candidate 3-momentum
c       goes to the light remnant: for a meson the pseudoscalar (q2 q1bar),
c       for a baryon the lightest ground-state baryon (q1 q2 qv) - baryon
c       number, charm and 3-momentum are conserved per candidate.
c     (2026-10-19: baryons added; before, charm baryons always coalesced
c     while D mesons fragmented, which inflated Lambda_c/D0 3-4x.)
c     c-cbar states always keep coalescence.  ihmode(n): 0 none,
c     1 coalescence, 2 fragmentation.
      PARAMETER (MAXSTR=150001)
      DOUBLE PRECISION  PXSGS,PYSGS,PZSGS,PESGS,PMSGS,
     1     GXSGS,GYSGS,GZSGS,FTSGS
      double precision tmx,xq(3),yq(3),zq(3),dd,xmi,emc,pc(3),pr(3),
     &     pcm,pplus,pd,ed,xmd,et,qx,qy,qz
      double precision r4(4),p4(4),rr(3,4),vv(3,4)
      COMMON/HMAIN1/EATT,JATT,NATT,NT,NP,N0,N01,N10,N11
      COMMON/HMAIN2/KATT(MAXSTR,4),PATT(MAXSTR,4)
      COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &     K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &     PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)
      COMMON/SOFT/PXSGS(MAXSTR,3),PYSGS(MAXSTR,3),PZSGS(MAXSTR,3),
     &     PESGS(MAXSTR,3),PMSGS(MAXSTR,3),GXSGS(MAXSTR,3),
     &     GYSGS(MAXSTR,3),GZSGS(MAXSTR,3),FTSGS(MAXSTR,3),
     &     K1SGS(MAXSTR,3),K2SGS(MAXSTR,3),NJSGS(MAXSTR)
      COMMON /ARPRC/ ITYPAR(MAXSTR),
     &     GXAR(MAXSTR), GYAR(MAXSTR), GZAR(MAXSTR), FTAR(MAXSTR),
     &     PXAR(MAXSTR), PYAR(MAXSTR), PZAR(MAXSTR), PEAR(MAXSTR),
     &     XMAR(MAXSTR)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      common /hqmod/ ihmode(MAXSTR)
      common /hqfst/ nfcoal, nffrag
      include 'hqctl.inc'
      dimension kq(3), kl(3)
      logical lcrit, lcoal, ldbg0
      integer icdbg
      character*8 cdbg
      data ldbg0/.true./, icdbg/0/
      SAVE
c     src_v15: AMPT_HQCDBG=1 prints every criterion evaluation (diagnostic)
      if(ldbg0) then
         call getenv('AMPT_HQCDBG',cdbg)
         if(cdbg(1:1).eq.'1') icdbg=1
         ldbg0=.false.
      endif
      nfcoal=0
      nffrag=0
      n=0
      do 5 i=1,NATT
         ihmode(i)=0
 5    continue
      do 50 isg=1,NSG
         if(NJSGS(isg).eq.0) goto 50
         n=n+1
         nq=NJSGS(isg)
         if(nq.ne.2.and.nq.ne.3) goto 50
         nc=0
         ic=0
         do 6 l=1,nq
            kq(l)=K2SGS(isg,l)
            if(iabs(kq(l)).eq.4) then
               nc=nc+1
               ic=l
            endif
 6       continue
         if(nc.ne.1) goto 50
         ihmode(n)=1
         if(ihqfrg.ne.1) goto 50
         kc=kq(ic)
c     criterion 1: largest separation at the latest freeze-out time
         tmx=FTSGS(isg,1)
         do 9 l=2,nq
            tmx=max(tmx,FTSGS(isg,l))
 9       continue
         do 10 l=1,nq
            xq(l)=GXSGS(isg,l)+PXSGS(isg,l)/PESGS(isg,l)*
     &           (tmx-FTSGS(isg,l))
            yq(l)=GYSGS(isg,l)+PYSGS(isg,l)/PESGS(isg,l)*
     &           (tmx-FTSGS(isg,l))
            zq(l)=GZSGS(isg,l)+PZSGS(isg,l)/PESGS(isg,l)*
     &           (tmx-FTSGS(isg,l))
 10      continue
c     src_v14 AMPT_HQDREST=1: d in the rest frame of the quarks, as Fan et
c     al. (arXiv:2601.12287 Eq. 1) define it: freeze-out points and
c     momenta boosted to the candidate rest frame, each quark propagated
c     to the latest rest-frame freeze-out time
         if(ihqdrf.eq.1) then
            et=0d0
            qx=0d0
            qy=0d0
            qz=0d0
            do 13 l=1,nq
               et=et+PESGS(isg,l)
               qx=qx+PXSGS(isg,l)
               qy=qy+PYSGS(isg,l)
               qz=qz+PZSGS(isg,l)
 13         continue
            do 14 l=1,nq
               r4(1)=GXSGS(isg,l)
               r4(2)=GYSGS(isg,l)
               r4(3)=GZSGS(isg,l)
               r4(4)=FTSGS(isg,l)
               p4(1)=PXSGS(isg,l)
               p4(2)=PYSGS(isg,l)
               p4(3)=PZSGS(isg,l)
               p4(4)=PESGS(isg,l)
               call hqlbst(r4,-qx/et,-qy/et,-qz/et)
               call hqlbst(p4,-qx/et,-qy/et,-qz/et)
               do 18 m=1,4
                  rr(l,m)=r4(m)
                  vv(l,m)=p4(m)
 18            continue
 14         continue
            tmx=rr(1,4)
            do 16 l=2,nq
               tmx=max(tmx,rr(l,4))
 16         continue
            do 17 l=1,nq
               xq(l)=rr(l,1)+vv(l,1)/vv(l,4)*(tmx-rr(l,4))
               yq(l)=rr(l,2)+vv(l,2)/vv(l,4)*(tmx-rr(l,4))
               zq(l)=rr(l,3)+vv(l,3)/vv(l,4)*(tmx-rr(l,4))
 17         continue
         endif
c     meson: c-qbar distance; baryon (src_v5): distance from the c to the
c     centre of the two light quarks (two-body analogue; the max pairwise
c     separation of dm3 let almost no Lambda_c coalesce)
         if(nq.eq.2) then
            dd=dsqrt((xq(1)-xq(2))**2+(yq(1)-yq(2))**2+(zq(1)-zq(2))**2)
         else
            xc=0.
            yc=0.
            zc=0.
            do 11 l=1,3
               if(l.ne.ic) then
                  xc=xc+0.5*sngl(xq(l))
                  yc=yc+0.5*sngl(yq(l))
                  zc=zc+0.5*sngl(zq(l))
               endif
 11         continue
            dd=dsqrt((xq(ic)-dble(xc))**2+(yq(ic)-dble(yc))**2
     &           +(zq(ic)-dble(zc))**2)
c     src_v16 AMPT_HQBDIST=1: equal-mass Jacobi lambda coordinate
            if(ihqbd.eq.1) dd=dd*dsqrt(2d0/3d0)
         endif
c     criterion 2: invariant mass
         et=0d0
         qx=0d0
         qy=0d0
         qz=0d0
         smq=0.
         do 12 l=1,nq
            et=et+PESGS(isg,l)
            qx=qx+PXSGS(isg,l)
            qy=qy+PYSGS(isg,l)
            qz=qz+PZSGS(isg,l)
            if(ihqcon.eq.1.and.iabs(kq(l)).le.4) then
c     src_v16 AMPT_HQCONST=1: constituent mass (PARF(101-104) = d 0.325,
c     u 0.325, s 0.5, c 1.6 GeV), as Fan et al. Eq. 2 define m_Q
               smq=smq+PARF(100+iabs(kq(l)))
            else
               smq=smq+PMSGS(isg,l)
            endif
 12      continue
         xmi=dsqrt(max(0d0,et**2-qx**2-qy**2-qz**2))
         xmh=XMAR(n)
         xmt=sqrt(PXAR(n)**2+PYAR(n)**2+XMAR(n)**2)
c     src_v7 (2026-10-20): Fan et al. arXiv:2601.12287 assign hadrons with
c     mT below the threshold to coalescence DIRECTLY and test the d and
c     m_inv criteria only above it (AMPT_HQMTMODE=1, default).  dm2-dm6
c     used mT < mT_max AND criteria (AMPT_HQMTMODE=0), which forced almost
c     every charm hadron to fragment.
         lcrit=dd.lt.hqpr.and.xmi.lt.smq+hqpm*(xmh-smq)
         if(ihqmtm.eq.1) then
            lcoal=(xmt.lt.hqmtmx).or.lcrit
         else
            lcoal=lcrit.and.(xmt.lt.hqmtmx)
         endif
         if(icdbg.eq.1) write(6,209) nq,ITYPAR(n),
     &        sqrt(PXAR(n)**2+PYAR(n)**2),xmt,dd,xmi,
     &        smq+hqpm*(xmh-smq),smq,xmh,lcoal
 209     format(' # hqcdb',i2,i6,7f9.4,l2)
         if(lcoal) then
            nfcoal=nfcoal+1
            goto 50
         endif
c     --- fragmentation of the charm quark
         if(NATT.ge.MAXSTR-1) goto 50
         nffrag=nffrag+1
         ihmode(n)=2
         pc(1)=PXSGS(isg,ic)
         pc(2)=PYSGS(isg,ic)
         pc(3)=PZSGS(isg,ic)
c     remnant momentum = all light constituents
         pr(1)=qx-pc(1)
         pr(2)=qy-pc(2)
         pr(3)=qz-pc(3)
         pcm=dsqrt(pc(1)**2+pc(2)**2+pc(3)**2)
         emc=dsqrt(pcm**2+PMSGS(isg,ic)**2)
         pplus=emc+pcm
c     vacuum flavour u:d:s = 1:1:0.3
         r=hqrand()*2.3
         if(r.lt.1.) then
            kf2=2
         elseif(r.lt.2.) then
            kf2=1
         else
            kf2=3
         endif
         isc=isign(1,kc)
c     diquark channel (meson candidates only): c + vacuum diquark (f1 f2)
c     -> charm baryon, remnant qbar1 + anti-diquark -> light antibaryon
         ibar=0
         if(nq.eq.2.and.hqrand().lt.hqfbar) then
            ibar=1
            kf3=2
            r=hqrand()*2.3
            if(r.ge.1..and.r.lt.2.) kf3=1
            if(r.ge.2.) kf3=3
         endif
         ispin=0
         if(hqrand().lt.PARJ(13)) ispin=1
         if(ibar.eq.0) then
c     D: charm kc with antiflavour -isc*kf2 (AMPT/PDG sign rule)
            kfd=(400+10*kf2+2*ispin+1)*isign(1,kc-isc*kf2)
         else
            kfd=isc*kcbary(kf2,kf3)
         endif
         xmd=dble(ULMASS(iabs(kfd)))
c     Peterson z, with the D on shell: need z*pplus > xmd
         zb=0.5*(1.-hqeps)
         fmax=1./(zb*(1.-1./zb-hqeps/(1.-zb))**2)
         do 15 l=1,40
            zt=0.5+0.0125*l
            if(zt.lt.1.) fmax=max(fmax,1./(zt*(1.-1./zt
     &           -hqeps/(1.-zt))**2))
 15      continue
         nz=0
 20      nz=nz+1
         z=hqrand()
         fz=1./(z*(1.-1./z-hqeps/(1.-z))**2)
         if(hqrand()*fmax*1.2.gt.fz.and.nz.lt.1000) goto 20
         if(z*pplus.le.xmd) then
            if(nz.lt.1000) goto 20
            z=min(1.0,real(xmd/pplus)+1e-3)
         endif
         pd=0.5d0*(z*pplus-xmd**2/(z*pplus))
         pd=max(pd,0d0)
         pd=min(pd,pcm)
         ed=dsqrt(pd**2+xmd**2)
c     light remnant species
         if(ibar.eq.1) then
            kp=kq(3-ic)
            kfl=-isc*kbary(iabs(kp),kf2,kf3)
         elseif(nq.eq.2) then
            kp=kq(3-ic)
            if(kf2.eq.iabs(kp)) then
               if(kf2.le.2) then
                  kfl=111
               else
                  kfl=221
               endif
            else
               kmx=max(kf2,iabs(kp))
               kmn=min(kf2,iabs(kp))
               kfl=(100*kmx+10*kmn+1)*isign(1,isc*kf2+kp)*(-1)**kmx
            endif
         else
            m=0
            do 25 l=1,3
               if(l.ne.ic) then
                  m=m+1
                  kl(m)=iabs(kq(l))
               endif
 25         continue
            kl(3)=kf2
            kfl=isc*kbary(kl(1),kl(2),kl(3))
         endif
c     overwrite hadron n with the fragmented charm meson
         KATT(n,1)=kfd
         ITYPAR(n)=kfd
         XMAR(n)=sngl(xmd)
         PXAR(n)=sngl(pc(1)*pd/max(pcm,1d-9))
         PYAR(n)=sngl(pc(2)*pd/max(pcm,1d-9))
         PZAR(n)=sngl(pc(3)*pd/max(pcm,1d-9))
         EATT=EATT-PEAR(n)
         PEAR(n)=sngl(ed)
         EATT=EATT+PEAR(n)
         PATT(n,1)=PXAR(n)
         PATT(n,2)=PYAR(n)
         PATT(n,3)=PZAR(n)
         PATT(n,4)=PEAR(n)
c     light remnant gets the rest of the candidate 3-momentum
         NATT=NATT+1
         KATT(NATT,1)=kfl
         KATT(NATT,2)=0
         KATT(NATT,3)=0
         KATT(NATT,4)=0
         ITYPAR(NATT)=kfl
         XMAR(NATT)=ULMASS(kfl)
         PXAR(NATT)=sngl(pr(1)+pc(1))-PXAR(n)
         PYAR(NATT)=sngl(pr(2)+pc(2))-PYAR(n)
         PZAR(NATT)=sngl(pr(3)+pc(3))-PZAR(n)
         PEAR(NATT)=sqrt(PXAR(NATT)**2+PYAR(NATT)**2+PZAR(NATT)**2
     &        +XMAR(NATT)**2)
         PATT(NATT,1)=PXAR(NATT)
         PATT(NATT,2)=PYAR(NATT)
         PATT(NATT,3)=PZAR(NATT)
         PATT(NATT,4)=PEAR(NATT)
         EATT=EATT+PEAR(NATT)
         GXAR(NATT)=GXAR(n)
         GYAR(NATT)=GYAR(n)
         GZAR(NATT)=GZAR(n)
         FTAR(NATT)=FTAR(n)
         ihmode(NATT)=0
 50   continue
      return
      end

c-----------------------------------------------------------------------
      integer FUNCTION KBARY(k1,k2,k3)
c     lightest ground-state baryon (KF > 0) of three light flavours
c     (1=d, 2=u, 3=s); all species exist in ART.  uds -> Lambda.
      dimension k(3)
      k(1)=k1
      k(2)=k2
      k(3)=k3
      nu=0
      nd=0
      ns=0
      do 10 i=1,3
         if(k(i).eq.2) nu=nu+1
         if(k(i).eq.1) nd=nd+1
         if(k(i).eq.3) ns=ns+1
 10   continue
      kbary=2212
      if(nu.eq.3) kbary=2224
      if(nu.eq.2.and.nd.eq.1) kbary=2212
      if(nu.eq.1.and.nd.eq.2) kbary=2112
      if(nd.eq.3) kbary=1114
      if(nu.eq.2.and.ns.eq.1) kbary=3222
      if(nu.eq.1.and.nd.eq.1.and.ns.eq.1) kbary=3122
      if(nd.eq.2.and.ns.eq.1) kbary=3112
      if(nu.eq.1.and.ns.eq.2) kbary=3322
      if(nd.eq.1.and.ns.eq.2) kbary=3312
      if(ns.eq.3) kbary=3334
      return
      end

c-----------------------------------------------------------------------
      integer FUNCTION KCBARY(k1,k2)
c     ground-state singly-charmed baryon (KF > 0) c + light diquark
c     (1=d, 2=u, 3=s): ud Lambda_c+, uu Sigma_c++, dd Sigma_c0,
c     us Xi_c+, ds Xi_c0, ss Omega_c0
      ka=max(k1,k2)
      kb=min(k1,k2)
      kcbary=4122
      if(ka.eq.2.and.kb.eq.2) kcbary=4222
      if(ka.eq.1.and.kb.eq.1) kcbary=4112
      if(ka.eq.3.and.kb.eq.2) kcbary=4232
      if(ka.eq.3.and.kb.eq.1) kcbary=4132
      if(ka.eq.3.and.kb.eq.3) kcbary=4332
      return
      end

c=======================================================================
c     src_v14: modern PDFs for the triggered cc-bar (AMPT_HQPDF) and
c     Cronin broadening of the pair (AMPT_HQCRONIN).  docs/PLAN_v14.md
c=======================================================================
      SUBROUTINE HQXPDF(X,Q2,XPQ,JBT)
c     PYSTFU, triggered hard process only: x*f of a proton from CTEQ6.1M
c     (LHAPDF cteq61 member 0; Zhang et al. PRC 99 064906) at
c     muF^2 = hqmuf*Q2 in the hard cross section (ihqhrd=1, Zheng et
c     al. PRC 101 034905: muF^2 = 2 mT^2) and at Q2 in the initial-state
c     shower (frozen at Q2 = 4 GeV^2 there, see below).  AMPT_HQPDF=1 and
c     a nucleus (A>2): EPS09s NLO bound-proton
c     modification r_i = 1 + sum_j c_j(x,Q) T_A(s)^j at the nucleon's
c     transverse distance s from its own nucleus centre (Helenius et
c     al. JHEP 07 (2012) 073 Eq. 2.9), off during PYINIT (IHNT2(16)=1)
c     as HIJING's own shadowing.  Isospin and charge conjugation are
c     left to PYSTFU (label 230).  Scale floors: see below.
c     src_v15: the EPS09s factor is HQNMOD; with AMPT_HQSHW=1 it is not
c     applied here but carried as the event weight (HQSHWT).
      include 'hqctl.inc'
      COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      real X,Q2,XPQ(-6:6)
      integer JBT,k
      double precision xd,qd,f(-6:6)
      xd=dble(X)
      qd=dble(Q2)
c     hard cross section: muF^2 = hqmuf*Q2 down to the grid minimum
c     (1.3 GeV)^2.  Initial-state shower: Q2 frozen at 4 GeV^2, exactly as
c     the Duke-Owens branch of PYSTFU it replaces (SD uses MAX(Q2,4.)).
c     Below that the CTEQ charm density vanishes (threshold m_c) and the
c     PYSSPA backward evolution of an incoming charm quark divides by it:
c     the step in TEVB underflows and the clin-12/2018 retry loop (label
c     188) never ends (found at event 504 of a p+p test, 2026-10-22).
      if(ihqhrd.eq.1) then
         qd=dsqrt(max(qd*dble(hqmuf),1.69d0))
      else
         qd=dsqrt(max(qd,4d0))
      endif
      call evolvepdf(xd,qd,f)
      if(ihqpdf.eq.1.and.ihqshw.eq.0.and.(JBT.eq.1.or.JBT.eq.2).and.
     &     IHPR2(6).ne.0.and.IHNT2(16).ne.1) call hqnmod(JBT,xd,qd,f)
      do 10 k=-6,6
         XPQ(k)=real(max(f(k),0d0))
 10   continue
      return
      end

c-----------------------------------------------------------------------
      SUBROUTINE HQNMOD(JBT,xd,qd,f)
c     src_v15 (moved out of HQXPDF, same arithmetic): EPS09s NLO bound-
c     proton modification of the free x*f(-6:6) at (x, Q) for beam JBT,
c     r_i = 1 + sum_j c_j T_A(s)^j at the transverse distance s of the
c     nucleon of the current collision from its own nucleus centre.
      include 'hqctl.inc'
      COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      COMMON/hjcrdn/YP(3,300),YT(3,300)
      integer JBT,ia,io
      double precision xd,qd,f(-6:6),s,t,cuv(4),cdv(4),cu(4),cd(4),
     &     cs(4),cc(4),cb(4),cg(4),ruv,rdv,ru,rd,rs,rc,rb,rg,uv,dv,
     &     ta_eps09s,hqrnp
      external ta_eps09s,hqrnp
         ia=IHNT2(2*JBT-1)
         if(ia.gt.2) then
            if(JBT.eq.1) then
               s=dsqrt(dble(YP(1,IHNT2(11)))**2
     &              +dble(YP(2,IHNT2(11)))**2)
            else
               s=dsqrt(dble(YT(1,IHNT2(12)))**2
     &              +dble(YT(2,IHNT2(12)))**2)
            endif
            io=2
            call eps09s_c(io,ihqpst,ia,xd,qd,cuv,cdv,cu,cd,cs,cc,cb,cg)
            t=ta_eps09s(ia,s)
            ruv=hqrnp(cuv,t)
            rdv=hqrnp(cdv,t)
            ru=hqrnp(cu,t)
            rd=hqrnp(cd,t)
            rs=hqrnp(cs,t)
            rc=hqrnp(cc,t)
            rb=hqrnp(cb,t)
            rg=hqrnp(cg,t)
            uv=f(2)-f(-2)
            dv=f(1)-f(-1)
            f(-2)=ru*f(-2)
            f(2)=ruv*uv+f(-2)
            f(-1)=rd*f(-1)
            f(1)=rdv*dv+f(-1)
            f(3)=rs*f(3)
            f(-3)=rs*f(-3)
            f(4)=rc*f(4)
            f(-4)=rc*f(-4)
            f(5)=rb*f(5)
            f(-5)=rb*f(-5)
            f(0)=rg*f(0)
         endif
      return
      end

c-----------------------------------------------------------------------
      SUBROUTINE HQSHWT
c     src_v15 AMPT_HQSHW=1, PYRAND after the channel choice of the
c     triggered hard process.  The pair was sampled with free-nucleon
c     PDFs; its nuclear modification is carried as the event weight
c     hqwsh = R_1(x1) R_2(x2) of the two incoming partons (flavours
c     MINT(15), MINT(16), x = VINT(41), VINT(42), Q^2 = VINT(52)).
c     HQPDF=1: EPS09s, ratio of the bound to the free nucleon density
c     of that flavour (isospin conjugate for a neutron, PYSTFU label
c     230); HQPDF=0: HIJING's shadowing as PYSTFU applies it (HQRRX);
c     HQPDF=2 (free CTEQ6.1M): 1.  A triggered event carries exactly one
c     pair, so without the weight the per-collision yield modification
c     is normalised away (R_AA per c quark integrates to 1).
      include 'hqctl.inc'
      COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      COMMON/PYINT1/MINT(400),VINT(400)
      integer jbt,kfl,k,ia
      real hqrrx
      double precision xd,qd,f(-6:6),g(-6:6),r(2)
      external hqrrx
      do 20 jbt=1,2
         r(jbt)=1d0
         ia=IHNT2(2*jbt-1)
         if(ia.le.2.or.IHPR2(6).eq.0.or.IHNT2(16).eq.1) goto 20
         kfl=MINT(14+jbt)
         if(kfl.eq.21) kfl=0
         if(iabs(kfl).gt.6) goto 20
         if(ihqpdf.eq.1) then
            xd=dble(VINT(40+jbt))
            qd=dsqrt(max(dble(VINT(52))*dble(hqmuf),1.69d0))
            call evolvepdf(xd,qd,f)
            do 10 k=-6,6
               g(k)=f(k)
 10         continue
            call hqnmod(jbt,xd,qd,g)
            k=kfl
            if(MINT(10+jbt).eq.2112.and.k.ne.0.and.iabs(k).le.2)
     &           k=isign(3-iabs(k),k)
            if(f(k).gt.0d0) r(jbt)=max(g(k),0d0)/f(k)
         elseif(ihqpdf.eq.0) then
            r(jbt)=dble(hqrrx(jbt,VINT(40+jbt)))
         endif
 20   continue
      hqwsh=real(r(1)*r(2))
      khqw(1)=MINT(15)
      khqw(2)=MINT(16)
      hqwx(1)=VINT(41)
      hqwx(2)=VINT(42)
      hqwq2=VINT(52)
      hqwr(1)=real(r(1))
      hqwr(2)=real(r(2))
      return
      end

c-----------------------------------------------------------------------
      REAL FUNCTION HQRRX(JBT,X)
c     src_v15: HIJING 1.383 nuclear shadowing factor exactly as PYSTFU
c     applies it to the triggered hard process (AMPT_HQSHAD-scaled when
c     HQSHAD >= 0), at the nucleon of the current collision
      include 'hqctl.inc'
      COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      COMMON/hjcrdn/YP(3,300),YT(3,300)
      integer JBT
      real X,ATNM,BBR2,ABX,APX,AAX,RRX
      ATNM=IHNT2(2*JBT-1)
      IF(JBT.EQ.1) THEN
         BBR2=(YP(1,IHNT2(11))**2+YP(2,IHNT2(11))**2)/1.44/
     1        ATNM**0.66666
      ELSE
         BBR2=(YT(1,IHNT2(12))**2+YT(2,IHNT2(12))**2)/1.44/
     1        ATNM**0.66666
      ENDIF
      BBR2=MIN(1.0,BBR2)
      ABX=(ATNM**0.33333333-1.0)
      APX=HIPR1(6)*4.0/3.0*ABX*SQRT(1.0-BBR2)
      AAX=1.192*ALOG(ATNM)**0.1666666
      RRX=AAX*(X**3-1.2*X**2+0.21*X)+1.0
     &     -(APX-1.079*ABX*SQRT(X)/ALOG(ATNM+1.0))
     1     *EXP(-X**2.0/0.01)
      IF(hqshad.GE.0.) RRX=1.0+hqshad*(RRX-1.0)
      HQRRX=RRX
      return
      end

c-----------------------------------------------------------------------
      double precision FUNCTION hqrnp(c,t)
c     EPS09s power series 1 + c1 T + c2 T^2 + c3 T^3 + c4 T^4
      double precision c(4),t
      hqrnp=1d0+t*(c(1)+t*(c(2)+t*(c(3)+t*c(4))))
      return
      end

c-----------------------------------------------------------------------
      SUBROUTINE HQTRGG(np,nt,xp,yp,xt,yt)
c     HIJING, at the choice of the triggered NN collision: number of
c     primary NN collisions of its projectile (np) and target (nt)
c     nucleon, and their transverse distances from their own nucleus
c     centres (Cronin width, Zhang et al. PLB 846 138219 Eqs. 2-3)
      include 'hqctl.inc'
      integer np,nt
      real xp,yp,xt,yt
      nhqcp=np+nt
      hqwsh=1.
      hqwr(1)=1.
      hqwr(2)=1.
      hqsab(1)=sqrt(xp**2+yp**2)
      hqsab(2)=sqrt(xt**2+yt**2)
      return
      end

c-----------------------------------------------------------------------
      SUBROUTINE HQCKICK(npar)
c     AMPT_HQCRONIN=delta, end of HTOP after the HQDIRECT replacement:
c     Zhang, Zheng, Shi & Lin, PLB 846 138219 (2023), Eqs. 1-3.  Each
c     c-cbar pair gets a transverse kick kT from exp(-kT^2/w^2)/(pi w^2),
c     w = w0 sqrt(1 + (ncoll - 2) delta), ncoll = primary NN collisions
c     of the two nucleons of the triggered collision; w0 = 0.35 GeV
c     sqrt(b0 (2+a0) / (bL (2+aL))), a0 = 0.5, b0 = 0.9 GeV^-2, aL = 0.8,
c     bL from the local scaling of Zhang et al. PRC 104 014908 Eq. 3
c     (b_pp = 0.7 GeV^-2, T_p = 0.22 fm^-2, beta(200 GeV) = 0.620; used
c     here only for w0, the bulk Lund b is unchanged).  The kick is a
c     transverse boost of the pair at fixed pair rapidity; the light
c     partons share -sum kT equally (momentum conservation).
      PARAMETER (MAXSTR=150001, MAXPTN=400001)
      DOUBLE PRECISION  GX0, GY0, GZ0, FT0, PX0, PY0, PZ0, E0, XMASS0
      DOUBLE PRECISION vxp0,vyp0,vzp0,xstrg0,ystrg0,xstrg,ystrg
      COMMON /prec1/GX0(MAXPTN),GY0(MAXPTN),GZ0(MAXPTN),FT0(MAXPTN),
     &     PX0(MAXPTN), PY0(MAXPTN), PZ0(MAXPTN), E0(MAXPTN),
     &     XMASS0(MAXPTN), ITYP0(MAXPTN)
      COMMON /ilist7/ LSTRG0(MAXPTN), LPART0(MAXPTN)
      common /precpa/vxp0(MAXPTN),vyp0(MAXPTN),vzp0(MAXPTN),
     1     xstrg0(MAXPTN),ystrg0(MAXPTN),
     2     xstrg(MAXPTN),ystrg(MAXPTN),istrg0(MAXPTN),istrg(MAXPTN)
      COMMON/HPARNT/HIPR1(100),IHPR2(50),HINT1(100),IHNT2(50)
      COMMON/RNDF77/NSEED
      COMMON /AREVT/ IAEVT, IARUN, MISS
      include 'hqctl.inc'
      integer npar,i,j,k,nc,nl,npr,ia,ib,ipc(50),jpc(50),iuse(50)
      double precision ta,tb,tp,bl,beta,w0,w,akt,phi,sx,sy,
     &     pe,ppx,ppy,ppz,am2,yy,amt,pe2,ppx2,ppy2,ppz2,
     &     bx,by,bz,bx2,by2,bz2,q(4),ta_eps09s,dcx,dcy,dkx,dky,ek
      external ta_eps09s
      real RANART
      external RANART
      save
      if(hqcron.lt.0.) return
c     the width of this event
      tp=0.22d0
      ia=IHNT2(1)
      ib=IHNT2(3)
      ta=tp
      tb=tp
      if(ia.gt.2) ta=max(tp,ta_eps09s(ia,dble(hqsab(1))))
      if(ib.gt.2) tb=max(tp,ta_eps09s(ib,dble(hqsab(2))))
      beta=0.620d0
      if(HINT1(1).gt.200.) beta=beta+0.112d0*dlog(dble(HINT1(1))/200d0)
      bl=0.7d0/(dsqrt(ta*tb)/tp)**beta
      w0=0.35d0*dsqrt(0.9d0*2.5d0/(bl*2.8d0))
      w=w0*dsqrt(max(0d0,1d0+dble(nhqcp-2)*dble(hqcron)))
c     pair the charm partons: c and cbar on the same string, else the
c     only c with the only cbar of the event
      nc=0
      do 10 i=1,npar
         if(ITYP0(i).ne.4) goto 10
         if(nc.ge.50) goto 10
         nc=nc+1
         ipc(nc)=i
         jpc(nc)=0
 10   continue
      do 12 k=1,50
         iuse(k)=0
 12   continue
      npr=0
      do 20 k=1,nc
         i=ipc(k)
         do 15 j=1,npar
            if(ITYP0(j).ne.-4) goto 15
            if(istrg0(LSTRG0(j)).ne.istrg0(LSTRG0(i))) goto 15
            do 14 ia=1,k-1
               if(jpc(ia).eq.j) goto 15
 14         continue
            jpc(k)=j
            npr=npr+1
            goto 20
 15      continue
 20   continue
      if(npr.eq.0.and.nc.eq.1) then
         nb=0
         do 22 j=1,npar
            if(ITYP0(j).eq.-4) then
               nb=nb+1
               jpc(1)=j
            endif
 22      continue
         if(nb.eq.1) npr=1
         if(nb.ne.1) jpc(1)=0
      endif
      sx=0d0
      sy=0d0
      do 40 k=1,nc
         if(jpc(k).eq.0) goto 40
         i=ipc(k)
         j=jpc(k)
         akt=w*dsqrt(-dlog(max(1d-300,dble(RANART(NSEED)))))
         phi=6.283185307179586d0*dble(RANART(NSEED))
         dkx=akt*dcos(phi)
         dky=akt*dsin(phi)
c     old pair: velocity and rapidity
         pe=E0(i)+E0(j)
         ppx=PX0(i)+PX0(j)
         ppy=PY0(i)+PY0(j)
         ppz=PZ0(i)+PZ0(j)
         am2=pe**2-ppx**2-ppy**2-ppz**2
         yy=0.5d0*dlog((pe+ppz)/(pe-ppz))
c     new pair: pT + kT at the same rapidity and invariant mass
         ppx2=ppx+dkx
         ppy2=ppy+dky
         amt=dsqrt(am2+ppx2**2+ppy2**2)
         pe2=amt*dcosh(yy)
         ppz2=amt*dsinh(yy)
         bx=ppx/pe
         by=ppy/pe
         bz=ppz/pe
         bx2=ppx2/pe2
         by2=ppy2/pe2
         bz2=ppz2/pe2
         do 30 ia=1,2
            ib=i
            if(ia.eq.2) ib=j
            q(1)=PX0(ib)
            q(2)=PY0(ib)
            q(3)=PZ0(ib)
            q(4)=E0(ib)
            call hqlbst(q,-bx,-by,-bz)
            call hqlbst(q,bx2,by2,bz2)
            PX0(ib)=q(1)
            PY0(ib)=q(2)
            PZ0(ib)=q(3)
            E0(ib)=dsqrt(q(1)**2+q(2)**2+q(3)**2+XMASS0(ib)**2)
            if(ihqdir.eq.1) FT0(ib)=0.197d0*E0(ib)
     &           /(PX0(ib)**2+PY0(ib)**2+XMASS0(ib)**2)
 30      continue
         sx=sx+(PX0(i)+PX0(j))-ppx
         sy=sy+(PY0(i)+PY0(j))-ppy
 40   continue
c     recoil shared by the light partons
      nl=0
      do 50 i=1,npar
         if(iabs(ITYP0(i)).ne.4) nl=nl+1
 50   continue
      ek=0d0
      if(nl.gt.0.and.npr.gt.0) then
         dcx=sx/dble(nl)
         dcy=sy/dble(nl)
         do 55 i=1,npar
            if(iabs(ITYP0(i)).eq.4) goto 55
            ek=ek-E0(i)
            PX0(i)=PX0(i)-dcx
            PY0(i)=PY0(i)-dcy
            E0(i)=dsqrt(PX0(i)**2+PY0(i)**2+PZ0(i)**2+XMASS0(i)**2)
            ek=ek+E0(i)
 55      continue
      endif
      write(6,201) IAEVT,nhqcp,hqsab,bl,w0,w,nc,npr,sx,sy,ek
 201  format(' # hqcron ev',i6,' ncoll',i4,' sA,sB',2f7.3,' bL',f6.3,
     &     ' w0',f6.3,' w',f7.3,' nc',i3,' npair',i3,' kT',2f8.4,
     &     ' dE_light',e11.3)
      return
      end

c-----------------------------------------------------------------------
      SUBROUTINE HQLBST(q,bx,by,bz)
c     Lorentz boost of the four-vector q (px,py,pz,E) by velocity b
      double precision q(4),bx,by,bz,b2,g,bp,g2
      b2=bx**2+by**2+bz**2
      if(b2.le.1d-20) return
      g=1d0/dsqrt(1d0-b2)
      bp=bx*q(1)+by*q(2)+bz*q(3)
      g2=(g-1d0)/b2
      q(1)=q(1)+g2*bp*bx+g*bx*q(4)
      q(2)=q(2)+g2*bp*by+g*by*q(4)
      q(3)=q(3)+g2*bp*bz+g*bz*q(4)
      q(4)=g*(q(4)+bp)
      return
      end

c-----------------------------------------------------------------------
      SUBROUTINE HQVMAX
c     once per event: report a new maximum of SIGS/XSEC_max inside the
c     triggered hard process (hqvio, recorded in PYRAND); > 1 means the
c     maximum was violated there and AMPT_HQXSF must be raised
      include 'hqctl.inc'
      real vmx
      data vmx/0./
      save vmx
      if(ihqpdf.le.0) return
      if(hqvio.gt.vmx) then
         vmx=hqvio
         write(6,'(a,f9.4)') ' # hqpdf: trigger SIGS/max ratio ',vmx
      endif
      return
      end

c-----------------------------------------------------------------------
      SUBROUTINE HQHADR(nt,ntmax)
c     src_v17 AMPT_HQHAD=1: hadronic rescattering of charm hadrons in
c     ART.  In stock AMPT charm hadrons enter ART as LB = PDG + 10000 and
c     free-stream.  Here, once per ART time step (after RELCOL, before
c     the position update), each charm hadron may scatter elastically
c     and isotropically (in the pair c.m.) off one formed light hadron,
c     with the He-Fries-Rapp effective cross sections (PLB 701 (2011)
c     445: isotropic 7 mb with mesons and 10 mb with baryons reproduce
c     their empirical-amplitude D relaxation rate; AMPT_HQHADM/B change
c     them).  Default species: D, D*, (anti); Ds/Ds* with AMPT_HQHADDS=1,
c     charm baryons with AMPT_HQHADLC=1; charmonia never.  Collision
c     criterion as in ART/ZPC: closest approach in the pair c.m. frame
c     with b < sqrt(sigma/pi), reached (lab time) within this step; the
c     earliest such partner is taken, at most one collision per charm
c     hadron per step, not the same partner twice in a row.  The charm
c     hadron gets the full two-body elastic kinematics; by default the
c     light partner keeps its momentum (probe approximation, as in
c     Langevin/Boltzmann heavy-flavour transport: the bulk evolves
c     exactly as without rescattering, so every event stays paired with
c     the same seeds without HQHAD); AMPT_HQHADREC=1 also recoils the
c     partner (exact energy-momentum conservation).  Angles come from a
c     private generator (HQRNH), never RANART, so the off/on runs share
c     every other random number.  nhqhc(1..3):
c     collisions with mesons, with baryons, eligible charm hadrons present
c     at the last active step (printed as '# hqhad' at nt = ntmax-1).
      include 'hqctl.inc'
      PARAMETER (MAXSTR=150001,MAXR=1)
      COMMON /AA/ R(3,MAXSTR)
      COMMON /BB/ P(3,MAXSTR)
      COMMON /CC/ E(MAXSTR)
      COMMON /EE/ ID(MAXSTR),LB(MAXSTR)
      COMMON /RR/ MASSR(0:MAXR)
      COMMON /RUN/ NUM
      common/input1/ MASSPR,MASSTA,ISEED,IAVOID,DT
      COMMON/RNDF77/NSEED
      integer nt,ntmax,i,j,k,irun,iso,kf,ka,ilb,jbest,ibar,lastp
      real hqrnh
      external hqrnh
      double precision e1,e2,b(3),bb,g,q1(4),q2(4),x1(4),x2(4),t0,
     &     v1(3),v2(3),dr(3),dv(3),dv2,tau,b2,bmax2,tl,tbest,pst,cth,
     &     sth,phi,xm(4),sg,bp,pin(4),dmx
      dimension lastp(MAXSTR)
      integer iseedh
      save lastp,dmx,iseedh
      data dmx/0d0/, iseedh/0/
      if(ihqhad.ne.1) return
c     seed the private generator once per job from ART's current seed
c     (job specific; identical in the runs with and without HQHAD)
      if(iseedh.eq.0) then
         ihqhsd=1+mod(iabs(NSEED),2147483646)
         iseedh=1
      endif
      if(nt.eq.1) then
         do 5 i=1,MAXSTR
            lastp(i)=0
 5       continue
      endif
c     at nt = ntmax-1 ART reads in all unformed hadrons: stop before
      if(nt.ge.ntmax-1) then
         if(nt.eq.ntmax-1) then
            write(6,210) nhqhc(3),nhqhc(1),nhqhc(2),dmx
 210        format(' # hqhad charm',i5,' coll mesons',i6,' baryons',i6,
     &           ' max dP',1p,e10.2)
            dmx=0d0
            nhqhc(1)=0
            nhqhc(2)=0
            nhqhc(3)=0
         endif
         return
      endif
      iso=0
      do 100 irun=1,NUM
      iso=iso+MASSR(irun-1)
      do 90 k=1,MASSR(irun)
         i=k+iso
         if(E(i).le.0.) goto 90
         kf=LB(i)-10000
         ka=iabs(kf)
c     open charm only: mesons 4xy (not 44y), baryons 4xyz
         if(ka.ge.400.and.ka.lt.500) then
            if(ka/10.eq.44) goto 90
            if(ka/10.eq.43.and.ihqhds.ne.1) goto 90
         elseif(ka.ge.4000.and.ka.lt.5000) then
            if(ihqhlc.ne.1) goto 90
         else
            goto 90
         endif
         if(nt.eq.ntmax-2) nhqhc(3)=nhqhc(3)+1
         e1=dsqrt(dble(E(i))**2+dble(P(1,i))**2+dble(P(2,i))**2
     &        +dble(P(3,i))**2)
         jbest=0
         tbest=dble(DT)
         do 80 j=1+iso,MASSR(irun)+iso
            if(j.eq.i.or.E(j).le.0.) goto 80
            ilb=iabs(LB(j))
            if(ilb.gt.45.or.ilb.eq.42) goto 80
            if(j.eq.lastp(i)) goto 80
c     quick spatial preselection (|dr| < 3 fm)
            if((R(1,i)-R(1,j))**2+(R(2,i)-R(2,j))**2
     &           +(R(3,i)-R(3,j))**2.gt.9.) goto 80
            ibar=0
            if(ilb.eq.1.or.ilb.eq.2.or.(ilb.ge.6.and.ilb.le.17)
     &           .or.ilb.eq.40.or.ilb.eq.41.or.ilb.eq.45) ibar=1
            sg=hqsgm
            if(ibar.eq.1) sg=hqsgb
            if(sg.le.0.) goto 80
            bmax2=0.1d0*dble(sg)/3.141592653589793d0
            e2=dsqrt(dble(E(j))**2+dble(P(1,j))**2+dble(P(2,j))**2
     &           +dble(P(3,j))**2)
c     boost to the pair c.m.; both particles are at lab time 0
            do 10 m=1,3
               b(m)=(dble(P(m,i))+dble(P(m,j)))/(e1+e2)
 10         continue
            bb=b(1)**2+b(2)**2+b(3)**2
            if(bb.ge.1d0) goto 80
            g=1d0/dsqrt(1d0-bb)
            q1(4)=e1
            q2(4)=e2
            x1(4)=0d0
            x2(4)=0d0
            do 12 m=1,3
               q1(m)=dble(P(m,i))
               q2(m)=dble(P(m,j))
               x1(m)=dble(R(m,i))
               x2(m)=dble(R(m,j))
 12         continue
            call hqbst4(q1,b,g,bb)
            call hqbst4(q2,b,g,bb)
            call hqbst4(x1,b,g,bb)
            call hqbst4(x2,b,g,bb)
            t0=max(x1(4),x2(4))
            do 14 m=1,3
               v1(m)=q1(m)/q1(4)
               v2(m)=q2(m)/q2(4)
               dr(m)=(x1(m)+v1(m)*(t0-x1(4)))-(x2(m)+v2(m)*(t0-x2(4)))
               dv(m)=v1(m)-v2(m)
 14         continue
            dv2=dv(1)**2+dv(2)**2+dv(3)**2
            if(dv2.lt.1d-12) goto 80
            tau=-(dr(1)*dv(1)+dr(2)*dv(2)+dr(3)*dv(3))/dv2
c     tau may be < 0: the closest approach can lie between the c.m. times
c     of the two lab-time-0 events; the lab-time window below decides
            b2=(dr(1)+dv(1)*tau)**2+(dr(2)+dv(2)*tau)**2
     &           +(dr(3)+dv(3)*tau)**2
            if(b2.gt.bmax2) goto 80
c     lab time of the collision point (pair midpoint at c.m. time t0+tau)
            do 16 m=1,3
               xm(m)=0.5d0*((x1(m)+v1(m)*(t0+tau-x1(4)))
     &              +(x2(m)+v2(m)*(t0+tau-x2(4))))
 16         continue
            bp=b(1)*xm(1)+b(2)*xm(2)+b(3)*xm(3)
            tl=g*((t0+tau)+bp)
            if(tl.lt.0d0.or.tl.ge.tbest) goto 80
            tbest=tl
            jbest=j
 80      continue
         if(jbest.eq.0) goto 90
c     elastic isotropic scattering in the pair c.m.
         j=jbest
         e2=dsqrt(dble(E(j))**2+dble(P(1,j))**2+dble(P(2,j))**2
     &        +dble(P(3,j))**2)
         do 20 m=1,3
            b(m)=(dble(P(m,i))+dble(P(m,j)))/(e1+e2)
 20      continue
         bb=b(1)**2+b(2)**2+b(3)**2
         g=1d0/dsqrt(1d0-bb)
         q1(4)=e1
         q2(4)=e2
         do 22 m=1,3
            q1(m)=dble(P(m,i))
            q2(m)=dble(P(m,j))
 22      continue
         do 21 m=1,3
            pin(m)=q1(m)+q2(m)
 21      continue
         pin(4)=e1+e2
         call hqbst4(q1,b,g,bb)
         call hqbst4(q2,b,g,bb)
         pst=dsqrt(q1(1)**2+q1(2)**2+q1(3)**2)
         cth=1d0-2d0*dble(hqrnh())
         sth=dsqrt(max(0d0,1d0-cth*cth))
         phi=6.283185307179586d0*dble(hqrnh())
         q1(1)=pst*sth*dcos(phi)
         q1(2)=pst*sth*dsin(phi)
         q1(3)=pst*cth
         q2(1)=-q1(1)
         q2(2)=-q1(2)
         q2(3)=-q1(3)
         do 24 m=1,3
            b(m)=-b(m)
 24      continue
         call hqbst4(q1,b,g,bb)
         call hqbst4(q2,b,g,bb)
         do 26 m=1,3
            P(m,i)=sngl(q1(m))
            if(ihqhrc.eq.1) P(m,j)=sngl(q2(m))
            dmx=max(dmx,dabs(q1(m)+q2(m)-pin(m)))
 26      continue
         dmx=max(dmx,dabs(q1(4)+q2(4)-pin(4)))
         lastp(i)=j
         ilb=iabs(LB(j))
         if(ilb.eq.1.or.ilb.eq.2.or.(ilb.ge.6.and.ilb.le.17)
     &        .or.ilb.eq.40.or.ilb.eq.41.or.ilb.eq.45) then
            nhqhc(2)=nhqhc(2)+1
         else
            nhqhc(1)=nhqhc(1)+1
         endif
 90   continue
 100  continue
      return
      end

c-----------------------------------------------------------------------
      SUBROUTINE HQBST4(x,b,g,bb)
c     src_v17: Lorentz boost of the 4-vector x = (x1,x2,x3,x0) into the
c     frame moving with velocity b (g = gamma, bb = b.b):
c     x0' = g (x0 - b.x),  x' = x + b [(g-1) b.x / bb - g x0]
      double precision x(4),b(3),g,bb,bx,x0,f
      bx=b(1)*x(1)+b(2)*x(2)+b(3)*x(3)
      x0=x(4)
      x(4)=g*(x0-bx)
      if(bb.gt.1d-20) then
         f=(g-1d0)*bx/bb-g*x0
      else
         f=-g*x0
      endif
      x(1)=x(1)+b(1)*f
      x(2)=x(2)+b(2)*f
      x(3)=x(3)+b(3)*f
      return
      end

c-----------------------------------------------------------------------
      real FUNCTION hqrnh()
c     src_v17: private Park-Miller generator of the charm rescattering
c     angles (state ihqhsd), so AMPT_HQHAD never changes RANART/NSEED
      include 'hqctl.inc'
      integer*8 ls
      ls=mod(int(ihqhsd,8)*48271_8,2147483647_8)
      if(ls.le.0) ls=ls+2147483647_8
      ihqhsd=int(ls)
      hqrnh=real(ihqhsd)/2147483647.
      return
      end
