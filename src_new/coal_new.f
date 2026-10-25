      SUBROUTINE CNPAR(ION, RBM)
*
*     He & Lin's improved quark coalescence, PRC 96 014910 (arXiv:1703.02673),
*     reimplemented because it is published but not distributed: the newest
*     public AMPT (v1.26t9b/v2.26t9b, 12/2018) carries the old `coales`, and
*     Lin & Zheng's review arXiv:2110.02989 states the improved coalescence
*     "will be released in the future".
*
*     Switch: AMPT_COALNEW=1 (off by default, so an unset environment is
*     byte-identical to the old algorithm). Strength: AMPT_RBM, default 0.61,
*     which is the value used throughout He & Lin's paper.
*
      INTEGER ION, ILEN, ISTAT, IERR, ITMP, IFIRST
      DOUBLE PRECISION RBM, RTMP
      CHARACTER*32 IENV
      COMMON /cnord/ CNORD
      DOUBLE PRECISION RBMHQ
      COMMON /cnhq/ RBMHQ
      SAVE IFIRST
      DATA IFIRST /0/

      ION = 0
      RBM = 0.61d0
      CALL GET_ENVIRONMENT_VARIABLE('AMPT_COALNEW', IENV, ILEN, ISTAT)
      IF (ISTAT .EQ. 0 .AND. ILEN .GT. 0) THEN
         READ (IENV(1:ILEN), *, IOSTAT=IERR) ITMP
         IF (IERR .EQ. 0 .AND. ITMP .EQ. 1) ION = 1
      END IF
      IF (ION .NE. 1) RETURN

*     The paper does not specify the order in which partons are processed, and a
*     greedy algorithm depends on it. AMPT_CNORDER=1 reverses the sweep, as a
*     sensitivity probe on that choice.
*     D-meson framework src_v9 (2026-10-21): AMPT_CNORDER=2 sweeps the seeds in
*     a random order, drawn per event from a private LCG seeded by the event
*     number (RANART untouched).  The stored order is projectile strings, then
*     target strings, then independent strings, so the default sweep gives the
*     projectile partons first pick of partners: a projectile/target asymmetry
*     (rapidity-even v1 relative to the reaction plane, VALIDATION G15).
      CALL GET_ENVIRONMENT_VARIABLE('AMPT_CNORDER', IENV, ILEN, ISTAT)
      CNORD = 0
      IF (ISTAT .EQ. 0 .AND. ILEN .GT. 0) THEN
         READ (IENV(1:ILEN), *, IOSTAT=IERR) ITMP
         IF (IERR .EQ. 0 .AND. (ITMP .EQ. 1 .OR. ITMP .EQ. 2)) CNORD = ITMP
      END IF

      CALL GET_ENVIRONMENT_VARIABLE('AMPT_RBM', IENV, ILEN, ISTAT)
      IF (ISTAT .EQ. 0 .AND. ILEN .GT. 0) THEN
         READ (IENV(1:ILEN), *, IOSTAT=IERR) RTMP
         IF (IERR .EQ. 0 .AND. RTMP .GT. 0d0) RBM = RTMP
      END IF

*     D-meson framework (2026-10-20): separate r_BM for groups containing a
*     charm quark, as in the improved AMPT (Zheng et al. PRC 101 034905:
*     r_BM(HQ) = 1.0; Fan et al. arXiv:2601.12287: 2).  Default = light r_BM.
      RBMHQ = RBM
      CALL GET_ENVIRONMENT_VARIABLE('AMPT_RBMHQ', IENV, ILEN, ISTAT)
      IF (ISTAT .EQ. 0 .AND. ILEN .GT. 0) THEN
         READ (IENV(1:ILEN), *, IOSTAT=IERR) RTMP
         IF (IERR .EQ. 0 .AND. RTMP .GT. 0d0) RBMHQ = RTMP
      END IF

      IF (IFIRST .EQ. 0) THEN
         IFIRST = 1
         WRITE (6, *) 'improved quark coalescence ON (He & Lin), r_BM =',
     1        RBM, ' r_BM(charm) =', RBMHQ, ' sweep order =', CNORD
      END IF
      RETURN
      END

c=======================================================================

      SUBROUTINE CNDIS2(I, J, D)
*
*     Relative distance between two flat-list partons in AMPT's own convention:
*     boost to the pair rest frame, propagate the earlier one to the later
*     formation time, take the separation. This is `locldr(2,...)`, reused so the
*     number means the same thing as everywhere else in the code.
*
      PARAMETER (MAXPTN=400001)
      DOUBLE PRECISION gxp, gyp, gzp, ftp, pxp, pyp, pzp, pep, pmp
      common /loclco/gxp(3),gyp(3),gzp(3),ftp(3),
     1     pxp(3),pyp(3),pzp(3),pep(3),pmp(3)
      DOUBLE PRECISION CNX, CNY, CNZ, CNT, CNPX, CNPY, CNPZ, CNE, CNM
      COMMON /cnlist/ CNX(MAXPTN), CNY(MAXPTN), CNZ(MAXPTN), CNT(MAXPTN),
     1     CNPX(MAXPTN), CNPY(MAXPTN), CNPZ(MAXPTN), CNE(MAXPTN),
     2     CNM(MAXPTN), CNK(MAXPTN), CNUSE(MAXPTN), CNSTR(MAXPTN), CNP
      INTEGER I, J
      DOUBLE PRECISION D

      ftp(1) = CNT(I)
      gxp(1) = CNX(I)
      gyp(1) = CNY(I)
      gzp(1) = CNZ(I)
      pxp(1) = CNPX(I)
      pyp(1) = CNPY(I)
      pzp(1) = CNPZ(I)
      pep(1) = CNE(I)
      pmp(1) = CNM(I)
      ftp(2) = CNT(J)
      gxp(2) = CNX(J)
      gyp(2) = CNY(J)
      gzp(2) = CNZ(J)
      pxp(2) = CNPX(J)
      pyp(2) = CNPY(J)
      pzp(2) = CNPZ(J)
      pep(2) = CNE(J)
      pmp(2) = CNM(J)
      D = 0d0
      CALL LOCLDR(2, D)
      RETURN
      END

c=======================================================================

      SUBROUTINE CNDIS3(I, J, K, D)
*
*     He & Lin's d_B: the average of the three relative distances among three
*     partons in their common rest frame. locldr's 3-body branch boosts and
*     propagates but does NOT return a distance (it leaves the positions in
*     /prtn23/), so the average is formed here from those positions.
*
      PARAMETER (MAXPTN=400001)
      DOUBLE PRECISION gxp, gyp, gzp, ftp, pxp, pyp, pzp, pep, pmp
      common /loclco/gxp(3),gyp(3),gzp(3),ftp(3),
     1     pxp(3),pyp(3),pzp(3),pep(3),pmp(3)
      DOUBLE PRECISION gxp0, gyp0, gzp0, ft0fom
      common /prtn23/ gxp0(3),gyp0(3),gzp0(3),ft0fom
      DOUBLE PRECISION CNX, CNY, CNZ, CNT, CNPX, CNPY, CNPZ, CNE, CNM
      COMMON /cnlist/ CNX(MAXPTN), CNY(MAXPTN), CNZ(MAXPTN), CNT(MAXPTN),
     1     CNPX(MAXPTN), CNPY(MAXPTN), CNPZ(MAXPTN), CNE(MAXPTN),
     2     CNM(MAXPTN), CNK(MAXPTN), CNUSE(MAXPTN), CNSTR(MAXPTN), CNP
      INTEGER I, J, K, IDX(3), L, M
      DOUBLE PRECISION D, DUM

      IDX(1) = I
      IDX(2) = J
      IDX(3) = K
      DO 10 L = 1, 3
         M = IDX(L)
         ftp(L) = CNT(M)
         gxp(L) = CNX(M)
         gyp(L) = CNY(M)
         gzp(L) = CNZ(M)
         pxp(L) = CNPX(M)
         pyp(L) = CNPY(M)
         pzp(L) = CNPZ(M)
         pep(L) = CNE(M)
         pmp(L) = CNM(M)
   10 CONTINUE
      DUM = 0d0
      CALL LOCLDR(3, DUM)
      D = ( DSQRT((gxp0(1)-gxp0(2))**2 + (gyp0(1)-gyp0(2))**2
     1          + (gzp0(1)-gzp0(2))**2)
     2    + DSQRT((gxp0(1)-gxp0(3))**2 + (gyp0(1)-gyp0(3))**2
     3          + (gzp0(1)-gzp0(3))**2)
     4    + DSQRT((gxp0(2)-gxp0(3))**2 + (gyp0(2)-gyp0(3))**2
     5          + (gzp0(2)-gzp0(3))**2) ) / 3d0
      RETURN
      END

c=======================================================================

      SUBROUTINE COALNEW
*
*     Replaces `coales` when AMPT_COALNEW=1.
*
*     Structure: PTOH loops ISG = 1, NSG and skips NJSGS(ISG) = 0, reading
*     K2SGS(ISG,J) for J = 1, NJSGS(ISG). So the group layout is already dynamic
*     and the new algorithm can run as a re-partitioning: flatten the SGS arrays
*     into a parton list, decide the groupings, then rewrite the groups and NSG.
*     PTOH's species assignment by flavour then works unchanged, provided every
*     group is well formed (2 = one quark and one antiquark, 3 = three of one
*     sign).
*
*     Bookkeeping, which is where this goes wrong if done carelessly. With B
*     baryons, BB antibaryons and M mesons from NQ quarks and NA antiquarks, using
*     every parton exactly once:
*
*         3*B + M = NQ,   3*BB + M = NA   =>   B - BB = (NQ - NA)/3
*
*     so net-baryon number is conserved automatically, but ONLY if every parton is
*     used, and that requires M = NQ (mod 3). A greedy pass does not respect that
*     residue, so a repair step is mandatory. CNCHK enforces it and aborts rather
*     than letting a parton be dropped or double-counted.
*
      PARAMETER (MAXSTR=150001)
      PARAMETER (MAXPTN=400001)
      DOUBLE PRECISION PXSGS,PYSGS,PZSGS,PESGS,PMSGS,
     &     GXSGS,GYSGS,GZSGS,FTSGS
      COMMON/SOFT/PXSGS(MAXSTR,3),PYSGS(MAXSTR,3),PZSGS(MAXSTR,3),
     &     PESGS(MAXSTR,3),PMSGS(MAXSTR,3),GXSGS(MAXSTR,3),
     &     GYSGS(MAXSTR,3),GZSGS(MAXSTR,3),FTSGS(MAXSTR,3),
     &     K1SGS(MAXSTR,3),K2SGS(MAXSTR,3),NJSGS(MAXSTR)
      COMMON /orgsgs/ ISTRGS(MAXSTR,3)
      COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &     K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &     PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)
      COMMON /AREVT/ IAEVT, IARUN, MISS
      DOUBLE PRECISION CNX, CNY, CNZ, CNT, CNPX, CNPY, CNPZ, CNE, CNM
      COMMON /cnlist/ CNX(MAXPTN), CNY(MAXPTN), CNZ(MAXPTN), CNT(MAXPTN),
     1     CNPX(MAXPTN), CNPY(MAXPTN), CNPZ(MAXPTN), CNE(MAXPTN),
     2     CNM(MAXPTN), CNK(MAXPTN), CNUSE(MAXPTN), CNSTR(MAXPTN), CNP
*     the chosen groups: CNG(n,1..3) parton indices, CNGN(n) size
      COMMON /cngrp/ CNG(MAXSTR,3), CNGN(MAXSTR), CNNG
      COMMON /cnord/ CNORD
      DOUBLE PRECISION RBMHQ, RUSE
      COMMON /cnhq/ RBMHQ
      INTEGER ION, I, IL, J, K, L, NQ, NA, NB, NBB, NM, ISG, JJ
      DOUBLE PRECISION RBM, DM, DB, DTMP, DBEST
      INTEGER JM, JB1, JB2, ISIGN, NUSEDQ, NUSEDA
      INTEGER CNPERM(MAXPTN), ISW
      INTEGER*8 LSD

      CALL CNPAR(ION, RBM)
      IF (ION .NE. 1) THEN
         CALL COALES
         RETURN
      END IF

*     ---- 1. flatten the SGS groups into one parton list ----
      CNP = 0
      DO 20 ISG = 1, NSG
         IF (NJSGS(ISG) .EQ. 0) GOTO 20
         DO 15 J = 1, NJSGS(ISG)
            CNP = CNP + 1
            CNK(CNP)   = K2SGS(ISG, J)
            CNPX(CNP)  = PXSGS(ISG, J)
            CNPY(CNP)  = PYSGS(ISG, J)
            CNPZ(CNP)  = PZSGS(ISG, J)
            CNE(CNP)   = PESGS(ISG, J)
            CNM(CNP)   = PMSGS(ISG, J)
            CNX(CNP)   = GXSGS(ISG, J)
            CNY(CNP)   = GYSGS(ISG, J)
            CNZ(CNP)   = GZSGS(ISG, J)
            CNT(CNP)   = FTSGS(ISG, J)
            CNSTR(CNP) = ISTRGS(ISG, J)
            CNUSE(CNP) = 0
   15    CONTINUE
   20 CONTINUE

      NQ = 0
      NA = 0
      DO 25 I = 1, CNP
         IF (CNK(I) .GT. 0) THEN
            NQ = NQ + 1
         ELSE
            NA = NA + 1
         END IF
   25 CONTINUE

*     random sweep order (CNORD = 2): Fisher-Yates with a Park-Miller LCG
      IF (CNORD .EQ. 2) THEN
         DO 26 I = 1, CNP
            CNPERM(I) = I
   26    CONTINUE
         LSD = MOD(INT(IAEVT, 8) * 7919_8 + 104729_8, 2147483647_8)
         IF (LSD .LE. 0) LSD = LSD + 2147483646_8
         DO 27 IL = CNP, 2, -1
            LSD = MOD(LSD * 16807_8, 2147483647_8)
            J = 1 + INT(DBLE(LSD) / 2147483647d0 * DBLE(IL))
            IF (J .GT. IL) J = IL
            ISW = CNPERM(IL)
            CNPERM(IL) = CNPERM(J)
            CNPERM(J) = ISW
   27    CONTINUE
      END IF

*     ---- 2. greedy pass, He & Lin eq. 1 ----
      CNNG = 0
      DO 100 IL = 1, CNP
         I = IL
         IF (CNORD .EQ. 1) I = CNP + 1 - IL
         IF (CNORD .EQ. 2) I = CNPERM(IL)
         IF (CNUSE(I) .NE. 0) GOTO 100
         ISIGN = 1
         IF (CNK(I) .LT. 0) ISIGN = -1

*        closest opposite-sign partner -> d_M
         JM = 0
         DM = 1d30
         DO 30 J = 1, CNP
            IF (J .EQ. I .OR. CNUSE(J) .NE. 0) GOTO 30
            IF (CNK(J) * ISIGN .GT. 0) GOTO 30
            CALL CNDIS2(I, J, DTMP)
            IF (DTMP .LT. DM) THEN
               DM = DTMP
               JM = J
            END IF
   30    CONTINUE

*        closest same-sign partner, then the third giving the smallest d_B
         JB1 = 0
         DTMP = 1d30
         DO 40 J = 1, CNP
            IF (J .EQ. I .OR. CNUSE(J) .NE. 0) GOTO 40
            IF (CNK(J) * ISIGN .LT. 0) GOTO 40
            CALL CNDIS2(I, J, DBEST)
            IF (DBEST .LT. DTMP) THEN
               DTMP = DBEST
               JB1 = J
            END IF
   40    CONTINUE
         JB2 = 0
         DB = 1d30
         IF (JB1 .NE. 0) THEN
            DO 50 K = 1, CNP
               IF (K .EQ. I .OR. K .EQ. JB1 .OR. CNUSE(K) .NE. 0) GOTO 50
               IF (CNK(K) * ISIGN .LT. 0) GOTO 50
               CALL CNDIS3(I, JB1, K, DTMP)
               IF (DTMP .LT. DB) THEN
                  DB = DTMP
                  JB2 = K
               END IF
   50       CONTINUE
         END IF

*        decide. The criteria are not needed when only one option exists.
*        r_BM(charm) when the seed or either baryon partner is a c/cbar
         RUSE = RBM
         IF (IABS(NINT(CNK(I))) .EQ. 4) RUSE = RBMHQ
         IF (JB1 .NE. 0) THEN
            IF (IABS(NINT(CNK(JB1))) .EQ. 4) RUSE = RBMHQ
         END IF
         IF (JB2 .NE. 0) THEN
            IF (IABS(NINT(CNK(JB2))) .EQ. 4) RUSE = RBMHQ
         END IF
         IF (JB2 .NE. 0 .AND. (JM .EQ. 0 .OR. DB .LT. DM * RUSE)) THEN
            CNNG = CNNG + 1
            CNG(CNNG, 1) = I
            CNG(CNNG, 2) = JB1
            CNG(CNNG, 3) = JB2
            CNGN(CNNG) = 3
            CNUSE(I) = 1
            CNUSE(JB1) = 1
            CNUSE(JB2) = 1
         ELSE IF (JM .NE. 0) THEN
            CNNG = CNNG + 1
            CNG(CNNG, 1) = I
            CNG(CNNG, 2) = JM
            CNG(CNNG, 3) = 0
            CNGN(CNNG) = 2
            CNUSE(I) = 1
            CNUSE(JM) = 1
         END IF
*        if neither a meson nor a baryon partner exists the parton is left for
*        the repair step below
  100 CONTINUE

*     No repair step is needed, and that is a result rather than an assumption:
*     the partons come from melted hadrons, each meson giving q+qbar and each
*     baryon three quarks of one sign, so NQ - NA = 3*(net baryon number) is
*     divisible by 3. If the greedy pass leaves u unused quarks and consumes every
*     antiquark then 3B + M = NQ - u and 3Bbar + M = NA, whose difference forces
*     u = 0 (mod 3); and u <= 2, because three same-sign leftovers would have
*     formed a baryon, and the leftovers are all of one sign, because whichever of
*     two opposite-sign leftovers came later in the loop would have found the
*     earlier one free. Hence u = 0. CNCHK asserts it instead of trusting it.
      CALL CNCHK(NQ, NA)
      CALL CNWRIT

      RETURN
      END

c=======================================================================

      SUBROUTINE CNCHK(NQ, NA)
*
*     Gate 1: parton bookkeeping, as an assertion. A hadronisation rewrite that
*     silently drops or duplicates a parton would poison every observable
*     downstream, so this aborts rather than warns.
*
      PARAMETER (MAXSTR=150001)
      PARAMETER (MAXPTN=400001)
      DOUBLE PRECISION CNX, CNY, CNZ, CNT, CNPX, CNPY, CNPZ, CNE, CNM
      COMMON /cnlist/ CNX(MAXPTN), CNY(MAXPTN), CNZ(MAXPTN), CNT(MAXPTN),
     1     CNPX(MAXPTN), CNPY(MAXPTN), CNPZ(MAXPTN), CNE(MAXPTN),
     2     CNM(MAXPTN), CNK(MAXPTN), CNUSE(MAXPTN), CNSTR(MAXPTN), CNP
      COMMON /cngrp/ CNG(MAXSTR,3), CNGN(MAXSTR), CNNG
      COMMON /AREVT/ IAEVT, IARUN, MISS
      INTEGER NQ, NA, I, J, G, NB, NBB, NM, NUSED, IBAD
      INTEGER CNT2(MAXPTN)

      IBAD = 0
      DO 5 I = 1, CNP
         CNT2(I) = 0
    5 CONTINUE
      NB = 0
      NBB = 0
      NM = 0
      DO 20 G = 1, CNNG
         DO 10 J = 1, CNGN(G)
            I = CNG(G, J)
            IF (I .LT. 1 .OR. I .GT. CNP) THEN
               WRITE (6,*) 'COALNEW: group', G, 'has bad parton index', I
               IBAD = 1
               GOTO 10
            END IF
            CNT2(I) = CNT2(I) + 1
   10    CONTINUE
         IF (CNGN(G) .EQ. 2) THEN
            NM = NM + 1
            IF (CNK(CNG(G,1)) * CNK(CNG(G,2)) .GT. 0) THEN
               WRITE (6,*) 'COALNEW: meson group', G, 'is not q qbar'
               IBAD = 1
            END IF
         ELSE IF (CNGN(G) .EQ. 3) THEN
            IF (CNK(CNG(G,1)) .GT. 0) THEN
               NB = NB + 1
            ELSE
               NBB = NBB + 1
            END IF
            IF (CNK(CNG(G,1)) * CNK(CNG(G,2)) .LT. 0 .OR.
     1          CNK(CNG(G,1)) * CNK(CNG(G,3)) .LT. 0) THEN
               WRITE (6,*) 'COALNEW: baryon group', G, 'mixes signs'
               IBAD = 1
            END IF
         ELSE
            WRITE (6,*) 'COALNEW: group', G, 'has size', CNGN(G)
            IBAD = 1
         END IF
   20 CONTINUE

      NUSED = 0
      DO 30 I = 1, CNP
         IF (CNT2(I) .EQ. 1) THEN
            NUSED = NUSED + 1
         ELSE IF (CNT2(I) .EQ. 0) THEN
            WRITE (6,*) 'COALNEW: parton', I, 'flavour', CNK(I),
     1           'was never used (the u = 0 argument failed)'
            IBAD = 1
         ELSE
            WRITE (6,*) 'COALNEW: parton', I, 'used', CNT2(I), 'times'
            IBAD = 1
         END IF
   30 CONTINUE

      IF (NUSED .NE. CNP) IBAD = 1
      IF (3 * NB + NM .NE. NQ) THEN
         WRITE (6,*) 'COALNEW: 3B + M =', 3*NB+NM, ' but NQ =', NQ
         IBAD = 1
      END IF
      IF (3 * NBB + NM .NE. NA) THEN
         WRITE (6,*) 'COALNEW: 3Bbar + M =', 3*NBB+NM, ' but NA =', NA
         IBAD = 1
      END IF
      IF (3 * (NB - NBB) .NE. NQ - NA) THEN
         WRITE (6,*) 'COALNEW: net baryon number not conserved'
         IBAD = 1
      END IF

      IF (IBAD .NE. 0) THEN
         WRITE (6,*) 'COALNEW: bookkeeping failed in event', IAEVT,
     1        ' partons', CNP, ' NQ', NQ, ' NA', NA,
     2        ' groups', CNNG, ' B', NB, ' Bbar', NBB, ' M', NM
         STOP 'COALNEW bookkeeping failure'
      END IF
      RETURN
      END

c=======================================================================

      SUBROUTINE CNWRIT
*
*     Write the chosen groups back into the SOFT arrays and reset NSG, so that
*     PTOH sees exactly the layout it expects, only with different group sizes.
*
      PARAMETER (MAXSTR=150001)
      PARAMETER (MAXPTN=400001)
      DOUBLE PRECISION PXSGS,PYSGS,PZSGS,PESGS,PMSGS,
     &     GXSGS,GYSGS,GZSGS,FTSGS
      COMMON/SOFT/PXSGS(MAXSTR,3),PYSGS(MAXSTR,3),PZSGS(MAXSTR,3),
     &     PESGS(MAXSTR,3),PMSGS(MAXSTR,3),GXSGS(MAXSTR,3),
     &     GYSGS(MAXSTR,3),GZSGS(MAXSTR,3),FTSGS(MAXSTR,3),
     &     K1SGS(MAXSTR,3),K2SGS(MAXSTR,3),NJSGS(MAXSTR)
      COMMON /orgsgs/ ISTRGS(MAXSTR,3)
      COMMON/HJJET2/NSG,NJSG(MAXSTR),IASG(MAXSTR,3),K1SG(MAXSTR,100),
     &     K2SG(MAXSTR,100),PXSG(MAXSTR,100),PYSG(MAXSTR,100),
     &     PZSG(MAXSTR,100),PESG(MAXSTR,100),PMSG(MAXSTR,100)
      DOUBLE PRECISION CNX, CNY, CNZ, CNT, CNPX, CNPY, CNPZ, CNE, CNM
      COMMON /cnlist/ CNX(MAXPTN), CNY(MAXPTN), CNZ(MAXPTN), CNT(MAXPTN),
     1     CNPX(MAXPTN), CNPY(MAXPTN), CNPZ(MAXPTN), CNE(MAXPTN),
     2     CNM(MAXPTN), CNK(MAXPTN), CNUSE(MAXPTN), CNSTR(MAXPTN), CNP
      COMMON /cngrp/ CNG(MAXSTR,3), CNGN(MAXSTR), CNNG
      INTEGER G, J, I, NCLR

      NCLR = NSG
      IF (CNNG .GT. NCLR) NCLR = CNNG
      IF (NCLR .GT. MAXSTR) STOP 'COALNEW: more groups than MAXSTR'
      DO 20 G = 1, NCLR
         NJSGS(G) = 0
         DO 10 J = 1, 3
            K1SGS(G, J) = 0
            K2SGS(G, J) = 0
            PXSGS(G, J) = 0d0
            PYSGS(G, J) = 0d0
            PZSGS(G, J) = 0d0
            PESGS(G, J) = 0d0
            PMSGS(G, J) = 0d0
            GXSGS(G, J) = 0d0
            GYSGS(G, J) = 0d0
            GZSGS(G, J) = 0d0
            FTSGS(G, J) = 0d0
            ISTRGS(G, J) = 0
   10    CONTINUE
   20 CONTINUE

      DO 40 G = 1, CNNG
         NJSGS(G) = CNGN(G)
         DO 30 J = 1, CNGN(G)
            I = CNG(G, J)
            K2SGS(G, J) = CNK(I)
            PXSGS(G, J) = CNPX(I)
            PYSGS(G, J) = CNPY(I)
            PZSGS(G, J) = CNPZ(I)
            PESGS(G, J) = CNE(I)
            PMSGS(G, J) = CNM(I)
            GXSGS(G, J) = CNX(I)
            GYSGS(G, J) = CNY(I)
            GZSGS(G, J) = CNZ(I)
            FTSGS(G, J) = CNT(I)
            ISTRGS(G, J) = CNSTR(I)
   30    CONTINUE
   40 CONTINUE
      NSG = CNNG
      RETURN
      END
