c     D-meson framework switches (charm_mod.f); explicit types so the
c     layout is identical under any IMPLICIT rule of the including unit
      double precision sighqr
      real hqptmn, hqlqf
      integer ihqtrg, ihqcln, ihqdir, ichout, lqseed, idstar, ncolhj,
     &     ihqiso, ihqfrg, ihqmtm
      real hqpr, hqpm, hqmtmx, hqeps, hqpcut, hqfbar, hqtilt
      real hqecr, hqtstp
      integer ihqsout
      real hqlate
      real hqshad
      integer ihqtrx
c     src_v14: modern PDFs (AMPT_HQPDF) and Cronin broadening
      integer ihqpdf, ihqpst, ihqhrd, nhqcp, ihqdrf
      real hqmuf, hqcron, hqxsf, hqsab, hqvio
c     src_v15: nuclear modification of the triggered pair as an event
c     weight (AMPT_HQSHW) instead of inside the sampling
      integer ihqshw, khqw
c     src_v16: coalescence criterion variants (AMPT_HQCONST, AMPT_HQBDIST)
      integer ihqcon, ihqbd
c     src_v17: hadronic rescattering of charm hadrons in ART (AMPT_HQHAD)
      integer ihqhad, ihqhds, ihqhlc, nhqhc, ihqhrc, ihqhsd
      real hqsgm, hqsgb
      real hqwsh, hqwx, hqwq2, hqwr
      common /hqctl/ sighqr, hqptmn, hqlqf, ihqtrg, ihqcln, ihqdir,
     &     ichout, lqseed, idstar, ncolhj, ihqiso, ihqfrg, ihqmtm,
     &     hqpr, hqpm, hqmtmx, hqeps, hqpcut, hqfbar, hqtilt,
     &     hqecr, hqtstp, ihqsout, hqlate, hqshad, ihqtrx,
     &     ihqpdf, ihqpst, ihqhrd, nhqcp, hqmuf, hqcron, hqxsf,
     &     hqsab(2), hqvio, ihqdrf, ihqshw, hqwsh, khqw(2), hqwx(2),
     &     hqwq2, hqwr(2), ihqcon, ihqbd, ihqhad, ihqhds, ihqhlc,
     &     nhqhc(3), hqsgm, hqsgb, ihqhrc, ihqhsd
