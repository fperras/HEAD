;HEAD: Hahn-echo assisted deconvolution
;
;2D experiment to measure the T2' lineshapes in a 1H fast-MAS experiment
;For secondary HEAD processing, the digital resolution in F1 and F2 needs to be identical
;
;
;
;parameters:
;pl1 : RF power level for 90/180
;p1 : 90 degree pulse
;p2 : =p1*2, 180 degree pulse
;d1 : recycle delay
;d6 : echo delay (calculated)
;d7 : echo delay (calculated)
;d0 : 0
;cnst31 : =MAS spin rate (or =1e6 for static)
;ns : 4*n
;
;FnMODE: States
;
;$CLASS=Solids
;$DIM=2D
;$TYPE=direct excitation
;$SUBTYPE=2D
;$COMMENT=Hahn-Echo experiment, MAS.

;"p2=p1*2"
"d6=(1s/cnst31)-(p1/2)-(p2/2)"
"d7=(1s/cnst31)-(p2/2)"
;cnst11 : to adjust t=0 for acquisition, if digmod = baseopt
"acqt0=0"

"in0=inf1/2"
"d0=0"

1 ze
2 10m 
  d1
  (p1 pl1 ph1):f1
  d6
  d0
  (p2 ph2):f1
  d7
  d0
  go=2 ph31
  10m wr #0 if #0 zd 
  1m if #0            ; increase FID number, that is, skip antiecho signal,
  1m id0              ; increment delay d0 by in0 (t1 increment),
  lo to 2 times td1
exit

ph0=0
ph1=0 1 2 3
ph2=0 0 0 0 1 1 1 1 2 2 2 2 3 3 3 3
ph30=0
ph31=0 3 2 1 2 1 0 3
