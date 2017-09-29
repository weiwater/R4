rem RF module
rem
    fl32/c rf_v7a.for
    fl32/c subrs.for
    fl32 rf_v7a subrs
rem 
rem RO module
rem 
    fl32/c ro_v6s.for
    fl32/c rdoutp4.for
    fl32/c lnsetup4.for
    fl32/c rv_sim4.for
    fl32 ro_v6s rdoutp4 lnsetup4 rv_sim4 subrs
rem
rem RZ module
rem
    fl32 rz_v3s2.for subrs
rem
