rem RF model
fl32/c qrfv7a20.for
fl32/c subrsv2.for
fl32 qrfv7a20 subrsv2

rem RO model
fl32/c qro_v6a4.for
fl32/c lnsetup5q.for
fl32/c rv_sim2.for
fl32/c rdoutp6.for
link qro_v6a4 lnsetup5q rv_sim2 subrsv2 rdoutp6

rem RZ model
fl32 qrz_v3s7b.for subrsv2

