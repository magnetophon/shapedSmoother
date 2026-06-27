declare name "lookahead-gr-smoother";
declare version "0.5";
declare author "Bart Brouns";
declare license "AGPL-3.0-only";
declare copyright "2026, Bart Brouns";
import("stdfaust.lib");

//============================ compile-time constants ============================
// These size par/seq bounds, so they MUST be compile-time literals. ma.SR is
// runtime — do NOT use it here.
SR      = 48000;
minFreq = 23.5;
maxN    = 1 << int(ceil(log2(SR / minFreq)));   // -> 2048   (2^11)
nBits   = int(floor(log2(maxN))) + 1;           // -> 12     (W = 1,2,...,2048)
look    = maxN;                                  // -> 2048   bank delay base; latency = look-1 = 2047
log2(x) = log(x) / log(2);                       // no ma.log2 dependency; keep int() on shifts
tau(k)  = pow2(k) / ma.SR;                       // W/SR seconds; pole = exp(-2*PI/W) at design rate
pow2(k) = 1 << k;

//============================ GUI groups + built-in test signal ==================
// MainGroup (hgroup) holds two columns: the Test signal controls and the Smoother
// controls. testSignal is generated here and fed straight into the bank as the input
// GR (so process has no audio input -- it is a self-contained generator for auditioning).
MainGroup(x)     = hgroup("[0]shapedSmoother", x);
TestGroup(x)     = vgroup("[0]Test signal", x);
SmootherGroup(x) = vgroup("[1]Smoother", x);
// Endpoint subgroups: the big block (largest window) and small block (smallest window)
// each get their own mult/order/self-correct. Every window in between interpolates.
BigGroup(x)      = SmootherGroup(hgroup("[1]Big block (largest window)",   x));
SmallGroup(x)    = SmootherGroup(hgroup("[2]Small block (smallest window)", x));

// --- Test signal ---
testNoiseLevel = TestGroup(hslider("[0]noise level", 0, 0, 1, 0.001));
testNoiseRate  = TestGroup(hslider("[1]noise rate", 42, 1, 1000, 1));
testBlockscale = TestGroup(hslider("[2]blockscale", 1, 0.01, 10, 0.01));
testFreq       = TestGroup(hslider("[3]freq", 1, 0.001, 30, 0.001));
testStep1      = TestGroup(hslider("[4]step1", 0.75, -1, 1, 0.001));
testStep2      = TestGroup(hslider("[5]step2", 0.125, -1, 1, 0.001));
testSelect     = TestGroup(checkbox("[6]signal select"));

testSignal = select2(testSelect, testSignal1, testSignal2) * 0.5 + 0.5;

testSignal1 = it.interpolate_linear(testNoiseLevel,
    (loop ~ _),
    no.lfnoise(testNoiseRate))
with {
    loop(prev) = no.lfnoise0(testBlockscale * (abs(prev*69)%9 : pow(0.75)*5+1));
};
testSignal2 = os.lf_squarewave(testFreq) * 0.5;

//============================ process ===========================================
// in : raw GR (linear gain <= 1).  out : (delayed raw GR, smoothed GR), aligned.
//
// `n` is now a runtime gate driven by a frequency control. The largest ACTIVE
// window is ~SR/freq samples, so freq sets the smoothing / lookahead corner:
// lower freq -> larger windows active -> more lookahead. The dyadic ladder
// itself stays fixed at compile time (maxN); the slider only chooses which
// rungs participate (see useValSize / `pow2(i) <= n`). Bounds are compile-time
// literals (SR, maxN); the Hz->samples conversion uses ma.SR so the control is
// a true frequency at run time. The window *times* are fixed at the SR=48000
// design rate, so the freq<->window mapping is exact at 48k (swap ma.SR -> SR
// in nWin to make the gating sample-rate-invariant instead).
fMin = SR / maxN;                              // ~23.44 Hz : W=2048 corner, all windows on
fMax = SR / 2;                                 // Nyquist    : only W=1,2 on (least smoothing)
freq = SmootherGroup(hslider("[0]smoothing corner[unit:Hz][scale:log]", fMin, fMin, fMax, 0.01));
nWin = ma.SR / freq;                           // Hz -> samples

// DEBUG BUILD. in : built-in testSignal (above).  out : (delayed raw GR, smoothed GR),
// look-1=2047. The whole graph is wrapped in MainGroup so the Test signal and Smoother
// control columns both sit under the "shapedSmoother" group in the GUI.
//   ch1  Output A   : t0, the W=1 tap delayed by look-1 (linear, unchanged).
//   ch2  smoothed GR: the FULL self-correcting dB fold over t1..t11, back to linear.
process = MainGroup(testSignal : lin2dB : bank : combine);

// lin2dB now runs ONCE up front: min commutes with linear2db (monotonic), so a sliding-min
// in linear then per-tap dB == convert once then sliding-min in dB. disabledVal flips 1 -> 0
// so a gated-off window is unity in the NEW domain (0 dB instead of linear 1.0).
bank = slidingReducePar(min, nWin, maxN, 0);   // n = nWin gates which windows are active; 0 = dB-domain unity

//============================ self-correcting dB fold (ALL blocks) ===============
// Every window in the dyadic ladder is now the SAME kind of block: a self-correcting
// dB step with its own (mult, order, self-correct). The fold runs IN dB, threaded
// LARGEST-WINDOW-FIRST. Each block is fed (its dB residue - sc * its own prev output);
// because residue - prev_output == target_dB - combined_output_dB, each block drives
// the running COMBINED output toward THAT window's target rather than stacking blindly.
// Residue is clamped at 0 dB (reduce-only): a block may only duck, never boost, so the
// combined gain stays <= unity and the startup ring is killed.
//
// PER-BLOCK PARAMETERS. The largest window (fold stage 0) takes the Big-block sliders,
// the smallest window (stage nBits-2) takes the Small-block sliders, and every stage in
// between is a LINEAR interpolation of the two endpoints (LinArray, indexed by stage):
//   mult  : attack (duck) slowdown, x base grAttack. release is left at base.
//   order : smoother order. Compiled at 8 stages for every block (seq size is
//           compile-time); `order` gates how many are active poles via (i < order).
//           CPU is fixed at 8 poles/block regardless of the slider; interpolated
//           values are fractional (the pole count steps at integer boundaries).
//   sc    : self-correct amount [0,1]. 0 = plain dB fold; 1 = full self-correction
//           (loop DC gain 1/2, under-reaches at the peak); ~0.6 splits the difference.
// Defaults reproduce the previous build at the endpoints: big = 8 / 1.5x / 0 sc,
// small = 4 / 1.0x / 0.6 sc.
bigMult    = BigGroup  (hslider("[0]mult [x base att]", 1.5, 0.5, 4.0, 0.05));
bigOrder   = BigGroup  (hslider("[1]order",             8,   1,   8,   1   ));
bigSC      = BigGroup  (hslider("[2]self-correct",      0.0, 0.0, 1.0, 0.01));
smallMult  = SmallGroup(hslider("[0]mult [x base att]", 1.0, 0.5, 4.0, 0.05));
smallOrder = SmallGroup(hslider("[1]order",             4,   1,   8,   1   ));
smallSC    = SmallGroup(hslider("[2]self-correct",      0.6, 0.0, 1.0, 0.01));

// 12 -> 2 : (t0, G). t0 is now a dB tap -> db2linear back to linear for Output A.
combine     = (ba.db2linear, (foldChainDB : ba.db2linear));

// PARAMETER DELIVERY VIA ro.interleave (param-arrays-and-interleave pattern).
// Instead of each stage rebuilding the full ramp and selecting its own element, we
// build each parameter ONCE as an N-wide array (interpolated across stages), temporally
// smooth it (zipper-free slider moves), then transpose all of them — plus the incoming
// tap bus, treated as a 4th data column — from parameter-major to element-major so each
// stage is handed exactly its (mult, order, sc, tap) tuple.
//
//   element axis  : LinArray spaces small-endpoint -> big-endpoint across the stages.
//   time axis     : si.smooth(paramPole) per element, so live slider moves don't zip.
//   routing       : ro.interleave(N, K) transposes K arrays of N into N tuples of K.
//
// Column order in == (mult[0..N-1], order[0..N-1], sc[0..N-1], tap[0..N-1]).
// Element-major out == bundle r = (mult_r, order_r, sc_r, t_{r+1}).  The fold consumes
// bundles from the BACK, so bundle N-1 (= big endpoint, paired with t11) is folded first.
N           = nBits - 1;                               // 11 stages / array width
paramPole   = ba.tau2pole(0.02);                       // 20 ms control smoothing
smoothBus   = par(i, N, si.smooth(paramPole));
multArray   = LinArray(smallMult,  bigMult,  N) : smoothBus;   // elem 0 = small, elem N-1 = big
orderArray  = LinArray(smallOrder, bigOrder, N) : smoothBus;
scArray     = LinArray(smallSC,    bigSC,    N) : smoothBus;

// in: t1..t11 (the bus that becomes the tap column).  out: gdB.
foldChainDB = (multArray, orderArray, scArray, si.bus(N))  // K=4 columns, each N wide
            : ro.interleave(N, 4)                          // -> N bundles of (mult,ord,sc,tap)
            : (si.bus(N * 4), 0.0)                         // append 0 dB backstop at the back
            : foldDB;                                      // 4N+1 -> 1

// N = 11 stages. stage(0) folds t11 (largest), stage(N-1) folds t1 (smallest). Each stage
// peels the last 4-tuple (its interleaved params + tap) plus the running gain off the back.
foldDB   = seq(j, N, stage(j));
stage(j) = (si.bus(4 * (N - 1 - j)), step(N - j));

// One self-correcting dB block. k = window exponent (timing), compile-time (0 = largest).
//   signal inputs (mult, ord, sc, t, g): mult/ord/sc = this stage's interleaved params,
//   t = this window's dB min tap, g = running combined gain in dB.
step(k, mult, ord, sc, t, g) = g + sres
with {
    resDB = min(0.0, t - g);                           // t already in dB; residue (<= 0 dB)
    sres  = resDB : (loop ~ _);                        // self-correcting block output
    loop(fb, x) = min(0.0, smootherARorder(8, ord, ord,
                          mult * grAttack * tau(k), grRelease * tau(k),
                          x - sc * fb));
};

// dB fold primitive. lin2dB floors the gain before the log so a 0 (full mute) tap maps
// to a finite, very low dB instead of -inf.
lin2dB(x) = ba.linear2db(max(ma.EPSILON, x));

// GR ATTACK / RELEASE timing. smootherARorder is wired for gain reduction (falling
// input = attack), so att = the duck and rel = the recovery. Both scale tau(k) per
// window, so one value sets the timing across every window (lookahead and smoother
// time both ~ 2^k). mult (per block) further scales the attack only.
//   grAttack : the duck (falling edge). 1.0 = frequency-matched smoother, which
//              completes ~0.3*lookahead too early; 1.45 lands the duck right at
//              Output A's drop (order 4). Order-coupled via cutoffCorrection.
//   grRelease: the recovery (rising edge). 1.0 = frequency-matched; raise for gentler.
grAttack  = 1.45;
grRelease = 1.45;

//============================ library functions (verbatim) ======================
// Linear ramp of nrElements points from bottom to top inclusive (user-provided).
LinArray(bottom,top,0) =   0:! ;
LinArray(bottom,top,nrElements) =     par(i,nrElements,   ((top-bottom)*(i/(nrElements-1)))+bottom);

// Custom variant of basics.lib slidingReduce — see LOCKED #11. Exposes aligned
// parallel taps (parTaps) instead of reducing; delay scheme is look - pow2(i).
// Intentional; do NOT reconcile with the stock definition.
slidingReducePar(op, n, maxN, gainIsLinear) =
    sequentialOperatorParOut(nBits - 1, op) : parTaps
with {
    nBits       = maxNrBits(maxN);
    look        = maxN;
    disabledVal = gainIsLinear;
    parTaps     = par(i, nBits, _@(look - pow2(i)) : useValSize(i));
    useValSize(i) = select2(pow2(i) <= n, disabledVal, _);
    sequentialOperatorParOut(N, op) = seq(i, N, operator(i));
    operator(i) = si.bus(i), (_ <: _, op(_, _@pow2(i)));
    maxNrBits(x) = int(floor(log(x) / log(2)) + 1);
    pow2(i) = 1 << i;
};

smootherARorder(maxOrder, orderAtt, orderRel, att, rel, xx) =
    xx : seq(i, maxOrder, loop(i) ~ _)
with {
    loop(i, fb, x) = coeff(i) * fb + (1.0 - coeff(i)) * x
    with {
        cutoffCorrection(order) = 1.0 / sqrt(pow(2.0, 1.0 / order) - 1.0);
        // Wired for GAIN REDUCTION: GR ducks DOWNWARD on attack, so a FALLING
        // input (x < fb) selects attCoeff (the att param) and a RISING input
        // selects relCoeff (rel param). This is the swap that makes att = GR
        // attack and rel = GR release. (A standard amplitude-envelope follower
        // would use x > fb here instead.)
        coeff(i) = ba.if(x < fb, attCoeff(i), relCoeff(i));
        attCoeff(i) = exp(-TWOPIT * cutoffCorrection(orderAtt) / max(ma.EPSILON, att)) * (i < orderAtt);
        relCoeff(i) = exp(-TWOPIT * cutoffCorrection(orderRel) / max(ma.EPSILON, rel)) * (i < orderRel);
        TWOPIT = 2 * ma.PI * ma.T;
    };
};
