declare name "lookahead-gr-smoother";
declare version "0.1";
declare author "Bart Brouns";
declare license "AGPL-3.0-only";
declare copyright "2026, Bart Brouns";
import("stdfaust.lib");

//============================ compile-time constants ============================
// These size par/seq bounds, so they MUST be compile-time literals. ma.SR is
// runtime — do NOT use it here.
SR = 48000;
minFreq = 23.5;
maxN = 1<<int(ceil(log2(SR/minFreq)));
// -> 2048   (2^11)
nBits = int(floor(log2(maxN)))+1;
// -> 12     (W = 1,2,...,2048)
look = maxN;
// -> 2048   bank delay base; latency = look-1 = 2047
log2(x) = log(x)/log(2);
// no ma.log2 dependency; keep int() on shifts
tau(k) = pow2(k)/ma.SR;
// W/SR seconds; pole = exp(-2*PI/W) at design rate
pow2(k) = 1<<k;

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
fMin = SR/maxN;
// ~23.44 Hz : W=2048 corner, all windows on
fMax = SR/2;
// Nyquist    : only W=1,2 on (least smoothing)
freq = hslider("[0]smoothing corner[unit:Hz][scale:log]", fMin, fMin, fMax, 0.01);
nWin = ma.SR/freq;
// Hz -> samples

// DEBUG BUILD. out = (delayed raw GR, fold of 1 biggest, fold of 2 biggest), look-1=2047.
// The fold now runs IN dB: residue = tap_dB - product_dB, combine = add in dB, residue
// clamped at 0 dB (a stage may only reduce, never boost -> output gain stays <= unity).
//   ch1  Output A   : t0, the W=1 tap delayed by look-1 (linear, unchanged).
//   ch2  fold of 1  : LARGEST window t11 (W=2048) vs the 0 dB base -> G1 (back to linear).
//   ch3  fold of 2  : fold W=1024 into the product -> G2 (back to linear). Threads
//        largest-window-first (LOCKED #9). The 0 dB clamp is the dB form of res<=1: it
//        keeps gain <= unity and removes the startup blowup an unclamped fold would show.
// t10/t11 are gated by the freq slider (active at the default fMin). For always-on taps
// regardless of the slider, drive bank with fixed n=maxN instead of nWin.
process = testSignal:bank:pickTaps:(_, foldTwo);

bank = slidingReducePar(min, nWin, maxN, 1);
// n = nWin gates which windows are active

// keep t0, t11, t10 from the 12-tap bank (drop the 9 middle taps). The route swaps the
// last two so the fold sees them largest-first: (t0, t11, t10).
pickTaps = _, par(i, nBits-3, !), route(2, 2, 1, 2, 2, 1);
// 12 -> 3 : (t0, t11, t10)

// ---- self-correcting dB fold (two biggest blocks) ----
// The big block runs slow + high-order so it deliberately under-ducks, leaving the small block
// real work. The small block is SELF-CORRECTING: it is fed (residue - scAmt*prev_output).
// Note residue - prev_output == target_dB - G2dB, i.e. the COMBINED tracking error, so the loop
// drives the combined output toward the target instead of blindly stacking onto the big block.
//   bigOrder  : big-block smoother order (the small block stays order 4).
//   bigAttMul : big-block attack (the duck) slowdown, x the base grAttack.
//   scAmt     : self-correction strength, in [0,1].
//               0   = plain fold (the small block can overshoot the floor under a slow big block);
//               1   = full self-correction (no overshoot, but the loop DC gain is 1/2 so the
//                     approach under-reaches the floor at the peak);
//               ~0.6 splits the difference (lands the approach near the floor, little overshoot).
bigOrder = 8;
bigAttMul = 1.5;
scAmt = 0.6;

foldTwo(t11, t10) = ba.db2linear(G1dB), ba.db2linear(G2dB)
    with {
        G1dB = stepDBbig(nBits-1, lin2dB(t11), 0.0);
        // big block: slow, high order
        G2dB = stepDBsc(nBits-2, lin2dB(t10), G1dB);
        // small block: self-correcting
    };

// big block: plain dB step at its own order and attack multiplier (release stays x1).
stepDBbig(k, tdB, gdB) = gdB+smootherARorder(bigOrder,
    bigOrder,
    bigOrder,
    bigAttMul*grAttack*tau(k),
    grRelease*tau(k),
    min(0.0, tdB-gdB));

// small block: self-correcting smoother. loop(fb, x) feeds (x - scAmt*fb) to the order-4
// smoother and clamps the result at 0 dB (reduce-only; this also tames the loop's startup ring).
// fb is the previous loop output -- the ~ delays it one sample.
stepDBsc(k, tdB, gdB) = gdB+sres
    with {
        resDB = min(0.0, tdB-gdB);
        // this block's residue (<= 0 dB)
        sres = resDB:(loop~_);
        // self-correcting small-smoother output
        loop(fb, x) = min(0.0,
            smootherARorder(4,
                4,
                4,
                grAttack*tau(k),
                grRelease*tau(k),
                x-scAmt*fb));
    };

// dB fold primitive. lin2dB floors the gain before the log so a 0 (full mute) tap maps to a
// finite, very low dB instead of -inf.
lin2dB(x) = ba.linear2db(max(ma.EPSILON, x));

// Bus widths traced: 12 -> 2. t0 is Output A only; the smoothed product runs over
// t1..t11 with a constant-1 backstop (see #7, #9). The `,` partition splits inputs
// left-to-right by operand arity. faust -wall confirms the arities.
combine = (_, foldChain);
// 12 -> 2 : (t0, G)

// Inject the constant-1 base at the BACK, then fold largest-window-first.
foldChain = (si.bus(nBits-1), 1.0):fold;
// 11 -> 12 -> 1 : (t1..t11, 1) -> G

fold = seq(j, nBits-1, stage(j));
// 12 -> 1 ; largest window first

// Pass the front (smaller-window) taps; fold the back tap into the running product.
// stage(0) folds t11 (slowest smoother) with g=1; stage(10) folds t1 last.
stage(j) = (si.bus(nBits-2-j), step(nBits-1-j));

// step folds one window into the running product. Inputs are (target, gain):
//   target = this window's min t_k ;  gain = running product G, threaded from the
//   previous (larger) window. Do NOT swap to (gain, target): G lives at the BACK
//   of the fold bus, so it is step's 2nd input. (LOCKED #9.)
//
// GR ATTACK / RELEASE timing. smootherARorder above is wired for gain reduction
// (falling input = attack), so here att = the duck and rel = the recovery, named
// correctly at last. Both scale tau(k) per window, so one value sets the timing
// across every window (lookahead and smoother time both ~ 2^k).
//   grAttack : the duck (falling edge). 1.0 = frequency-matched smoother, which
//              completes ~0.3*lookahead too early; 1.45 lands the duck right at
//              Output A's drop (order 4). Order-coupled: cutoffCorrection shifts
//              the on-time value if you change the smoother order.
//   grRelease: the recovery (rising edge). 1.0 = frequency-matched; raise for a
//              gentler, longer recovery.
grAttack = 1.45;
grRelease = 1.45;
// min(1.0, ...) clamps the residue t_k/G <= 1 so a fold stage only REDUCES the running
// product, never lifts it. Resolves the open clamp decision (ON): it kills the intermediate-
// fold startup blowup (residue t_k/G with G~0) and preserves the larger window's anticipation.
step(k, t, g) = g*smootherARorder(4, 4, 4, grAttack*tau(k), grRelease*tau(k), min(1.0, t/divGuard(g)))
    with {
        divGuard(x) = max(ma.EPSILON, x);
        // raw GR can be 0 -> guard the divide (NaN otherwise)
    };

//============================ library functions (verbatim) ======================
// Custom variant of basics.lib slidingReduce — see LOCKED #11. Exposes aligned
// parallel taps (parTaps) instead of reducing; delay scheme is look - pow2(i).
// Intentional; do NOT reconcile with the stock definition.
slidingReducePar(op, n, maxN, gainIsLinear) = sequentialOperatorParOut(nBits-1, op):parTaps
    with {
        nBits = maxNrBits(maxN);
        look = maxN;
        disabledVal = gainIsLinear;
        parTaps = par(i, nBits, _@(look-pow2(i)):useValSize(i));
        useValSize(i) = select2(pow2(i)<=n, disabledVal, _);
        sequentialOperatorParOut(N, op) = seq(i, N, operator(i));
        operator(i) = si.bus(i), (_<:_, op(_, _@pow2(i)));
        maxNrBits(x) = int(floor(log(x)/log(2))+1);
        pow2(i) = 1<<i;
    };

smootherARorder(maxOrder, orderAtt, orderRel, att, rel, xx) = xx:seq(i, maxOrder, loop(i)~_)
    with {
        loop(i, fb, x) = coeff(i)*fb+(1.0-coeff(i))*x
            with {
                cutoffCorrection(order) = 1.0/sqrt(pow(2.0, 1.0/order)-1.0);
                // Wired for GAIN REDUCTION: GR ducks DOWNWARD on attack, so a FALLING
                // input (x < fb) selects attCoeff (the att param) and a RISING input
                // selects relCoeff (rel param). This is the swap that makes att = GR
                // attack and rel = GR release. (A standard amplitude-envelope follower
                // would use x > fb here instead.)
                coeff(i) = ba.if(x<fb, attCoeff(i), relCoeff(i));
                attCoeff(i) = exp(-TWOPIT*cutoffCorrection(orderAtt)/max(ma.EPSILON, att))*(i<orderAtt);
                relCoeff(i) = exp(-TWOPIT*cutoffCorrection(orderRel)/max(ma.EPSILON, rel))*(i<orderRel);
                TWOPIT = 2*ma.PI*ma.T;
            };
    };

MainGroup(x) = hgroup("[0]shapedSmoother", x);
TestGroup(x) = vgroup("[0]Test signal", x);
SmootherGroup(x) = vgroup("[1]Smoother", x);

// --- Test signal ---
testNoiseLevel = TestGroup(hslider("[0]noise level", 0, 0, 1, 0.001));
testNoiseRate = TestGroup(hslider("[1]noise rate", 42, 1, 1000, 1));
testBlockscale = TestGroup(hslider("[2]blockscale", 1, 0.01, 10, 0.01));
testFreq = TestGroup(hslider("[3]freq", 1, 0.001, 30, 0.001));
testStep1 = TestGroup(hslider("[4]step1", 0.75, -1, 1, 0.001));
testStep2 = TestGroup(hslider("[5]step2", 0.125, -1, 1, 0.001));
testSelect = TestGroup(checkbox("[6]signal select"));

testSignal = select2(testSelect, testSignal1, testSignal2)*0.5+0.5;

testSignal1 = it.interpolate_linear(testNoiseLevel,
    (loop~_),
    no.lfnoise(testNoiseRate))
    with {
        loop(prev) = no.lfnoise0(testBlockscale*(abs(prev*69)%9:pow(0.75)*5+1));
    };
testSignal2 = os.lf_squarewave(testFreq)*0.5;
