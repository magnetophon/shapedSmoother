import("stdfaust.lib");

// ========================================================================
//  multidetector : parallel sliding-min bank + combined + deepWin
// ========================================================================

// ---- user-defined controls ---------------------------------------------
SR = 48000;
// assumed runtime rate (compile-time, NOT ma.SR)
minFreq = 23.5;
// lowest freq the largest *single* block must span
gainIsLinear = 1;
// 1 = gains are linear, 0 = gains are in dB

// maxN = smallest power of two >= one cycle of minFreq (the -1e-9 stops an
// exact power of two rounding up a needless octave).
//   48000 / 23.5 = 2042.55 -> maxN = 2048  (single block reaches 23.4375 Hz)
maxN = 1<<int(ceil(log(SR/minFreq)/log(2)-1e-9));
look = 2*maxN;
// full lookahead window (= deepWin at max size)

// Runtime size control, in Hz. n = samples for one full cycle of freq, ceil()'d
// so the window is always *at least* a full cycle. The slider spans from
// "half the lookahead" (n = maxN, ~SR/maxN Hz) up to Nyquist (n = 2). deepWin
// then fills the other half -> reaches SR/look (~11.72 Hz), an octave lower.
//   high freq -> short window -> few samples ;  low freq -> long window.
freq = hslider("freq[unit:Hz][scale:log]", SR/maxN, SR/maxN, SR/2, 0.01);
n = int(min(maxN, max(1, ceil(SR/freq))));

process = slidingMinPar(n, maxN, gainIsLinear);
slidingMinPar(n, maxN, lin) = slidingReducePar(min, n, maxN, lin);

// ========================================================================
// Parallel sliding-reduce bank + combined + deepWin.
//
// nBits power-of-two blocks (1,2,4,…,maxN) are computed once and fanned out to:
//   * parTaps  : nBits single-window minima (sizes 1..maxN)
//   * combined : one variable-size (1..maxN) reduce, slidingReduce-style
//   * deepWin  : op(combined, combined@n) -> a 2n window (combined flows in),
//                the longest window, reaching one octave below combined.
//
// look = 2*maxN is deepWin's max length, so deepWin almost DOUBLES combined's
// reach -> ~one octave lower:
//   combined max maxN=2048 -> 23.44 Hz   ;   deepWin max 4096 -> 11.72 Hz.
//
// disabledVal = unity gain. We min() gains together, and a limiter only ever
// attenuates, so unity is the *largest* value a gain can take. Feeding unity
// into a min therefore can never pull the result down -> a disabled block has
// no effect. Unity is 1 in linear, 0 in dB, hence select2(lin, 0, 1).
//
// Alignment: every output shares its trailing edge at t-(look-1) and reaches
// forward by its size, so the *max* deepWin size (n=maxN -> 2n=look) lands at
// 0 delay (reaching the present sample t). A window of size S gets delay
// look-S: taps -> look-pow2(i) (constant), combined -> look-n, deepWin -> look-2n.
//
// ---- outputs @ SR=48000, minFreq=23.5 -> maxN=2048, look=4096 -----------
//  out  size        delay     win(ms)  covers >= Hz
//   1      1         4095      0.0208   48000.0000
//   2      2         4094      0.0417   24000.0000
//   3      4         4092      0.0833   12000.0000
//   4      8         4088      0.1667    6000.0000
//   5     16         4080      0.3333    3000.0000
//   6     32         4064      0.6667    1500.0000
//   7     64         4032      1.3333     750.0000
//   8    128         3968      2.6667     375.0000
//   9    256         3840      5.3333     187.5000
//  10    512         3584     10.6667      93.7500
//  11   1024         3072     21.3333      46.8750
//  12   2048         2048     42.6667      23.4375
//  13  n   (2..2048) look-n   <=42.6667  >= 23.4375   <- combined
//  14  2n  (4..4096) look-2n  <=85.3333  >= 11.7188   <- deepWin (octave below)
// shared trailing edge t-4095 ; system latency = look-1 = 4095 samp.
// A tap is disabled (-> unity) when its window is longer than n.
// ========================================================================
slidingReducePar(op, n, maxN, gainIsLinear) = sequentialOperatorParOut(nBits-1, op)<:parTaps, (combined<:alignedCombined, deepWin)
    with {
        nBits = maxNrBits(maxN);
        // # of blocks: 1,2,4,…,maxN
        look = 2*maxN;
        // deepWin max length / alignment ref
        disabledVal = gainIsLinear;
        // unity: 1 lin / 0 dB (see header)

        // single-window taps, aligned to the shared trailing edge; a tap longer
        // than the current n is switched off to unity (the size-based useVal).
        parTaps = par(i, nBits, _@(look-pow2(i)):useValSize(i));
        useValSize(i) = select2(pow2(i)<=n, disabledVal, _);

        // combined reduce of variable size n, leading edge still at t (it fans out
        // to its own aligned output AND into deepWin).
        combined = par(i, nBits, _@sumOfPrevBlockSizes(i):useValBit(i)):parallelOp(op, nBits);
        useValBit(i) = select2(isUsed(i), disabledVal, _);
        // the bit-based useVal

        // align combined to the shared trailing edge.
        alignedCombined = _@(look-n);

        // deepWin: double combined to a 2n window, then align so n=maxN -> 0 delay.
        deepWin = (_<:op(_, _@n)):_@(look-2*n);

        // ---- shared helpers (from slidingReduce) ----
        sequentialOperatorParOut(N, op) = seq(i, N, operator(i));
        operator(i) = si.bus(i), (_<:_, op(_, _@pow2(i)));

        sumOfPrevBlockSizes(0) = 0;
        sumOfPrevBlockSizes(i) = ba.subseq(allBlockSizes, 0, i):>_;
        allBlockSizes = par(i, maxNrBits(maxN-1), pow2(i)*isUsed(i));
        isUsed(i) = ba.take(i+1, int2bin(n));

        parallelOp(op, 1) = _;
        parallelOp(op, N) = op(parallelOp(op, N-1), _);

        int2bin(x) = par(j, nBits, int(floor(x/pow2(j)))%2);
        maxNrBits(x) = int(floor(log(x)/log(2))+1);
        pow2(i) = 1<<i;
    };
