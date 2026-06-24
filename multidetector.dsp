import("stdfaust.lib");

// ========================================================================
//  multidetector : parallel sliding-min bank + combined variable reduce
// ========================================================================

// ---- user-defined controls ---------------------------------------------
SR           = 48000;    // assumed runtime rate (compile-time, NOT ma.SR)
minFreq      = 23.5;     // lowest freq the largest *single* block must span
gainIsLinear = 1;        // 1 = gains are linear, 0 = gains are in dB

// maxN = smallest power of two >= one cycle of minFreq (the -1e-9 stops an
// exact power of two rounding up a needless octave).
//   48000 / 23.5 = 2042.55 -> maxN = 2048  (single block reaches 23.4375 Hz)
maxN   = 1 << int(ceil(log(SR/minFreq)/log(2) - 1e-9));
maxLen = 2*maxN - 1;     // longest combined window (samples)

// Runtime window control, in Hz. n = samples for one full cycle of freq,
// ceil()'d so the window is always *at least* a full cycle, then clamped to
// [1, maxLen]. The slider floor SR/maxLen (~11.72 Hz) is the combined output's
// deepest reach -- a full octave below minFreq, available only at runtime.
//   high freq -> short window -> few samples ;  low freq -> long window.
freq = hslider("freq[unit:Hz][scale:log]", minFreq, SR/maxLen, SR/2, 0.01);
n    = int(min(maxLen, max(1, ceil(SR/freq))));

process = slidingMinPar(n, maxN, gainIsLinear);
slidingMinPar(n, maxN, lin) = slidingReducePar(min, n, maxN, lin);

// ========================================================================
// Parallel sliding-reduce bank + combined reduce.
//
// nBits power-of-two blocks (1,2,4,…,maxN) are computed once and fanned out to:
//   * parTaps  : nBits single-window minima (sizes 1..maxN)
//   * combined : one variable-size (1..2*maxN-1) reduce, slidingReduce-style
//
// maxLen = 1+2+…+maxN = 2*maxN-1 is the longest combined window, so combining
// almost DOUBLES the reach of a single block -> ~one octave lower:
//   single maxN=2048 -> 23.44 Hz   ;   combined 4095 -> 11.72 Hz.
//
// disabledVal = unity gain. We min() gains together, and a limiter only ever
// attenuates, so unity is the *largest* value a gain can take. Feeding unity
// into a min therefore can never pull the result down -> a disabled block has
// no effect. Unity is 1 in linear, 0 in dB, hence select2(lin, 0, 1).
//
// Alignment: every output shares its trailing edge at t-(maxLen-1) and reaches
// forward by its size, so the *max* combined size (n=maxLen) lands at 0 delay
// (reaching the present sample t). Taps get the constant delay maxLen-size;
// the combined gets the variable delay maxLen-n.
//
// ---- outputs @ SR=48000, minFreq=23.5 -> maxN=2048, maxLen=4095 ---------
//  out  size     delay    win(ms)  covers >= Hz
//   1      1      4094     0.0208   48000.0000
//   2      2      4093     0.0417   24000.0000
//   3      4      4091     0.0833   12000.0000
//   4      8      4087     0.1667    6000.0000
//   5     16      4079     0.3333    3000.0000
//   6     32      4063     0.6667    1500.0000
//   7     64      4031     1.3333     750.0000
//   8    128      3967     2.6667     375.0000
//   9    256      3839     5.3333     187.5000
//  10    512      3583    10.6667      93.7500
//  11   1024      3071    21.3333      46.8750
//  12   2048      2047    42.6667      23.4375
//  13   n(1..4095) maxLen-n  <=85.3125  >= 11.7216   <- combined, doubles range
// shared trailing edge t-4094 ; system latency = maxLen-1 = 4094 samp.
// A tap is disabled (-> unity) when its window is longer than n.
// ========================================================================
slidingReducePar(op, n, maxN, gainIsLinear) =
    sequentialOperatorParOut(nBits-1, op) <: parTaps , combined
with {
    nBits       = maxNrBits(maxN);          // # of blocks: 1,2,4,…,maxN
    maxLen      = (1 << nBits) - 1;          // 2*maxN-1, longest combined window
    disabledVal = select2(gainIsLinear, 0, 1);  // unity: 1 lin / 0 dB (see header)

    // single-window taps, aligned to the shared trailing edge; a tap longer
    // than the current n is switched off to unity (the size-based useVal).
    parTaps        = par(i, nBits, _ @ (maxLen - pow2(i)) : useValSize(i));
    useValSize(i)  = select2(pow2(i) <= n, disabledVal, _);

    // combined reduce: tile the blocks selected by the bits of n, then delay
    // the whole result by (maxLen-n) so n=maxLen lands at 0 delay.
    combined       = par(i, nBits, _ @ sumOfPrevBlockSizes(i) : useValBit(i))
                   : parallelOp(op, nBits)
                   : _ @ (maxLen - n);
    useValBit(i)   = select2(isUsed(i), disabledVal, _);   // the original (bit) useVal

    // ---- shared helpers (from slidingReduce) ----
    sequentialOperatorParOut(N, op) = seq(i, N, operator(i));
    operator(i) = si.bus(i), (_ <: _ , op(_, _ @ pow2(i)));

    sumOfPrevBlockSizes(0) = 0;
    sumOfPrevBlockSizes(i) = ba.subseq(allBlockSizes, 0, i) :> _;
    allBlockSizes          = par(i, maxNrBits(maxN-1), pow2(i) * isUsed(i));
    isUsed(i)              = ba.take(i+1, int2bin(n, maxLen));

    parallelOp(op, 1) = _;
    parallelOp(op, N) = op(parallelOp(op, N-1), _);

    int2bin(x, m)  = par(j, maxNrBits(m-1), int(floor(x / pow2(j))) % 2);
    maxNrBits(x)   = int(floor(log(x)/log(2)) + 1);
    pow2(i)        = 1 << i;
};
