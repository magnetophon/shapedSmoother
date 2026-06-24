import("stdfaust.lib");

// ---- user-defined lookahead sizing -------------------------------------
// One cycle of minFreq lasts SR/minFreq samples; the largest window must be at
// least that long, and the bank needs a power-of-two maxN, so:
//     maxN = 2 ^ ceil( log2(SR / minFreq) )
// (the -1e-9 keeps an exact power of two from rounding up a needless octave).
//
// NOTE: SR here is the *assumed runtime rate* used to turn minFreq into a fixed
// sample count at compile time — it is not ma.SR. It must match the rate you
// actually run at, or the Hz figures below shift.
//
// Defaults: SR = 48000, minFreq = 23.5 Hz
//     48000 / 23.5 = 2042.55 samples -> ceil(log2) = 11 -> maxN = 2048
//     2048 actually reaches 48000/2048 = 23.4375 Hz, so 23.5 Hz is fully covered.
//
// ---- bank @ SR=48000, minFreq=23.5 -> maxN=2048 : 12 outputs ------------
// Every output is a sliding-min over N samples, all sharing the trailing edge
// at t-(maxN-1); each reaches forward by N. Outputs come out smallest-window
// first. The window-1 branch IS the raw signal at that edge (min over 1 sample
// = identity), so there is no separate dry tap.
//   #    N  delay   win(ms)  covers >= Hz
//   1     1  2047    0.0208  48000.0000  <- raw signal at the shared edge
//   2     2  2046    0.0417  24000.0000
//   3     4  2044    0.0833  12000.0000
//   4     8  2040    0.1667   6000.0000
//   5    16  2032    0.3333   3000.0000
//   6    32  2016    0.6667   1500.0000
//   7    64  1984    1.3333    750.0000
//   8   128  1920    2.6667    375.0000
//   9   256  1792    5.3333    187.5000
//  10   512  1536   10.6667     93.7500
//  11  1024  1024   21.3333     46.8750
//  12  2048     0   42.6667     23.4375  <- minFreq window, full lookahead
// Edge sits at t-(maxN-1); system latency = maxN-1 = 2047 samp = 42.6458 ms.
// ------------------------------------------------------------------------
SR = 48000;
minFreq = 23.5;
maxN = 16;
// maxN = 1<<int(ceil(log(SR/minFreq)/log(2)-1e-9));

// Parallel sliding-reduce bank.
// For power-of-two n, emits the windows 1,2,4,…,n (ascending), each delayed by
// (n-windowSize) so every window shares its trailing edge at the output point
// and reaches forward by its size: 1@(n-1), 2@(n-2), … , (n/2)@(n/2), n@0.
// The 1@(n-1) branch is the raw signal at that edge.
slidingReducePar(op, n) = sequentialOperatorParOut(nBits-1, op)// w_0..w_(nBits-1): sizes 1,2,4,…,n
:par(i, nBits, _@(n-pow2(i)))// size pow2(i) -> delay n-pow2(i)
    with {
        nBits = int(floor(log(n)/log(2))+1);
        sequentialOperatorParOut(N, op) = seq(i, N, operator(i));
        operator(i) = si.bus(i), (_<:_, op(_, _@pow2(i)));
        pow2(i) = 1<<i;
    };

slidingMinPar(n) = slidingReducePar(min, n);

process = slidingMinPar(maxN);
// 12 outs at the defaults above
