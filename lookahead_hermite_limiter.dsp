//============================================================================
// lookahead_hermite_limiter.dsp  —  ~50 ms lookahead
// Lookahead comp-limiter, parallel multi-scale Hermite. Linear gain domain.
//
// maxN is now ~50 ms. Faust requires structural sizes (delay lengths, the
// parallel window bank) to be COMPILE-TIME constants, so the lookahead is fixed
// in samples at a design SR: L = round(maxMs*1e-3*SRdesign) = 2400 @ 48 kHz.
// Lookahead-in-ms = L / runtimeSR, so at 44.1 kHz this is ~54 ms, at 96 kHz
// ~25 ms. Change `maxMs` (and rebuild) to retune; keep L NOT a power of two
// (2400 is fine) so the exact-maxN top window is distinct from the dyadic bank.
//
// SCALING NOTES (what changed vs the 16-sample version):
//  1. Window bank is GENERATED, not hand-written: dyadic reaches 2^0..2^nTop
//     plus the exact maxN=L on top (13 windows at L=2400).
//  2. The nested windows are computed as ONE prefix scan over DISJOINT blocks
//     carrying (min,pos) together, so total cost is O(L), not O(sum of W).
//     This is what makes it compile (~35 s, ~7.5k lines of C++) instead of
//     blowing up. `dist` (nearest tap, tie -> nearest) falls out of the scan.
//  3. Hermite tangents are now MONOTONICITY-CLAMPED (Fritsch-Carlson:
//     |m| <= 3|p1-p0|, same sign as the secant). At maxN=16 the raw per-span
//     tangent m = vel*d was bounded; at d up to 2400 it is not — an unclamped
//     far window commits with m0 = v*d ~ hundreds and the cubic runs away.
//     Clamping keeps slope-matching where it is geometrically sane, forbids
//     overshoot/undershoot, guarantees exact landing, and kills the runaway.
//     (To recover the raw spec behavior at short L, make clampM the identity.)
//
// BEHAVIOUR (verified, Faust 2.70.3 @ 48 kHz):
//  * Lone 0.95 spike (thr 0.5): gain eases monotonically from unity starting the
//    instant the peak enters the 50 ms window and lands EXACTLY on 0.5/0.95 =
//    0.5263 the sample the peak hits the output, then releases. Smooth ease-in.
//  * 0.8 burst: settles EXACTLY on 0.625 (out = 0.5). DC step: lands on 0.625.
//  * 80k-sample dense random bipolar peaks up to |1.0|, many seeds, 2 ms
//    release: max|out| == threshold, zero overshoot, no NaN/Inf. Hard-safe.
//  * Startup: state primed to UNITY for the first L samples (~50 ms) while the
//    lookahead fills; the limiter opens at unity gain rather than fading in.
//  * Release holds until a dip is within the 50 ms window (relCeil), i.e. gain
//    does not recover until the loud part is ~50 ms gone — correct for long LA.
//============================================================================
declare name "lookahead_hermite_limiter";
import("stdfaust.lib");

maxMs    = 50.0;                                  // lookahead in ms (edit me)
SRdesign = 48000.0;                               // design sample rate
L        = int(maxMs * 0.001 * SRdesign + 0.5);   // 2400 smp ~= 50 ms @ 48 kHz
BIG      = 1.0e9;

gAt(k, g) = g @ (L - k);

hermite(t, p0, m0, p1, m1) = h00*p0 + h10*m0 + h01*p1 + h11*m1
  with {
    t2 = t*t; t3 = t2*t;
    h00 =  2*t3 - 3*t2 + 1;  h10 = t3 - 2*t2 + t;
    h01 = -2*t3 + 3*t2;      h11 = t3 -   t2;
  };

// per-span tangent clamped to the monotone region: same sign as the secant
// d=p1-p0, magnitude <= 3|d|. Keeps slope-matching where sane; forbids
// overshoot/undershoot/runaway at any lookahead length.
clampM(m, d) = r * d
  with { z = (d == 0.0); dd = d + z; r = min(max(m/dd, 0.0), 3.0) * (1.0 - z); };
hermiteM(t, p0, m0, p1, m1) = hermite(t, p0, clampM(m0, dd), p1, clampM(m1, dd))
  with { dd = p1 - p0; };

limCore(relTime, g) = (transition ~ si.bus(9)) : (_,!,!,!,!,!,!,!,!)
  with {
    gd   = g - g';
    g0   = gAt(0, g);
    relC = 1.0 - exp(-1.0/(max(relTime,1.0e-6)*ma.SR));

    // ---- nested prefix (min,pos) over the lookahead, sampled at dyadic reaches
    // combine two (val,pos): smaller val wins; tie -> smaller pos (the nearer)
    cmb(v1,p1,v2,p2) = select2(lt,v2,v1), select2(lt,p2,p1) with { lt = v1 <= v2; };
    // balanced reduce over n (val,pos) pairs; left subtree = nearer taps
    redP(1) = _,_;
    redP(n) = (redP(n-h), redP(h)) : cmb  with { h = int(n/2); };
    // (val,pos) over taps a+1 .. a+n
    blockVP(a, n) = par(k, n, (gAt(a+1+k, g), float(a+1+k))) : redP(n);

    nTop  = int(floor(log(L)/log(2.0)));   // 11 for L=2400
    Wr(j) = int(pow(2.0, j));              // 1,2,4,...,2^nTop
    // prefix (min,pos) for dyadic window 2^j, folded from disjoint blocks
    prefVP(0) = blockVP(0, 1);
    prefVP(j) = (prefVP(j-1), blockVP(Wr(j-1), Wr(j)-Wr(j-1))) : cmb;
    // top window reaches the exact maxN = L
    prefTop   = (prefVP(nTop), blockVP(Wr(nTop), L - Wr(nTop))) : cmb;

    // ---- one candidate from a (min,pos) prefix
    candVP(vp, G, v) = prop, mm, dd, vv
      with {
        mm   = vp:(_,!);
        dd   = vp:(!,_);
        vv   = de.delay(L, int(L - dd), gd);   // raw-GR diff across the min tap
        prop = hermiteM(1.0/dd, G, v*dd, mm, vv*dd);
      };
    pick2(pa,ma,da,va, pb,mb,db,vb) =
        select2(c,pb,pa), select2(c,mb,ma), select2(c,db,da), select2(c,vb,va)
      with { c = pa <= pb; };
    bank(G,v,0) = candVP(prefVP(0), G, v);
    bank(G,v,n) = (bank(G,v,n-1), candVP(prefVP(n), G, v)) : pick2;
    allWin(G,v) = (bank(G,v,nTop), candVP(prefTop, G, v)) : pick2;

    transition(G,v,k,inc,P0,V0,P1,V1,c) =
        Gn, vn, kn, incn, P0n, V0n, P1n, V1n, cn
      with {
        cn    = min(c + 1.0, float(L));
        live  = c >= float(L);

        incG = max(inc, 1.0e-6);
        win  = allWin(G, v);
        Pwin = win:(_,!,!,!); minW = win:(!,_,!,!);
        dW   = win:(!,!,_,!); velW = win:(!,!,!,_);

        kc    = k + 1.0;
        tc    = min(kc*inc, 1.0);
        Gcont = hermiteM(tc, P0, V0/incG, P1, V1/incG);

        attacking = min(Pwin, Gcont) < G;
        overtake  = Pwin < Gcont;

        relCeil = min(g0, prefTop:(_,!));   // deepest dip still in lookahead
        Grel    = G + (relCeil - G)*relC;

        Gn0 = select2(attacking, Grel, select2(overtake, Gcont, Pwin));
        GnL = min(Gn0, g0);
        vnL = GnL - G;

        pick3(a,b,c2) = select2(attacking, c2, select2(overtake, b, a));
        Gn   = select2(live, 1.0, GnL);
        vn   = select2(live, 0.0, vnL);
        kn   = select2(live, 1.0, pick3(1.0,    kc,  1.0));
        incn = select2(live, 1.0, pick3(1.0/dW, inc, 1.0));
        P0n  = select2(live, 1.0, pick3(G,      P0,  GnL));
        V0n  = select2(live, 0.0, pick3(v,      V0,  0.0));
        P1n  = select2(live, 1.0, pick3(minW,   P1,  GnL));
        V1n  = select2(live, 0.0, pick3(velW,   V1,  0.0));
      };
  };

thresh = hslider("[0]threshold[unit:lin]", 0.5, 0.01, 1.0, 0.001);
relT   = hslider("[1]release[unit:s]",     0.050, 0.001, 1.0, 0.001);
rawGain(pk) = min(1.0, thresh / max(pk, 1.0e-9));
limiterStereo(l, r) = (l @ L)*gn, (r @ L)*gn
  with { gn = limCore(relT, rawGain(max(abs(l), abs(r)))); };
process = limiterStereo;
